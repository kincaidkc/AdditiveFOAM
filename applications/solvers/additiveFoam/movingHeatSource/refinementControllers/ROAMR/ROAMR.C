/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2023 Oak Ridge National Laboratory                
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "ROAMR.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace refinementControllers
{
    defineTypeNameAndDebug(ROAMR, 0);
    addToRunTimeSelectionTable(refinementController, ROAMR, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementControllers::ROAMR::ROAMR
(
    const PtrList<heatSourceModel>& sources,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    uniformIntervals(sources, dict, mesh, true),

    coeffs_(refinementDict_.optionalSubDict(typeName + "Coeffs")),
    cellsPerProc_(coeffs_.lookupOrDefault<int>("cellsPerProc", 20000))
{
    Info << "Calculating initial AMR interval size..." << endl;

    //- Get average cell volume and cross-sectional area
    label totalCells = mesh_.nCells();
    reduce(totalCells, sumOp<label>());
    scalar vAvg = gSum(mesh_.V()) / totalCells;
    
    //- Find longest path and estimate total scan area
    scalar maxLen = 0.0;
    scalar maxDim = 0.0;
    scalar scanArea = 0.0;

    forAll(sources_, i)
    {
        maxLen = max(sources_[i].beam().length(), maxLen);
        
        treeBoundBox beamBb
        (
            min(vector::zero, boundingBox_.min()),
            max(1.5 * sources_[i].dimensions(), boundingBox_.max())
        );
        
        point bbMin = beamBb.min();
        point bbMax = beamBb.max();
        
        scalar bbMaxDim = max(bbMax[0] - bbMin[0], bbMax[1] - bbMin[1]);

        maxDim = max(maxDim, bbMaxDim);
        
        scanArea += 4.0 * bbMaxDim * Foam::pow(Foam::pow(bbMax[2] - bbMin[2], 2.0), 0.5);
    }

    //- Calculate maximum number of intervals or shortest interval
    //  size to maintain at least 1 bounding box between updates
    scalar maxIntervals = maxLen / maxDim;
    minIntervalTime_ = endTime_ / maxIntervals;

    //- Calculate number of intervals to optimize cells per processor
    scalar targetCells = Pstream::nProcs() * cellsPerProc_;

    if (targetCells > totalCells)
    {
        intervals_ = maxLen * scanArea / vAvg
                     / (targetCells - totalCells)
                     * (Foam::pow(2.0, 3.0 * nLevels_) - 1.0);
    }

    //- Bound number of intervals between 1 and maxIntervals
    intervals_ = max(min(intervals_, maxIntervals), 1.0);

    Info << "Setting initial number of intervals to: " << intervals_ << endl;

    //- Set number of intervals in uniformIntervals class    
    intervalTime_ = endTime_ / intervals_;

    Info << "Setting initial interval time to: " << intervalTime_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::refinementControllers::ROAMR::update(const bool& force)
{
    //- Update if mesh time equals update time
    //  OR if time index is equal to the max refinement level.
    //  This second condition adjusts the mesh after the guess at the first
    //  refinement interval size to prevent an overly long first interval.
    if ((updateTime_ - mesh_.time().value() < small)
        ||
        (mesh_.time().timeIndex() == nLevels_ + 1))
    {
        //- Guard against rescaling until full refinement is reached
        if (mesh_.time().timeIndex() >= nLevels_ + 1)
        {
            //- Scale interval size based on current cells/proc
            label totalCells = mesh_.nCells();
            reduce(totalCells, sumOp<label>());
            scalar currCellsPerProc = totalCells / Pstream::nProcs();

            Info << "Current cells per processor: " << currCellsPerProc << endl;
            Info << "Current interval time: " << intervalTime_ << endl;

            //- Rescale interval time
            scalar scale = min(2.0, max(0.5, cellsPerProc_ / currCellsPerProc));
            
            intervalTime_ *= scale;
            
            //- Ensure interval time is above minimum time
            intervalTime_ = max(intervalTime_, minIntervalTime_);
            
            Info << "New interval time: " << intervalTime_ << endl;
        }

        //- Force update refinement field using uniform intervals functions
        uniformIntervals::update(true);
    }

    return true;
}


bool Foam::refinementControllers::ROAMR::read()
{
    if (uniformIntervals::read())
    {
        refinementDict_ = optionalSubDict(type() + "Coeffs");

        //- Mandatory entries
        refinementDict_.lookup("cellsPerProc") >> cellsPerProc_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
