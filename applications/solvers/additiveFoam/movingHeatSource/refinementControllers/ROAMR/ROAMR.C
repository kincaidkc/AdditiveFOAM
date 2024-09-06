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
    cellsPerProc_(coeffs_.lookupOrDefault<int>("cellsPerProc", 10000))
{
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
        treeBoundBox beamBb
        (
            min(vector::zero, boundingBox_.min()),
            max(1.5 * sources_[i].dimensions(), boundingBox_.max())
        );
        
        point bbMin = beamBb.min();
        point bbMax = beamBb.max();
        
        scalar bbMaxDim = max(bbMax[0] - bbMin[0], bbMax[1] - bbMin[1]);

        maxDim = max(maxDim, bbMaxDim);
        
        scanArea +=
            4.0 * bbMaxDim
          * Foam::pow(Foam::pow(bbMax[2] - bbMin[2], 2.0), 0.5);
          
        //- Calculate beam scan length + effective spot melt length
        scalar beamLen =
            sources_[i].beam().totalLength()
          + sources_[i].beam().totalSpots() * bbMaxDim;
        
        maxLen = max(beamLen, maxLen);
    }

    //- Calculate maximum number of intervals or shortest interval
    //  size to maintain at least 1 bounding box between updates
    scalar maxIntervals = maxLen / maxDim;
    minIntervalTime_ = endTime_ / maxIntervals;

    //- Calculate number of intervals to optimize cells per processor
    scalar targetCells = Pstream::nProcs() * cellsPerProc_;

    if (targetCells > totalCells)
    {
        intervals_ =
            maxLen * scanArea / vAvg
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
    //- Calculate ratio of actual cells in mesh to target mesh size
    const scalar nCells = returnReduce(mesh_.nCells(), sumOp<label>());
    const scalar ratio = nCells / (cellsPerProc_ * Pstream::nProcs());
   
    //- Force an update if ratio is off by +/- 20% or more 
    bool forceUpdate = false;

    if ((ratio > 1.2) || (ratio < 0.8))
    {
        //- Only force update if the current interval is larger than the
        //  calculated minimum interval size
        if (intervalTime_ > (1.05 * minIntervalTime_))
        {
            //- Don't force update as the mesh size decreases at the end
            //  of the scan path
            if ((ratio < 0.8) && (updateTime_ >= endTime_))
            {
                Info << "Mesh contains fewer cells than optimal, "
                     << "but update time is larger than path end time. "
                     << "Not forcing refinement." << endl;
            }
            else
            {
                Info << "Mesh contains " << nCells << " cells, which is "
                     << ratio << " times the desired number of cells. "
                     << "Forcing mesh update..." << endl;

                forceUpdate = true;
            }
        }
    }
    
    //- Update if mesh time equals update time or if the forceUpdate flag is
    //  set due to an improperly sized mesh
    if ((updateTime_ - mesh_.time().value() < small) || (forceUpdate))
    {
        //- Guard against rescaling until full refinement is reached
        if (mesh_.time().timeIndex() >= nLevels_ + 1)
        {
            //- Scale interval size based on current cells/proc
            scalar totalCells = mesh_.nCells();
            reduce(totalCells, sumOp<scalar>());
            scalar currCellsPerProc = totalCells / Pstream::nProcs();

            Info << "Current cells per processor: "
                 << currCellsPerProc << endl;
            Info << "Current interval time: " << intervalTime_ << endl;

            //- Rescale interval time
            scalar scale =
                min(2.0, max(0.5, cellsPerProc_ / currCellsPerProc));

            Info << "Scale factor: " << scale << endl;
            
            intervalTime_ *= scale;
            
            //- Ensure interval time is above minimum time
            intervalTime_ = min(max(intervalTime_, minIntervalTime_), endTime_ - mesh_.time().value());
            
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
