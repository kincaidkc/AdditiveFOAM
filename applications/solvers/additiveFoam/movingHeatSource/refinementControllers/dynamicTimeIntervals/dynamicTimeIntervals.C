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

#include "dynamicTimeIntervals.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace refinementControllers
{
    defineTypeNameAndDebug(dynamicTimeIntervals, 0);
    addToRunTimeSelectionTable
    (
        refinementController,
        dynamicTimeIntervals,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementControllers::dynamicTimeIntervals::dynamicTimeIntervals
(
    const PtrList<heatSourceModel>& sources,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    refinementController(typeName, sources, dict, mesh),
    coeffs_(refinementDict_.optionalSubDict(typeName + "Coeffs")),
    cellsPerProc_(coeffs_.lookupOrDefault<int>("cellsPerProc", 10000)),
    relax_(coeffs_.lookupOrDefault<scalar>("relax", 0.9)),
    minIntervalTime_(0.0),
    intervalLength_(0.0),
    updateTime_(0.0)
{
    //- Estimate the length of the first refinement interval by using the swept
    //  volume of the beam(s) and refined cell volume to guess the mesh size
    //  associated with a given volume.
    
    //- Get average un-refined cell volume and cross-sectional area
    label totalCells = mesh_.nCells();
    reduce(totalCells, sumOp<label>());
    scalar vAvg = gSum(mesh_.V()) / totalCells;
    
    //- Find longest path and estimate total scan area
    scalar maxLen = 0.0;
    scalar maxDim = 0.0;
    scalar scanArea = 0.0;

    forAll(sources_, i)
    {
        maxLen = max(sources_[i].beam().totalLength(), maxLen);
        
        treeBoundBox beamBb
        (
            min(-1.5 * sources_[i].dimensions(), -buffer_),
            max(1.5 * sources_[i].dimensions(), buffer_)
        );
        
        point bbMin = beamBb.min();
        point bbMax = beamBb.max();
        
        scalar bbMaxDim = max(bbMax[0] - bbMin[0], bbMax[1] - bbMin[1]);

        maxDim = max(maxDim, bbMaxDim);
        
        scanArea +=
            4.0 * bbMaxDim
          * Foam::pow(Foam::pow(bbMax[2] - bbMin[2], 2.0), 0.5);
    }
    
    //- Calculate maximum number of intervals or shortest interval size so
    //  that each AMR interval will refine a distance of at least the beam
    //  bounding box dimension.
    scalar maxIntervals = maxLen / maxDim;
    minIntervalTime_ = endTime_ / maxIntervals;

    //- Calculate number of intervals to reach target cells per processor
    scalar targetCells = Pstream::nProcs() * cellsPerProc_;
    
    scalar intervals = 0.0;

    if (targetCells > totalCells)
    {
        intervals =
            maxLen * scanArea / vAvg
          / (targetCells - totalCells)
          * (Foam::pow(2.0, 3.0 * nLevels_) - 1.0);
    }
    
    //- Bound number of intervals between 1 and maxIntervals
    intervals = max(min(intervals, maxIntervals), 1.0);
    
    //- Set number of intervals in uniformIntervals class    
    intervalLength_ = endTime_ / intervals;
    
    Info << "dynamicTimeIntervals: set first interval to " << intervalLength_
         << " s, which corresponds to approx. " << intervals
         << " total intervals." << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::refinementControllers::dynamicTimeIntervals::update()
{
    //- Update if mesh time equals update time
    //  OR if time index is equal to the max refinement level.
    //  This second condition adjusts the mesh after the guess at the first
    //  refinement interval size to prevent an overly long first interval.
    if ((updateTime_ - mesh_.time().value() < small)
      ||
        (mesh_.time().timeIndex() == nLevels_ + 1))
    {
        //- Refine in regions above specified temperature
        refinementController::refineUsingTemperature();

        //- Don't perform additional refinements if scan path is completed
        if ((endTime_ - mesh_.time().value()) < small)
        {
            Info << "dynamicTimeIntervals: Scan path completed. Continuing AMR"
                 << " checks for possible mesh coarsening" << endl;
                 
            updateTime_ = mesh_.time().value() + intervalLength_;
                 
            return true;
        }
        
        //- Calculate current CPU load (cells per processor)
        label nCells = mesh_.nCells();
        reduce(nCells, sumOp<label>());
        scalar currCellsPerProc = nCells / Pstream::nProcs();
        
        Info << "dynamicTimeIntervals: Current CPU load is "
             << currCellsPerProc
             << " cells per processor, with an interval of length "
             << intervalLength_ << " s." << endl;
        
        //- Rescale interval length
        intervalLength_
            = relax_ * cellsPerProc_ / currCellsPerProc * intervalLength_
              + (1.0 - relax_) * intervalLength_;
        
        //- Update next refinement time
        updateTime_ = mesh_.time().value() + intervalLength_;
        
        Info << "dynamicTimeIntervals: rescaled interval to "
             << intervalLength_ << " s. Next update will occur at " 
             << updateTime_ << "s. Updating AMR marker field." << endl;
        
        //- Update marker field using calculated update time
        refinementController::refineUsingTime(updateTime_);
    }

    return true;
}


bool Foam::refinementControllers::dynamicTimeIntervals::read()
{
    if (refinementController::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
