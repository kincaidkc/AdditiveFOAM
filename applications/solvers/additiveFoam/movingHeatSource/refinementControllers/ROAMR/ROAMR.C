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
    cellsPerProc_(coeffs_.lookupOrDefault<int>("cellsPerProc", 10000)),
    lastRefi_(0)
{
    //- Calculate total mesh volume
    scalar Vtot = 0.0;
    forAll(mesh_.C(), celli)
    {
        Vtot += mesh_.V()[celli];
    }
    reduce(Vtot, sumOp<scalar>());
    Info << "Total mesh volume: " << Vtot << endl;
    
    //- Find average cell dimension
    scalar dx = Foam::pow(Vtot / mesh_.nCells(), 1.0 / 3.0);
    
    //- Calculate the bounding box for each cell
    List<treeBoundBox> cellBbs(mesh_.nCells());    
    const pointField& points = mesh_.points();
    const vector extend = 1e-10 * vector::one;
    
    forAll(mesh_.cells(), celli)
    {
        treeBoundBox cellBb(point::max, point::min);

        const labelList& vertices = mesh_.cellPoints()[celli];

        forAll(vertices, j)
        {
            cellBb.min()
                = min(cellBb.min(), points[vertices[j]] - extend);
            cellBb.max()
                = max(cellBb.max(), points[vertices[j]] + extend);
        }

        cellBbs[celli] = cellBb;
    }
    
    //- Calculate scan path volume
    forAll(sources_, i)
    {
        const movingBeam& beam_ = sources_[i].beam();
        
        scalar time_ = 0.0;

        vector offset_ = 1.5*sources_[i].dimensions();

        while (time_ < beam_.endTime())
        {
            vector position_ = beam_.position(time_);

            treeBoundBox beamBb
            (
                position_ - min(offset_, boundingBox_.min()),
                position_ + max(offset_, boundingBox_.max())
            );
            
            forAll(mesh_.cells(), celli)
            {
                if (refinementField_[celli] > 0)
                {
                    // Do nothing, cell already marked for refiment
                }
                else if (cellBbs[celli].overlaps(beamBb))
                {
                    refinementField_[celli] = 1;
                }
            }
            
            //- Calculate time step required to resolve beam motion on mesh
            label index_ = beam_.findIndex(time_);
            segment path_ = beam_.getSegment(index_);
            scalar timeToNextPath_ = path_.time() - time_;

            //- If the path end time is directly hit, step to next path
            while (mag(timeToNextPath_) < small)
            {
                index_ = index_ + 1;
                path_ = beam_.getSegment(index_);
                timeToNextPath_ = path_.time() - time_;
            }

            scalar dt_ = timeToNextPath_;

            if (path_.mode() == 0)
            {
                const scalar scanTime_ =
                    sources_[i].D2sigma() / path_.parameter();

                dt_ = min(timeToNextPath_, scanTime_);
            }

            time_ += dt_;
        }
    }
    
    scalar Vref = fvc::domainIntegrate(refinementField_).value();
    Info << "Total volume of refined scan path: " << Vref << endl;
    
    //- Calculate number of cells required to fully resolve entire scan path
    scalar nTot1 = (Vtot + Vref * (Foam::pow(2.0, 3.0 * nLevels_) - 1.0))
                   / Foam::pow(dx, 3.0);
            
    //- Calculate target mesh size
    scalar targetCells = Pstream::nProcs() * cellsPerProc_;
    
    //- Estimate number of intervals
    intervals_ = max(nTot1 / targetCells, 1.0);
    
    //- TODO: set minimum interval size?
    
    Info << "Estimated that " << intervals_
         << " intervals are required." << endl;
         
    //- Set number of intervals in uniformIntervals class    
    intervalTime_ = endTime_ / intervals_;

    Info << "Setting initial interval time to: " << intervalTime_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::refinementControllers::ROAMR::update(const bool& force)
{
    //- Check number of cells in mesh against desired number of cells
    bool forceUpdate = false;
    const label nCells = returnReduce(mesh_.nCells(), sumOp<label>());
    const scalar ratio = nCells / (cellsPerProc_ * Pstream::nProcs());

    if (((ratio > 1.5) || (ratio < 0.5))
        && 
        (mesh_.time().timeIndex() >= (lastRefi_ + nLevels_ + 1)))
    {
        Info << "Mesh contains " << nCells << " cells, which is " << ratio
             << " desired number of cells. Forcing mesh update..." << endl;
        forceUpdate = true;
    }

    //- Update if mesh time equals update time or if the forceUpdate flag is
    //  set due to an improperly sized mesh
    if ((updateTime_ - mesh_.time().value() < small) || (forceUpdate))
    {
        //- Guard against rescaling until full refinement is reached
        if (mesh_.time().timeIndex() >= (lastRefi_ + nLevels_ + 1))
        {
            //- Scale interval size based on current cells/proc
            label totalCells = mesh_.nCells();
            reduce(totalCells, sumOp<label>());
            scalar currCellsPerProc = totalCells / Pstream::nProcs();

            Info << "Current cells per processor: "
                 << currCellsPerProc << endl;
            Info << "Current interval time: " << intervalTime_ << endl;

            //- Rescale interval time
            scalar scale =
                min(2.0, max(0.5, cellsPerProc_ / currCellsPerProc));
            
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
