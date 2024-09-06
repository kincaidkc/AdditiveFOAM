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

#include "dynamicIntervals.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace refinementControllers
{
    defineTypeNameAndDebug(dynamicIntervals, 0);
    addToRunTimeSelectionTable
    (
        refinementController,
        dynamicIntervals,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementControllers::dynamicIntervals::dynamicIntervals
(
    const PtrList<heatSourceModel>& sources,
    const dictionary& dict,
    const fvMesh& mesh,
    const bool& roamr
)
:
    refinementController(typeName, sources, dict, mesh),
    coeffs_(roamr ? refinementDict_.optionalSubDict("ROAMRCoeffs")
                  : refinementDict_.optionalSubDict(typeName + "Coeffs")),
    boundingBox_
    (
        coeffs_.lookupOrDefault<boundBox>
        (
            "boundingBox",
            boundBox(point::max, point::min)
        )
    ),
    cellsPerProc_(coeffs_.lookupOrDefault<label>("cellsPerProc", 5000)),
    hCoarse_(coeffs_.lookup<scalar>("hCoarse")),
    totalCells_(cellsPerProc_ * Pstream::nProcs()),
    nCells0_(0),
    updateTime_(0.0),
    endTime_(0.0)
{
    //- Set AMR update end time to minimum of solution time and max beam time
    forAll(sources_, i)
    {
        endTime_ = max(sources_[i].beam().endTime(), endTime_);
    }

    endTime_ = min(endTime_, mesh.time().endTime().value());
    
    nCells0_ = mesh_.nCells();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::refinementControllers::dynamicIntervals::update(const bool& force)
{
    if (updateTime_ - mesh_.time().value() < small)
    {
        //- Mark cells above refinement temperature
        refinementController::setRefinementField();
        
        if ((endTime_ - mesh_.time().value()) < small)
        {
            Info << "Continuing AMR checks for possible mesh coarsening" << endl;
            return true;
        }
        
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
        
        //- Number of cells added per cell refined
        scalar cellsAdded = Foam::pow(2.0, 3.0 * nLevels_) - 1.0;
        
        //- Mark current positions of all beams for refinement
        forAll(sources_, i)
        {
            const movingBeam& beam_ = sources_[i].beam();
            
            scalar time_ = mesh_.time().value();

            vector offset_ = 1.5*sources_[i].dimensions();

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
            
            refinementField_.correctBoundaryConditions();
        }
        
        //- Estimate number of cells required to refine regions marked 
        //  based on temperature and current beam position
        scalar vRef = fvc::domainIntegrate(refinementField_).value();
        scalar nRef = vRef / Foam::pow(hCoarse_, 3.0) * cellsAdded;
                        
        scalar nTot = nCells0_ + nRef;
        
        scalar time_ = mesh_.time().value();
        
        //- Mark cells ahead of beam for refinement until estimated mesh size
        //  is equal to desired mesh size
        while (nTot < totalCells_)
        {
            //- Find minimum time step for all beams
            scalar dt_ = 0.0;
            
            forAll(sources_, i)
            {
                const movingBeam& beam_ = sources_[i].beam();
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

                dt_ = timeToNextPath_;

                if (path_.mode() == 0)
                {
                    const scalar scanTime_ =
                        sources_[i].D2sigma() / path_.parameter();

                    dt_ = min(dt_, scanTime_);
                }
            }
            
            Info << "Pseudo time step size: " << dt_ << endl;
            
            //- Increment psuedo time
            time_ += dt_;
            
            Info << "Pseudo time: " << time_ << endl;
            
            if (time_ > endTime_)
            {
                break;
            }
            
            //- Update refinement marker field for all beams at new time
            forAll(sources_, i)
            {
                const movingBeam& beam_ = sources_[i].beam();
                
                if (time_ <= beam_.endTime())
                {
                    vector offset_ = 1.5*sources_[i].dimensions();

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
                            vRef += mesh_.V()[celli];
                        }
                    }
                }
            }
            
            nRef = vRef / Foam::pow(hCoarse_, 3.0) * cellsAdded;
            nTot = nCells0_ + nRef;
            
            Info << "Estimated mesh size: " << nTot << endl;
            
            refinementField_.correctBoundaryConditions(); 
        }
        
        //- Set update time equal to beam pseudo time at last check
        updateTime_ = time_;
    }

    return true;
}


bool Foam::refinementControllers::dynamicIntervals::read()
{
    if (refinementController::read())
    {
        refinementDict_ = optionalSubDict(type() + "Coeffs");

        //- Mandatory entries
        refinementDict_.lookup("boundingBox") >> boundingBox_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
