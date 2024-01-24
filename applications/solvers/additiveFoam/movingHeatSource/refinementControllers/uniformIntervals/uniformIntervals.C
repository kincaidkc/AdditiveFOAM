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

#include "uniformIntervals.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace refinementControllers
{
    defineTypeNameAndDebug(uniformIntervals, 0);
    addToRunTimeSelectionTable(refinementController, uniformIntervals, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementControllers::uniformIntervals::uniformIntervals
(
    const PtrList<heatSourceModel>& sources,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    refinementController(typeName, sources, dict, mesh)
{
    intervals_ = refinementDict_.lookup<int>("intervals");
    dt_ = mesh_.time().endTime().value() / intervals_;
    updateTime_ = 0.0;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::refinementControllers::uniformIntervals::update()
{
    if (mesh_.time().value() >= updateTime_)
    {
        //- Initialize refinement marker field in base class
        refinementController::initializeRefinementField();
    
        //- Update next refinement time
        const scalar currTime = mesh_.time().value();
        updateTime_ = currTime + dt_;
        
        //- Buffer for bounding box calculations
        const vector extend = 1e-10 * vector::one;
        
        const pointField& points = mesh_.points();
        
        //- Update refinement marker field
        forAll(sources_, i)
        {        
            //- Create dummy time for forward march scheme
            scalar beamTime = currTime;
            
            //- March along beam path until next update time is reached
            while (beamTime <= updateTime_)
            {
                //- Initialize time step size and velocity
                scalar beamDt = 0.0;
                
                //- Calculate max beam timestep based on velocity for
                //  moving beams, or time until next path for stationary
                //scalar beamMode = sources_[i].beam().mode(beamTime);
                scalar beamMode = sources_[i].beam().mode(beamTime);
                
                if (beamMode == 1.0)
                {
                    beamDt = sources_[i].D2sigma()
                             / sources_[i].beam().parameter(beamTime);
                }
                
                else
                {
                    beamDt = sources_[i].beam().timeToNextPath(beamTime);
                }
                
                //- Ensure final time step refines at next update time
                beamTime = min(beamTime + beamDt, updateTime_);
                
                //- Mark cells near beam at current time/location for
                //  refinement
                vector beamDims = sources_[i].dimensions();
                vector currPos = sources_[i].beam().position(beamTime);
                treeBoundBox beamBb(currPos - beamDims, currPos + beamDims);
                
                forAll(mesh_.C(), celli)
                {
                    //- Don't redo checks if cell is already marked
                    if (refinementField_[celli] > 0.0)
                    {
                        continue;
                    }
                    
                    treeBoundBox cellBb(point::max, point::min);
                    
                    const labelList& vertices = mesh_.cellPoints()[celli];
                    
                    forAll(vertices, j)
                    {
                        cellBb.min()
                            = min(cellBb.min(), points[vertices[j]] - extend);
                        cellBb.max()
                            = max(cellBb.max(), points[vertices[j]] + extend);
                    }
                    
                    if (cellBb.overlaps(beamBb))
                    {
                        refinementField_[celli] = 1.0;
                    }
                }
            }            
        }
        
        refinementField_ = pos(refinementField_);
        
        //- Return true to call mesh.update()
        return true;
    }
    
    else
    {
        return false;
    }
}


bool Foam::refinementControllers::uniformIntervals::read()
{
    if (refinementController::read())
    {
        refinementDict_ = optionalSubDict(type() + "Coeffs");

        //- Mandatory entries
        refinementDict_.lookup("intervals") >> intervals_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
