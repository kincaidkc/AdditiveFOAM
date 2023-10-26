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

Application
    additiveFoam

Description
    Transient solver for additive manufacturing simulations
    
\*---------------------------------------------------------------------------*/

#include "Timer.H"

#include "pimpleControl.H"

#include "Polynomial.H"

#include "interpolateXY/interpolateXY.H"
#include "graph/graph.H"

#include "movingHeatSourceModel.H"
#include "foamToExaCA/foamToExaCA.H"

#include "argList.H"
#include "timeSelector.H"
#include "zeroGradientFvPatchFields.H"
#include "IFstream.H"
#include "uniformDimensionedFields.H"
#include "pressureReference.H"
#include "findRefCell.H"

#include "fvmDiv.H"
#include "fvmDdt.H"
#include "fvmLaplacian.H"
#include "constrainPressure.H"
#include "constrainHbyA.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    using namespace Foam;

    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createDyMControls.H"
    //#include "createControl.H"
    #include "createFields.H"
    //#include "createTimeControls.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Timers timer(runTime);
    
    // initialize time-stepping controls
    scalar FoNum = 0.0;

    scalar alphaCoNum = 0.0;

    foamToExaCA ExaCA(T);

    movingHeatSourceModel sources(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
    
        #include "updateProperties.H"

        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"
        
        // Update the mesh for topology change, mesh to mesh mapping
        timer.start("Mesh Update");
        mesh.update();
        timer.stop("Mesh Update");

        timer.start("Heat Source");
        sources.update();
        timer.stop("Heat Source");
        
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "solutionControls.H"
        
        while (pimple.loop() && fluidInDomain)
        {
            if (pimple.firstPimpleIter() || pimple.moveMeshOuterCorrectors())
            {
                // Move the mesh
                mesh.move();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        //#include "correctPhi.H"
                    }

                    if (checkMeshCourantNo)
                    {
                        //#include "meshCourantNo.H"
                    }
                }
            }
            
            timer.start("pUEqn");
            #include "pU/UEqn.H"
            #include "pU/pEqn.H"
            timer.stop("pUEqn");
        }

        timer.start("TEqn");
        #include "thermo/TEqn.H"
        timer.stop("TEqn");

        timer.start("ExaCA Update");
        ExaCA.update();
        timer.stop("ExaCA Update");

        runTime.write();

        if (runTime.writeTime())
        {
            sources.qDot().write();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    timer.start("ExaCA Write");
    ExaCA.write();
    timer.stop("ExaCA Write");
    
    timer.write();

    return 0;
}

// ************************************************************************* //
