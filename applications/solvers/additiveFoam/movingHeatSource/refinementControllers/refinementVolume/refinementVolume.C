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

#include "refinementVolume.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace refinementControllers
{
    defineTypeNameAndDebug(refinementVolume, 0);
    addToRunTimeSelectionTable
    (
        refinementController,
        refinementVolume,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementControllers::refinementVolume::refinementVolume
(
    const PtrList<heatSourceModel>& sources,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    refinementController(typeName, sources, dict, mesh),
    coeffs_(refinementDict_.optionalSubDict(typeName + "Coeffs")),
    cellsPerProc_(coeffs_.lookupOrDefault<int>("cellsPerProc", 10000)),
    unrefinedSize_(dimLength, coeffs_.lookup<scalar>("unrefinedSize")),
    refVol_(dimVolume, 0.0),
    updateTime_(dimTime, 0.0)
{
    //- Estimate the volume of the refined region corresponding to the target
    //  mesh size, calculated from the number of CPUs and target cells per CPU
    scalar targetCells = cellsPerProc_ * Pstream::nProcs();
    
    //- Find initial mesh size
    scalar nCells0 = mesh_.nCells();
    reduce(nCells0, sumOp<scalar>());
    
    //- Provide warning if target mesh size is smaller than initial mesh size
    if (nCells0 > targetCells)
    {
        Info << "refinementVolume: WARNING - initial mesh size larger than "
                "target mesh size." << endl;
    }
    //- Otherwise, estimate refined volume required to hit target mesh size
    else
    {
        refVol_ = Foam::pow(unrefinedSize_, 3.0) * (targetCells - nCells0)
                  / (Foam::pow(2.0, 3.0 * nLevels_) - 1.0);
                  
        Info << "refinementVolume: estimated refinement volume is "
             << refVol_.value() << " m^3" << endl;
    }    
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::refinementControllers::refinementVolume::update()
{
    //- Update if past update time
    if (updateTime_.value() - mesh_.time().value() < small)
    {
        //- Refine in regions above specified temperature
        refinementController::refineUsingTemperature();

        //- Don't perform additional refinements if scan path is completed
        if ((endTime_ - mesh_.time().value()) < small)
        {
            Info << "refinementVolume: Scan path completed. Continuing AMR"
                 << " checks for possible mesh coarsening" << endl;
                 
            //- TODO: increase updateTime_ by some increment
                 
            return true;
        }
        
        //- Calculate current CPU load (cells per processor)
        label nCells = mesh_.nCells();
        reduce(nCells, sumOp<label>());
        scalar currCellsPerProc = nCells / Pstream::nProcs();
        
        Info << "refinementVolume: Current CPU load is "
             << currCellsPerProc << " cells per processor." << endl;
        
        //- Update marker field using refinement volume strategy
        updateTime_ = refinementController::refineUsingVolume(refVol_);
    }

    return true;
}


bool Foam::refinementControllers::refinementVolume::read()
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
