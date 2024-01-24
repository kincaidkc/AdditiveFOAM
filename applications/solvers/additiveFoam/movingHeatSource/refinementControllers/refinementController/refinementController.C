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

#include "refinementController.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(refinementController, 0);
    defineRunTimeSelectionTable(refinementController, dictionary);
}

const Foam::word Foam::refinementController::refinementControllerDictName
(
    "refinementControllerDict"
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::refinementController::createIOobject
(
    const dictionary& dict,
    const fvMesh& mesh
) const
{
    typeIOobject<IOdictionary> io
    (
        dict.name(),
        mesh.time().constant(),
        mesh.thisDb(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.headerOk())
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementController::refinementController
(
    const word& type,
    const PtrList<heatSourceModel>& sources,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    IOdictionary(createIOobject(dict, mesh)),

    sources_(sources),
    heatSourceDict_(dict),
    mesh_(mesh),
    refinementDict_(heatSourceDict_.optionalSubDict(type + "Coeffs")),
    refine_(refinementDict_.lookup<bool>("refine")),
    refinementField_
    (
        IOobject
        (
            "refinementField",
            mesh_.time().timeName(),
            mesh_,
            refine_ ? IOobject::READ_IF_PRESENT : IOobject::NO_READ,
            refine_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, 0.0)
    ),
    resolveTail_(refine_ ? refinementDict_.lookupOrDefault<bool>("resolveTail", false) : false),
    persistence_(refine_ ? refinementDict_.lookupOrDefault<scalar>("persistence", 0.0) : 0.0),
    solidificationTime_
    (
        IOobject
        (
            "solidificationTime",
            mesh_.time().timeName(),
            mesh_,
            persistence_ > 0.0
                ? IOobject::READ_IF_PRESENT : IOobject::NO_READ,
            persistence_ > 0.0
                ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTime, GREAT)
    )
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::refinementController::initializeRefinementField()
{
    //- Initialize refinement field to capture melt pool
    if (resolveTail_)
    {
        const volScalarField alphaSol
            = mesh_.lookupObject<volScalarField>("alpha1");
            
        const volScalarField alphaSol0 = alphaSol.oldTime();
        
        refinementField_ = pos0(1.0 - alphaSol);
        
        if (persistence_ > 0.0)
        {
            const scalar currTime = mesh_.time().value();
        
            //- Update solidification time field
            forAll(mesh_.C(), celli)
            {
                if ((alphaSol[celli] > 0.99) && (alphaSol0[celli] < 0.99))
                {
                    //- Update solidification time field
                    solidificationTime_[celli] = currTime;
                }
                
                //- Ensure recently solidified regions stay refined
                //  for persistence time
                if (solidificationTime_[celli] + persistence_ < currTime)
                {
                    refinementField_[celli] = 1.0;
                }
            }
        }        
    }
    
    //- Reset refinement field if no tail capturing is desired
    else
    {
        refinementField_ = dimensionedScalar(dimless, 0.0);
    }
}

bool Foam::refinementController::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
