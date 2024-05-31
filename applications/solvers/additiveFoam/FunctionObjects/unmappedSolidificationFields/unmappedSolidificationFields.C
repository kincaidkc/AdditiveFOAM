/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "unmappedSolidificationFields.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvc.H"
#include "OSspecific.H"
#include "zeroField.H"

#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "polyDistributionMap.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(unmappedSolidificationFields, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        unmappedSolidificationFields,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::unmappedSolidificationFields::correct()
{
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::unmappedSolidificationFields::unmappedSolidificationFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    dict_(dict),
    AMR_(dict_.lookup<bool>("AMR")),
    refinedSize_(0.0),
    Tl_(dict_.lookup<scalar>("Tl")),
    T_(mesh_.lookupObject<VolField<scalar>>("T")),
    R_
    (
        IOobject
        (
            "R",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTemperature/dimTime, 0.0)
    )
{
    read(dict);
    
    if (AMR_)
    {
        refinedSize_ = dict_.lookup<scalar>("refinedSize");
    }

    R_ = fvc::ddt(T_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::unmappedSolidificationFields::~unmappedSolidificationFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::unmappedSolidificationFields::read(const dictionary& dict)
{
    return true;
}


Foam::wordList Foam::functionObjects::unmappedSolidificationFields::fields() const
{
    return wordList::null();
}


bool Foam::functionObjects::unmappedSolidificationFields::execute()
{
    //- Get current time
    const scalar& time = mesh_.time().value();
    
    //- Get old temperature and time derivative of temperature
    const volScalarField& T0 = T_.oldTime();
    
    //- Update cooling rate
    R_.oldTime();
    R_ = fvc::ddt(T_);

    //- Calculate the thermal gradient
    const volScalarField G = mag(fvc::grad(T_));
    
    //- Calculate refined mesh volume
    const scalar Vr = Foam::pow(1.1 * refinedSize_, 3.0);
    
    label nEvents = 0;
    
    forAll(mesh_.C(), celli)
    {
        //- Check that cells are in max refinement level for AMR cases
        bool maxLevel = true;
        
        if (AMR_)
        {
            if (mesh_.V()[celli] > Vr)
            {
                maxLevel = false;
            }
        }
        
        //- Check for solidification events
        if ((T0[celli] > Tl_) && (T_[celli] <= Tl_) && (maxLevel))
        {
            vector C = mesh_.C()[celli];

            List<scalar> eventi(6);

            eventi[0] = C[0];
            eventi[1] = C[1];
            eventi[2] = C[2];
            eventi[3] = time;
            eventi[4] = R_.oldTime()[celli];
            eventi[5] = G[celli];

            events_.append(eventi);
            
            ++nEvents;
        }
    }
    
    reduce(nEvents, sumOp<label>());
    
    Info << "unmappedSolidificationFields recorded " << nEvents 
         << " events." << endl;
    
    //- Clean cells which have been remelted (not yet implemented)
    //removeRemelts();
    
    return true;
}


bool Foam::functionObjects::unmappedSolidificationFields::end()
{
    return true;
}


bool Foam::functionObjects::unmappedSolidificationFields::write()
{
    //- Get current time
    const fileName currTime = Foam::name(mesh_.time().value());
    
    //- Create file path for current time data
    const fileName currTimePath(mesh_.time().rootPath()
                                /mesh_.time().globalCaseName()
                                /"unmappedSolidificationFields"/currTime);
    
    mkDir(currTimePath);
    
    //- Open file for each proc
    OFstream os(currTimePath + "/" + "data_" 
                + Foam::name(Pstream::myProcNo()) + ".csv");
    
    //- Write header
    os << "x,y,z,t,R,G\n";
    
    //- Write each event in series to file
    for (int i = 0; i < events_.size(); ++i)
    {
        int n = events_[i].size() - 1;
        
        for (int j = 0; j < n; ++j)
        {
            os << events_[i][j] << ",";
        }
        
        os << events_[i][n] << "\n";
    }
    
    return true;
}


void Foam::functionObjects::unmappedSolidificationFields::movePoints(const polyMesh& mesh)
{
    if (&mesh == &mesh_)
    {
        correct();
    }
}


void Foam::functionObjects::unmappedSolidificationFields::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        correct();
    }
}


void Foam::functionObjects::unmappedSolidificationFields::mapMesh(const polyMeshMap& map)
{
    if (&map.mesh() == &mesh_)
    {
        correct();
    }
}


void Foam::functionObjects::unmappedSolidificationFields::distribute
(
    const polyDistributionMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        correct();
    }
}

// ************************************************************************* //
