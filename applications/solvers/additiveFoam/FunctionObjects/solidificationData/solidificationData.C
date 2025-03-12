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

#include "solidificationData.H"
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
    defineTypeNameAndDebug(solidificationData, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        solidificationData,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::solidificationData::correct()
{
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::solidificationData::solidificationData
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    dict_(dict),
    AMR_(dict_.lookupOrDefault<bool>("AMR", false)),
    refinedSize_(0.0),
    trackMelting_(dict_.lookupOrDefault<bool>("trackMelting", false)),
    Tl_(dict_.lookup<scalar>("Tl")),
    T_(mesh_.lookupObject<VolField<scalar>>("T")),
    dTdt_
    (
        IOobject
        (
            "dTdt",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTemperature/dimTime, 0.0)
    ),
    R_
    (
        IOobject
        (
            "R",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            AMR_ ? IOobject::NO_WRITE : IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTemperature/dimTime, 0.0)
    ),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            AMR_ ? IOobject::NO_WRITE : IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTemperature/dimLength, 0.0)
    ),
    tm_
    (
        IOobject
        (
            "tm",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            AMR_ ? IOobject::NO_WRITE
                 : (trackMelting_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE)
        ),
        mesh_,
        dimensionedScalar(dimTime, 0.0)
    ),
    ts_
    (
        IOobject
        (
            "ts",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            AMR_ ? IOobject::NO_WRITE
                 : (trackMelting_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE)
        ),
        mesh_,
        dimensionedScalar(dimTime, 0.0)
    )
{
    read(dict);

    if (AMR_)
    {
        refinedSize_ = dict_.lookup<scalar>("refinedSize");
    }

    dTdt_ = fvc::ddt(T_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::solidificationData::~solidificationData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::solidificationData::read(const dictionary& dict)
{
    return true;
}


Foam::wordList Foam::functionObjects::solidificationData::fields() const
{
    return wordList::null();
}


bool Foam::functionObjects::solidificationData::execute()
{
    //- Get current time
    const scalar& time = mesh_.time().value();
    const scalar& deltaT = mesh_.time().deltaTValue();

    //- Get old temperature and time derivative of temperature
    const volScalarField& T0_ = T_.oldTime();

    //- Update cooling rate
    dTdt_.oldTime();
    dTdt_ = fvc::ddt(T_);

    //- Calculate the thermal gradient
    const volScalarField gradT("Gtmp", mag(fvc::grad(T_)));

    //- Calculate refined mesh volume
    const scalar Vr = Foam::pow(1.1 * refinedSize_, 3.0);

    label nSolidificationEvents = 0;
    label nMeltingEvents = 0;

    forAll(mesh_.C(), celli)
    {
        //- Check that cells are in max refinement level for AMR cases
        if (AMR_)
        {
            if (mesh_.V()[celli] > Vr)
            {
                continue;
            }
        }
        
        const scalar T = T_[celli];
        const scalar T0 = T0_[celli];

        //- Check for solidification events
        if ((T0 > Tl_) && (T <= Tl_))
        {
            if (AMR_)
            {
                vector C = mesh_.C()[celli];

                List<scalar> eventi(6);

                eventi[0] = C[0];
                eventi[1] = C[1];
                eventi[2] = C[2];
                eventi[3] = time;
                eventi[4] = dTdt_.oldTime()[celli];
                eventi[5] = gradT[celli];

                solidificationEvents_.append(eventi);
            }
            else
            {
                R_[celli] = dTdt_.oldTime()[celli];
                G_[celli] = gradT[celli];
                
                if (trackMelting_)
                {
                    ts_[celli] = time - deltaT * min(max((Tl_ - T) / (T0 - T), 0), 1);
                }
            }

            ++nSolidificationEvents;
        }
        
        //- Check for melting events if trackMelting is set
        if (trackMelting_ && (T0 < Tl_) && (T >= Tl_))
        {
            if (AMR_)
            {
                vector C = mesh_.C()[celli];
                
                List<scalar> eventi(4);
                
                eventi[0] = C[0];
                eventi[1] = C[1];
                eventi[2] = C[2];
                eventi[3] = time;
                
                meltingEvents_.append(eventi);
            }
            else
            {
                tm_[celli] = time - deltaT * min(max((T - Tl_) / (T - T0), 0), 1);
            }
            
            ++nMeltingEvents;
        }
    }

    reduce(nSolidificationEvents, sumOp<label>());
    reduce(nMeltingEvents, sumOp<label>());

    Info << "solidificationData recorded " << nSolidificationEvents 
         << " solidification events." << endl;
         
    if (trackMelting_)
    {
        Info << "solidificationData recorded " << nMeltingEvents
             << " melting events." << endl;
    }

    //- Clean cells which have been remelted (not yet implemented)
    //removeRemelts();

    return true;
}


bool Foam::functionObjects::solidificationData::end()
{
    return true;
}


bool Foam::functionObjects::solidificationData::write()
{
    if (AMR_)
    {
        //- Get current time
        const fileName currTime = Foam::name(mesh_.time().value());

        //- Create file path for current time data
        const fileName currTimePath(mesh_.time().rootPath()
                                    /mesh_.time().globalCaseName()
                                    /"solidificationData"/currTime);
                                    
        mkDir(currTimePath);

        //- Open file for each proc
        OFstream os(currTimePath + "/" + "data_" 
                    + Foam::name(Pstream::myProcNo()) + ".csv");

        //- Write header
        os << "x,y,z,t,R,G\n";

        //- Write each event in series to file
        for (int i = 0; i < solidificationEvents_.size(); ++i)
        {
            int n = solidificationEvents_[i].size() - 1;

            for (int j = 0; j < n; ++j)
            {
                os << solidificationEvents_[i][j] << ",";
            }

            os << solidificationEvents_[i][n] << "\n";
        }

        solidificationEvents_.clear();
        
        if (trackMelting_)
        {
            const fileName meltingDir(currTimePath/"meltingEvents");
            
            mkDir(meltingDir);
            
            OFstream mos(meltingDir + "/" + "data_"
                         + Foam::name(Pstream::myProcNo()) + ".csv");
                         
            mos << "x,y,z,t\n";
            
            for (int i = 0; i < meltingEvents_.size(); ++i)
            {
                int n = meltingEvents_[i].size() - 1;

                for (int j = 0; j < n; ++j)
                {
                    os << meltingEvents_[i][j] << ",";
                }

                os << meltingEvents_[i][n] << "\n";
            }
        }
    }

    return true;
}


void Foam::functionObjects::solidificationData::movePoints(const polyMesh& mesh)
{
    if (&mesh == &mesh_)
    {
        correct();
    }
}


void Foam::functionObjects::solidificationData::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        correct();
    }
}


void Foam::functionObjects::solidificationData::mapMesh(const polyMeshMap& map)
{
    if (&map.mesh() == &mesh_)
    {
        correct();
    }
}


void Foam::functionObjects::solidificationData::distribute
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
