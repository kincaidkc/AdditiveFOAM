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

#include "timeAboveMelting.H"
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
    defineTypeNameAndDebug(timeAboveMelting, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        timeAboveMelting,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::timeAboveMelting::correct()
{
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::timeAboveMelting::timeAboveMelting
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
    T_(mesh_.lookupObject<VolField<scalar>>("T"))
{
    read(dict);

    if (AMR_)
    {
        refinedSize_ = dict_.lookup<scalar>("refinedSize");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::timeAboveMelting::~timeAboveMelting()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::timeAboveMelting::read(const dictionary& dict)
{
    return true;
}


Foam::wordList Foam::functionObjects::timeAboveMelting::fields() const
{
    return wordList::null();
}


bool Foam::functionObjects::timeAboveMelting::execute()
{
    //- Get current time
    const scalar& time = mesh_.time().value();
    const scalar& dt = mesh_.time().deltaTValue();

    //- Get old temperature and time derivative of temperature
    const volScalarField& T0_ = T_.oldTime();

    //- Calculate refined mesh volume
    const scalar Vr = Foam::pow(1.1 * refinedSize_, 3.0);

    forAll(mesh_.C(), celli)
    {
        //- Check that cells are in max refinement level for AMR cases
        if (AMR_)
        {
            //- Skip non-fully refined cells
            if (mesh_.V()[celli] > Vr)
            {
                continue;
            }
        }

        const scalar T = T_[celli];
        const scalar T0 = T0_[celli];

        //- Check for melted cells
        if (((T >= Tl_) && (T0 < Tl_)) || ((T <= Tl_) && (T0 > Tl_)))
        {
            if (mag(T - T0) > small)
            {
                vector C = mesh_.C()[celli];

                List<scalar> eventi(5);

                eventi[0] = C[0];
                eventi[1] = C[1];
                eventi[2] = C[2];
                eventi[3] = time;
                eventi[4] = dt * min(max((Tl_ - T0) / (T - T0), 0), 1);

                events_.append(eventi);
            }
        }
        else if (T >= Tl_)
        {
            vector C = mesh_.C()[celli];

            List<scalar> eventi(5);

            eventi[0] = C[0];
            eventi[1] = C[1];
            eventi[2] = C[2];
            eventi[3] = time;
            eventi[4] = dt;

            events_.append(eventi);
        }
    }

    //- Consolidate duplicate entries (not yet implemented)
    //consolidateTAMs();

    return true;
}


bool Foam::functionObjects::timeAboveMelting::end()
{
    return true;
}


bool Foam::functionObjects::timeAboveMelting::write()
{
    //- Get current time
    const fileName currTime = Foam::name(mesh_.time().value());

    //- Create file path for current time data
    const fileName currTimePath(mesh_.time().rootPath()
                                /mesh_.time().globalCaseName()
                                /"timeAboveMelting"/currTime);

    mkDir(currTimePath);

    //- Open file for each proc
    OFstream os(currTimePath + "/" + "data_" 
                + Foam::name(Pstream::myProcNo()) + ".csv");

    //- Write header
    os << "x,y,z,t,dtm\n";

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

    events_.clear();

    return true;
}


void Foam::functionObjects::timeAboveMelting::movePoints(const polyMesh& mesh)
{
    if (&mesh == &mesh_)
    {
        correct();
    }
}


void Foam::functionObjects::timeAboveMelting::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        correct();
    }
}


void Foam::functionObjects::timeAboveMelting::mapMesh(const polyMeshMap& map)
{
    if (&map.mesh() == &mesh_)
    {
        correct();
    }
}


void Foam::functionObjects::timeAboveMelting::distribute
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
