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
    T_(mesh_.lookupObject<VolField<scalar>>("T"))
{
    read(dict);
    
    if (AMR_)
    {
        refinedSize_ = dict_.lookup<scalar>("refinedSize");
    }
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
    Info << "executing unmappedSolidifaction fields" << endl;
    //- Get current time
    const scalar& time = mesh_.time().value();
    
    //- Get old temperature and time derivative of temperature
    const volScalarField& T0 = T_.oldTime();
    
    const volScalarField dTdt = fvc::ddt(T_);
    
    forAll(mesh_.C(), celli)
    {
        //- Check that cells are in max refinement level for AMR cases
        bool maxLevel = true;
        
        if (AMR_)
        {
            if (mesh_.V()[celli] > (pow(refinedSize_, 3.0) + VSMALL))
            {
                maxLevel = false;
            }
        }
        
        //- Check for solidification events
        if ((T0[celli] > Tl_) && (T_[celli] <= Tl_) && (maxLevel))
        {
            vector C = mesh_.C()[celli];
            
            List<scalar> eventi(5);
            
            eventi[0] = C[0];
            eventi[1] = C[1];
            eventi[2] = C[2];
            
            eventi[3] = time;
            
            eventi[4] = dTdt[celli];
            
            events_.append(eventi);
        }
    }
    
    //- Clean cells which have been remelted (not yet implemented)
    //removeRemelts();
    
    return true;
}


bool Foam::functionObjects::unmappedSolidificationFields::end()
{
    Info << "ending unmapped" << endl;
    const fileName writePath
    (
        mesh_.time().rootPath()/mesh_.time().globalCaseName()/"unmappedSolidificationFields"
    );
    
    Info << "unmapped write dir: " << writePath << endl;
    
    mkDir(writePath);
    
    OFstream os
    (
        writePath + "/" + "data_" + Foam::name(Pstream::myProcNo()) + ".csv"
    );
    
    os << "x,y,z,t,dTdt\n";
    
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


bool Foam::functionObjects::unmappedSolidificationFields::write()
{
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
