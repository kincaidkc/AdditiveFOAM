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
    writeCsv_(dict_.lookupOrDefault<bool>("writeCsv", false)),
    bb_(point::min, point::max),
    nBgPts_(0, 0, 0),
    refinedSize_(0.0),
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
            AMR_ ? IOobject::NO_READ : IOobject::READ_IF_PRESENT,
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
            AMR_ ? IOobject::NO_READ : IOobject::READ_IF_PRESENT,
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
            AMR_ ? IOobject::NO_READ : IOobject::READ_IF_PRESENT,
            AMR_ ? IOobject::NO_WRITE : IOobject::AUTO_WRITE
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
            AMR_ ? IOobject::NO_READ : IOobject::READ_IF_PRESENT,
            AMR_ ? IOobject::NO_WRITE : IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTime, 0.0)
    ),
    tam_
    (
        IOobject
        (
            "tam",
            mesh_.time().timeName(),
            mesh_,
            AMR_ ? IOobject::NO_READ : IOobject::READ_IF_PRESENT,
            AMR_ ? IOobject::NO_WRITE : IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTime, 0.0)
    )
{
    read(dict);

    if (AMR_)
    {
        refinedSize_ = dict_.lookup<scalar>("refinedSize");
        bb_ = dict_.lookup<treeBoundBox>("boundingBox");
        
        //- Resize solidification data to size of background mesh
        nBgPts_[0] = static_cast<int>((bb_.max()[0] - bb_.min()[0]) / refinedSize_ + 0.5);
        nBgPts_[1] = static_cast<int>((bb_.max()[1] - bb_.min()[1]) / refinedSize_ + 0.5);
        nBgPts_[2] = static_cast<int>((bb_.max()[2] - bb_.min()[2]) / refinedSize_ + 0.5);
        
        Info << "Background mesh measures " << nBgPts_[0] << " by "
             << nBgPts_[1] << " by " << nBgPts_[2] << endl;
        
        int numPts
            = static_cast<int>(((bb_.max()[0] - bb_.min()[0])
                              * (bb_.max()[1] - bb_.min()[1])
                              * (bb_.max()[2] - bb_.min()[2])
                              / Foam::pow(refinedSize_, 3.0) + 0.5));
            
        solidData_.resize(numPts);
        
        // Populate solidification data with cell center locations
        label index = 0;
        
        for (int z = 0; z < nBgPts_[2]; ++z)
        {
            for (int y = 0; y < nBgPts_[1]; ++y)
            {
                for (int x = 0; x < nBgPts_[0]; ++x)
                {
                    List<scalar> datai(8);
                    
                    datai[0] = bb_.min()[0] + x * refinedSize_ + refinedSize_ / 2.0;
                    datai[1] = bb_.min()[1] + y * refinedSize_ + refinedSize_ / 2.0;
                    datai[2] = bb_.min()[2] + z * refinedSize_ + refinedSize_ / 2.0;
                    datai[3] = 0.0; // Temperature gradient
                    datai[4] = 0.0; // Cooling rate
                    datai[5] = 0.0; // Melting time
                    datai[6] = 0.0; // Solidification time
                    datai[7] = 0.0; // Time above melting
                    
                    solidData_[index] = datai;
                    
                    ++index;
                }
            }
        }
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
    
    const auto& meshC = mesh_.C();
    const auto& meshV = mesh_.V();
    const pointField& points = mesh_.points();

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
    label ptNum = 0;

    forAll(meshC, celli)
    {
        //- Check that current cell is within bounding box
        treeBoundBox cellBb(point::max, point::min);
        
        const labelList& vertices = mesh_.cellPoints()[celli];
        
        forAll(vertices, i)
        {
            cellBb.min() = min(cellBb.min(), points[vertices[i]]);
            cellBb.max() = max(cellBb.max(), points[vertices[i]]);
        }
        
        if (!cellBb.overlaps(bb_))
        {
            continue;
        }
        
        //- Check that cells are in max refinement level for AMR cases
        if (AMR_)
        {
            if (meshV[celli] > Vr)
            {
                continue;
            }
            else
            {
                //- Find index of current location in structured mesh
                const vector& C = meshC[celli];
                int nx = static_cast<int>(((C[0] - bb_.min()[0] - refinedSize_ / 2.0)
                            / refinedSize_ + 0.5));
                int ny = static_cast<int>(((C[1] - bb_.min()[1] - refinedSize_ / 2.0)
                            / refinedSize_ + 0.5));
                int nz = static_cast<int>(((C[2] - bb_.min()[2] - refinedSize_ / 2.0)
                            / refinedSize_ + 0.5));
                
                ptNum = static_cast<int>((nz * nBgPts_[0] * nBgPts_[1]
                               + ny * nBgPts_[0] + nx + 0.5));
            }
        }
        
        const scalar T = T_[celli];
        const scalar T0 = T0_[celli];
        
        if (T >= Tl_)
        {
            if (AMR_)
            {
                solidData_[ptNum][7] += deltaT;
            }
            else
            {
                tam_[celli] += deltaT;
            }
        }

        //- Check for solidification events
        if ((T0 > Tl_) && (T <= Tl_))
        {
            if (AMR_)
            {                
                //- Update data for this point
                solidData_[ptNum][3] = gradT[celli];
                solidData_[ptNum][4] = dTdt_.oldTime()[celli];
                solidData_[ptNum][6]
                    = time - deltaT * min(max((Tl_ - T) / (T0 - T), 0), 1);
                solidData_[ptNum][7] += deltaT * min(max((Tl_ - T) / (T0 - T), 0), 1);
            }
            else
            {
                R_[celli] = dTdt_.oldTime()[celli];
                G_[celli] = gradT[celli];
                
                ts_[celli]
                    = time - deltaT * min(max((Tl_ - T) / (T0 - T), 0), 1);
                    
                tam_[celli] += deltaT - deltaT * min(max((Tl_ - T) / (T0 - T), 0), 1);
            }

            ++nSolidificationEvents;
        }
        
        //- Check for melting events
        if ((T0 < Tl_) && (T >= Tl_))
        {
            if (AMR_)
            {
                solidData_[ptNum][5]
                    = time - deltaT * min(max((T - Tl_) / (T - T0), 0), 1);
                    
                solidData_[ptNum][7]
                    += deltaT * min(max((T - Tl_) / (T - T0), 0), 1);
            }
            else
            {
                tm_[celli]
                    = time - deltaT * min(max((T - Tl_) / (T - T0), 0), 1);
                    
                tam_[celli]
                    += deltaT * min(max((T - Tl_) / (T - T0), 0), 1);
            }
            
            ++nMeltingEvents;
        }
    }

    reduce(nSolidificationEvents, sumOp<label>());
    reduce(nMeltingEvents, sumOp<label>());

    Info << "solidificationData recorded " << nSolidificationEvents 
         << " solidification events." << endl;

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
        os << "x,y,z,G,R,tm,ts,tam,melted?\n";

        //- Write each event in series to file
        for (int i = 0; i < solidData_.size(); ++i)
        {
            //Info << "writing line " << i << " of " << solidData_.size() << endl;
            int n = solidData_[i].size() - 1;

            for (int j = 0; j < n; ++j)
            {
                os << solidData_[i][j] << ",";
            }

            os << solidData_[i][n] << ","
               << (solidData_[i][n] > 0.0 ? 1.0 : 0.0) << "\n";
        }
    }
    
    else if (writeCsv_)
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
        os << "x,y,z,G,R,tm,ts,tam,melted?\n";
        
        const auto& meshC = mesh_.C();
        const pointField& points = mesh_.points();
        
        for (int i = 0; i < mesh_.nCells(); ++i)
        {
            //- Check that current cell is within bounding box
            treeBoundBox cellBb(point::max, point::min);
            
            const labelList& vertices = mesh_.cellPoints()[i];
            
            forAll(vertices, j)
            {
                cellBb.min() = min(cellBb.min(), points[vertices[j]]);
                cellBb.max() = max(cellBb.max(), points[vertices[j]]);
            }
            
            if (!cellBb.overlaps(bb_))
            {
                continue;
            }
        
            const auto& C = meshC[i];
            
            os << C[0] << "," << C[1] << "," << C[2] << "," << G_[i] << ","
               << R_[i] << "," << tm_[i] << "," << ts_[i] << "," << tam_[i]
               << "," << (tam_[i] > 0.0 ? 1.0 : 0.0) << endl;
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
