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

#include "stork.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "interpolation.H"

#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "polyDistributionMap.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "labelVector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(stork, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        stork,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::stork::correct()
{
    // Update the stork cell model list
    if (cellToSortedVertices_.size() != mesh_.nCells())
    {
        cellToSortedVertices_.setSize(mesh_.nCells());
    }
   
    forAll(cellToSortedVertices_, celli)
    {
        cellToSortedVertices_[celli].clear();
    }
    
    // Update indexing for stork cell model
    const pointField& points = mesh_.points();
    
    // Fully refined cells will have a "hex" cell shape
    const cellModel& hex = *(cellModeller::lookup("hex"));
    
    forAll(mesh_.cells(), celli)
    {
        if ( mesh_.cellShapes()[celli].model() == hex )
        {            
            boundBox cellBb(point::max, point::min);
            
            const labelList& vertices = mesh_.cellPoints()[celli];

            for (const label& pointi : vertices)
            {
                cellBb.min() = min(cellBb.min(), points[pointi]);
                cellBb.max() = max(cellBb.max(), points[pointi]);
            }

            if (!cellBb.overlaps(box_))
            {
                continue;
            }
          
            const vector span = cellBb.span();

            labelList sortedVrtList = vertices;

            for (const label& pointi : vertices)
            {                
                const vector pt =
                    cmptDivide(mesh_.points()[pointi] -  cellBb.min(), span);
            
                const label index = 
                        4 * static_cast<label>(pt[0] + 0.5)
                      + 2 * static_cast<label>(pt[1] + 0.5)
                      + 1 * static_cast<label>(pt[2] + 0.5);
                      
               sortedVrtList[index] = pointi;
            }
            
            cellToSortedVertices_[celli] = sortedVrtList;
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::stork::stork
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    dict_(dict),
    outputPath_
    (
        mesh_.time().globalPath()
       /writeFile::outputPrefix
       /(mesh_.name() != polyMesh::defaultRegion ? mesh_.name() : word())
       /name
    ),
    T_(mesh_.lookupObject<VolField<scalar>>("T")),
    vpi_(volPointInterpolation::New(mesh_)),
    Tp_
    (
        IOobject
        (
            "Tp_",
            mesh_.time().name(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        vpi_.interpolate(T_)
    ),
    referenceBb
    (
        point::max,
        point::min
    ),
    cellToSortedVertices_(mesh_.nCells())
{
    read(dict);
    
    initialize();
    
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::stork::~stork()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::stork::read(const dictionary& dict)
{
    box_ = dict.lookup("box");
    
    isoValue_ = dict.lookup<scalar>("isoValue");
    
    spacing_  = dict.lookup<scalar>("spacing");
    
    return true;
}


void Foam::functionObjects::stork::initialize()
{
    const pointField& points = mesh_.points();

    forAll(mesh_.cells(), celli)
    {
        boundBox cellBb(point::max, point::min);

        const labelList& vertices = mesh_.cellPoints()[celli];

        forAll(vertices, i)
        {
            cellBb.min() = min(cellBb.min(), points[vertices[i]]);
            cellBb.max() = max(cellBb.max(), points[vertices[i]]);
        }

        if (cellBb.overlaps(box_))
        {            
            referenceBb.min() = min(referenceBb.min(), cellBb.min());
            referenceBb.max() = max(referenceBb.max(), cellBb.max());
        }
    }
    
    reduce(referenceBb.min(), minOp<vector>());
    reduce(referenceBb.max(), maxOp<vector>());
    
    Info << "target bounding box for interpolation: " << box_ << endl;
    Info << "reference bounding box for interpolation: " << referenceBb << endl;
}


Foam::wordList Foam::functionObjects::stork::fields() const
{    
    return wordList::null();
}


bool Foam::functionObjects::stork::execute()
{
    const pointScalarField Tp0_("Tp0_", Tp_);

    Tp_ = vpi_.interpolate(T_);
          
    // Capture stork events: {i,j,k,t0,t,Tp0[8],Tp[8]}
    forAll(mesh_.cells(), celli)
    {
        const labelList& sortedVrtList = cellToSortedVertices_[celli];

        if (sortedVrtList.size() == 0)
        {
            continue;   // not a refined hex cell
        }
        
        label c0 = 0;
        label c1 = 0;

        for (const label pointi : sortedVrtList)
        {
            if (Tp0_[pointi] >= isoValue_)
            {
                c0++;
            }
            if (Tp_[pointi]  >= isoValue_)
            {
                c1++;
            }
        }

        // Format event and append        
        if ( (c0 % 8) || (c1 % 8) || (c0 != c1) )
        {
            const vector pt =
                (mesh_.points()[sortedVrtList[0]] - referenceBb.min())
              / spacing_;

            List<scalar> event(21);
            
            event[0] = static_cast<label>(pt[0] + 0.5);
            event[1] = static_cast<label>(pt[1] + 0.5);
            event[2] = static_cast<label>(pt[2] + 0.5);
            event[3] = mesh_.time().value() - mesh_.time().deltaTValue();
            event[4] = mesh_.time().value();
            
            label index = 5;
            
            for (const label pointi : sortedVrtList)
            {
                event[index] = Tp0_[pointi];
                event[index + 8] = Tp_[pointi];
                index++;
            }
     
            events.append(event);
        }
    }
    
    return true;
}


bool Foam::functionObjects::stork::end()
{
    if (!mesh_.time().writeTime())
    {
        write();
    }
     
    return true;
}


bool Foam::functionObjects::stork::write()
{
    if (mesh_.time().writeTime())
    {
        events.shrink();
                
        const uint32_t nEvents = returnReduce(events.size(), sumOp<scalar>());
        
        if (nEvents == 0)
        {
            return true;    // no solidification events: skip time step
        }

        labelVector localMin = labelVector::max;
        labelVector localMax = labelVector::min;
                   
        for (const auto& event : events)
        {
            localMin.x() = min(localMin.x(), event[0]);
            localMin.y() = min(localMin.y(), event[1]);
            localMin.z() = min(localMin.z(), event[2]);

            localMax.x() = max(localMax.x(), event[0]);
            localMax.y() = max(localMax.y(), event[1]);
            localMax.z() = max(localMax.z(), event[2]);           
            
        }
        labelVector globalMin = returnReduce(localMin, minOp<labelVector>());
        labelVector globalMax = returnReduce(localMax, maxOp<labelVector>());
                
        Info<< "Solidification Events:" << endl;
        Info<< "    Global Size:\t" << nEvents << endl;    
        Info<< "    Index Span:\t" << globalMax - globalMin << endl;
        Info<< "    Reference Point:\t" << referenceBb.min() << endl;
        
        //- Stork binary writer
        const fileName path_
        (
            mesh_.time().rootPath()/mesh_.time().globalCaseName()/"Stork"
           /mesh_.time().name()
        );

        mkDir(path_);
                
        //- write global header file   
        if (Pstream::master())
        {
            std::ofstream os
            (
                path_ + "/info.bin", 
                std::ios::binary
            );

            const uint32_t uint32Header[7] =
                {
                    static_cast<uint32_t>(nEvents),
                    static_cast<uint32_t>(globalMin[0]),
                    static_cast<uint32_t>(globalMin[1]),
                    static_cast<uint32_t>(globalMin[2]),
                    static_cast<uint32_t>(globalMax[0]),
                    static_cast<uint32_t>(globalMax[1]),
                    static_cast<uint32_t>(globalMax[2])
                };
            
            const float floatHeader[5] =
                {
                    static_cast<float>(referenceBb.min().x()),
                    static_cast<float>(referenceBb.min().y()),
                    static_cast<float>(referenceBb.min().z()),
                    static_cast<float>(spacing_),
                    static_cast<float>(isoValue_)
                };

            os.write
            (
                reinterpret_cast<const char*> (&uint32Header),
                7 * sizeof(uint32_t)
            );

            os.write
            (
                reinterpret_cast<const char*> (&floatHeader),
                5 * sizeof(float)
            );
            
            os.close();
        }
        
        //- write data
        if (events.size() > 0)
        {
            std::ofstream os
            (
                path_ + "/processor"
              + Foam::name(Pstream::myProcNo()) + ".bin",
                std::ios::binary
            );
            
            //- local header
            const uint32_t uint32Header[7] =
                {
                    static_cast<uint32_t>(events.size()),
                    static_cast<uint32_t>(localMin[0]),
                    static_cast<uint32_t>(localMin[1]),
                    static_cast<uint32_t>(localMin[2]),
                    static_cast<uint32_t>(localMax[0]),
                    static_cast<uint32_t>(localMax[1]),
                    static_cast<uint32_t>(localMax[2])
                };
            
            const float floatHeader[5] =
                {
                    static_cast<float>(referenceBb.min().x()),
                    static_cast<float>(referenceBb.min().y()),
                    static_cast<float>(referenceBb.min().z()),
                    static_cast<float>(spacing_),
                    static_cast<float>(isoValue_)
                };

            os.write
            (
                reinterpret_cast<const char*> (&uint32Header),
                7 * sizeof(uint32_t)
            );

            os.write
            (
                reinterpret_cast<const char*> (&floatHeader),
                5 * sizeof(float)
            );
            
            //- local data
            for (const auto& event : events)
            {
                for (int i=0; i < 3; i++)
                {
                    uint32_t e = static_cast<uint32_t>(event[i]);
                    os.write
                    (
                        reinterpret_cast<const char*> (&e),
                        sizeof(uint32_t)
                    );
                }
            }

            for (const auto& event : events)
            {
                for (int i=3; i < 5; i++)
                {
                    float ef = static_cast<float>(event[i]);
                    os.write
                    (
                        reinterpret_cast<const char*> (&ef),
                        sizeof(float)
                    );
                }
            }

            for (const auto& event : events)
            {
                for (int i=5; i < 21; i++)
                {
                    float ef = static_cast<float>(event[i]);
                    os.write
                    (
                        reinterpret_cast<const char*> (&ef),
                        sizeof(float)
                    );
                }
            }

            os.close();
        }
        
        events.setSize(0);
    }

    return true;
}


void Foam::functionObjects::stork::movePoints(const polyMesh& mesh)
{
    if (&mesh == &mesh_)
    {
        correct();
    }
}


void Foam::functionObjects::stork::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        correct();
    }
}


void Foam::functionObjects::stork::mapMesh(const polyMeshMap& map)
{
    if (&map.mesh() == &mesh_)
    {
        correct();
    }
}


void Foam::functionObjects::stork::distribute
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
