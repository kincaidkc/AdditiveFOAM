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

Class
    Foam::functionObjects::stork

Description
    This function object...

    Example of function object specification:
    \verbatim
    stork1
    {
        type           stork;
        libs ("libstorkFunctionObject.so");
        ...
    }
    \endverbatim

Usage
    \table
        Property     | Description               | Required | Default value
        type         | type name: stork | yes      |
    \endtable

SourceFiles
    stork.C

\*---------------------------------------------------------------------------*/

#ifndef stork_H
#define stork_H

#include "fvMeshFunctionObject.H"
#include "meshSearch.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "writeFile.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class stork Declaration
\*---------------------------------------------------------------------------*/

class stork
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Input dictionary
        dictionary dict_;

        //- Output path
        fileName outputPath_;

        const volScalarField& T_;

        const volPointInterpolation& vpi_;
        
        mutable pointScalarField Tp_;

        boundBox  referenceBb;
        
        List<labelList> cellToSortedVertices_;

        DynamicList<List<scalar>> events;
                
        //- Input for stork reference mesh
        boundBox  box_;
        
        scalar    spacing_;
        
        scalar    isoValue_;
        
    // Private Member Functions

        //- Correct for mesh changes
        void correct();

public:

    //- Runtime type information
    TypeName("stork");


    // Constructors

        //- Construct from Time and dictionary
        stork
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        stork(const stork&) = delete;


    //- Destructor
    virtual ~stork();


    // Member Functions

        //- Read the stork data
        virtual bool read(const dictionary&);

        //- Set up stork reference mesh
        virtual void initialize();

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Calculate and write stork
        virtual bool write();

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh&);

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const stork&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
