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
    Foam::functionObjects::unmappedSolidificationFields
Description
    This function object records the solidification cooling rate and thermal
    gradient for additive manufacturing simulations. Each time a cell crosses
    the specified liquidus/melting temperature, a new event is created with
    the information on the location, time, and solidification data. Events
    are appended to a single list so that multiple remelting events can be
    accounted for. This requires some consolidation in post-processing.
    
    The option is provided to track the melting time of each point in a
    separate list for use with solidifcation codes. This feature also requires
    post-processing to match each melting event with a corresponding
    solidification event.
    
    Example of function object specification:
    \verbatim
    unmappedSolidificationFields1
    {
        type            unmappedSolidificationFields;
        
        libs            ("libadditiveFoamFunctionObjects.so");
        
        writeControl    writeTime;
        
        AMR             true;
        refinedSize     10e-6;
        Tl              1620;
        trackMelting    true;
    }
    \endverbatim
Usage
    \table
        Property     | Description               | Required | Default value
        type         | type name: unmappedSolidificationFields | yes      |
        AMR          | is this an AMR simulation?  | no       | no
        refinedSize  | length scale of fully refined cell | no | 0
        Tl           | liquidus/melting temperature of material | yes |
        trackMelting | record melting events       | no         | false
    \endtable
SourceFiles
    unmappedSolidificationFields.C
\*---------------------------------------------------------------------------*/

#ifndef unmappedSolidificationFields_H
#define unmappedSolidificationFields_H

#include "fvMeshFunctionObject.H"
#include "meshSearch.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "writeFile.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class unmappedSolidificationFields Declaration
\*---------------------------------------------------------------------------*/

class unmappedSolidificationFields
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Input dictionary
        dictionary dict_;

        //- AMR switch
        const bool AMR_;

        //- Side length of fully refined cells
        scalar refinedSize_;
        
        //- Switch to track melting time
        const bool trackMelting_;

        //- Liquidus temperature
        const scalar Tl_;

        //- Pointer to temperature field
        const volScalarField& T_;

        //- Cooling rate field
        volScalarField R_;

        //- List of solidification events
        DynamicList<List<scalar>> solidificationEvents_;
        
        //- List of melting events
        DynamicList<List<scalar>> meltingEvents_;


    // Private Member Functions

        //- Correct for mesh changes
        void correct();

public:

    //- Runtime type information
    TypeName("unmappedSolidificationFields");


    // Constructors

        //- Construct from Time and dictionary
        unmappedSolidificationFields
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        unmappedSolidificationFields(const unmappedSolidificationFields&) = delete;


    //- Destructor
    virtual ~unmappedSolidificationFields();


    // Member Functions

        //- Read the unmappedSolidificationFields data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Calculate and write unmappedSolidificationFields
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
        void operator=(const unmappedSolidificationFields&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
