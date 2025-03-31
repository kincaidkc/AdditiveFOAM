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
    Foam::functionObjects::solidificationData
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
    solidificationData1
    {
        type            solidificationData;
        
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
        type         | type name: solidificationData | yes      |
        AMR          | is this an AMR simulation?  | no       | no
        refinedSize  | length scale of fully refined cell | no | 0
        Tl           | liquidus/melting temperature of material | yes |
        trackMelting | record melting events       | no         | false
    \endtable
SourceFiles
    solidificationData.C
\*---------------------------------------------------------------------------*/

#ifndef solidificationData_H
#define solidificationData_H

#include "fvMeshFunctionObject.H"
//#include "meshSearch.H"
#include "treeBoundBox.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "writeFile.H"
#include "OFstream.H"
#include "labelVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class solidificationData Declaration
\*---------------------------------------------------------------------------*/

class solidificationData
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Input dictionary
        dictionary dict_;

        //- AMR switch
        const bool AMR_;
        
        //- Write csv file for non-AMR case
        const bool writeCsv_;
        
        //- Bounding box for AMR background mesh
        treeBoundBox bb_;
        
        //- Number of cells in each direction of background mesh
        labelVector nBgPts_;

        //- Side length of fully refined cells
        scalar refinedSize_;

        //- Liquidus temperature
        const scalar Tl_;

        //- Pointer to temperature field
        const volScalarField& T_;
        
        //- Time derivative of temperature
        volScalarField dTdt_;
        
        //- Solidification cooling rate (for non-AMR case)
        volScalarField R_;
        
        //- Solidification temperature gradient (for non-AMR case)
        volScalarField G_;
        
        //- Melt time (for non-AMR case)
        volScalarField tm_;
        
        //- Solidification time (for non-AMR case)
        volScalarField ts_;
        
        //- Time above melting (for non-AMR case)
        volScalarField tam_;

        //- List of solidification events (for AMR case)
        DynamicList<List<scalar>> solidData_;


    // Private Member Functions

        //- Correct for mesh changes
        void correct();

public:

    //- Runtime type information
    TypeName("solidificationData");


    // Constructors

        //- Construct from Time and dictionary
        solidificationData
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        solidificationData(const solidificationData&) = delete;


    //- Destructor
    virtual ~solidificationData();


    // Member Functions

        //- Read the solidificationData data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Calculate and write solidificationData
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
        void operator=(const solidificationData&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
