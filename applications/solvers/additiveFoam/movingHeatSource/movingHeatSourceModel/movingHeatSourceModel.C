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

#include "movingHeatSourceModel.H"
#include "DynamicList.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingHeatSourceModel::movingHeatSourceModel
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    dict_
    (    
        IOobject
        (
            "heatSourceDict",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    sourceNames_(dict_.lookup("sources")),
    refine_(dict_.lookupOrDefault<bool>("refine", false)),
    nSteps_(0),
    currStep_(0),
    qDot_
    (
        IOobject
        (
            "qDot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimPower/dimVolume, 0.0)
    ),
    refinementField_
    (
    	IOobject
    	(
    	    "refinementField",
    	    mesh_.time().timeName(),
    	    mesh_,
    	    IOobject::READ_IF_PRESENT,
    	    IOobject::NO_WRITE
    	),
    	mesh_,
    	dimensionedScalar(dimless, 0.0)
    )  
{
    sources_.resize(sourceNames_.size());
        
    //- Create new instance of movingBeam for each beam
    forAll(sources_, i)
    {
        Info << "Adding heatSourceModel for " << sourceNames_[i] << endl;
        sources_.set
        (
            i,
            heatSourceModel::New
            (
                sourceNames_[i],
                dict_,
                mesh_
            ).ptr()
        );
    }
    
    if(refine_)
        nSteps_ = dict_.lookup<int>("nSteps");
}

// * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * * //

Foam::movingHeatSourceModel::~movingHeatSourceModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::movingHeatSourceModel::adjustDeltaT(scalar& deltaT)
{
    forAll(sources_, i)
    {
        sources_[i].beam().adjustDeltaT(deltaT);
    }
}

void Foam::movingHeatSourceModel::update()
{
    //- Reset qDot field at every time step
    qDot_ = dimensionedScalar("Zero", qDot_.dimensions(), 0.0);
    
    ++currStep_;
    
    //- Reset refinement field after desired interval
    if (currStep_ == nSteps_)
        refinementField_ = dimensionedScalar(refinementField_.dimensions(), 0.0);
    
    //- Subcycle each moving heat source in time and combine into a single field
    forAll(sources_, i)
    {
        if (sources_[i].beam().activePath())
        {
            sources_[i].updateDimensions();

            // integrate volumetric heat source over desired time step
            scalar pathTime = mesh_.time().value();

            const scalar nextTime = pathTime + mesh_.time().deltaTValue();

            const scalar beam_dt = sources_[i].beam().deltaT();

            //- Update individual beam heat source contribution
            volScalarField qDoti
            (
                IOobject
                (
                    "qDoti",
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar("Zero", qDot_.dimensions(), 0.0)
            );
            
            scalar sumWeights = 0.0;

            while ((nextTime - pathTime) > small)
            {
                scalar dt = min(beam_dt, max(0, nextTime - pathTime));

                pathTime += dt;

                sources_[i].beam().move(pathTime);
                                
                qDoti += dt*sources_[i].qDot();

                sumWeights += dt;
            }
            
            qDoti /= sumWeights;
            
            qDot_ += qDoti;
            
            if((refine_) && (currStep_ == nSteps_))
            {
                //- Set refinement around melt pool
                const volScalarField& alpha1
                    = mesh_.lookupObject<volScalarField>("alpha.solid");
                
                refinementField_ = pos(1.0 - alpha1);
                
                //- Add refinement around individual beam position for next n timesteps
                scalar currTime = mesh_.time().value();
                scalar currDt = mesh_.time().deltaT().value();
                vector currPos = Zero;
                scalar currPow = 0.0;
                
                for(int i = 0; i < nSteps_; ++i)
                {
                    currTime += i * currDt;
                    
                    sources_[i].beam().move(currPos, currPow, currTime);
                    
                    Info << "Current position: " << currPos << endl;
                    Info << "Current power: " << currPow << endl;
                    
                    refinementField_
                        += pos0(sources_[i].D2sigma().value() - mag(mesh_.C() - currPos));
                }
            }
        }
    }
    
    if (currStep_ == nSteps_)
        currStep_ = 0;
}

// ************************************************************************* //
