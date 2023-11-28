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
    resolveTail_(dict_.lookupOrDefault<bool>("resolveTail", false)),
    nSteps_(0),
    nRevSteps_(0),
    nRefineSteps_(0),
    currStep_(0),
    nextRefinement_(0),
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
    
    Info << "reading refine fields" << endl;
    
    if(refine_)
    {
        nSteps_ = dict_.lookup<label>("nSteps");
        nRevSteps_ = dict_.lookup<label>("nRevSteps");
        currStep_ = nSteps_ - 1;
        nRefineSteps_ = dict_.lookupOrDefault<label>("nRefine", 1);
    }
}

// * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * * //

Foam::movingHeatSourceModel::~movingHeatSourceModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::movingHeatSourceModel::refineTime()
{
    if(mesh_.time().timeIndex() == nextRefinement_)
    {
        Info << "IT'S REFINE TIME" << endl;
        
        nextRefinement_ += nSteps_;
        
        return true;
    }

    else if (mesh_.time().timeIndex() <= nextRefinement_ - nSteps_ + nRefineSteps_)
    {
        Info << "Performing additional refinement operations." << endl;
        
        return true;
    }
    
    else
    {
        return false;
    }
}

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
    
    const volScalarField& alpha1
        = mesh_.lookupObject<volScalarField>("alpha.solid");
    
    //- Reset refinement field after desired interval
    if (currStep_ == nSteps_)
    {
        if (resolveTail_)
        {
            //- Reset to capture melt pool tail
            refinementField_ = pos(1.0 - alpha1);
        }
        
        else
        {
            //- Reset to zero
            refinementField_
                = dimensionedScalar(refinementField_.dimensions(), 0.0);
        }
    }
    
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
            
            Info << "Current beam position: " << sources_[i].beam().position() << endl;
            
            qDoti /= sumWeights;
            
            qDot_ += qDoti;
            
            if((refine_) && (currStep_ == nSteps_))
            {
                //- Add refinement around individual beam position for next n timesteps
                scalar currTime = mesh_.time().value();
                scalar currDt = mesh_.time().deltaT().value();
                vector currPos = Zero;
                scalar currPow = 0.0;
    
                dimensionedScalar charLen_ = fvc::domainIntegrate(Foam::pow(mesh_.V(), 1.0/3.0))
                                             / fvc::domainIntegrate(mesh_.V() / mesh_.V());
                                             
                Info << "Mesh characteristic length: " << charLen_ << endl;

                const int ratio = 2;
                
                dimensionedScalar searchRadius(dimLength, sources_[i].D2sigma());
                
                //- Repeat search until refinement field is non-zero
                //  This is required for cases with extremely coarse starting mesh
                while (max(refinementField_).value() < 1.0)
                {
                    Info << "Searching for cells to refine in radius " << searchRadius.value() << endl;
                
                    //- Integrate forward ahead of the beam
                    for(int j = 0; j < ratio * nSteps_ + 5; ++j)
                    {
                        scalar nowTime = currTime + j * currDt;
                        
                        sources_[i].beam().move(currPos, currPow, nowTime);
                        
                        //dimensionedVector Cxy(mesh_.C().x(), mesh_.C().y(), 0.0);
                        dimensionedVector currentPosition(dimLength, currPos);
                        
                        refinementField_
                            += pos0(searchRadius - mag(mesh_.C() - currentPosition));
                    }
                    
                    //- Integrate behind the beam to capture the tail
                    for(int j = 0; j < nRevSteps_; ++j)
                    {
                        scalar nowTime = currTime - j * currDt;
                        
                        if (nowTime < 0.0)
                            break;
                        
                        sources_[i].beam().move(currPos, currPow, nowTime);
                        
                        dimensionedVector currentPosition(dimLength, currPos);
                        
                        refinementField_
                            += pos0(searchRadius - mag(mesh_.C() - currentPosition));
                    }
                    
                    //- Increase search radius for next loop
                    //searchRadius *= 4.0;
                    searchRadius += charLen_;
                }
            }
        }
    }
    
    //- Rescale refinement field to be 1 or 0
    refinementField_ = pos(refinementField_);
    
    if (currStep_ == nSteps_)
        currStep_ = 0;
}

// ************************************************************************* //
