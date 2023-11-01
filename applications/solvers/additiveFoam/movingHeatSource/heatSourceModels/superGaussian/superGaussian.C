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

#include "superGaussian.H"
#include "addToRunTimeSelectionTable.H"
#include "multiDirRefinement.H"
#include "labelVector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatSourceModels
{
    defineTypeNameAndDebug(superGaussian, 0);
    addToRunTimeSelectionTable(heatSourceModel, superGaussian, dictionary);
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatSourceModels::superGaussian::superGaussian
(
    const word& sourceName,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heatSourceModel(typeName, sourceName, dict, mesh),
    mesh_(mesh)
{
    k_ = heatSourceModelCoeffs_.lookup<scalar>("k");
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar
Foam::heatSourceModels::superGaussian::factor(const vector& d)
{
    vector s = dimensions_ / Foam::pow(2.0, 1.0/k_);

    scalar x = Foam::pow(magSqr(cmptDivide(d, s)), k_/2.0);

    return Foam::exp(-x);
}

inline Foam::dimensionedScalar
Foam::heatSourceModels::superGaussian::I0()
{
    const scalar power_ = movingBeam_->power();

    const scalar AR = 
        dimensions_.z() / min(dimensions_.x(), dimensions_.y());
    
    const scalar eta_ = absorptionModel_->eta(AR);

    vector s = dimensions_ / Foam::pow(2.0, 1.0/k_);

    const dimensionedScalar I0
    (
        "I0",
        dimPower / dimVolume,
        eta_*power_*k_
      / (s.x()*s.y()*s.z()*2.0*pi*Foam::tgamma(3.0/k_))
    );

    return I0;
}



/*
Foam::tmp<Foam::volScalarField>
Foam::heatSourceModels::superGaussian::qDot()
const
{
    vector s = dimensions_ / Foam::pow(2.0, 1.0/k_);

    scalar x = Foam::pow(magSqr(cmptDivide(d, s)), k_/2.0);

    return Foam::exp(-x);
}

inline Foam::dimensionedScalar
Foam::heatSourceModels::superGaussian::V0()
{
    vector s = dimensions_ / Foam::pow(2.0, 1.0/k_);

    const dimensionedScalar V0
    (
        "V0",
        dimVolume,
        (2.0 / 3.0)*s.x()*s.y()*s.z()*pi*Foam::tgamma(1.0 + 3.0/k_)
    );

    const scalar power_ = movingBeam_->power();
    
    if (power_ > small)
    {        
        const scalar AR = 
            dimensions_.z() / min(dimensions_.x(), dimensions_.y());
        
        const scalar eta_ = absorptionModel_->eta(AR);
        
        //- Calculate volumetric intesity
        const vector s = dimensions_ / Foam::pow(2.0, 1.0/k_);

        const dimensionedScalar I0
        (
            "I0",
            dimPower / dimVolume,
            eta_*power_*k_
          / (s.x()*s.y()*s.z()*2.0*pi*Foam::tgamma(3.0/k_))
        );

        // change sampling location for gaussian
        volScalarField factor
        (
            IOobject
            (
                "factor",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("Zero", dimless, 0.0)          
        );

        vector dx = vector(5e-6,  5e-6, 5e-6);

        const cellList& cells = mesh_.cells();
        const faceList& faces = mesh_.faces();
        const pointField& points = mesh_.points();

        forAll(mesh_.cells(), celli)
        {
            treeBoundBox beamBb
            (
                movingBeam_->position() - 3.0*dimensions_,
                movingBeam_->position() + 3.0*dimensions_
            );

            const cell& c = cells[celli];

            treeBoundBox cellBb(point::max, point::min);
            forAll(c, facei)
            {
                const face& f = faces[c[facei]];
                forAll(f, fp)
                {
                    cellBb.min() = min(cellBb.min(), points[f[fp]]);
                    cellBb.max() = max(cellBb.max(), points[f[fp]]);
                }
            }

            if (cellBb.overlaps(beamBb))
            {
                labelVector nPoints
                (
                    cmptDivide
                    (
                        (cellBb.span() + small*vector::one),
                        dx
                    )
                );

                nPoints = max(nPoints, labelVector(1, 1, 1));

                vector dxi = cmptDivide(cellBb.span(), vector(nPoints));

                scalar dVi = dxi.x() * dxi.y() * dxi.z();

                scalar Vi = 0;

                for (label k=0; k < nPoints.z(); ++k)
                {
                    for (label j=0; j < nPoints.y(); ++j)
                    {
                        for (label i=0; i < nPoints.x(); ++i)
                        {
                            const point pt
                            (
                                cellBb.max()
                              - cmptMultiply
                                (
                                    vector(i + 0.5, j + 0.5, k + 0.5),
                                    dxi
                                )
                            );

                            point d = cmptMag(pt - movingBeam_->position());

                            scalar f = magSqr(cmptDivide(d, s));

                            factor[celli] += 
                                Foam::exp(-Foam::pow(f, k_/2.0)) * dVi;

                            Vi += dVi;
                        }
                    }
                }

                factor[celli] /= Vi;
            }
        }

        qDot_ = I0 * factor;
    }

    return tqDot;
}
*/

bool Foam::heatSourceModels::superGaussian::read()
{
    if (heatSourceModel::read())
    {
        heatSourceModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        //- Mandatory entries
        heatSourceModelCoeffs_.lookup("k") >> k_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
