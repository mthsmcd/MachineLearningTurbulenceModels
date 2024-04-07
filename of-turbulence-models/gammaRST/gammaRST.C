/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 M. S. S. Macedo, Coppe/UFRJ
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

#include "gammaRST.H"
#include "fvOptions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gammaRST, 0);
addToRunTimeSelectionTable(RASModel, gammaRST, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void gammaRST::correctNut()
{
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

gammaRST::gammaRST
(
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<incompressible::RASModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    //Dummy declarations required by OpenFOAM's baseline ShihQuadraticKE
    //turbulence model used as a template
    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	    mesh_,
        dimensionedScalar("k",dimensionSet(0,2,-2,0,0,0,0),1.0)
    ),

    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	    mesh_,
        dimensionedScalar("epsilon",dimensionSet(0,2,-3,0,0,0,0),1.0)
    ),

    //model coefficients
    couplingFactor_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "couplingFactor",
            this->coeffDict_,
            0.0
        )
    ),

    //fields
    nut_
    (
        IOobject
        (
            IOobject::groupName("nut", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    R_
    (
        IOobject
        (
            IOobject::groupName("R", alphaRhoPhi.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    Gamma_
    (
        IOobject
        (
            IOobject::groupName("Gamma", alphaRhoPhi.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )

{
    if (type == typeName)
    {
        
        if (couplingFactor_.value() < 0.0 || couplingFactor_.value() > 1.0)
        {
            FatalErrorInFunction
                << "couplingFactor = " << couplingFactor_
                << " is not in range 0 - 1" << nl
                << exit(FatalError);
        }
        printCoeffs(type);
        this->boundNormalStress(this->R_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool gammaRST::read()
{
    if (eddyViscosity<incompressible::RASModel>::read())
    {   
        couplingFactor_.readIfPresent(this->coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}

void gammaRST::boundNormalStress
(
    volSymmTensorField& R
) const
{
    scalar kMin = this->kMin_.value();

    R.max
    (
        dimensionedSymmTensor
        (
            "zero",
            R.dimensions(),
            symmTensor
            (
                kMin, -GREAT, -GREAT,
                kMin, -GREAT,
                kMin
            )
        )
    );
}

void gammaRST::correct()
{
    if (!turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volSymmTensorField& Gamma = this->Gamma_;
    volSymmTensorField& R = this->R_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<incompressible::RASModel>::correct();
    
    // Reynolds stress equation
    tmp<fvSymmTensorMatrix> REqn
    (
        fvm::ddt(alpha, rho, R)
      + fvm::div(alphaRhoPhi, R)
      - fvm::laplacian(alpha*rho*nuEff(), R)
      ==
        alpha*rho*Gamma
      + fvOptions(alpha, rho, R)
    );

    REqn.ref().relax();
    fvOptions.constrain(REqn.ref());
    solve(REqn);
    fvOptions.correct(R);

    this->boundNormalStress(R);
}

Foam::tmp<Foam::fvVectorMatrix>
gammaRST::divDevRhoReff
(
    volVectorField& U
) const
{

    if (couplingFactor_.value() > 0.0)
    {
        return
        (
            fvc::laplacian
            (
                (1.0 - couplingFactor_)*this->alpha_*this->rho_*this->nut(),
                U,
                "laplacian(nuEff,U)"
            )
            + fvc::div
            (
                this->alpha_*this->rho_*dev(R_)
                + couplingFactor_
                *this->alpha_*this->rho_*this->nut()*fvc::grad(U),
                "div(devRhoReff)"
            )
            - fvc::div(this->alpha_*this->rho_*this->nu()*dev2(T(fvc::grad(U))))
            - fvm::laplacian(this->alpha_*this->rho_*this->nuEff(), U)
        );
    }
    else
    {
        return
        (
            fvc::laplacian
            (
                this->alpha_*this->rho_*this->nut(),
                U,
                "laplacian(nuEff,U)"
            )
            + fvc::div(this->alpha_*this->rho_*dev(R_))
            - fvc::div(this->alpha_*this->rho_*this->nu()*dev2(T(fvc::grad(U))))
            - fvm::laplacian(this->alpha_*this->rho_*this->nuEff(), U)
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
