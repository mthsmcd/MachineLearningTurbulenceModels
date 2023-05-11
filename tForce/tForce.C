/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "tForce.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(tForce, 0);
addToRunTimeSelectionTable(RASModel, tForce, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void tForce::correctNut()
{
    this->nut_ = dimensionedScalar("nut", dimensionSet(0,2,-1,0,0,0,0), 0.0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tForce::tForce
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
    
    //Reynolds force vector
    t_
    (
        IOobject
        (
            IOobject::groupName("t", alphaRhoPhi.group()),
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
        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool tForce::read()
{
    if (eddyViscosity<incompressible::RASModel>::read())
    {   
        return true;
    }
    else
    {
        return false;
    }
}


void tForce::correct()
{
    if (!turbulence_)
    {
        return;
    }

    eddyViscosity<incompressible::RASModel>::correct();
}

Foam::tmp<Foam::fvVectorMatrix>
tForce::divDevRhoReff
(
    volVectorField& U
) const
{
    return
    (
      - fvc::div((this->alpha_*this->rho_*this->nu())*dev2(T(fvc::grad(U))))
      - fvm::laplacian(this->alpha_*this->rho_*this->nu(), U)
      + t_	
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
