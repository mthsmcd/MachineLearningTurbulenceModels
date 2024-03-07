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

Class
    Foam::incompressible::RASModels::gammaRST

Description
    Data-driven Reynolds Stress Model (RSM).

    Reynolds Stress Tensor R is calculated through an adapted RSM equation
    using a RANS turbulent viscosity nut.

    The adpated RSM equation does not contain the production terms involving
    the velocity gradient tensor. These terms are absorbed by the source-term
    Gamma. 
    Gamma is also composed by all terms of the canonical RSM equation that
    require modelling.

    Proposed in the paper:

    Macedo, Matheus S. S., et al. (2024) - A data-driven turbulence modeling 
    for the Reynolds stress tensor transport equation.

SourceFiles
    gammaRST.C

\*---------------------------------------------------------------------------*/

#ifndef gammaRST_H
#define gammaRST_H

#include "turbulentTransportModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class gammaRST Declaration
\*---------------------------------------------------------------------------*/

class gammaRST
:
    public eddyViscosity<incompressible::RASModel>
{

protected:

    // Protected data
    
        //Dummy declarations required by OpenFOAM's baseline ShihQuadraticKE
        //turbulence model used as a template
        volScalarField k_;
        volScalarField epsilon_;

        // Model coefficients
        dimensionedScalar couplingFactor_;

        // Fields
        volScalarField nut_;
	    volSymmTensorField R_;
        volSymmTensorField Gamma_;



    // Protected Member Functions

        virtual void correctNut();


public:

    //- Runtime type information
    TypeName("gammaRST");


    // Constructors

        //- Construct from components
        gammaRST
        (
            const geometricOneField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~gammaRST()
    {}


    // Member Functions

        //- Read RASProperties dictionary
        virtual bool read();

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Bound normal stresses from R tensor
        virtual void boundNormalStress(volSymmTensorField& R) const;

	    //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
