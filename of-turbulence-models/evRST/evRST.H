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

Class
    Foam::incompressible::RASModels::evRST

Description
    Data-driven Reynolds Stress tensor model.

    R is decomposed into a proportional-to S part and a orthogonal part "Rperp".
    Where S is the mean strain rate tensor.

    The proportional-to S segment is explicitly calculated every iteration
    using the previous solution for the velocity field U, and a optimal 
    turbulent viscosity nut.

    Both segments are injected into the mean momentum equation.

    Implicit diffusive term containing nut was proposed in the article:

    Wu et al. (2018) - Physics-informed machine learning approach for augmenting 
    turbulence models: A comprehensive framework.
    DOI: 10.1103/PhysRevFluids.3.074602

    Explicit diffusive term containing nut was proposed in the article:
    
    Brener et al. (2021) - Conditioning and accurate solutions of Reynolds 
    average Navier–Stokes equations with data-driven turbulence closures.
    DOI: 10.1017/jfm.2021.148

    To use the implicit method define implicitFactor = 1 in turbulenceProperties
    dictionary.

    To use the explicit method define implicitFactor = 0 in turbulenceProperties
    dictionary.

SourceFiles
    evRST.C

\*---------------------------------------------------------------------------*/

#ifndef evRST_H
#define evRST_H

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
                      Class RST Declaration
\*---------------------------------------------------------------------------*/

class evRST
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
        dimensionedScalar implicitFactor_;

        // Fields
        volScalarField nut_;
	    volSymmTensorField R_;


    // Protected Member Functions

        virtual void correctNut();


public:

    //- Runtime type information
    TypeName("evRST");


    // Constructors

        //- Construct from components
        evRST
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
    virtual ~evRST()
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
