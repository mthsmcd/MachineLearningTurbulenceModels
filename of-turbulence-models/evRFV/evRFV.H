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
    Foam::incompressible::RASModels::RFV

Description
    A variation of the method proposed by Cruz et al. (2019).
    
    A variation of the Reynolds Force Vector tStar is introduced into the mean
    momentum equations.
    An optimal turbulent viscosity is also introduced in the mean momentum balance.

    This method was proposed by:

    Brener et al. (2021) - Conditioning and accurate solutions of Reynolds 
    average Navier–Stokes equations with data-driven turbulence closures.
    DOI:10.1017/jfm.2021.148

    Based on:
    
    Cruz, Matheus A., et al. (2019) - The use of the Reynolds force vector in a
    physics informed machine learning approach for predictive turbulence modeling 
    DOI: 10.1016/j.compfluid.2019.104258

    To use the implicit method define implicitFactor = 1 in turbulenceProperties
    dictionary.

    To use the explicit method define implicitFactor = 0 in turbulenceProperties
    dictionary.

SourceFiles
    evRFV.C

\*---------------------------------------------------------------------------*/

#ifndef evRFV_H
#define evRFV_H

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
                      Class RFV Declaration
\*---------------------------------------------------------------------------*/

class evRFV
:
    public eddyViscosity<incompressible::RASModel>
{

protected:

    // Protected data

        // Fields

        //Dummy declarations required by OpenFOAM's baseline ShihQuadraticKE
        //turbulence model used as a template
        volScalarField k_;
        volScalarField epsilon_;

        // Model coefficients
        dimensionedScalar implicitFactor_;

        // Fields
        volScalarField nut_;
	    volVectorField t_;


    // Protected Member Functions

        virtual void correctNut();


public:

    //- Runtime type information
    TypeName("evRFV");


    // Constructors

        //- Construct from components
        evRFV
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
    virtual ~evRFV()
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
