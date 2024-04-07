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

Application
    calculateGamma

Description
    Calculates the symmetric tensor Gamma, which is the source term of the 
    modified Reynolds stress tensor transport equation presented in 
    Macedo et al (2024).

    You can find more details in:

    Macedo et al (2024) - A data-driven turbulence modeling for the Reynolds 
    stress tensor transport equation (2024)
    DOI: 10.1002/fld.5284

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"  

    Info<< "Reading Transport Properties\n" << endl;
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading kinematic viscosity nu\n" << endl;
    dimensionedScalar nu ("nu", transportProperties);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // return set of times based on arglist
    const instantList& timeDirs = timeSelector::select0(runTime, args);

    forAll(timeDirs, timeI)
    {
	    runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << nl << endl;
         	
	    //Reading and loading fields
	    Info<< "Reading field U\n" << endl;
        volVectorField U
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        Info<< "Reading the eddy-viscosity nut\n" << endl;
        volScalarField nut
        (
            IOobject
            (
                "nut",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        #include "createPhi.H"

        Info<< "Reading the Reynolds stress tensor R\n" << endl;
        volSymmTensorField R
        (
            IOobject
            (
                "R",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
	
	    volScalarField nuEff("nuEff", nu + nut);
    
	    Info<< "Calculating the source-term Gamma\n" << endl;
	    volSymmTensorField Gamma
	    (
		    IOobject
		    (
			    "Gamma",
			    runTime.timeName(),
			    mesh,
			    IOobject::NO_READ,
			    IOobject::AUTO_WRITE
		    ),
		    fvc::div(phi,R)
            - fvc::div(U)*R
		    - fvc::laplacian(nuEff,R)
	    );
        Gamma.write();		
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
