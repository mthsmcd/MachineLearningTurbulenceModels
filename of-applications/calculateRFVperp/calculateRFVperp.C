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
    calculateRFVperp

Description
    Calculates the nonlinear part of the Reynolds force vector tStar.

    You can find more details in:

    Brener et al. (2021) - Conditioning and accurate solutions of Reynolds 
    average Navier–Stokes equations with data-driven turbulence closures.
    DOI: 10.1017/jfm.2021.148

    Brener et al. (2024) - A highly accurate strategy for data-driven turbulence 
    modeling. 
    DOI: 10.1007/s40314-023-02547-9

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // OF headers
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
        // Set the time to the value of current time directory.
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

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

        Info<< "Reading/calculating the modified Reynolds force vector t"<< endl;
        volVectorField t
        (
            IOobject
            (
                "t",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            - fvc::div(phi,U) + fvc::laplacian(nu,U)
            // non-continuity correction
            + fvc::div(U)*U + (1/3)*nu*fvc::div(T(fvc::grad(U)))
        );

        Info<< "Calculating the strain-rate tensor S"<< endl;
        volTensorField gradU("gradU", fvc::grad(U));
        volSymmTensorField S("S",symm(gradU));

        volVectorField tStar
        (   "tStar",
            t + 2*fvc::div(nut*S)
        );
        tStar.write();
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
