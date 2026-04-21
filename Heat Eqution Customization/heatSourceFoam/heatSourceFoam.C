/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

Application
    scalarTransportFoam

Description
    Solves the steady or transient transport equation for a passive scalar.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Function to define the time-dependent source term
dimensionedScalar computeQ(const Time& runTime)
{
    scalar Q_value = 10.0 * std::sin(2.0 * constant::mathematical::pi * static_cast<double>(runTime.value()));
    return dimensionedScalar("Q", dimEnergy*dimTemperature*dimTime/(dimMass*dimLength*dimLength), Q_value);
}


//dimensionedScalar rho("rho", dimDensity, transportProperties.lookup("rho"));
//dimensionedScalar cp("cp", dimEnergy/(dimMass*dimTemperature), transportProperties.lookup("cp"));


int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        fvModels.correct();

        while (simple.correctNonOrthogonal())
        {
			
			// Compute the time-dependent Q
			dimensionedScalar Q = computeQ(runTime);
        
            fvScalarMatrix TEqn
            (
                fvm::ddt(T)
              - fvm::laplacian(DT, T) + Q
             ==
                fvModels.source(T)
            );

            TEqn.relax();
            fvConstraints.constrain(TEqn);
            TEqn.solve();
            fvConstraints.constrain(T);
        }

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
