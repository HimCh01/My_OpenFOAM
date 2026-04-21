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

#include "cellSet.H" // Include the cellSet library
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Function to calculate time-dependent swirling velocity field
void calculateSwirlingVelocity(volVectorField& U, const Time& runTime, const scalar omega)
{
    forAll(U, celli)
    {
        const vector& pos = U.mesh().C()[celli]; // Cell center position
        scalar x = pos.x();
        scalar y = pos.y();
        scalar t = static_cast<double>(runTime.value());

        U[celli].x() = -omega * y * std::cos(omega * t);
        U[celli].y() =  omega * x * std::cos(omega * t);
        U[celli].z() = 0; // 2D swirling flow
    }
}
int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    
    // Time-varying angular velocity (swirling strength)
    scalar omega = 0.5;
    
     // Load the cellSet for the centerZone
    cellSet centerZone
    (
        mesh,
        "centerZone" // Name of the cellSet created by topoSet
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        fvModels.correct();

         // Update the swirling velocity field
       // calculateSwirlingVelocity(U, runTime, omega);

        // Compute volume flux from velocity field
        //surfaceScalarField phi = fvc::interpolate(U) & mesh.Sf();

		#include "createPhi.H"

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T)
              + fvm::div(phi, T)
              - fvm::laplacian(DT, T)
             ==
                fvModels.source(T)
            );

/*
            // Apply custom source term in the centerZone
            forAll(T.internalField(), celli)
			{
				if (centerZone.found(celli)) // Check if the cell is in centerZone
				{
					TEqn.source()[celli] += 10; // Source term value
				}
			}
*/
            
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
