/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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

#include "parabolicVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parabolicVelocityFvPatchVectorField::
parabolicVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    origin_(Zero),
    axis_(Zero),
    umax_(Zero),
    radius_(Zero)
{}


Foam::parabolicVelocityFvPatchVectorField::
parabolicVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    origin_(dict.lookup("origin")),
    axis_(dict.lookup("axis")),
    umax_(dict.lookup<scalar>("umax")),
    radius_(dict.lookup<scalar>("radius"))    
   {}


Foam::parabolicVelocityFvPatchVectorField::
parabolicVelocityFvPatchVectorField
(
    const parabolicVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    umax_(ptf.umax_),
    radius_(ptf.radius_)
{}


Foam::parabolicVelocityFvPatchVectorField::
parabolicVelocityFvPatchVectorField
(
    const parabolicVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    umax_(ptf.umax_),
    radius_(ptf.radius_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::parabolicVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Compute the radial distance of each patch face from the origin
    vectorField r(patch().Cf() - origin_); // Vector from origin to patch center
    scalarField magr(mag(r));             // Magnitude (radial distance)

    // Compute the parabolic profile
    scalarField parabolicVelocity = umax_ * (1 - sqr(magr / radius_));

    // Normalize the axis direction
    vector normalizedAxis = axis_ / mag(axis_);

    // Assign the velocity to the patch
    operator==(parabolicVelocity * normalizedAxis);


    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::parabolicVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    writeEntry(os, "origin", origin_);
    writeEntry(os, "axis", axis_);
    writeEntry(os, "umax", umax_);
    writeEntry(os, "radius", radius_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       parabolicVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
