/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2022 OpenFOAM Foundation
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
    Foam::parabolicVelocityFvPatchVectorField

Description
    This boundary condition describes an inlet vector boundary condition in
    swirl co-ordinates given a central axis, central point, axial, radial and
    tangential velocity profiles.

Usage
    \table
        Property       | Description             | Required    | Default value
        axis           | Axis of rotation        | yes         |
        origin         | Origin of rotation      | yes         |
        axialVelocity  | Axial velocity profile [m/s] | yes    |
        radialVelocity | Radial velocity profile [m/s] | yes  |
        tangentialVelocity | tangential velocity profile [m/s] | yes  |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            parabolicVelocity;
        axis            (0 0 1);
        origin          (0 0 0);
        axialVelocity   constant 30;
        radialVelocity  constant 10;
        tangentialVelocity  constant 100;
    }
    \endverbatim

    Note:
        The \c axialVelocity, \c radialVelocity and \c tangentialVelocity
        entries are Function1 types, able to describe time varying functions.
        The example above gives the usage for supplying constant values.

See also
    Foam::fixedValueFvPatchField
    Foam::Function1s

SourceFiles
    parabolicVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef parabolicVelocityFvPatchVectorField_H
#define parabolicVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
            Class parabolicVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class parabolicVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data
        const vector origin_;
        const vector axis_;
		scalar umax_;    // Maximum velocity at the center of the parabolic profile
		scalar radius_;  // Radius defining the extent of the parabolic profile		
		
public:

   //- Runtime type information
   TypeName("parabolicVelocity");


   // Constructors

        //- Construct from patch and internal field
        parabolicVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        parabolicVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateInletVelocityFvPatchVectorField
        //  onto a new patch
        parabolicVelocityFvPatchVectorField
        (
            const parabolicVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        parabolicVelocityFvPatchVectorField
        (
            const parabolicVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        parabolicVelocityFvPatchVectorField
        (
            const parabolicVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new parabolicVelocityFvPatchVectorField(*this, iF)
            );
        }

    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
