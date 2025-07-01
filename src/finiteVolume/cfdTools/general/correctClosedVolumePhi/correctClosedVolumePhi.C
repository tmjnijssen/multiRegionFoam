/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "correctClosedVolumePhi.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "processorFvsPatchFields.H"
#include "inletOutletFvPatchFields.H"
#include "fvc.H"
#include"regionCoupledMassTransferVelocityValue.H"
#include"regionCoupledVelocityValue.H"
// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::correctClosedVolumePhi
(
    surfaceScalarField& phi,
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField& rAU
)
{
    label intPatchID_ = -1;

    forAll (U.boundaryField(), patchi)
    {
        const fvPatchVectorField& phip = U.boundaryField()[patchi];
        if
        (
            isA<regionCoupledMassTransferVelocityValue>(phip) ||
            isA<regionCoupledVelocityValue>(phip)

        )
        {
            intPatchID_ = patchi;
        }
    }

    phi.boundaryField()[intPatchID_] =
    (
        U.boundaryField()[intPatchID_]
        & phi.mesh().Sf().boundaryField()[intPatchID_]
    );

    scalarField weights =
        mag(phi.boundaryField()[intPatchID_]);

    if(mag(gSum(weights)) > VSMALL)
    {
        weights /= gSum(weights);
    }

    scalar uncorrectPhi = gSum(phi.boundaryField()[intPatchID_]);

    phi.boundaryField()[intPatchID_] -=
        weights*gSum(phi.boundaryField()[intPatchID_]);

    // phi.boundaryField()[intPatchID_] +=
    //     p.boundaryField()[intPatchID_].snGrad()
    //    *phi.mesh().magSf().boundaryField()[intPatchID_]
    //    *rAU.boundaryField()[intPatchID_];

    scalar correctPhi = gSum(phi.boundaryField()[intPatchID_]);

    if (fvMesh::debug)
    {
    Info<< "bool Foam::correctClosedVolumePhi(...) integral uncorrectPhi: " << uncorrectPhi
        << " integral correctPhi: " << correctPhi
        << endl;
    }
}


// ************************************************************************* //
