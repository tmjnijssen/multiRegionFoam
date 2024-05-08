/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "massTransferInterface.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionInterfaces
{
    defineTypeNameAndDebug(massTransferInterface, 0);

    addToRunTimeSelectionTable
    (
        regionInterfaceType,
        massTransferInterface,
        IOdictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionInterfaces::massTransferInterface::massTransferInterface
(
    const word& type,
    const dictionary& dict,
    const Time& runTime,
    const fvPatch& patchA,
    const fvPatch& patchB
)
:
    regionInterfaceType(type, dict, runTime, patchA, patchB),

    dict_(dict.subDict(type + "Coeffs")),

    hlv0_
    (
        dict_.lookup("hlv")
    ),
    hlv_
    (
        areaScalarField
        (
            IOobject
            (
                "hlvInterface",
                runTime.timeName(),
                aMesh().thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aMesh(),
            hlv0_,
            zeroGradientFaPatchScalarField::typeName
        )
    ),
    TSat0_
    (
        dict_.lookup("TSat")
    ),
    mDotsPtr_()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionInterfaces::massTransferInterface::clearOut() const
{
    mDotsPtr_.clear();

    regionInterfaceType::clearOut();
}

void Foam::regionInterfaces::massTransferInterface::makeMDotS() const
{
    if (!mDotsPtr_.empty())
    {
        FatalErrorIn("regionInterfaceType::makemDotS()")
            << "surface fluid flux already exists"
            << abort(FatalError);
    }
    // Set patch field types for Us
    wordList patchFieldTypes
    (
        aMesh().boundary().size(),
        zeroGradientFaPatchVectorField::typeName
    );

    forAll(aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == wedgeFaPatch::typeName
        )
        {
            patchFieldTypes[patchI] =
                wedgeFaPatchVectorField::typeName;
        }
        else
        {
            label ngbPolyPatchID =
                aMesh().boundary()[patchI].ngbPolyPatchIndex();

            if (ngbPolyPatchID != -1)
            {
                if
                (
                    meshA().boundary()[ngbPolyPatchID].type()
                 == wallFvPatch::typeName
                )
                {
                    WarningIn("regionInterfaceType::makeUs() const")
                        << "Patch neighbouring to interface is wall" << nl
                        << "Not appropriate for inlets/outlets" << nl
                        << endl;

                    patchFieldTypes[patchI] =
                        slipFaPatchVectorField::typeName;
                }
            }
        }
    }

    mDotsPtr_.reset
    (
        new areaScalarField
        (
            IOobject
            (
                patchA().name() + "MDotS",
                runTime().timeName(),
                meshA(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aMesh(),
            dimensioned<scalar>("mDot", dimMass/dimTime, 0.0),
            patchFieldTypes
        )
    );
}

void Foam::regionInterfaces::massTransferInterface::updateMDotS()
{
   
    scalarField sF = saturatedFlux();

    mDotS().internalField() = sF;
}

void Foam::regionInterfaces::massTransferInterface::correct()
{
    // Update transport properties
    Info << "CORRECT MASS TRANSFER" << endl;
    updateMDotS();

    // Update interface physics
    // TODO: call function to calculate new sigma_ with
    // interface equation of state for contaminated surfaces
}

Foam::scalarField Foam::regionInterfaces::massTransferInterface::saturatedFlux()
{
    volScalarField TA = meshA().lookupObject<volScalarField>("T");
    volScalarField TB = meshB().lookupObject<volScalarField>("T");

    dimensionedScalar kA
    (
        meshA().lookupObject<IOdictionary>("transportProperties")
        .lookup("k")
    );

    dimensionedScalar kB
    (
        meshB().lookupObject<IOdictionary>("transportProperties")
        .lookup("k")
    );
    
    TA.boundaryField().set
    (
        patchA().index() , fvPatchField<scalar>::New("fixedValue", meshA().boundary()[patchA().index() ],TA)
    );

    forAll (meshA().boundaryMesh()[patchA().index()],facei) 
    {
         TA.boundaryField()[patchA().index()][facei] = TSat0_.value();
    }

    TB.boundaryField().set
    (
        patchB().index() , fvPatchField<scalar>::New("fixedValue", meshA().boundary()[patchA().index() ],TB)
    );

    forAll (meshB().boundaryMesh()[patchA().index()],facei) 
    {
         TA.boundaryField()[patchB().index()][facei] = TSat0_.value();
    }

    surfaceScalarField snGradTA = fvc::snGrad(TA);
    surfaceScalarField snGradTB = fvc::snGrad(TB);
    
    scalarField saturateFlux = kA.value()*snGradTA.boundaryField()[patchA().index()] + kB.value()*snGradTB.boundaryField()[patchB().index()]; 
    return saturateFlux;
}

Foam::scalar Foam::regionInterfaces::massTransferInterface::getMinDeltaT()
{
    return GREAT;
}


// ************************************************************************* //
