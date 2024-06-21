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
    mDots0_
    (
        dict_.lookup("mDots")
    ),
    mDotInterface
    (
        areaScalarField
        (
            IOobject
            (
                "mDotInterface",
                runTime.timeName(),
                aMesh().thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aMesh(),
            mDots0_,
            zeroGradientFaPatchScalarField::typeName
        )
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
    mDotS().internalField() = sF/hlv_;
    mDotInterface = mDotS();
    Info << meshA().time().value() << " mDot Inteface: " << sum(mDotInterface*aMesh().S()*meshA().time().deltaT().value()).value() << " " << sum(mDotInterface*aMesh().S()).value() << " " << gSum(meshA().V()) << " " << gSum(meshB().V()) << endl;
}

void Foam::regionInterfaces::massTransferInterface::correct()
{
    // Update transport properties
    updateMDotS();

    // Update interface physics

}

Foam::scalarField Foam::regionInterfaces::massTransferInterface::saturatedFlux()
{
    // Define temperature fields
    volScalarField TA = meshA().lookupObject<volScalarField>("T");
    volScalarField TB = meshB().lookupObject<volScalarField>("T");

    // recall Thermal Conductivities
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

    // Set new BCs (Maybe unecessary)
    TA.boundaryField().set
    (
        patchA().index() , fvPatchField<scalar>::New("fixedValue", meshA().boundary()[patchA().index() ],TA)
    );

    TB.boundaryField().set
    (
        patchB().index() , fvPatchField<scalar>::New("fixedValue", meshB().boundary()[patchB().index() ],TB)
    );

    // Change BC values to saturation temperature
    forAll (meshA().boundaryMesh()[patchA().index()],facei) 
    {
         TA.boundaryField()[patchA().index()][facei] = TSat0_.value();
    }

    forAll (meshB().boundaryMesh()[patchB().index()],facei) 
    {
         TB.boundaryField()[patchB().index()][facei] = TSat0_.value();
    }

    // Compute gradients
    surfaceScalarField snGradTA = fvc::snGrad(TA);
    surfaceScalarField snGradTB = fvc::snGrad(TB);
   
    // Compute Saturated Heat Flux
    scalarField saturateFlux = kA.value()*snGradTA.boundaryField()[patchA().index()] + kB.value()*snGradTB.boundaryField()[patchB().index()]; 
    return saturateFlux;
}

Foam::scalar Foam::regionInterfaces::massTransferInterface::getMinDeltaT()
{
    return GREAT;
}


// ************************************************************************* //
