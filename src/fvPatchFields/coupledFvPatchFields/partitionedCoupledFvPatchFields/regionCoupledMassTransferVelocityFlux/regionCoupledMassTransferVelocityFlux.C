/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "regionCoupledMassTransferVelocityFlux.H"
#include "regionCoupledVelocityValue.H"
#include "addToRunTimeSelectionTable.H"
#include "primitiveFieldsFwd.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCoupledMassTransferVelocityFlux::
regionCoupledMassTransferVelocityFlux
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    genericRegionCoupledFluxFvPatchField<vector>(p, iF)
{}


Foam::regionCoupledMassTransferVelocityFlux::
regionCoupledMassTransferVelocityFlux
(
    const regionCoupledMassTransferVelocityFlux& icvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    genericRegionCoupledFluxFvPatchField<vector>(icvf, p, iF, mapper)
{}


Foam::regionCoupledMassTransferVelocityFlux::
regionCoupledMassTransferVelocityFlux
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    genericRegionCoupledFluxFvPatchField<vector>(p, iF, dict)
{}


Foam::regionCoupledMassTransferVelocityFlux::
regionCoupledMassTransferVelocityFlux
(
    const regionCoupledMassTransferVelocityFlux& icvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    genericRegionCoupledFluxFvPatchField<vector>(icvf, iF)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> regionCoupledMassTransferVelocityFlux::fluxJump() const
{
    const vectorField& nf = capInterface().aMesh().faceAreaNormals();

    vectorField nEps = nf;

    if (capInterface().meshA() != refMesh())
    {
        nEps *= -1;
    }

    // Lookup neighbouring patch field
    const volVectorField& nbrUField =
        nbrMesh().lookupObject<volVectorField>
        (
            // same field name as on this side
            this->dimensionedInternalField().name()
        );

    // Get velocity face values from neighbour patch
    tmp<vectorField> tnbrU =
        refCast<const genericRegionCoupledJumpFvPatchField<vector>>
        (
            nbrPatch()
            .patchField<volVectorField, vector>(nbrUField)
        );

    const vectorField& nbrU = tnbrU();

    // Lookup neighbouring patch field
    const volVectorField& UField =
        refMesh().lookupObject<volVectorField>
        (
            // same field name as on this side
            this->dimensionedInternalField().name()
        );

    // Get velocity face values from this patch
    const vectorField& refU = *this;

    // Get flux face values from neighbour patch
    tmp<vectorField> tnbrUFlux =
        refCast<const genericRegionCoupledJumpFvPatchField<vector>>
        (
            nbrPatch()
            .patchField<volVectorField, vector>(nbrUField)
        ).flux();

    const vectorField& nbrUFlux = tnbrUFlux();

    // Calculate interpolated patch field
    vectorField UfluxNbrToOwn = interpolateFromNbrField<vector>(nbrUFlux);

    // Enforce flux matching
    UfluxNbrToOwn *= -1.0;

    // surface tension
    const areaScalarField& sigma = capInterface().sigma();

    // surface velocity terms
    const areaVectorField& Us = capInterface().Us();

    areaScalarField divSU = fac::div(Us);

    // [MP] Already done in fac::div
    //divSU.correctBoundaryConditions();

    areaTensorField gradSU = fac::grad(Us);

    vectorField surfaceTensionForce =
        sigma
       *fac::edgeIntegrate
        (
            capInterface().aMesh().Le()
            *capInterface().aMesh().edgeLengthCorrection()
        )().internalField();

    vectorField tangentialSurfaceTensionForce =
        surfaceTensionForce
      - sigma
       *capInterface().aMesh().faceCurvatures().internalField()*nf;
       
   dimensionedScalar muFluidNbr
    (
        nbrMesh().lookupObject<IOdictionary>("transportProperties")
        .lookup("mu")
    );

    dimensionedScalar muFluid
    (
        refMesh().lookupObject<IOdictionary>("transportProperties")
        .lookup("mu")
    );

    const scalarField mDots = interpolateFromNbrField<scalar>(massTrInterface().mDotS());

    dimensionedScalar rhoFluidNbr
    (
        nbrMesh().lookupObject<IOdictionary>("transportProperties")
        .lookup("rho")
    );

    dimensionedScalar rhoFluid
    (
        refMesh().lookupObject<IOdictionary>("transportProperties")
        .lookup("rho")
    );

    if (refMesh().foundObject<volScalarField>("pKin"))
    {
        dimensionedScalar nuFluid
        (
            refMesh().lookupObject<IOdictionary>("transportProperties")
            .lookup("nu")
        );   
        
        dimensionedScalar rhoFluid
        (
            refMesh().lookupObject<IOdictionary>("transportProperties")
            .lookup("rho")
        );  
        
        dimensionedScalar muFluidCalc = nuFluid*rhoFluid; 

        if (muFluidCalc.value() != muFluid.value())
        {
            FatalErrorInFunction  << this->typeName 
                << " Dynamic viscosity entry of region "
                << refMesh().name() << " " << muFluid.value() << " is different from the product of rho and nu "
                << muFluidCalc.value()
                << endl
                << exit(FatalError);        
        } 
    }  

    if (nbrMesh().foundObject<volScalarField>("pKin"))
    {
        dimensionedScalar nuFluidNbr
        (
            nbrMesh().lookupObject<IOdictionary>("transportProperties")
            .lookup("nu")
        );   
        
        dimensionedScalar rhoFluidNbr
        (
            nbrMesh().lookupObject<IOdictionary>("transportProperties")
            .lookup("rho")
        );  
        
        dimensionedScalar muFluidCalcNbr = nuFluidNbr*rhoFluidNbr; 

        if (muFluidCalcNbr.value() != muFluidNbr.value())
        {
            FatalErrorInFunction << this->typeName 
                << " Dynamic viscosity entry of region "
                << refMesh().name() << " " << muFluidNbr.value() << " is different from the product of rho and nu "
                << muFluidCalcNbr.value()
                << endl
                << exit(FatalError);        
        } 
    } 

    return
    (
      - nf*(nf & UfluxNbrToOwn)
      + tangentialSurfaceTensionForce
      - muFluid.value()*nf*divSU.internalField()
      + (muFluidNbr.value() - muFluid.value())*(gradSU.internalField()&nf)
      + ((mDots/rhoFluid.value()*nEps - refU) - (mDots/rhoFluidNbr.value()*nEps - nbrU))*mDots
    );
}

const regionInterfaces::capillaryInterface&
regionCoupledMassTransferVelocityFlux::capInterface() const
{
    if(   rgInterface().type()
       != regionInterfaces::capillaryInterface::typeName )
    {
        FatalErrorInFunction
            << this->typeName << " BC can only "
            << "be used in combination with a "
            << regionInterfaces::capillaryInterface::typeName
            << endl
            << exit(FatalError);
    }

    return refCast<const regionInterfaces::capillaryInterface>
        (
            rgInterface()
        );
}

const regionInterfaces::massTransferInterface&
regionCoupledMassTransferVelocityFlux::massTrInterface() const
{
    return refCast<const regionInterfaces::massTransferInterface>
        (
            rgInterface(regionInterfaces::massTransferInterface::typeName)
        );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchVectorField,
    regionCoupledMassTransferVelocityFlux
);

} // End namespace Foam

// ************************************************************************* //
