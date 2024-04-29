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

#include "regionCoupledMassTransferVelocityValue.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCoupledMassTransferVelocityValue::
regionCoupledMassTransferVelocityValue
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    genericRegionCoupledJumpFvPatchField<vector>(p, iF)
{}


Foam::regionCoupledMassTransferVelocityValue::
regionCoupledMassTransferVelocityValue
(
    const regionCoupledMassTransferVelocityValue& icvv,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    genericRegionCoupledJumpFvPatchField<vector>(icvv, p, iF, mapper)
{}


Foam::regionCoupledMassTransferVelocityValue::
regionCoupledMassTransferVelocityValue
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    genericRegionCoupledJumpFvPatchField<vector>(p, iF, dict)
{}


Foam::regionCoupledMassTransferVelocityValue::
regionCoupledMassTransferVelocityValue
(
    const regionCoupledMassTransferVelocityValue& icvv,
    const DimensionedField<vector, volMesh>& iF
)
:
    genericRegionCoupledJumpFvPatchField<vector>(icvv, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionCoupledMassTransferVelocityValue::updatePhi()
{
    //- Non-const access to flux on patch
    fvsPatchField<scalar>& patchPhiField = const_cast<fvsPatchField<scalar>& >
    (
        this->db().lookupObject<surfaceScalarField>("phi")
        .boundaryField()[this->patch().index()]
    );

    //- Get the flux on neighboring patch
    surfaceScalarField nbrPhi =
        nbrMesh().lookupObject<surfaceScalarField>("phi");

    //- Impose interpolated flux field
    patchPhiField = interpolateFromNbrField<scalar>
        (
            nbrPatch().patchField<surfaceScalarField, scalar>(nbrPhi)
        )*(-1.); // consider outer normals pointing in opposite directions
}


//- Zero velocity jump
tmp<vectorField> Foam::regionCoupledMassTransferVelocityValue::valueJump() const
{
    const vectorField nf = refMesh().boundary()[refPatchID()].nf();

    vectorField UsNbrToOwn = interpolateFromNbrField<vector>(capInterface().Us());

    const volVectorField& U =
        refMesh().objectRegistry::lookupObject<volVectorField>("U");
        
    // MP-Start
    //const volScalarField& rho =
        //refMesh().objectRegistry::lookupObject<volScalarField>("rho");


    scalarField meshPhi = 0.0*fvc::meshPhi(U)().boundaryField()[refPatchID()];
    
    if (refMesh().objectRegistry::foundObject<volScalarField>("rho"))
    {
        const volScalarField& rho =
            refMesh().objectRegistry::lookupObject<volScalarField>("rho");
            
        meshPhi = fvc::meshPhi(rho, U)().boundaryField()[refPatchID()];
    }
    else
    {
        meshPhi = fvc::meshPhi(U)().boundaryField()[refPatchID()];      
    }

    dimensionedScalar rhoFluid
    (
        refMesh().lookupObject<IOdictionary>("transportProperties")
        .lookup("rho")
    );

    dimensionedScalar massTransferRate
    (
        refMesh().lookupObject<IOdictionary>("transportProperties")
        .lookup("mDot")
    );

    Pout << "BOUNDARY" << endl;
    // MP-End
    return
    (
        nf*(-(nf & UsNbrToOwn)
         + meshPhi/
           refMesh().boundary()[refPatchID()].magSf()
         + massTransferRate.value()/rhoFluid.value())//*UCoeff.value()
    );
}

const regionInterfaces::capillaryInterface&
regionCoupledMassTransferVelocityValue::capInterface() const
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchVectorField,
    regionCoupledMassTransferVelocityValue
);

} // End namespace Foam

// ************************************************************************* //

