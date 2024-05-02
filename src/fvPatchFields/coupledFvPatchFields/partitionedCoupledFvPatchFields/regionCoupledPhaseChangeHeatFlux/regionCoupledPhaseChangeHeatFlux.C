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

#include "regionCoupledPhaseChangeHeatFlux.H"
#include "regionCoupledScalarJump.H"
#include "addToRunTimeSelectionTable.H"
#include "heatTransferInterface.H"
#include "primitiveFieldsFwd.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCoupledPhaseChangeHeatFlux::
regionCoupledPhaseChangeHeatFlux
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    genericRegionCoupledFluxFvPatchField<scalar>(p, iF)
{}


Foam::regionCoupledPhaseChangeHeatFlux::
regionCoupledPhaseChangeHeatFlux
(
    const regionCoupledPhaseChangeHeatFlux& icpf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    genericRegionCoupledFluxFvPatchField<scalar>(icpf, p, iF, mapper)
{}


Foam::regionCoupledPhaseChangeHeatFlux::
regionCoupledPhaseChangeHeatFlux
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    genericRegionCoupledFluxFvPatchField<scalar>(p, iF, dict)
{}


Foam::regionCoupledPhaseChangeHeatFlux::
regionCoupledPhaseChangeHeatFlux
(
    const regionCoupledPhaseChangeHeatFlux& icpf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    genericRegionCoupledFluxFvPatchField<scalar>(icpf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> regionCoupledPhaseChangeHeatFlux::fluxJump() const
{

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    // Lookup neighbouring patch field
    volScalarField nbrTField =
        nbrMesh().lookupObject<volScalarField>
        (
            // same field name as on this side
            this->dimensionedInternalField().name()
        );

    // Interpolate flux face values from neighbour patch
    tmp<scalarField> tnbrTFlux =
        refCast<const genericRegionCoupledJumpFvPatchField<scalar>>
        (
            nbrPatch()
            .patchField<volScalarField, scalar>(nbrTField)
        ).flux();

    const scalarField& nbrTFlux = tnbrTFlux();

    // Calculate interpolated patch field
    scalarField TfluxNbrToOwn = interpolateFromNbrField<scalar>(nbrTFlux);

    // Enforce flux matching
    TfluxNbrToOwn *= -1.0;

    dimensionedScalar TSat
    (
        nbrMesh().lookupObject<IOdictionary>("transportProperties")
        .lookup("TSat")
    );

    dimensionedScalar k
    (
        mesh.lookupObject<IOdictionary>("transportProperties")
        .lookup("k")
    );

    dimensionedScalar nbrk
    (
        nbrMesh().lookupObject<IOdictionary>("transportProperties")
        .lookup("k")
    );

    // Lookup patch field
    volScalarField TSatField =
        mesh.lookupObject<volScalarField>
        (
            // same field name as on this side
            this->dimensionedInternalField().name()
        );

    
    
    TSatField.boundaryField().set(patch().index() , fvPatchField<scalar>::New("fixedValue", mesh.boundary()[patch().index() ],TSatField));

    forAll (mesh.boundaryMesh()[patch().index()],facei) 
    {
         TSatField.boundaryField()[patch().index()][facei] = TSat.value();
    }
    
    // Lookup neighbouring patch field
    volScalarField nbrTSatField =
        nbrMesh().lookupObject<volScalarField>
        (
            // same field name as on this side
            this->dimensionedInternalField().name()
        );

    forAll (nbrMesh().boundaryMesh()[nbrPatch().index()],facei) 
    {
         nbrTSatField.boundaryField()[nbrPatch().index()][facei] = TSat.value();
    }

    Info << patch().name() << endl;
    surfaceScalarField snGradT = fvc::snGrad(TSatField);
    surfaceScalarField nbrSnGradT = fvc::snGrad(nbrTSatField);
    scalarField nbrSnGradTInterpolated = interpolateFromNbrField<scalar>(nbrSnGradT.boundaryField()[nbrPatch().index()]);

    scalarField saturateFlux = k.value()*snGradT.boundaryField()[patch().index()] + nbrk.value()*nbrSnGradTInterpolated; //nbrSnGradT.boundaryField()[nbrPatch().index()];

    Info << k.value() << endl;
    return ( -TfluxNbrToOwn + saturateFlux);


}

const regionInterfaces::heatTransferInterface&
regionCoupledPhaseChangeHeatFlux::heatTransInterface() const
{
    if(   rgInterface().type()
       != regionInterfaces::heatTransferInterface::typeName )
    {
        FatalErrorInFunction
            << this->typeName << " BC can only "
            << "be used in combination with a "
            << regionInterfaces::heatTransferInterface::typeName
            << endl
            << exit(FatalError);
    }

    return refCast<const regionInterfaces::heatTransferInterface>
        (
            rgInterface()
        );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    regionCoupledPhaseChangeHeatFlux
);

} // End namespace Foam

// ************************************************************************* //
