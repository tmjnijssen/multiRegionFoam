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

#include "fvCFD.H"
#include "transportSpecie.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(transportSpecie, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        transportSpecie,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::transportSpecie::transportSpecie
(
    const Time& runTime,
    const word& regionName
)
:
    regionType(runTime, regionName),

    regionName_(regionName),

    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            mesh().time().constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    DCO2_(transportProperties_.lookup("DCO2")),
    DH2O_(transportProperties_.lookup("DH2O")),

    U_(nullptr),
    DeCO2_(nullptr),
    DeH2O_(nullptr),
    phi_(nullptr),
    CO2_(nullptr),
    H2O_(nullptr)
{
    // set velocity field
    // Postponing field creation since U is probably provided by
    // another regionType, e.g. icoFluid, and thus to be re-used.
    U_ = lookupOrRead<volVectorField>(mesh(), "U");

    // set flux field
    phi_ = lookupOrRead<surfaceScalarField>
    (
        mesh(),
        "phi",
        false,
        true,
        linearInterpolate(U_()) & mesh().Sf()
    );

    // set diffusivity fields
    DeCO2_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "DCO2", 
        DCO2_,
        true
    );

    DeH2O_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "DH2O", 
        DH2O_,
        true
    );

    // set concentration fields
    CO2_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "CO2", 
        true,
        true
    );

    H2O_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "H2O", 
        true,
        true
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::transportSpecie::~transportSpecie()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::transportSpecie::correct()
{
    DeCO2_().correctBoundaryConditions();
    DeH2O_().correctBoundaryConditions();
}


Foam::scalar Foam::regionTypes::transportSpecie::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::transportSpecie::setCoupledEqns(word fieldName)
{
    CO2Eqn =
    (  
       (
            fvm::ddt(CO2_())
          + fvm::div(phi_(), CO2_())
        )
     ==
        fvm::laplacian(DeCO2_(), CO2_())
    );

    fvScalarMatrices.set
    (
        CO2_().name()
      + mesh().name() + "Mesh"
      + transportSpecie::typeName + "Type"
      + "Eqn",
        &CO2Eqn()
    );

    H2OEqn =
    (  
       (
            fvm::ddt(H2O_())
          + fvm::div(phi_(), H2O_())
        )
     ==
        fvm::laplacian(DeH2O_(), H2O_())
    );

    fvScalarMatrices.insert
    (
        H2O_().name()
      + mesh().name() + "Mesh"
      + transportSpecie::typeName + "Type"
      + "Eqn",
        &H2OEqn()
    );
}

void Foam::regionTypes::transportSpecie::postSolve()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportSpecie::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportSpecie::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportSpecie::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportSpecie::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportSpecie::meshMotionCorrector()
{
    // do nothing, add as required
}

// ************************************************************************* //

