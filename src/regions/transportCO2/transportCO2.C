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
#include "transportCO2.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(transportCO2, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        transportCO2,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::transportCO2::transportCO2
(
    const Time& runTime,
    const word& regionName
)
:
    regionType(runTime, regionName),

    regionName_(regionName),

    U_(nullptr),
    CO2_(nullptr),
    De_(nullptr),
    phi_(nullptr),
    rho_(nullptr)
{
    // set velocity field
    // Postponing field creation since U is probably provided by
    // another regionType, e.g. icoFluid, and thus to be re-used.
    U_ = lookupOrRead<volVectorField>(mesh(), "U");
    rho_ = lookupOrRead<volScalarField>(mesh(), "rho");

    // set flux field
    phi_ = lookupOrRead<surfaceScalarField>
    (
        mesh(),
        "phi",
        false,
        true,
        linearInterpolate(rho_()*U_()) & mesh().Sf()
    );

    // set diffusivity field
    De_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "De", 
        true,
        true
    );

    CO2_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "CO2", 
        true,
        true
    );

    // set temperature field
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::transportCO2::~transportCO2()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::transportCO2::correct()
{
    kappa_().correctBoundaryConditions();
}


Foam::scalar Foam::regionTypes::transportCO2::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::transportCO2::setCoupledEqns()
{
    CO2Eqn =
    (
        
       (
            fvm::ddt(CO2_())
          + fvm::div(phi_(), CO2())
        )
     ==
        fvm::laplacian(kappa_(), CO2())
    );

    fvScalarMatrices.set
    (
        CO2_().name()
      + mesh().name() + "Mesh"
      + transportCO2::typeName + "Type"
      + "Eqn",
        &CO2Eqn()
    );
}

void Foam::regionTypes::transportCO2::postSolve()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportCO2::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportCO2::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportCO2::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportCO2::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportCO2::meshMotionCorrector()
{
    // do nothing, add as required
}

// ************************************************************************* //

