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
#include "transportH2O.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(transportH2O, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        transportH2O,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::transportH2O::transportH2O
(
    const Time& runTime,
    const word& regionName
)
:
    regionType(runTime, regionName),

    regionName_(regionName),

    U_(nullptr),
    De_(nullptr),
    rho_(nullptr),
    phi_(nullptr),
    H2O_(nullptr)
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

    H2O_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "H2O", 
        true,
        true
    );

    // set temperature field
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::transportH2O::~transportH2O()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::transportH2O::correct()
{
    De_().correctBoundaryConditions();
}


Foam::scalar Foam::regionTypes::transportH2O::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::transportH2O::setCoupledEqns(word fieldName)
{
    H2OEqn =
    (
        
       (
            fvm::ddt(H2O_())
          + fvm::div(phi_(), H2O())
        )
     ==
        fvm::laplacian(De_(), H2O())
    );

    fvScalarMatrices.set
    (
        H2O_().name()
      + mesh().name() + "Mesh"
      + transportH2O::typeName + "Type"
      + "Eqn",
        &H2OEqn()
    );
}

void Foam::regionTypes::transportH2O::postSolve()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportH2O::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportH2O::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportH2O::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportH2O::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportH2O::meshMotionCorrector()
{
    // do nothing, add as required
}

// ************************************************************************* //

