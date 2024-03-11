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
#include "transportTemperature.H"
#include "rhoPimpleFluid.H"
#include "basicPsiThermo.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(transportTemperature, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        transportTemperature,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::transportTemperature::transportTemperature
(
    const Time& runTime,
    const word& regionName
)
:
    regionType(runTime, regionName),

    regionName_(regionName),

    U_(nullptr),
    rho_(nullptr),
    cp_(nullptr),
    k_(nullptr),
    phi_(nullptr),
    T_(nullptr)
{
    // set velocity field
    // Postponing field creation since U is probably provided by
    // another regionType, e.g. icoFluid, and thus to be re-used.
    U_ = lookupOrRead<volVectorField>(mesh(), "U");
    rho_ = lookupOrRead<volScalarField>(mesh(), "rho");
    cp_ = lookupOrRead<volScalarField>(mesh(), "cp");
    k_ = lookupOrRead<volScalarField>(mesh(), "k");

    // set flux field
    phi_ = lookupOrRead<surfaceScalarField>
    (
        mesh(),
        "phi",
        false,
        true,
        linearInterpolate(rho_()*U_()) & mesh().Sf()
    );

    // set temperature field
    T_ = lookupOrRead<volScalarField>(mesh(), "T");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::transportTemperature::~transportTemperature()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::transportTemperature::correct()
{
    k_().correctBoundaryConditions();
}


Foam::scalar Foam::regionTypes::transportTemperature::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::transportTemperature::setCoupledEqns()
{
    TEqn =
    (
        rho_()*cp_()
       *(
            fvm::ddt(T())
          + fvm::div(phi_(), T())
        )
     ==
        fvm::laplacian(k_(), T())
    );

    fvScalarMatrices.set
    (
        T_().name()
      + mesh().name() + "Mesh"
      + transportTemperature::typeName + "Type"
      + "Eqn",
        &TEqn()
    );
}

void Foam::regionTypes::transportTemperature::postSolve()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportTemperature::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportTemperature::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportTemperature::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportTemperature::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::transportTemperature::meshMotionCorrector()
{
    // do nothing, add as required
}

// ************************************************************************* //
