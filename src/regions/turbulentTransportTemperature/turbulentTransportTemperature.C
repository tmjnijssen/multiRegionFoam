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
#include "turbulentTransportTemperature.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(turbulentTransportTemperature, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        turbulentTransportTemperature,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::turbulentTransportTemperature::turbulentTransportTemperature
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

    Pr_(transportProperties_.lookup("Pr")),
    Prt_(transportProperties_.lookup("Prt")),

    U_(nullptr),
    kappa_(nullptr),
    phi_(nullptr),
    nu_(nullptr),
    nut_(nullptr),
    T_(nullptr)
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

    nu_ = lookupOrRead<volScalarField>(mesh(), "nu");
    nut_ = lookupOrRead<volScalarField>(mesh(), "nut");

    // set thermal conductivity field
    kappa_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "kappa",
        false,
        false,
        nu_()/Pr_
    );

    // set temperature field
    T_ = lookupOrRead<volScalarField>(mesh(), "T");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::turbulentTransportTemperature::~turbulentTransportTemperature()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::turbulentTransportTemperature::correct()
{
    kappa_().correctBoundaryConditions();
}


Foam::scalar Foam::regionTypes::turbulentTransportTemperature::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::turbulentTransportTemperature::setCoupledEqns()
{
    kappa_() = nu_()/Pr_ + nut_()/Prt_;

    TEqn =
    (
        fvm::ddt(T())
      + fvm::div(phi_(), T())
     ==
        fvm::laplacian(kappa_(), T())
    );

    fvScalarMatrices.set
    (
        T_().name()
      + mesh().name() + "Mesh"
      + turbulentTransportTemperature::typeName + "Type"
      + "Eqn",
        &TEqn()
    );
}

void Foam::regionTypes::turbulentTransportTemperature::postSolve()
{
    // do nothing, add as required
}

void Foam::regionTypes::turbulentTransportTemperature::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::turbulentTransportTemperature::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::turbulentTransportTemperature::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::turbulentTransportTemperature::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::turbulentTransportTemperature::meshMotionCorrector()
{
    // do nothing, add as required
}

// ************************************************************************* //
