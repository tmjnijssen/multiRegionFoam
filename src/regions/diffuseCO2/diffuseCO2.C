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
#include "diffuseCO2.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(diffuseCO2, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        diffuseCO2,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::diffuseCO2::diffuseCO2
(
    const Time& runTime,
    const word& regionName
)
:
    regionType(runTime, regionName),

    regionName_(regionName),

    sorbentProperties_
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
    t_
    (
        IOobject
        (
            "t",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(transportProperties_.lookup("t"))
    ),
    eps_
    (
        IOobject
        (
            "eps",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(transportProperties_.lookup("eps"))
    ),

    Dpore_
    (
        IOobject
        (
            "Dpore",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(transportProperties_.lookup("Dpore"))
    ),

    Dp_(nullptr),
    CO2_(nullptr),
    q_(nullptr)
{
    // set thermal diffusivity field
    Dp_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "Dp",
        dimensionedScalar(transportProperties_.lookup("Dp")),
        true
    );

    // set CO2 concentration field
    CO2_ = lookupOrRead<volScalarField>(mesh(), "CO2");
    q_ = lookupOrRead<volScalarField>(mesh(), "CO2");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::diffuseCO2::~diffuseCO2()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::diffuseCO2::correct()
{
    Dp_().correctBoundaryConditions();
}


Foam::scalar Foam::regionTypes::diffuseCO2::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::diffuseCO2::setCoupledEqns()
{
    CO2Eqn =
    (
        fvm::ddt(CO2())
     ==
        fvm::laplacian(Dp_(), CO2(), "laplacian(Dp,CO2)") - (1-eps/eps)*(fvm::ddt(q))
    );

    fvScalarMatrices.set
    (
        CO2_().name()
      + mesh().name() + "Mesh"
      + diffuseCO2::typeName + "Type"
      + "Eqn",
        &CO2Eqn()
    );
}

void Foam::regionTypes::diffuseCO2::postSolve()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseCO2::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseCO2::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseCO2::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseCO2::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseCO2::meshMotionCorrector()
{
    // do nothing, add as required
}

// ************************************************************************* //
