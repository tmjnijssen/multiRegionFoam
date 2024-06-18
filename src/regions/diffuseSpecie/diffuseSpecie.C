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
#include "diffuseSpecie.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(diffuseSpecie, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        diffuseSpecie,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::diffuseSpecie::diffuseSpecie
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
            "sorbentProperties",
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
        dimensionedScalar(sorbentProperties_.lookup("t"))
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
        dimensionedScalar(sorbentProperties_.lookup("eps"))
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
        dimensionedScalar(sorbentProperties_.lookup("Dpore"))
    ),

    Dp_(nullptr),
    CO2_(nullptr),
    H2O_(nullptr),
    dqdtCO2_(nullptr),
    dqdtH2O_(nullptr),
    qHCO3_(nullptr),
    qR2NCO2_(nullptr)
{
    // set thermal diffusivity field
    Dp_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "Dpore",
        dimensionedScalar(sorbentProperties_.lookup("Dpore")),
        true
    );

    

    // set specie concentration fields and loadings in sorbent
    CO2_ = lookupOrRead<volScalarField>(mesh(), "CO2");
    H2O_ = lookupOrRead<volScalarField>(mesh(), "H2O");
    dqdtH2O_ = lookupOrRead<volScalarField>(mesh(), "dqdtH2O");
    dqdtCO2_ = lookupOrRead<volScalarField>(mesh(), "dqdtCO2");
    qR2NCO2_ = lookupOrRead<volScalarField>(mesh(), "R2NCO2");
    qHCO3_ = lookupOrRead<volScalarField>(mesh(), "HCO3");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::diffuseSpecie::~diffuseSpecie()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::diffuseSpecie::correct()
{
    Dp_().correctBoundaryConditions();
}


Foam::scalar Foam::regionTypes::diffuseSpecie::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::diffuseSpecie::setCoupledEqns()
{
    CO2Eqn =
    (
        fvm::ddt(CO2_())
     ==
        fvm::laplacian(Dp_(), CO2_(), "laplacian(Dp,CO2)") - (1-eps_/eps_)*(fvm::ddt(qR2NCO2_()))
    );

    H2OEqn =
    (
        fvm::ddt(H2O_())
     ==
        fvm::laplacian(Dp_(), H2O_(), "laplacian(Dp,H2O)") - (1-eps_/eps_)*(fvm::ddt(qHCO3_()))
    );

    fvScalarMatrices.set
    (
        (CO2_().name(), H2O_().name())
      + mesh().name() + "Mesh"
      + diffuseSpecie::typeName + "Type"
      + "Eqn",
        (&CO2Eqn(), &H2OEqn())
    );
}

void Foam::regionTypes::diffuseSpecie::postSolve()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseSpecie::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseSpecie::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseSpecie::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseSpecie::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseSpecie::meshMotionCorrector()
{
    // do nothing, add as required
}

// ************************************************************************* //
