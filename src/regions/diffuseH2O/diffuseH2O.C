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
#include "diffuseH2O.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(diffuseH2O, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        diffuseH2O,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::diffuseH2O::diffuseH2O
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
    H2O_(nullptr),
    q_(nullptr)
{
    // set thermal diffusivity field
    Dp_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "Dpore",
        dimensionedScalar(sorbentProperties_.lookup("Dpore")),
        true
    );

    

    // set H2O concentration field and rection rate
    H2O_ = lookupOrRead<volScalarField>(mesh(), "H2O");
    q_ = lookupOrRead<volScalarField>(mesh(), "q");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::diffuseH2O::~diffuseH2O()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::diffuseH2O::correct()
{
    Dp_().correctBoundaryConditions();
}


Foam::scalar Foam::regionTypes::diffuseH2O::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::diffuseH2O::setCoupledEqns()
{
    H2OEqn =
    (
        fvm::ddt(H2O())
     ==
        fvm::laplacian(Dp_(), H2O(), "laplacian(Dp,H2O)") - (1-eps_/eps_)*(fvm::ddt(q_()))
    );

    fvScalarMatrices.set
    (
        H2O_().name()
      + mesh().name() + "Mesh"
      + diffuseH2O::typeName + "Type"
      + "Eqn",
        &H2OEqn()
    );
}

void Foam::regionTypes::diffuseH2O::postSolve()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseH2O::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseH2O::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseH2O::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseH2O::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseH2O::meshMotionCorrector()
{
    // do nothing, add as required
}

// ************************************************************************* //
