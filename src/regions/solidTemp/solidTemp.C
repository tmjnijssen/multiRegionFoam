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
#include "solidTemp.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(solidTemp, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        solidTemp,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::solidTemp::solidTemp
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
    cv_
    (
        IOobject
        (
            "cv",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(transportProperties_.lookup("cv"))
    ),
    rho_
    (
        IOobject
        (
            "rho",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(transportProperties_.lookup("rho"))
    ),

    kappa_(nullptr),
    source_(nullptr),
    T_(nullptr)
{
    // set thermal diffusivity field
    kappa_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "k",
        dimensionedScalar(transportProperties_.lookup("k")),
        true
    );

    // read source field
    source_ = lookupOrRead<volScalarField>
    (
        mesh(),
        transportProperties_.lookupOrDefault<word>("heatSourceName", "heatSource"),
        dimensionedScalar("heatSource", dimEnergy/dimTime/dimVolume, 0.0),
        true
    );

    // set temperature field
    T_ = lookupOrRead<volScalarField>(mesh(), "T");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::solidTemp::~solidTemp()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::solidTemp::correct()
{
    kappa_().correctBoundaryConditions();
}


Foam::scalar Foam::regionTypes::solidTemp::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::solidTemp::setCoupledEqns(word fieldName)
{

}

void Foam::regionTypes::solidTemp::postSolve()
{
    // do nothing, add as required
}

void Foam::regionTypes::solidTemp::solveRegion()
{
    TEqn =
    (
        fvm::ddt(rho_*cv_, T())
     ==
        fvm::laplacian(kappa_(), T(), "laplacian(k,T)")
        + source_()
    );

    fvScalarMatrices.set
    (
        T_().name()
      + mesh().name() + "Mesh"
      + solidTemp::typeName + "Type"
      + "Eqn",
        &TEqn()
    );
}

void Foam::regionTypes::solidTemp::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::solidTemp::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::solidTemp::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::solidTemp::meshMotionCorrector()
{
    // do nothing, add as required
}

// ************************************************************************* //
