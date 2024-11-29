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
#include "conductPosPotentialTemperature.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "thermoPhysicalConstants.H"

using namespace thermoPhysicalConstant;


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(conductPosPotentialTemperature, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        conductPosPotentialTemperature,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::regionTypes::conductPosPotentialTemperature::calculateJouleHeating()
{

    volScalarField Qohm = sigmaPos_*(fvc::grad(faiPos_())&fvc::grad(faiPos_()));

    ST_() = Qohm;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::conductPosPotentialTemperature::conductPosPotentialTemperature
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

    sigmaPos_(transportProperties_.lookup("sigmaPos")),
    rho_(transportProperties_.lookup("rho")),
    cp_(transportProperties_.lookup("cp")),
    k_(transportProperties_.lookup("k")),
    C_(dimensionedScalar("C", dimensionSet(-1, -5, 4, 0, 0, 2, 0), 1)),
    ST_(nullptr),
    faiPos_(nullptr),
    T_(nullptr)
{

    // set summarized heat source terms field
    ST_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "ST",
        dimensionedScalar("STinit", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0),
        true
    );

    // set positive electrode potential field
    faiPos_ = lookupOrRead<volScalarField>(mesh(), "faiPos");

    // set temperature field
    T_ = lookupOrRead<volScalarField>(mesh(), "T");

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::conductPosPotentialTemperature::~conductPosPotentialTemperature()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::conductPosPotentialTemperature::correct()
{
    calculateJouleHeating();
}


Foam::scalar Foam::regionTypes::conductPosPotentialTemperature::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::conductPosPotentialTemperature::setCoupledEqns()
{
	faiPosEqn =
    (
       fvm::laplacian(sigmaPos_, faiPos(), "laplacian(sigma,fai)")
    );

    TEqn =
    (
         fvm::ddt(rho_*cp_, T(), "ddt(rho*cp,T)")
       - fvm::laplacian(k_, T(), "laplacian(k,T)")
       ==
         ST_()
    );

    fvScalarMatrices.set
    (
        faiPos_().name()
      + mesh().name() + "Mesh"
      + conductPosPotentialTemperature::typeName + "Type"
      + "Eqn",
        &faiPosEqn()
    );

    fvScalarMatrices.set
    (
        T_().name()
      + mesh().name() + "Mesh"
      + conductPosPotentialTemperature::typeName + "Type"
      + "Eqn",
        &TEqn()
    );

}

void Foam::regionTypes::conductPosPotentialTemperature::postSolve()
{
    this->correct();
}

void Foam::regionTypes::conductPosPotentialTemperature::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::conductPosPotentialTemperature::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::conductPosPotentialTemperature::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::conductPosPotentialTemperature::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::conductPosPotentialTemperature::meshMotionCorrector()
{
    // do nothing, add as required
}

// ************************************************************************* //
