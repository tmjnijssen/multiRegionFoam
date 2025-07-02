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
#include "conductNegPotentialTemperature.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "thermoPhysicalConstants.H"

using namespace thermoPhysicalConstant;


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(conductNegPotentialTemperature, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        conductNegPotentialTemperature,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::regionTypes::conductNegPotentialTemperature::calculateJouleHeating()
{

    volScalarField Qohm = sigmaNeg_*(fvc::grad(phiNeg_())&fvc::grad(phiNeg_()));

    ST_() = Qohm;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::conductNegPotentialTemperature::conductNegPotentialTemperature
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

    sigmaNeg_(transportProperties_.lookup("sigmaNeg")),
    rho_(transportProperties_.lookup("rho")),
    cp_(transportProperties_.lookup("cp")),
    k_(transportProperties_.lookup("k")),
    C_(dimensionedScalar("C", dimensionSet(-1, -5, 4, 0, 0, 2, 0), 1)),
    ST_(nullptr),
    phiNeg_(nullptr),
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

    // set negative electrode potential field
    phiNeg_ = lookupOrRead<volScalarField>(mesh(), "phiNeg");

    // set temperature field
    T_ = lookupOrRead<volScalarField>(mesh(), "T");

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::conductNegPotentialTemperature::~conductNegPotentialTemperature()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::conductNegPotentialTemperature::correct()
{
    calculateJouleHeating();
}


Foam::scalar Foam::regionTypes::conductNegPotentialTemperature::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::conductNegPotentialTemperature::setCoupledEqns()
{
	phiNegEqn =
    (
         C_*fvm::ddt(phiNeg(), "fai")
       - fvm::laplacian(sigmaNeg_, phiNeg(), "laplacian(sigma,fai)")
    );

    TEqn =
    (
         fvm::ddt(rho_*cp_, T())
       - fvm::laplacian(k_, T(), "laplacian(k,T)")
       ==
         ST_()
    );

    fvScalarMatrices.set
    (
        phiNeg_().name()
      + mesh().name() + "Mesh"
      + conductNegPotentialTemperature::typeName + "Type"
      + "Eqn",
        &phiNegEqn()
    );

    fvScalarMatrices.set
    (
        T_().name()
      + mesh().name() + "Mesh"
      + conductNegPotentialTemperature::typeName + "Type"
      + "Eqn",
        &TEqn()
    );

}

void Foam::regionTypes::conductNegPotentialTemperature::postSolve()
{
    this->correct();
}

void Foam::regionTypes::conductNegPotentialTemperature::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::conductNegPotentialTemperature::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::conductNegPotentialTemperature::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::conductNegPotentialTemperature::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::conductNegPotentialTemperature::meshMotionCorrector()
{
    // do nothing, add as required
}

// ************************************************************************* //
