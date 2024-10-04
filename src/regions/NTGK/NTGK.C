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
#include "NTGK.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "thermoPhysicalConstants.H"

using namespace thermoPhysicalConstant;


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(NTGK, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        NTGK,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::regionTypes::NTGK::calculateElectrochemicalParameters()
{
    dimensionedScalar dimVolt = 
        dimensionedScalar("dimVolt", dimensionSet(1, 2, -3, 0, 0, -1, 0), 1);
    
    volScalarField U0 = a0_
                      + a1_*DOD_()
                      + a2_*Foam::pow(DOD_(), 2) 
                      + a3_*Foam::pow(DOD_(), 3);
    
    volScalarField Y0 = a4_
                      + a5_*DOD_()
                      + a6_*Foam::pow(DOD_(), 2) 
                      + a7_*Foam::pow(DOD_(), 3)
                      + a8_*Foam::pow(DOD_(), 4)
                      + a9_*Foam::pow(DOD_(), 5);
    
    U_() = U0 - C2_*(T_() - TRef_);
    
    volScalarField Y = Y0*Foam::exp(-C1_*(1/T_() - 1/TRef_));
    
    j_() = spArea_*Y*(faiPos_() - faiNeg_() - U_())/dimVolt;
    
}


void Foam::regionTypes::NTGK::calculateThermalBehavior()
{    
    
    volScalarField QEch = j_()*(U_() - (faiPos_() - faiNeg_()));
    
    volScalarField Qohm = sigmaPos_*(fvc::grad(faiPos_())&fvc::grad(faiPos_()))
                        + sigmaNeg_*(fvc::grad(faiNeg_())&fvc::grad(faiNeg_()));
    
    ST_() = QEch + Qohm;
}


void Foam::regionTypes::NTGK::calculateThermalAbuse()
{
	
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::NTGK::NTGK
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
    
    electrochemicalProperties_
    (
        IOobject
        (
            "electrochemicalProperties",
            mesh().time().constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    
    sigmaPos_(transportProperties_.lookup("sigmaPos")),
    sigmaNeg_(transportProperties_.lookup("sigmaNeg")),
    rho_(transportProperties_.lookup("rho")),
    cp_(transportProperties_.lookup("cp")),
    k_(transportProperties_.lookup("k")),
    TRef_(transportProperties_.lookup("TRef")),
    spArea_(transportProperties_.lookup("spArea")),
    a0_(electrochemicalProperties_.lookup("a0")),
    a1_(electrochemicalProperties_.lookup("a1")),
    a2_(electrochemicalProperties_.lookup("a2")),
    a3_(electrochemicalProperties_.lookup("a3")),
    a4_(electrochemicalProperties_.lookup("a4")),
    a5_(electrochemicalProperties_.lookup("a5")),
    a6_(electrochemicalProperties_.lookup("a6")),
    a7_(electrochemicalProperties_.lookup("a7")),
    a8_(electrochemicalProperties_.lookup("a8")),
    a9_(electrochemicalProperties_.lookup("a9")),
    C1_(electrochemicalProperties_.lookup("C1")),
    C2_(electrochemicalProperties_.lookup("C2")),
    DOD_(nullptr),
    U_(nullptr),
    j_(nullptr),
    ST_(nullptr),
    faiPos_(nullptr),
    faiNeg_(nullptr),
    T_(nullptr)
{
    
    // set depth of discharge field
    DOD_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "DOD",
        dimensionedScalar("DODinit", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.7),
        true
    );

    // set open circuit potential field
    U_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "U",
        dimensionedScalar("Uinit", dimensionSet(1, 2, -3, 0, 0, -1, 0), 1),
        true
    );

    // set volumetric exchange current density field
    j_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "j",
        dimensionedScalar("j_init", dimensionSet(0, -3, 0, 0, 0, 1, 0), 0),
        true
    );

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

    // set negative electrode potential field
    faiNeg_ = lookupOrRead<volScalarField>(mesh(), "faiNeg");
    
    // set temperature field
    T_ = lookupOrRead<volScalarField>(mesh(), "T");

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::NTGK::~NTGK()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::NTGK::correct()
{
    
}


Foam::scalar Foam::regionTypes::NTGK::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::NTGK::setCoupledEqns()
{
	calculateElectrochemicalParameters();
    calculateThermalBehavior();
    calculateThermalAbuse();
     	
	faiPosEqn =
    (
        fvm::laplacian(sigmaPos_, faiPos(), "laplacian(sigma,fai)")
      ==
       - j_()
    );

    faiNegEqn =
    (
        fvm::laplacian(sigmaNeg_, faiNeg(), "laplacian(sigma,fai)")
      ==
        j_()
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
        faiPos_().name()
      + mesh().name() + "Mesh"
      + NTGK::typeName + "Type"
      + "Eqn",
        &faiPosEqn()
    );

    fvScalarMatrices.set
    (
        faiNeg_().name()
      + mesh().name() + "Mesh"
      + NTGK::typeName + "Type"
      + "Eqn",
        &faiNegEqn()
    );

    fvScalarMatrices.set
    (
        T_().name()
      + mesh().name() + "Mesh"
      + NTGK::typeName + "Type"
      + "Eqn",
        &TEqn()
    );
    
}

void Foam::regionTypes::NTGK::postSolve()
{
    // do nothing, add as required
}

void Foam::regionTypes::NTGK::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::NTGK::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::NTGK::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::NTGK::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::NTGK::meshMotionCorrector()
{
    // do nothing, add as required
}

// ************************************************************************* //
