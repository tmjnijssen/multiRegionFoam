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

    volScalarField QEch = Foam::mag(j_()*(U_() - (faiPos_() - faiNeg_())));

    volScalarField Qohm = Foam::mag(sigmaPos_*(fvc::grad(faiPos_())&fvc::grad(faiPos_()))
                        + sigmaNeg_*(fvc::grad(faiNeg_())&fvc::grad(faiNeg_())));

    ST_() = QEch /*+ Qohm*/;
}


void Foam::regionTypes::NTGK::calculateThermalAbuse()
{
	// thermal abuse model for li-ion cells
    // source: Kim, G. H., Pesaran, A., & Spotnitz, R. (2007).
    // A three-dimensional thermal abuse model for lithium-ion cells.
    // Journal of power sources, 170(2), 476-489.
    const dimensionedScalar TSEIScalar = dimensionedScalar("TSEI", dimensionSet(0, 0, 0, 1, 0, 0, 0), 363.15);
    const dimensionedScalar TNEScalar = dimensionedScalar("TNE", dimensionSet(0, 0, 0, 1, 0, 0, 0), 393.15);
    const dimensionedScalar TELEScalar = dimensionedScalar("TELE", dimensionSet(0, 0, 0, 1, 0, 0, 0), 473.15);


    volScalarField TSEI = Tdummy_() + TSEIScalar;
    volScalarField TNE = Tdummy_() + TNEScalar;
    volScalarField TELE = Tdummy_() + TELEScalar;

    if(T_() > TSEI && T_() <= TNE)
    {
        RSEI_() = ASEI_*Foam::exp(-EASEI_/R/T_())*Foam::pow(cSEI_(), mSEI_);

        volScalarField QSEI = HSEI_*WC_*RSEI_();

        ST_() += QSEI;

    }
    if(T_() > TNE && T_() <= TELE)
    {
        RSEI_() = ASEI_*Foam::exp(-EASEI_/R/T_())*Foam::pow(cSEI_(), mSEI_);

        RNE_() = ANE_*Foam::exp(-tSEI_()/tSEIRef_)*Foam::pow(cNE_(), mNE_)*Foam::exp(-EANE_/R/T_());

        RPE_() = APE_*Foam::pow(alpha_(), mPE1_)*Foam::pow((1-alpha_()), mPE2_)*Foam::exp(-EAPE_/R/T_());

        volScalarField QSEI = HSEI_*WC_*RSEI_();

        volScalarField QNE = HNE_*WC_*RNE_();

        volScalarField QPE = HPE_*WP_*RPE_();

        ST_() += (QSEI + QNE + QPE);

    }
    if(T_() > TELE)
    {
        RSEI_() = ASEI_*Foam::exp(-EASEI_/R/T_())*Foam::pow(cSEI_(), mSEI_);

        RNE_() = ANE_*Foam::exp(-tSEI_()/tSEIRef_)*Foam::pow(cNE_(), mNE_)*Foam::exp(-EANE_/R/T_());

        RPE_() = APE_*Foam::pow(alpha_(), mPE1_)*Foam::pow((1-alpha_()), mPE2_)*Foam::exp(-EAPE_/R/T_());

        RELE_() = AELE_*Foam::exp(-EAELE_/R/T_())*Foam::pow(cELE_(), mELE_);

        volScalarField QSEI = HSEI_*WC_*RSEI_();

        volScalarField QNE = HNE_*WC_*RNE_();

        volScalarField QPE = HPE_*WP_*RPE_();

        volScalarField QELE = HELE_*WELE_*RELE_();

        ST_() += (QSEI + QNE + QPE + QELE);

    }
    else
    {
        ST_() = ST_()*1;
    }



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

    thermalAbuseProperties_
    (
        IOobject
        (
            "thermalAbuseProperties",
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
    ASEI_(thermalAbuseProperties_.lookup("ASEI")),
    ANE_(thermalAbuseProperties_.lookup("ANE")),
    APE_(thermalAbuseProperties_.lookup("APE")),
    AELE_(thermalAbuseProperties_.lookup("AELE")),
    EASEI_(thermalAbuseProperties_.lookup("EASEI")),
    EANE_(thermalAbuseProperties_.lookup("EANE")),
    EAPE_(thermalAbuseProperties_.lookup("EAPE")),
    EAELE_(thermalAbuseProperties_.lookup("EAELE")),
    mSEI_(thermalAbuseProperties_.lookup("mSEI")),
    mNE_(thermalAbuseProperties_.lookup("mNE")),
    mPE1_(thermalAbuseProperties_.lookup("mPE1")),
    mPE2_(thermalAbuseProperties_.lookup("mPE2")),
    mELE_(thermalAbuseProperties_.lookup("mELE")),
    HSEI_(thermalAbuseProperties_.lookup("HSEI")),
    HNE_(thermalAbuseProperties_.lookup("HNE")),
    HPE_(thermalAbuseProperties_.lookup("HPE")),
    HELE_(thermalAbuseProperties_.lookup("HELE")),
    WC_(thermalAbuseProperties_.lookup("WC")),
    WP_(thermalAbuseProperties_.lookup("WP")),
    WELE_(thermalAbuseProperties_.lookup("WELE")),
    tSEIRef_(thermalAbuseProperties_.lookup("tSEIRef")),
    DOD_(nullptr),
    U_(nullptr),
    j_(nullptr),
    RSEI_(nullptr),
    RNE_(nullptr),
    RPE_(nullptr),
    RELE_(nullptr),
    Tdummy_(nullptr),
    ST_(nullptr),
    faiPos_(nullptr),
    faiNeg_(nullptr),
    cSEI_(nullptr),
    cNE_(nullptr),
    tSEI_(nullptr),
    alpha_(nullptr),
    cELE_(nullptr),
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

    // set reaction rate SEI decomposition field
    RSEI_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "RSEI",
        dimensionedScalar("RSEI0", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        true
    );

    // set reaction rate negative solvent reaction field
    RNE_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "RNE",
        dimensionedScalar("RNE0", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        true
    );

    // set reaction rate positive solvent reaction field
    RPE_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "RPE",
        dimensionedScalar("RPE0", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        true
    );

    // set reaction rate electrolyte decomposition field
    RELE_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "RELE",
        dimensionedScalar("RELE0", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        true
    );

    // set reaction rate electrolyte decomposition field
    Tdummy_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "Tdummy",
        dimensionedScalar("Tdummy0", dimensionSet(0, 0, 0, 1, 0, 0, 0), 0),
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

    // set dimensionless amount of Li-containg meta-stabel species in SEI field
    cSEI_ = lookupOrRead<volScalarField>(mesh(), "cSEI");

    // set dimensionless amount of Li amount intercalacted within the carbon field
    cNE_ = lookupOrRead<volScalarField>(mesh(), "cNE");

    // set dimensionless measure of SEI layer thickness field
    tSEI_ = lookupOrRead<volScalarField>(mesh(), "tSEI");

    // set degree of conversion field
    alpha_ = lookupOrRead<volScalarField>(mesh(), "alpha");

    // set dimensionless concentration of electrolyte field
    cELE_ = lookupOrRead<volScalarField>(mesh(), "cELE");

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
        -1.0 *
        fvc::average
        (
            fvc::interpolate(j_())
        )
    );

    faiNegEqn =
    (
        fvm::laplacian(sigmaNeg_, faiNeg(), "laplacian(sigma,fai)")
      ==
        fvc::average
        (
            fvc::interpolate(j_())
        )
    );

    cSEIEqn =
    (
        fvm::ddt(1, cSEI())
      ==
       - RSEI_()
    );

    cNEEqn =
    (
        fvm::ddt(1, cNE())
      ==
       - RNE_()
    );

    tSEIEqn =
    (
        fvm::ddt(1, tSEI())
      ==
        RNE_()
    );

    alphaEqn =
    (
        fvm::ddt(1, alpha())
      ==
        RPE_()
    );

    cELEEqn =
    (
        fvm::ddt(1, cELE())
      ==
       - RELE_()
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
        cSEI_().name()
      + mesh().name() + "Mesh"
      + NTGK::typeName + "Type"
      + "Eqn",
        &cSEIEqn()
    );

    fvScalarMatrices.set
    (
        cNE_().name()
      + mesh().name() + "Mesh"
      + NTGK::typeName + "Type"
      + "Eqn",
        &cNEEqn()
    );

    fvScalarMatrices.set
    (
        tSEI_().name()
      + mesh().name() + "Mesh"
      + NTGK::typeName + "Type"
      + "Eqn",
        &tSEIEqn()
    );

    fvScalarMatrices.set
    (
        alpha_().name()
      + mesh().name() + "Mesh"
      + NTGK::typeName + "Type"
      + "Eqn",
        &alphaEqn()
    );

    fvScalarMatrices.set
    (
        cELE_().name()
      + mesh().name() + "Mesh"
      + NTGK::typeName + "Type"
      + "Eqn",
        &cELEEqn()
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
