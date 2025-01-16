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
    // Store DOD from previous time step
    if (mesh().time().timeIndex() != timeIndex_)
    {
        DODoldTime_ = DOD_;

        DODFieldOldTime_ = DODField_;

        timeIndex_ = mesh().time().timeIndex();
    }

    // Calculate DOD
    DOD_ = DODoldTime_
         + dimensionedScalar
           (
               "DODNewTime",
                dimless,
                gSum(-j_().internalField() * mesh().V()/Dech_.value())
              * mesh().time().deltaT().value()
              / (QBat_.value())
            );

    DODField_ = DOD_;

    // DODField_ = DODFieldOldTime_
    //             + j_()
    //                *dimensionedScalar
    //                 (
    //                     "1byJdimensions",
    //                     dimensionSet(0, 2, 0, 0, 0, -1, 0),
    //                     1
    //                 )
    //                *(
    //                     gSum(mesh().V()) * mesh().time().deltaT().value()
    //                     /(QBat_.value() * Dech_.value())
    //                 );
    // boundMinMax
    // (
    //     DODField_,
    //     dimensionedScalar("minDOD", dimless, 0.0),
    //     dimensionedScalar("maxDOD", dimless, 1.0)
    // );

    Y_() = a0_*Foam::pow(DODField_, 0)
        + a1_*Foam::pow(DODField_, 1)
        + a2_*Foam::pow(DODField_, 2)
        + a3_*Foam::pow(DODField_, 3)
        + a4_*Foam::pow(DODField_, 4)
        + a5_*Foam::pow(DODField_, 5);

    Y_() *= Foam::exp(-C1_*(1/T_() - 1/TRef_));

    U_() = b0_*Foam::pow(DODField_, 0)
        + b1_*Foam::pow(DODField_, 1)
        + b2_*Foam::pow(DODField_, 2)
        + b3_*Foam::pow(DODField_, 3)
        + b4_*Foam::pow(DODField_, 4)
        + b5_*Foam::pow(DODField_, 5);

    U_() += C2_*(T_() - TRef_);

    j_() = Y_()*(faiPos_() - faiNeg_() - U_());

    Info << "DOD: " << DOD_.value() << endl;
    Info << "DOD field sum: " << gSum(DODField_.internalField()*mesh().V()/gSum(mesh().V())) << endl;

    Info << "minMax faiPos: " << gMin(faiPos_()) << " , " << gMax(faiPos_()) << endl;
    Info << "minMax faiNeg: " << gMin(faiNeg_()) << " , " << gMax(faiNeg_()) << endl;
    Info << "minMax U: " << gMin(U_()) << " , " << gMax(U_()) << endl;

    Info << "minMax Y: " << gMin(Y_()) << " , " << gMax(Y_()) << endl;

    Info << "minMax J: " << gMin(j_()) << " , " << gMax(j_()) << endl;

}


void Foam::regionTypes::NTGK::calculateThermalBehavior()
{

    volScalarField QEch = (1/Dech_)*j_()*(faiPos_() - faiNeg_() - U_() + C2_*T_());

    volScalarField Qohm = sigmaPos_*(fvc::grad(faiPos_())&fvc::grad(faiPos_()))
                        + sigmaNeg_*(fvc::grad(faiNeg_())&fvc::grad(faiNeg_()));

    ST_() = (QEch + Qohm);
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

Foam::tmp<fvScalarMatrix> Foam::regionTypes::NTGK::jPos()
{
    return
    (
      - fvm::Sp((1/Dp_)*Y_(), faiPos_())
      + (1/Dp_)*Y_()*(faiNeg_() + U_())
    );
}

Foam::tmp<fvScalarMatrix> Foam::regionTypes::NTGK::jNeg()
{
    return
    (
      - fvm::Sp((1/Dn_)*Y_(), faiNeg_())
      + (1/Dn_)*Y_()*(faiPos_() - U_())
    );
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

    timeIndex_(-1),

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
    TRef_(electrochemicalProperties_.lookup("TRef")),
    Dp_(electrochemicalProperties_.lookup("Dp")),
    Dn_(electrochemicalProperties_.lookup("Dn")),
    Dech_(electrochemicalProperties_.lookup("Dech")),
    a0_(electrochemicalProperties_.lookup("a0")),
    a1_(electrochemicalProperties_.lookup("a1")),
    a2_(electrochemicalProperties_.lookup("a2")),
    a3_(electrochemicalProperties_.lookup("a3")),
    a4_(electrochemicalProperties_.lookup("a4")),
    a5_(electrochemicalProperties_.lookup("a5")),
    b0_(electrochemicalProperties_.lookup("b0")),
    b1_(electrochemicalProperties_.lookup("b1")),
    b2_(electrochemicalProperties_.lookup("b2")),
    b3_(electrochemicalProperties_.lookup("b3")),
    b4_(electrochemicalProperties_.lookup("b4")),
    b5_(electrochemicalProperties_.lookup("b5")),
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
    QBat_(electrochemicalProperties_.lookup("QBat")),
    DOD_(electrochemicalProperties_.lookup("DOD")),
    DODoldTime_(DOD_),
    DODField_
    (
        IOobject
        (
            "DOD",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        DOD_
    ),
    DODFieldOldTime_
    (
        IOobject
        (
            "DODoldTime_",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        DOD_
    ),
    Y0_(a0_),
    U0_(b0_),
    Y_(nullptr),
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
    // set open circuit potential field
    Y_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "Y_ech",
        dimensionedScalar("Yinit", dimensionSet(-1, -4, 3, 0, 0, 2, 0), 1),
        true
    );

    // set open circuit potential field
    U_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "U",
        dimensionedScalar("Uinit", dimVoltage, 1),
        true
    );

    // set volumetric exchange current density field
    j_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "j",
        dimensionedScalar("jinit", dimensionSet(0, -2, 0, 0, 0, 1, 0), 0),
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

    // set dimensionless amount of Li-containing meta-stable species in SEI field
    cSEI_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "cSEI",
        dimensionedScalar
        (
            "cSEIinit",
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            thermalAbuseProperties_.lookup("cSEI")
        ),
        true
    );

    // set dimensionless amount of Li amount intercalacted within the carbon field
    cNE_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "cNE",
        dimensionedScalar
        (
            "cNEinit",
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            thermalAbuseProperties_.lookup("cNE")
        ),
        true
    );

    // set dimensionless measure of SEI layer thickness field
    tSEI_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "tSEI",
        dimensionedScalar
        (
            "tSEIinit",
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            thermalAbuseProperties_.lookup("tSEI")
        ),
        true
    );

    // set degree of conversion field
    alpha_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "alpha",
        dimensionedScalar
        (
            "alphainit",
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            thermalAbuseProperties_.lookup("alpha")
        ),
        true
    );

    // set dimensionless concentration of electrolyte field
    cELE_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "cELE",
        dimensionedScalar
        (
            "cELEinit",
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            thermalAbuseProperties_.lookup("cELE")
        ),
        true
    );

    // set temperature field
    T_ = lookupOrRead<volScalarField>(mesh(), "T");

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::NTGK::~NTGK()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::NTGK::correct()
{
    calculateElectrochemicalParameters();
    calculateThermalBehavior();
    calculateThermalAbuse();
}


Foam::scalar Foam::regionTypes::NTGK::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::NTGK::setCoupledEqns()
{
	faiPosEqn =
    (
        fvm::laplacian(sigmaPos_, faiPos(), "laplacian(sigma,fai)")
      ==
        // -(1/Dp_)*j_()
        jPos()
    );

    faiNegEqn =
    (
        fvm::laplacian(sigmaNeg_, faiNeg(), "laplacian(sigma,fai)")
      ==
        //(1/Dn_)*j_()
        jNeg()
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
    this->correct();
}

void Foam::regionTypes::NTGK::solveRegion()
{
    fvScalarMatrix cSEIEqn
    (
        fvm::ddt(cSEI())
      ==
       - RSEI_()
    );
    cSEIEqn.solve();

    fvScalarMatrix cNEEqn
    (
        fvm::ddt(cNE())
      ==
       - RNE_()
    );
    cNEEqn.solve();

    fvScalarMatrix tSEIEqn
    (
        fvm::ddt(tSEI())
      ==
        RNE_()
    );
    tSEIEqn.solve();

    fvScalarMatrix alphaEqn
    (
        fvm::ddt(alpha())
      ==
        RPE_()
    );
    alphaEqn.solve();

    fvScalarMatrix cELEEqn
    (
        fvm::ddt(cELE())
      ==
       - RELE_()
    );
    cELEEqn.solve();
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
