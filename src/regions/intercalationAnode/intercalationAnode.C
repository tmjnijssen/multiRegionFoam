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
#include "intercalationAnode.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "thermoPhysicalConstants.H"

using namespace thermoPhysicalConstant;


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(intercalationAnode, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        intercalationAnode,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::regionTypes::intercalationAnode::calculateOCV()
{
    // TO-DO: Exclude empirical models for reversible electrode potential calculation to sub classes (need: SoC, IApp as input variables!)

    SoC_() = cSE_() / cSMax_();
    dimensionedScalar dimensionVolt =
        dimensionedScalar("dimensionVolt", dimensionSet(1, 2, -3, 0, 0, -1, 0), 1);

    dimensionedScalar dimensionKelvin =
        dimensionedScalar("dimensionKelvin", dimensionSet(0, 0, 0, 1, 0, 0, 0), 1);


    if(electrodeMaterial_ == "graphite" || "Graphite")
    {
        //equilibrium potential of solid phase
        //data source: Mercer M P, Peng C, Soares C, et al. Voltage hysteresis during lithiation/delithiation of graphite associated with meta-stable carbon stackings[J]. Journal of Materials Chemistry A, 2021, 9(1): 492-504.
        scalar p1 = -1.031044745141316e+03;
        scalar p2 = 4.705957629441347e+03;
        scalar p3 = -8.994928528706332e+03;
        scalar p4 = 9.321254530787901e+03;
        scalar p5 = -5.664051502639276e+03;
        scalar p6 = 2.039265768656875e+03;
        scalar p7 = -4.203985227846820e+02;
        scalar p8 = 46.673532522233998;
        scalar p9 = -3.026252647513182;
        scalar p10 = 0.315992994754579;

        EEqRev_() = (p1*Foam::pow(SoC_(),9)
                + p2*Foam::pow(SoC_(),8)
                + p3*Foam::pow(SoC_(),7)
                + p4*Foam::pow(SoC_(),6)
                + p5*Foam::pow(SoC_(),5)
                + p6*Foam::pow(SoC_(),4)
                + p7*Foam::pow(SoC_(),3)
                + p8*Foam::pow(SoC_(),2)
                + p9*SoC_()
                + p10)*dimensionVolt;
    }
    if(electrodeMaterial_ == "silicium" || "Silicium")
    {
        // Calculation of solid equilibrium potential for silicium electrode
        // data source: Verbrugge M, Baker D, Xiao X. Formulation for the treatment of multiple electrochemical reactions and associated speciation for the Lithium-Silicon electrode[J]. Journal of The Electrochemical Society, 2015, 163(2): A262.
        scalar p1Lit = -96.63;
        scalar p2Lit = 372.6;
        scalar p3Lit = -587.6;
        scalar p4Lit = 489.9;
        scalar p5Lit = -232.8;
        scalar p6Lit = 62.99;
        scalar p7Lit = -9.286;
        scalar p8Lit = 0.8633;

        volScalarField EEqLit = (p1Lit * Foam::pow(SoC_(), 7)
                + p2Lit*Foam::pow(SoC_(), 6)
                + p3Lit*Foam::pow(SoC_(), 5)
                + p4Lit*Foam::pow(SoC_(), 4)
                + p5Lit*Foam::pow(SoC_(), 3)
                + p6Lit*Foam::pow(SoC_(), 2)
                + p7Lit*SoC_()
                + p8Lit)*dimensionVolt;

        scalar p1Delit = -51.02;
        scalar p2Delit = 161.3;
        scalar p3Delit = -205.7;
        scalar p4Delit = 140.2;
        scalar p5Delit = -58.76;
        scalar p6Delit =16.87;
        scalar p7Delit = -3.792;
        scalar p8Delit = 0.9937;

        volScalarField EEqDelit = (p1Delit * Foam::pow(SoC_(), 7)
                 + p2Delit*Foam::pow(SoC_(), 6)
                 + p3Delit*Foam::pow(SoC_(), 5)
                 + p4Delit*Foam::pow(SoC_(), 4)
                 + p5Delit*Foam::pow(SoC_(), 3)
                 + p6Delit*Foam::pow(SoC_(), 2)
                 + p7Delit*SoC_()
                 + p8Delit)*dimensionVolt;

        EEqRev_() = EEqLit * (IApp_.value() > 0) + EEqDelit * (IApp_.value() <= 0);
    }
    if(electrodeMaterial_ == "LiIonSimBa" || "liionsimba" || "lionsimba" || "Lionsimba" || "LionSimba")
    {
        //equilibrium potential of cathode material
        //data source: Torchio M, Magni L, Gopaluni R B, et al. Lionsimba: a matlab framework based on a finite volume model suitable for li-ion battery design, simulation, and control[J]. Journal of The Electrochemical Society, 2016, 163(7): A1192.
        volScalarField EEqRevRef = (0.7222
                  + 0.1387*SoC_()
                  + 0.029*Foam::pow(SoC_(), 0.5)
                  - 0.0172/SoC_()
                  + 0.0019/Foam::pow(SoC_(), 1.5)
                  + 0.2808*Foam::exp(0.9 - 15*SoC_())
                  - 0.7984*Foam::exp(0.4465*SoC_() - 0.4108))
                  *dimensionVolt;

        dEEqdT_() = (0.001*(0.005269056
                  + 3.299265709*SoC_()
                  - 91.79325798*Foam::pow(SoC_(), 2)
                  + 1004.911008*Foam::pow(SoC_(), 3)
                  - 5812.278127*Foam::pow(SoC_(), 4)
                  + 19329.75490*Foam::pow(SoC_(), 5)
                  - 37147.89470*Foam::pow(SoC_(), 6)
                  + 38379.18127*Foam::pow(SoC_(), 7)
                  - 16515.05308*Foam::pow(SoC_(), 8))
                  /(1 - 48.09287227*SoC_()
                  + 1017.234804*Foam::pow(SoC_(), 2)
                  - 10481.80419*Foam::pow(SoC_(), 3)
                  + 59431.30000*Foam::pow(SoC_(), 4)
                  - 195881.6488*Foam::pow(SoC_(), 5)
                  + 374577.3152*Foam::pow(SoC_(), 6)
                  - 385821.1607*Foam::pow(SoC_(), 7)
                  + 165708.8597*Foam::pow(SoC_(), 8)))
                  *dimensionVolt/dimensionKelvin;

        EEqRev_() = EEqRevRef + (T_() - TRef_)*dEEqdT_();
    }
    else
    {
        Info << "No valid material type input!" << endl;
        Info << "Check materialType entry in " << mesh().name() << " "
             << electrochemicalProperties_.name() << " dictionary!" << endl;
    }
}


void Foam::regionTypes::intercalationAnode::calculateButlerVolmer()
{
    // this version is from P2D Han et al.

    //dimensionedScalar spArea = epsS_/r_;

    eta_() = faiS_() - faiE_() - EEqRev_();

    // this version is from P2D Han et al.
    /*volScalarField iRef = kRct_*F*Foam::pow(cE_(), alphaA_)
                                 *Foam::pow(cSMax_() - cSE_(), alphaA_)
                                 *Foam::pow(cSE_(), alphaC_);

    volScalarField i = iRef*(Foam::exp(alphaA_*F/R/T_()*eta)
                           - Foam::exp(-alphaC_*F/R/T_()*eta));

    j_() = spArea*i;*/


    // this version is from LIONSIMBA - data source: Torchio M, Magni L, Gopaluni R B, et al. Lionsimba: a matlab framework based on a finite volume model suitable for li-ion battery design, simulation, and control[J]. Journal of The Electrochemical Society, 2016, 163(7): A1192.
    volScalarField keff = kRct_*Foam::exp(-EARct_/R*(1/T_()-1/TRef_));

    j_() = 2*keff*Foam::sqrt(cE_()*(cSMax_()- cSE_())*cSE_())*Foam::sinh(0.5*F/R/T_()*eta_());

    //Info << "max(j." << mesh().name() << " = " << max(j_()) << endl;
    //Info << "min(j." << mesh().name() << " = " << min(j_()) << endl;
}


void Foam::regionTypes::intercalationAnode::calculateTransportCoeffs()

{
	dimensionedScalar dimKappa =
        dimensionedScalar("dimKappa", dimensionSet(-1, -3, 3, 0, 0, 2, 0), 1);

    dimensionedScalar dimT =
        dimensionedScalar("dimT", dimensionSet(0, 0, 0, 1, 0, 0, 0), 1);

    dimensionedScalar dimD =
        dimensionedScalar("dimDE", dimensionSet(0, 2, -1, 0, 0, 0, 0), 1);


    volScalarField powKappa = -10.5
                    + 0.668e-3*cE_()/dimC_
                    + 0.494e-6*Foam::pow(cE_()/dimC_, 2)
                    + (0.074 - 1.78e-5*cE_()/dimC_
                    - 8.86e-10*Foam::pow(cE_()/dimC_, 2))*T_()/dimT
                    + (-6.96e-5 + 2.8e-8*cE_()/dimC_)*Foam::pow(T_()/dimT, 2);


    kappa_() = Foam::pow(epsE_, brugg_)*1.0e-4*cE_()/dimC_
              *Foam::pow(powKappa, 2)*dimKappa;

    Gamma_ = 2*(1-tNo_)*R/F;

    volScalarField powDE = -4.43 - (54/(T_()/dimT - 229 - 5e-3*cE_()/dimC_)) - 2.2e-4*cE_()/dimC_;
    DE_() = Foam::pow(epsE_, brugg_)*1.0e-4*Foam::pow(10, powDE)*dimD;

    DS_() = DSinit_*Foam::exp(-EADiff_/R*(1/T_()-1/TRef_));
}


void Foam::regionTypes::intercalationAnode::calculateInterfaceConcentration()
{
    //dimensionedScalar spArea = epsS_/r_;

    cSE_() = cS_() - j_()*r_/5/DS_();

    if(cSE_() <= 0*cSMin_())
    {
        cSE_() = cSMin_();
        Info << "Anode cSE too low! Value reset to 100" << endl;
    }
    if(cSE_() > cSMax_())
    {
        cSE_() = 0.99*cSMax_();
        Info << "Anode cSE too high! Value reset to 0.99*cSMax"
             << endl;
    }
}

void Foam::regionTypes::intercalationAnode::calculateHeatSourceTerms()
{
    volScalarField Qohm = sigma_()*(fvc::grad(faiS_())&fvc::grad(faiS_()))
        + kappa_()*(fvc::grad(faiE_())&fvc::grad(faiE_()))
        + kappa_()*Gamma_*T_()*(fvc::grad(Foam::log(cE_()/dimC_))&fvc::grad(faiE_()));

    volScalarField Qrct = F*spArea_*j_()*eta_();

    volScalarField Qrev = F*spArea_*j_()*T_()*dEEqdT_();

    ST_() = Qohm + Qrct + Qrev;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::intercalationAnode::intercalationAnode
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

    electrodeMaterial_(electrochemicalProperties_.lookup("electrodeMaterial")),
    epsE_(transportProperties_.lookup("epsE")),
    epsF_(transportProperties_.lookup("epsF")),
    epsS_(dimensionedScalar("epsS", dimless, 0)),
    r_(transportProperties_.lookup("r")),
    brugg_(transportProperties_.lookup("brugg")),
    spArea_(transportProperties_.lookup("spArea")),
    tNo_(electrochemicalProperties_.lookup("tNo")),
    Gamma_(dimensionedScalar("Gamma", dimensionSet(1, 2, -3, -1, 0, -1, 0), 0)),
    dimC_(dimensionedScalar("dimC", dimensionSet(0, -3, 0, 0, 1, 0, 0), 1)),
    rho_(transportProperties_.lookup("rho")),
    cp_(transportProperties_.lookup("cp")),
    kT_(transportProperties_.lookup("kT")),
    kappa_(nullptr),
    DE_(nullptr),
    DS_(nullptr),
    DSinit_(transportProperties_.lookup("DS")),
    EADiff_(transportProperties_.lookup("EADiff")),
    IApp_(electrochemicalProperties_.lookup("IApp")),
    alphaA_(electrochemicalProperties_.lookup("alphaA")),
    alphaC_(electrochemicalProperties_.lookup("alphaC")),
    EARct_(electrochemicalProperties_.lookup("EARct")),
    kRct_(electrochemicalProperties_.lookup("kRct")),
    TRef_(electrochemicalProperties_.lookup("TRef")),
    cSMax_(nullptr),
    cSMin_(nullptr),
    sigma_(nullptr),
    SoC_(nullptr),
    EEqRev_(nullptr),
    dEEqdT_(nullptr),
    eta_(nullptr),
    j_(nullptr),
    ST_(nullptr),
    faiS_(nullptr),
    faiE_(nullptr),
    cS_(nullptr),
    cE_(nullptr),
    cSE_(nullptr),
    T_(nullptr)
{


    // set electrolyte conductivity field
    kappa_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "kappa",
        dimensionedScalar("kappa0", dimensionSet(-1, -3, 3, 0, 0, 2, 0), 1),
        true
    );

    // set electrolyte diffusion coefficient field
    DE_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "DE",
        dimensionedScalar("DE0", dimensionSet(0, 2, -1, 0, 0, 0, 0), 1),
        true
    );

    // set solid diffusivity field
    DS_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "DS",
        dimensionedScalar(transportProperties_.lookup("DS")),
        true
    );

    // set maximum lithium solid surface concentration
    cSMax_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "cSMax",
        dimensionedScalar(transportProperties_.lookup("cSMax")),
        true
    );

    // set maximum lithium solid surface concentration
    cSMin_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "cSMin",
        dimensionedScalar(transportProperties_.lookup("cSMin")),
        true
    );

    // set effective electric diffusivity field
    sigma_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "sigma",
        dimensionedScalar(transportProperties_.lookup("sigmaInt"))*(1-epsE_-epsF_),
        true
    );

    // set state of charge field
    SoC_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "SoC",
        dimensionedScalar(electrochemicalProperties_.lookup("SoCInit")),
        true
    );

    // set reversible voltage field
    EEqRev_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "URev",
        dimensionedScalar("URev0", dimensionSet(1, 2, -3, 0, 0, -1, 0), 0),
        true
    );

    // set reversible voltage temperature dependency field
    dEEqdT_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "dEEqdT",
        dimensionedScalar("dEEqdT0", dimensionSet(1, 2, -3, -1, 0, -1, 0), 0),
        true
    );

    // set electrode overpotential field
    eta_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "eta",
        dimensionedScalar("eta0", dimensionSet(1, 2, -3, 0, 0, -1, 0), 0),
        true
    );

    // set volumetric exchange current density field
    j_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "j",
        dimensionedScalar("j_init", dimensionSet(0, -2, -1, 0, 1, 0, 0), 0),
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

    // set electric potential field
    faiS_ = lookupOrRead<volScalarField>(mesh(), "faiS");

    // set electrolyte potential field
    faiE_ = lookupOrRead<volScalarField>(mesh(), "faiE");

    // set solid surface lithium concentratio field
    cS_ = lookupOrRead<volScalarField>(mesh(), "cS");

    // set lithium electrolyte concentration field
    cE_ = lookupOrRead<volScalarField>(mesh(), "cE");

    // set lithium interface concentration field
    cSE_ = lookupOrRead<volScalarField>(mesh(), "cSE");

    // set temperature field
    T_ = lookupOrRead<volScalarField>(mesh(), "T");

    epsS_ = 1 - epsE_ - epsF_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::intercalationAnode::~intercalationAnode()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::intercalationAnode::correct()
{

}


Foam::scalar Foam::regionTypes::intercalationAnode::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::intercalationAnode::setCoupledEqns()
{
	calculateOCV();
    calculateButlerVolmer();
    calculateTransportCoeffs();
    calculateInterfaceConcentration();
    calculateHeatSourceTerms();

    faiSEqn =
    (
       - fvm::laplacian(sigma_(), faiS(), "laplacian(sigma,faiS)")
      ==
       - spArea_*F*j_()
    );

    faiEEqn =
    (
       - fvm::laplacian(kappa_(), faiE(), "laplacian(kappa,faiE)")
       + fvc::laplacian(kappa_()*Gamma_*T(), Foam::log(cE()/dimC_), "laplacian(DE,cE)")
       ==
         spArea_*F*j_()
    );

    cSEqn =
    (
        fvm::ddt(1, cS())
      ==
        - 3*j_()/r_
    );

    cEEqn =
    (
         epsE_*fvm::ddt(1, cE())
       - fvm::laplacian(DE_(), cE(), "laplacian(DE,cE)")
       ==
         spArea_*(1-tNo_)*j_()
    );

    TEqn =
    (
         fvm::ddt(rho_*cp_, T())
       - fvm::laplacian(kT_, T(), "laplacian(kT,T)")
       ==
         ST_()
    );

    fvScalarMatrices.set
    (
        faiS_().name()
      + mesh().name() + "Mesh"
      + intercalationAnode::typeName + "Type"
      + "Eqn",
        &faiSEqn()
    );

    fvScalarMatrices.set
    (
        faiE_().name()
      + mesh().name() + "Mesh"
      + intercalationAnode::typeName + "Type"
      + "Eqn",
        &faiEEqn()
    );

    fvScalarMatrices.set
    (
        cS_().name()
      + mesh().name() + "Mesh"
      + intercalationAnode::typeName + "Type"
      + "Eqn",
        &cSEqn()
    );

    fvScalarMatrices.set
    (
        cE_().name()
      + mesh().name() + "Mesh"
      + intercalationAnode::typeName + "Type"
      + "Eqn",
        &cEEqn()
    );

    fvScalarMatrices.set
    (
        T_().name()
      + mesh().name() + "Mesh"
      + intercalationAnode::typeName + "Type"
      + "Eqn",
        &TEqn()
    );

}

void Foam::regionTypes::intercalationAnode::postSolve()
{
    // do nothing, add as required
}

void Foam::regionTypes::intercalationAnode::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::intercalationAnode::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::intercalationAnode::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::intercalationAnode::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::intercalationAnode::meshMotionCorrector()
{
    // do nothing, add as required
}

// ************************************************************************* //
