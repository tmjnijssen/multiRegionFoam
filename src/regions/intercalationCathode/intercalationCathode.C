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
#include "intercalationCathode.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "thermoPhysicalConstants.H"

using namespace thermoPhysicalConstant;


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(intercalationCathode, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        intercalationCathode,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::regionTypes::intercalationCathode::calculateOCV()
{
    // TO-DO: Exclude empirical models for reversible electrode potential calculation to sub classes (need: SoC, IApp as input variables!)
    
    SoC_() = cSE_() / cSMax_();
    dimensionedScalar dimensionVolt = dimensionedScalar("dimensionVolt", dimensionSet(1, 2, -3, 0, 0, -1, 0), 1);
    dimensionedScalar dimensionKelvin = dimensionedScalar("dimensionKelvin", dimensionSet(0, 0, 0, 1, 0, 0, 0), 1);
    
    
    if(electrodeMaterial_ == "LFP" || "lfp")
    {
        //equilibrium potential of cathode material-LFP
        //data source: https://www.dandeliion.com/simulation
        EEqRev_() = (3.114559
                + 4.438792*Foam::atan(-71.7352*SoC_() + 70.85337) 
                - 4.240252*Foam::atan(-68.5605*SoC_() + 67.730082))
                *dimensionVolt;
    }
    if(electrodeMaterial_ == "NCA" || "nca")
    {
        //equilibrium potential of cathode material-NCA
        //data source: https://www.dandeliion.com/simulation
        EEqRev_() = (-1.30202*SoC_()
                - 0.214712*Foam::tanh(23.01*(SoC_() + 0.00350287)) 
                + 2.45808*Foam::tanh(2.90232*(SoC_() - 0.215657)) 
                - 1.26644*Foam::tanh(4.30574*(SoC_() - 0.329193)) 
                - 0.40112*Foam::tanh(9.07273*(SoC_() - 0.148644)) 
                - 0.0532656*Foam::exp(47.0417*(SoC_() - 0.95703)) 
                + 4.17032)*dimensionVolt;
    }
    if(electrodeMaterial_ == "LiIonSimBa" || "liionsimba" || "lionsimba" || "Lionsimba" || "LionSimba")
    {
        //equilibrium potential of cathode material
        //data source: Torchio M, Magni L, Gopaluni R B, et al. Lionsimba: a matlab framework based on a finite volume model suitable for li-ion battery design, simulation, and control[J]. Journal of The Electrochemical Society, 2016, 163(7): A1192.
        volScalarField EEqRevRef = (-4.656
                  + 88.669*Foam::pow(SoC_(),2)
                  - 401.119*Foam::pow(SoC_(),4)
                  + 342.909*Foam::pow(SoC_(),6)
                  - 462.471*Foam::pow(SoC_(),8)
                  + 433.434*Foam::pow(SoC_(),10))
                  /(-1 
                  + 18.933*Foam::pow(SoC_(),2)
                  - 79.532*Foam::pow(SoC_(),4)
                  + 37.311*Foam::pow(SoC_(),6)
                  - 73.083*Foam::pow(SoC_(),8) 
                  + 95.96*Foam::pow(SoC_(),10))
                  *dimensionVolt;

        dEEqdT_() = -0.001*(0.199521039
                  - 0.928373822*SoC_()
                  + 1.364550689000003*Foam::pow(SoC_(),2)
                  - 0.6115448939999998*Foam::pow(SoC_(),3))
                  /(1
                  - 5.661479886999997*SoC_()
                  + 11.47636191*Foam::pow(SoC_(),2)
                  - 9.82431213599998*Foam::pow(SoC_(),3)
                  + 3.048755063*Foam::pow(SoC_(),4))*dimensionVolt/dimensionKelvin;

        EEqRev_() = EEqRevRef + (T_() - TRef_)*dEEqdT_();
    }
    else
    {
        Info << "No valid material type input!" << endl;
        Info << "Check materialType entry in " << mesh().name() << " "
             << electrochemicalProperties_.name() << " dictionary!" << endl;
    }
}


void Foam::regionTypes::intercalationCathode::calculateButlerVolmer()
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


void Foam::regionTypes::intercalationCathode::calculateTransportCoeffs()

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


void Foam::regionTypes::intercalationCathode::calculateInterfaceConcentration()
{
    //dimensionedScalar spArea = epsS_/r_;
    
    cSE_() = cS_() - j_()*r_/5/DS_();
    
    if(cSE_() <= 0.427*cSMax_())
    {
        cSE_() = 0.427*cSMax_();
        Info << "Cathode cSE too low! Value reset to 0.427*cSMax" << endl;
    }
    if(cSE_() > cSMax_())
    {
        cSE_() = 0.99*cSMax_();
        Info << "Cathode cSE too high! Value reset to 0.99*cSMax"
             << endl;
    }
}

void Foam::regionTypes::intercalationCathode::calculateHeatSourceTerms()
{ 
    volScalarField Qohm = sigma_()*(fvc::grad(faiS_())&fvc::grad(faiS_()))
        + kappa_()*(fvc::grad(faiE_())&fvc::grad(faiE_()))
        + kappa_()*Gamma_*T_()*(fvc::grad(Foam::log(cE_()/dimC_))&fvc::grad(faiE_()));

    volScalarField Qrct = F*spArea_*j_()*eta_();

    volScalarField Qrev = F*spArea_*j_()*T_()*dEEqdT_();

    ST_() = Qohm + Qrct + Qrev;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::intercalationCathode::intercalationCathode
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
        dimensionedScalar(transportProperties_.lookup("sigmaInt"))*(1 - epsE_ - epsF_),
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

Foam::regionTypes::intercalationCathode::~intercalationCathode()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::intercalationCathode::correct()
{
    
}


Foam::scalar Foam::regionTypes::intercalationCathode::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::intercalationCathode::setCoupledEqns()
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
      + intercalationCathode::typeName + "Type"
      + "Eqn",
        &faiSEqn()
    );

    fvScalarMatrices.set
    (
        faiE_().name()
      + mesh().name() + "Mesh"
      + intercalationCathode::typeName + "Type"
      + "Eqn",
        &faiEEqn()
    );

    fvScalarMatrices.set
    (
        cS_().name()
      + mesh().name() + "Mesh"
      + intercalationCathode::typeName + "Type"
      + "Eqn",
        &cSEqn()
    );
    
    fvScalarMatrices.set
    (
        cE_().name()
      + mesh().name() + "Mesh"
      + intercalationCathode::typeName + "Type"
      + "Eqn",
        &cEEqn()
    );

    fvScalarMatrices.set
    (
        T_().name()
      + mesh().name() + "Mesh"
      + intercalationCathode::typeName + "Type"
      + "Eqn",
        &TEqn()
    );
    
}

void Foam::regionTypes::intercalationCathode::postSolve()
{
    // do nothing, add as required
}

void Foam::regionTypes::intercalationCathode::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::intercalationCathode::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::intercalationCathode::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::intercalationCathode::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::intercalationCathode::meshMotionCorrector()
{
    // do nothing, add as required
}

// ************************************************************************* //
