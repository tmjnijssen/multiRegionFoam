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
#include "diffuseAdsorbSpecie.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(diffuseAdsorbSpecie, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        diffuseAdsorbSpecie,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::diffuseAdsorbSpecie::diffuseAdsorbSpecie
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
        dimensionedScalar(transportProperties_.lookup("eps"))
    ),
    kCO2_(sorbentProperties_.lookup("kCO2")),
    kH2O_(sorbentProperties_.lookup("kH2O")),
    b0_(sorbentProperties_.lookup("b0")),
    tau_(sorbentProperties_.lookup("tau")),
    tau0_(sorbentProperties_.lookup("tau0")),
    qCO2inf_(sorbentProperties_.lookup("qCO2inf")),
    qCO2inf0_(sorbentProperties_.lookup("qCO2inf0")),
    qH2Oinf_(sorbentProperties_.lookup("qH2Oinf")),
    T0CO2_(sorbentProperties_.lookup("T0CO2")),
    T0H2O_(sorbentProperties_.lookup("T0H2O")),
    dH1_(sorbentProperties_.lookup("dH1")),
    dH2_(sorbentProperties_.lookup("dH2")),
    A_(sorbentProperties_.lookup("A")),
    B_(sorbentProperties_.lookup("B")),
    C_(sorbentProperties_.lookup("C")),
    D_(sorbentProperties_.lookup("D")),
    F_(sorbentProperties_.lookup("F")),
    G_(sorbentProperties_.lookup("G")),
    H_(sorbentProperties_.lookup("H")),
    J_(sorbentProperties_.lookup("J")),
    chi_(sorbentProperties_.lookup("chi")),
    alpha_(sorbentProperties_.lookup("alpha")),
    Ea1_(sorbentProperties_.lookup("Ea1")),
    
    DporeCO2_
    (
        IOobject
        (
            "DporeCO2",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(transportProperties_.lookup("DporeCO2"))
    ),
    DporeH2O_
    (
        IOobject
        (
            "DporeH2O",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(transportProperties_.lookup("DporeH2O"))
    ),

    dqdtCO2_
    (
        IOobject
        (
            "dqdtCO2",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("dqdtCO2", dimMoles/dimVolume/dimTime, 0.0)
    ),
    dqdtCO2ex_
    (
        IOobject
        (
            "dqdtCO2ex",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("dqdtCO2ex", dimMoles/dimVolume/dimTime, 0.0)
    ),
    dqdtCO2im_
    (
        IOobject
        (
            "dqdtCO2im",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("dqdtCO2im", dimless/dimTime, 0.0)
    ),
    dqdtH2O_
    (
        IOobject
        (
            "dqdtH2O",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("dqdtH2O", dimMoles/dimVolume/dimTime, 0.0)
    ),
    dqdtH2Oex_
    (
        IOobject  
        (
            "dqdtH2Oex",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("dqdtH2Oex", dimMoles/dimVolume/dimTime, 0.0)
    ),
    dqdtH2Oim_
    (
        IOobject  
        (
            "dqdtH2Oim",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("dqdtH2Oim", dimless/dimTime, 0.0)
    ), 
    qCO2ads_
    (
        IOobject
        (
            "qCO2ads",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("qCO2ads", dimMoles/dimMass, 0.0)
    ),    
    qH2Oads_
    (
        IOobject
        (
            "qH2Oads",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("qH2Oads", dimMoles/dimMass, 0.0)
    ),    
    
    T_(nullptr),
    heatSource_(nullptr),
    CO2_(nullptr),
    H2O_(nullptr)
{
    // temperature field
    T_ = lookupOrRead<volScalarField>(mesh(), "T");

    // heat source field
    heatSource_ = lookupOrRead<volScalarField>
    (
        mesh(), transportProperties_.lookupOrDefault<word>("heatSourceName", "heatSource"),
        dimensionedScalar("heatSource", dimEnergy/dimTime/dimVolume, 0.0),
        true
    );

    // set specie concentration fields and loadings in sorbent
    CO2_ = lookupOrRead<volScalarField>(mesh(), "CO2");
    H2O_ = lookupOrRead<volScalarField>(mesh(), "H2O");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::diffuseAdsorbSpecie::~diffuseAdsorbSpecie()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::diffuseAdsorbSpecie::correct()
{
   CO2_().correctBoundaryConditions();
   H2O_().correctBoundaryConditions();
}


Foam::scalar Foam::regionTypes::diffuseAdsorbSpecie::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::diffuseAdsorbSpecie::setCoupledEqns(word fieldName)
{
    if
    (
            fieldName == CO2_().name()
        ||  fieldName == word::null
    )
    {
        CO2Eqn =
        (
            eps_*fvm::ddt(CO2_())
         ==
            fvm::laplacian(DporeCO2_, CO2_())
          - fvm::SuSp
            (
                (1-eps_)/eps_* dqdtCO2im_, 
                CO2_()
            )
          - ((1-eps_)/eps_)*dqdtCO2ex_ 
        );

        fvScalarMatrices.set
        (
            CO2_().name()
          + mesh().name() + "Mesh"
          + diffuseAdsorbSpecie::typeName + "Type"
          + "Eqn",
            &CO2Eqn()
        );
    }

    if
    (
            fieldName == H2O_().name()
        ||  fieldName == word::null
    )
    {
        H2OEqn =
        (
            eps_*fvm::ddt(H2O_())
         ==
            fvm::laplacian(DporeH2O_, H2O_())
          - fvm::SuSp
            (
                ((1-eps_)/eps_)*dqdtH2Oim_,
                 H2O_()
            )
          - ((1-eps_)/eps_) * dqdtH2Oex_

        );

        fvScalarMatrices.set
        (
            H2O_().name()
          + mesh().name() + "Mesh"
          + diffuseAdsorbSpecie::typeName + "Type"
          + "Eqn",
            &H2OEqn()
        );
    }
}

void Foam::regionTypes::diffuseAdsorbSpecie::postSolve()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseAdsorbSpecie::solveRegion()
{
    dimensionedScalar R = dimensionedScalar("R", dimEnergy/dimMoles/dimTemperature, 8.314);
    volScalarField invRT = (1./(R*T_()));

    // reaction rates
    //volScalarField k2 = A2_ * exp(-Ea2_ * invRT);
    //volScalarField dG2 = dH2_ - T_() * dS2_;

    // Equilibrium constants
    // volScalarField K1 = exp(-dG1 * invRT);
    // volScalarField K2 = exp(-dG2 * invRT);

    // Calculate Relative Humidity and water equilibrium loading
    volScalarField Ctot = 101325*invRT;
    volScalarField pH2O = H2O_()/(Ctot*101325);
    volScalarField pH2Osat = 0.61078*exp((17.27*T_()/(T_()+237.3)));
    volScalarField RH = pH2O/pH2Osat;

    // CO2 rate implementation  --> Piscina 2024 
    volScalarField tauT = ((C_*RH*RH + D_*RH+F_)*qCO2inf_ + tau0_) + ((G_*RH*RH + H_*RH + J_ + alpha_))*(1-(T0CO2_/T_()));   // 
    volScalarField bT = (A_*RH*exp(B_*RH)+ b0_)*exp(-dH1_*invRT);
    volScalarField qCO2infT = qCO2inf0_*exp(chi_*(1-(T_()/T0CO2_)));
    


    // H2O rate implementation --> Low et al. 2025 https://doi.org/10.1016/j.cherd.2025.01.048 and young 2021 
    volScalarField Ea2 = cH2O_ - exp(dH2O_*T_());
    volScalarField Ea29 = fH2O_ + gH2O_*T_();
    volScalarField c = exp(Ea2 -(57220-44.38*T_())*invRT);
    volScalarField k = exp(Ea29);
    volScalarField qH2Oads = (qmH2O_*k*c*RH)/((1-k*RH)*(1+(k*RH)*(c-1)));

    // CO2 reaction rate
    volScalarField R1ex = -kCO2_*(qCO2ads_/bT_*qCO2infT); // explicit part
    volScalarField R1im = kCO2_*pow(pow(1-(qCO2ads_/qCO2infT),tauT_),1/tauT_)*(R*T_());            // implicit part
    volScalarField R1   = R1ex + R1im*CO2_();                // total carbamate rate

    // physical adsorption water
    volScalarField R2ex = -kH2O_ * (qH2Oads_);                // explicit part
    volScalarField R2im = kH2O_*(qH2Oinf_/H2O_());                                // implicit part
    volScalarField R2   = R2ex + R2im*H2O_();           // total water adsorption rate

    // solve adsorbed species
    solve(fvm::ddt(qCO2ads_ )  ==  R1 );
    solve(fvm::ddt(qH2Oads_ ) ==  R2 );

    // total adsorption rate
    dqdtCO2im_ = R1im;
    dqdtCO2ex_ = R1ex;;
    dqdtCO2_   = R1;

    dqdtH2Oim_ = R2im;
    dqdtH2Oex_ = R2ex;
    dqdtH2O_   = R2;

    // heat source
    heatSource_() = -(R1 * dH1_ + R2 * dH2_);

    CO2_().correctBoundaryConditions();
    H2O_().correctBoundaryConditions();
}

void Foam::regionTypes::diffuseAdsorbSpecie::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseAdsorbSpecie::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseAdsorbSpecie::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseAdsorbSpecie::meshMotionCorrector()
{
    // do nothing, add as required
}

// ************************************************************************* //
