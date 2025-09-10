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
    b0_(sorbentProperties_.lookup("b0")),
    tau0_(sorbentProperties_.lookup("tau0")),
    qCO2inf0_(sorbentProperties_.lookup("qCO2inf0")),
    dH0CO2_(sorbentProperties_.lookup("dH0CO2")),
    T0CO2_(sorbentProperties_.lookup("T0CO2")),
    kCO2_(sorbentProperties_.lookup("kCO2")),
    CH2O_(sorbentProperties_.lookup("CH2O")),
    DH2O_(sorbentProperties_.lookup("DH2O")),
    FH2O_(sorbentProperties_.lookup("FH2O")),
    GH2O_(sorbentProperties_.lookup("GH2O")),
    qmH2O_(sorbentProperties_.lookup("qmH2O")),
    dH0H2O_(sorbentProperties_.lookup("dH0H2O")),
    kH2O_(sorbentProperties_.lookup("kH2O")),   
    
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
        dimensionedScalar("dqdtCO2", dimMoles/dimMass/dimTime, 0.0)
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
        dimensionedScalar("dqdtCO2ex", dimMoles/dimMass/dimTime, 0.0)
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
        dimensionedScalar("dqdtCO2im", dimVolume/dimMass/dimTime, 0.0)
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
        dimensionedScalar("dqdtH2O", dimMoles/dimMass/dimTime, 0.0)
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
        dimensionedScalar("dqdtH2Oex", dimMoles/dimMass/dimTime, 0.0)
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
        dimensionedScalar("dqdtH2Oim", dimVolume/dimMass/dimTime, 0.0)
    ), 
    qCO2_
    (
        IOobject
        (
            "qCO2",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("qCO2", dimMoles/dimMass, 0.0)
    ),    
    qH2O_
    (
        IOobject
        (
            "qH2O",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("qH2O", dimMoles/dimMass, 0.0)
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
{   
    
}

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
                (1-eps_)/eps_*rho_*dqdtCO2im_, 
                CO2_()
            )
          - ((1-eps_)/eps_)*rho_*dqdtCO2ex_ 
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
                ((1-eps_)/eps_)*rho_*dqdtH2Oim_,
                 H2O_()
            )
          - ((1-eps_)/eps_)*rho_*dqdtH2Oex_

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
}
void Foam::regionTypes::diffuseAdsorbSpecie::solveRegion()
{
    dimensionedScalar R = dimensionedScalar("R", dimEnergy/dimMoles/dimTemperature, 8.314);
    dimensionedScalar Z = dimensionedScalar("Z", dimPressure, 611.21);
    dimensionedScalar Z1 = dimensionedScalar("Z1", dimTemperature, 1);
    volScalarField invRT = (1./(R*T_()));
    // Calculate relative humidity
    volScalarField pH2O = R*T_()*H2O_(); // water vapour pressure from ideal gas law, Pa
    volScalarField pH2Osat = Z*exp((18.678-((T_()-(273.15*Z1)))/(234.5*Z1))*((T_()-273.15*Z1)/(T_()-16.01*Z1))); // Buck equation, Pa
    volScalarField RH = pH2O/pH2Osat;
    // Humidity-dependent CO2 adsorption parameters, Piscina et al. 2024 http://doi.org/10.2139/ssrn.5068012
    volScalarField tau = ((C_*RH*RH + D_*RH + F_)*RH + tau0_) + ((G_*RH*RH + H_*RH + J_)*RH + alpha_)*(1-(T0CO2_/T_()));
    volScalarField b = (A_*RH*exp(B_*RH) + b0_)*exp(-dH0CO2_*invRT);
    volScalarField qCO2inf = qCO2inf0_*exp(chi_*(1-(T_()/T0CO2_)));
    volScalarField qCO22 = (qCO2inf*b*R*T_()*CO2_())/pow(1+pow(b*R*T_()*CO2_(),tau),1/tau);                         // CO2 loading piscina
    // CO2 adsorption rate, Driessen et al. 2020 https://doi.org/10.1021/acs.iecr.9b05503
    dqdtCO2ex_ = -kCO2_*qCO2_/(b*qCO2inf);                             // explicit part
    dqdtCO2im_ = kCO2_*pow(1-pow((qCO2_/qCO2inf),tau),1/tau)*(R*T_()); // implicit part
    dqdtCO2_   = dqdtCO2ex_ + dqdtCO2im_*CO2_();                       // total CO2 adsorption rate
    // H2O adsorption parameters --> Low et al. 2025 https://doi.org/10.1021/acs.jced.3c00401
    dimensionedScalar Jm = dimensionedScalar("Jm", dimEnergy/dimMoles, 1);
    dimensionedScalar Jmk = dimensionedScalar("Jmk", dimEnergy/dimMoles/dimTemperature, 1);
    volScalarField E1 = CH2O_ - Jm*exp(DH2O_*T_());
    volScalarField E29 = FH2O_ + GH2O_*T_();
    volScalarField E10 = Jm*57220-44.38*T_()*Jmk;
    volScalarField c = exp((E1-E10)*invRT);
    volScalarField k = exp((E29-E10)*invRT);
    volScalarField qH2Oeq = (qmH2O_*k*c*RH)/((1-k*RH)*(1+(k*RH)*(c-1)));
    // H2O adsorption rate
    dqdtH2Oex_ = -kH2O_ * (qH2O_);               // explicit part
    dimensionedScalar mm3 = dimensionedScalar("mm3", dimMoles/dimVolume, 1);
    dqdtH2Oim_ = kH2O_*(qH2Oeq/(H2O_()+SMALL*mm3));          // implicit part
    dqdtH2O_   = dqdtH2Oex_ + dqdtH2Oim_*H2O_(); // total water adsorption rate
    // solve adsorbed species
    solve(fvm::ddt(qCO2_ ) ==  dqdtCO2_);
    solve(fvm::ddt(qH2O_ ) ==  dqdtH2O_);
    // heat source
    heatSource_() = -rho_*(dqdtCO2_ * dH0CO2_ + dqdtH2O_ * dH0H2O_);
    CO2_().correctBoundaryConditions();
    H2O_().correctBoundaryConditions();
}

void Foam::regionTypes::diffuseAdsorbSpecie::prePredictor()
{

}

void Foam::regionTypes::diffuseAdsorbSpecie::momentumPredictor()
{
        
}

void Foam::regionTypes::diffuseAdsorbSpecie::pressureCorrector()
{
        
}

void Foam::regionTypes::diffuseAdsorbSpecie::meshMotionCorrector()
{
        
}

// ************************************************************************* //
