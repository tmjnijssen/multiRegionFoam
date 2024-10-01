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
#include "separator.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "thermoPhysicalConstants.H"

using namespace thermoPhysicalConstant;


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(separator, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        separator,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::regionTypes::separator::calculateTransportCoeffs()

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

}

void Foam::regionTypes::separator::calculateHeatSourceTerms()
{
    volScalarField Qohm = kappa_()*(fvc::grad(faiE_())&fvc::grad(faiE_()))
        + kappa_()*Gamma_*T_()*(fvc::grad(Foam::log(cE_()/dimC_))&fvc::grad(faiE_()));
 
    ST_() = Qohm;

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::separator::separator
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
    
    epsE_(transportProperties_.lookup("epsE")),
    brugg_(transportProperties_.lookup("brugg")),
    tNo_(electrochemicalProperties_.lookup("tNo")),
    Gamma_(dimensionedScalar("Gamma", dimensionSet(1, 2, -3, -1, 0, -1, 0), 0)),
    dimC_(dimensionedScalar("dimC", dimensionSet(0, -3, 0, 0, 1, 0, 0), 1)),
    rho_(transportProperties_.lookup("rho")),
    cp_(transportProperties_.lookup("cp")),
    kT_(transportProperties_.lookup("kT")),
    kappa_(nullptr),
    DE_(nullptr),
    ST_(nullptr),
    faiE_(nullptr),
    cE_(nullptr),
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
        dimensionedScalar(transportProperties_.lookup("DE")),
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

    // set electrolyte potential field
    faiE_ = lookupOrRead<volScalarField>(mesh(), "faiE");
    
    // set ion concentration field
    cE_ = lookupOrRead<volScalarField>(mesh(), "cE");
    
    // set temperature field
    T_ = lookupOrRead<volScalarField>(mesh(), "T");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::separator::~separator()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::separator::correct()
{
    
}


Foam::scalar Foam::regionTypes::separator::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::separator::setCoupledEqns()
{
	calculateTransportCoeffs();
    calculateHeatSourceTerms();
	
	faiEEqn =
    (
       - fvm::laplacian(kappa_(), faiE(), "laplacian(kappa,faiE)")
       - fvc::laplacian(kappa_()*Gamma_*T_(), Foam::log(cE()/dimC_), "laplacian(DE,cE)")
    );
    
    cEEqn =
    (
         epsE_*fvm::ddt(1, cE())
       - fvm::laplacian(DE_(), cE(), "laplacian(DE,cE)")
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
        faiE_().name()
      + mesh().name() + "Mesh"
      + separator::typeName + "Type"
      + "Eqn",
        &faiEEqn()
    );

    fvScalarMatrices.set
    (
        cE_().name()
      + mesh().name() + "Mesh"
      + separator::typeName + "Type"
      + "Eqn",
        &cEEqn()
    );

    fvScalarMatrices.set
    (
        T_().name()
      + mesh().name() + "Mesh"
      + separator::typeName + "Type"
      + "Eqn",
        &TEqn()
    );

}

void Foam::regionTypes::separator::postSolve()
{
    // do nothing, add as required
}

void Foam::regionTypes::separator::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::separator::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::separator::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::separator::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::separator::meshMotionCorrector()
{
    // do nothing, add as required
}

// ************************************************************************* //
