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
#include "diffuseSpecie.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(diffuseSpecie, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        diffuseSpecie,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::diffuseSpecie::diffuseSpecie
(
    const Time& runTime,
    const word& regionName
)
:
    regionType(runTime, regionName),

    regionName_(regionName),

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
    tau_
    (
        IOobject
        (
            "tau",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(sorbentProperties_.lookup("tau"))
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
        dimensionedScalar(sorbentProperties_.lookup("eps"))
    ),
    k1_
    (
        IOobject
        (
            "k1",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(sorbentProperties_.lookup("k1"))
    ),
    K1_
    (
        IOobject
        (
            "K1",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(sorbentProperties_.lookup("K1"))
    ),
    k2_
    (
        IOobject
        (
            "k2",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(sorbentProperties_.lookup("k2"))
    ),
    K2_
    (
        IOobject
        (
            "K2",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(sorbentProperties_.lookup("K2"))
    ),
    k3_
    (
        IOobject
        (
            "k2",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(sorbentProperties_.lookup("k3"))
    ),
    K3_
    (
        IOobject
        (
            "K3",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(sorbentProperties_.lookup("K3"))
    ),
    k4_
    (
        IOobject
        (
            "k4",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(sorbentProperties_.lookup("k4"))
    ),
    K4_
    (
        IOobject
        (
            "K4",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(sorbentProperties_.lookup("K4"))
    ),
    Dpore_
    (
        IOobject
        (
            "Dpore",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(sorbentProperties_.lookup("Dpore"))
    ),
    HrCO2_
    (
        IOobject
        (
            "HrCO2",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(sorbentProperties_.lookup("HrCO2"))
    ),
    HrH2O_
    (
        IOobject
        (
            "HrH2O",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(sorbentProperties_.lookup("HrH2O"))
    ),

    Dp_(nullptr),
    CO2_(nullptr),
    H2O_(nullptr),
    dqdtCO2_(nullptr),
    dqdtH2O_(nullptr),
    HCO3_(nullptr),
    R2NCO2_(nullptr),
    qHCO3_(nullptr),
    R2NH_(nullptr),
    R2NH2p_(nullptr)
{
    // set pore diffusivity field
    Dp_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "Dpore",
        dimensionedScalar(sorbentProperties_.lookup("Dpore")),
        true
    );

    

    // set specie concentration fields and loadings in sorbent
    CO2_ = lookupOrRead<volScalarField>(mesh(), "CO2");
    H2O_ = lookupOrRead<volScalarField>(mesh(), "H2O");
    R2NCO2_ = lookupOrRead<volScalarField>(mesh(), "R2NCO2");
    HCO3_ = lookupOrRead<volScalarField>(mesh(), "HCO3");
    R2NH_ = lookupOrRead<volScalarField>(mesh(), "R2NH");
    R2NH2p_ = lookupOrRead<volScalarField>(mesh(), "R2NH2p");
    dqdtH2O_ = lookupOrRead<volScalarField>(mesh(), "dqdtH2O");
    dqdtCO2_ = lookupOrRead<volScalarField>(mesh(), "dqdtCO2");
    qR2NCO2_ = lookupOrRead<volScalarField>(mesh(), "qR2NCO2");
    qHCO3_ = lookupOrRead<volScalarField>(mesh(), "qHCO3");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::diffuseSpecie::~diffuseSpecie()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::diffuseSpecie::correct()
{
    Dp_().correctBoundaryConditions();
}


Foam::scalar Foam::regionTypes::diffuseSpecie::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::diffuseSpecie::setCoupledEqns()
{
    dqdtCO2Eqn = (fvm::ddt(R2NCO2_()) + fvm::ddt(HCO3_()));
    dqdtH2OEqn = fvm::ddt(HCO3_());

    CO2Eqn =
    (
        fvm::ddt(CO2_())
     ==
        fvm::laplacian(Dp_(), CO2_(), "laplacian(Dp,CO2)") - ((1-eps_)/(eps_*eps_))*dqdtCO2Eqn
    );

    H2OEqn =
    (
        fvm::ddt(H2O_())
     ==
        fvm::laplacian(Dp_(), H2O_(), "laplacian(Dp,H2O)") - ((1-eps_)/(eps_*eps_))*dqdtH2OEqn
    );

    fvScalarMatrices.set
    (
        (CO2_().name(), H2O_().name())
      + mesh().name() + "Mesh"
      + diffuseSpecie::typeName + "Type"
      + "Eqn",
        (&CO2Eqn(), &H2OEqn())
    );
}

void Foam::regionTypes::diffuseSpecie::postSolve()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseSpecie::solveRegion()
{
    HCO2_() = HrCO2_*CO2_();
    HH2O_() = HrH2O_*H2O_();

    R2NCO2Eqn =
    (
        fvm::ddt(R2NCO2_())
     ==
        fvm::laplacian(Dp_(), HCO2_(), "laplacian(Dp,HCO2)") - ((1-eps_)/(eps_*eps_))*(k1_*(HCO2_()*R2NCO2_()*R2NCO2_())-(1/K1_)*R2NH2p_()*R2NCO2_())
    );

    
    HCO3Eqn =
    (
        fvm::ddt(HCO3_())
     ==
        fvm::laplacian(Dp_(), HH2O_(), "laplacian(Dp,HH2O)") - ((1-eps_)/(eps_*eps_))*(k2_*(HCO2_()*HH2O_()*R2NH_())-(1/K2_)*R2NH2p_()*HCO3_())
    );

    R2NH2pEqn =
    (
        fvm::ddt(R2NH2p_())
    ==
        fvm::ddt(R2NCO2_()) + fvm::ddt(HCO3_())
    );

    R2NHEqn =
    (
        fvm::ddt(R2NH_())
    ==
        (-2*fvm::ddt(R2NCO2_()) - fvm::ddt(HCO3_()))
    );

    fvScalarMatrices.set
    (
        (R2NCO2_().name(), HCO3_().name(), HCO2_().name(), HH2O_().name(), R2NH_().name(), R2NH2p_().name())
      + mesh().name() + "Mesh"
      + diffuseSpecie::typeName + "Type"
      + "Eqn",
        (&R2NCO2Eqn(), &HCO3Eqn(), &R2NH2pEqn(), &R2NHEqn())
    );
    
}

void Foam::regionTypes::diffuseSpecie::prePredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseSpecie::momentumPredictor()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseSpecie::pressureCorrector()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseSpecie::meshMotionCorrector()
{
    // do nothing, add as required
}

// ************************************************************************* //
