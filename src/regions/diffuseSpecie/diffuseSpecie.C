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
    R2NH_
    (
        IOobject
        (
            "R2NH",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("R2NH", dimMoles/dimVolume, 0.0)
    ),
    R2NH2p_
    (
        IOobject
        (
            "R2NH2p",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("R2NH2p", dimMoles/dimVolume, 0.0)
    ),
    R2NCO2m_
    (
        IOobject
        (
            "R2NHCO2m",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("R2NCO2m", dimMoles/dimVolume, 0.0)
    ),
    HCO3m_
    (
        IOobject
        (
            "HCO3m",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("HCO3m", dimMoles/dimVolume, 0.0)
    ),
    CO2_(nullptr),
    H2O_(nullptr)
{
    // set specie concentration fields and loadings in sorbent
    CO2_ = lookupOrRead<volScalarField>(mesh(), "CO2");
    H2O_ = lookupOrRead<volScalarField>(mesh(), "H2O");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::diffuseSpecie::~diffuseSpecie()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::diffuseSpecie::correct()
{
    // do nothing, add as required
}


Foam::scalar Foam::regionTypes::diffuseSpecie::getMinDeltaT()
{
    return GREAT;
}


void Foam::regionTypes::diffuseSpecie::setCoupledEqns(word fieldName)
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
            - ((1-eps_)/(eps_))*dqdtCO2_
        );

        fvScalarMatrices.set
        (
            CO2_().name()
          + mesh().name() + "Mesh"
          + diffuseSpecie::typeName + "Type"
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
            - ((1-eps_)/(eps_))*dqdtH2O_
        );

        fvScalarMatrices.set
        (
            H2O_().name()
          + mesh().name() + "Mesh"
          + diffuseSpecie::typeName + "Type"
          + "Eqn",
            &H2OEqn()
        );
    }
}

void Foam::regionTypes::diffuseSpecie::postSolve()
{
    // do nothing, add as required
}

void Foam::regionTypes::diffuseSpecie::solveRegion()
{
    // physisorbed species
    volScalarField HCO2 = HrCO2_*CO2_();
    volScalarField HH2O = HrH2O_*H2O_();

    // carbamate reaction rate
    volScalarField R1 = k1_ * (HCO2*R2NH_*R2NH_ - (1/K1_)*R2NH2p_*R2NCO2m_);

    // bicarbonate reaction rate
    volScalarField R2 = k2_ * (HCO2*HH2O*R2NH_ - (1/K2_)*R2NH2p_*HCO3m_);

    solve(fvm::ddt(R2NH_   ) == -2*R1 - R2);
    solve(fvm::ddt(R2NH2p_ ) ==    R1 + R2);
    solve(fvm::ddt(R2NCO2m_) ==    R1     );
    solve(fvm::ddt(HCO3m_  ) ==         R2);

    dqdtCO2_ = R1 + R2;
    dqdtH2O_ = R2;
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
