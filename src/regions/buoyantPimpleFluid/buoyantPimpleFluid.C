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

#include "label.H"
#include "bouyantPimpleFluid.H"

#include "fvCFD.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(bouyantPimpleFluid, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        bouyantPimpleFluid,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::bouyantPimpleFluid::bouyantPimpleFluid
(
    const Time& runTime,
    const word& regionName
)
:
    regionType(runTime, regionName),

    regionName_(regionName),

    pimple_(mesh()),

    laminarTransport_(nullptr),
    g_(
        IOobject
        (
            "g",
            runTime.constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    rho_(nullptr),
    beta_(nullptr),
    TRef_(nullptr),

    turbulence_(nullptr),

    U_(nullptr),
    phi_(nullptr),
    pKin_(nullptr),
    p_(nullptr),
    T_(nullptr),

    sigma_(nullptr),

    rAU_(nullptr),

    pRefCell_
    (
        pimple_.dict().lookupOrDefault<label>("pRefCell", 0)
    ),
    pRefValue_
    (
        pimple_.dict().lookupOrDefault<scalar>("pRefValue", 0.0)
    ),

    mrfZones_(mesh()),
    myTimeIndex_(mesh().time().timeIndex()),

    maxCo_
    (
        mesh().time().controlDict().lookupOrDefault<scalar>("maxCo", 1.0)
    ),
    maxDeltaT_
    (
        mesh().time().controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT)
    ),

    correctPhi_
    (
        pimple_.dict().lookupOrDefault<Switch>("correctPhi", false)
    ),
    checkMeshCourantNo_
    (
        pimple_.dict().lookupOrDefault<Switch>("checkMeshCourantNo", false)
    ),

    sumLocalContErr_(0),
    globalContErr_(0),
    cumulativeContErr_(0)
{
    U_ = lookupOrRead<volVectorField>
    (
        mesh(),
        "U",
        true,
        true
    );

    pKin_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "pKin",
        true,
        true
    );

    phi_ = lookupOrRead<surfaceScalarField>
    (
        mesh(),
        "phi",
        false,
        true,
        linearInterpolate(U_()) & mesh().Sf()
    );

    laminarTransport_.set(new singlePhaseTransportModel(U_(), phi_()));

    rho_.set(new dimensionedScalar(laminarTransport_().lookup("rho")));
    beta_.set(new dimensionedScalar(laminarTransport_().lookup("beta")));
    TRef_.set(new dimensionedScalar(laminarTransport_().lookup("TRef")));

    turbulence_ = incompressible::turbulenceModel::New
    (
        U_(), phi_(), laminarTransport_()
    );

    p_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "p",
        false,
        true,
        rho_().value()*pKin_()
    );

    T_ = lookupOrRead<volScalarField>(mesh(), "T");

    sigma_ = lookupOrRead<volSymmTensorField>
    (
        mesh(),
        "sigma",
        false,
        true,
        rho_().value()
       *(
            - pKin_()*symmTensor(1,0,0,1,0,1)
            - turbulence_().devReff()
        )
    );

    wordList rAUPatchFieldTypes
    (
        U_().boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );
    rAU_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "rAU",
        mesh().time().deltaT(),
        rAUPatchFieldTypes,
        true
    );

    setRefCell(pKin_(), pimple_.dict(), pRefCell_, pRefValue_);
    mesh().schemesDict().setFluxRequired(pKin_().name());
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::bouyantPimpleFluid::~bouyantPimpleFluid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::bouyantPimpleFluid::correct()
{
#       include "pimpleFluidCourantNo.H"
}


Foam::scalar Foam::regionTypes::bouyantPimpleFluid::getMinDeltaT()
{
    //- TODO: implement deltaT based on CFL criteria
    return GREAT;
}


void Foam::regionTypes::bouyantPimpleFluid::setCoupledEqns()
{
    // do nothing, add as required
}


void Foam::regionTypes::bouyantPimpleFluid::postSolve()
{
    if (myTimeIndex_ < mesh().time().timeIndex())
    {
        mrfZones_.translationalMRFs().correctMRF();

        mrfZones_.translationalMRFs().correctBoundaryVelocity(U_(), phi_());

        myTimeIndex_ = mesh().time().timeIndex();
    }
}


void Foam::regionTypes::bouyantPimpleFluid::solveRegion()
{
    // do nothing, add as required
}

void Foam::regionTypes::bouyantPimpleFluid::prePredictor()
{
    Info<< nl << "Pre-predictor for " << this->typeName
        << " in region " << mesh().name()
        << nl << endl;

    if (mesh().changing() && correctPhi_)
    {
#       include "pimpleFluidCorrectPhi.H"
    }
}

void Foam::regionTypes::bouyantPimpleFluid::momentumPredictor()
{
    Info<< nl << "Momentum predictor for " << this->typeName
        << " in region " << mesh().name()
        << nl << endl;

    tUEqn =
        (
            fvm::ddt(U_())
          + fvm::div(phi_(), U_())
          + turbulence_->divDevReff()
        );
    fvVectorMatrix& UEqn = tUEqn();

    mrfZones_.translationalMRFs().addFrameAcceleration(UEqn);

    UEqn.relax();

    if (pimple_.momentumPredictor())
    {
        solve
        (
            UEqn
          ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(1.0 - beta_()*(T_() - TRef_()))*(g_ & mesh().Sf())
                  - fvc::snGrad(p_())*mesh().magSf()
                )
            )
        );
    }
}

void Foam::regionTypes::bouyantPimpleFluid::pressureCorrector()
{
    Info<< nl << "Pressure corrector for " << this->typeName
        << " in region " << mesh().name()
        << nl << endl;

    // Get cached matricies from momentum predictor
    fvVectorMatrix& UEqn = tUEqn();

    // --- PISO loop
    while (pimple_.correct())
    {
        // Update pressure BCs
        pKin_().boundaryField().updateCoeffs();

        // Prepare clean 1/a_p without time derivative and under-relaxation
        // contribution
        rAU_() = 1.0/UEqn.A();
        surfaceScalarField rAUf("(1|A(U))", fvc::interpolate(rAU_()));

        // Calculate U from convection-diffusion matrix
        U_() = rAU_()*UEqn.H();

        // Calculate flux
        surfaceScalarField phiU
        (
            (fvc::interpolate(U_()) & mesh().Sf())
          + fvc::ddtPhiCorr(rAU_(), U_(), phi_())
        );

        phi_() = phiU
                + rAUf*fvc::interpolate(1.0 - beta_()*(T_() - TRef_()))
                  *(g_ & mesh().Sf());

        // Global flux balance
        adjustPhi(phi_(), U_(), pKin_());

        while (pimple_.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(rAUf, pKin_()) == fvc::div(phi_())
            );

            pEqn.setReference(pRefCell_, pRefValue_);
            pEqn.solve
            (
                mesh().solutionDict()
                .solver(pKin_().select(pimple_.finalInnerIter()))
            );

            if (pimple_.finalNonOrthogonalIter())
            {
                phi_() -= pEqn.flux();
            }
        }

        //- Pressure relaxation except for last corrector
        if (!pimple_.finalIter())
        {
            pKin_().relax();
        }

#       include "pimpleFluidMovingMeshContinuityErrs.H"

        // Correct the momentum source with the pressure gradient flux
        // calculated from the relaxed pressure
        U_() += rAU_()*fvc::reconstruct((phi_() - phiU)/rAUf);
        U_().correctBoundaryConditions();

        // Update pressure field
        p_() = rho_().value()*pKin_();

        // Update sigma field
        sigma_() = rho_().value()
           *(
                - pKin_()*symmTensor(1,0,0,1,0,1)
                - turbulence_().devReff()
            );
    }

    turbulence_().correct();

    Info<< nl
        << mesh().name() << " Pressure:" << nl
        << "  max: " << gMax(p_()) << nl
        << "  min: " << gMin(p_()) << nl
        << "  mean: " << gAverage(p_()) << nl
        << mesh().name() << " Velocity:" << nl
        << "  max: " << gMax(U_()) << nl
        << "  min: "<< gMax(U_()) << nl
        << "  mean: " << gAverage(U_()) << nl
        << mesh().name() << " Volume: "
        << gSum(mesh().V()) << nl
        << endl;
}

void Foam::regionTypes::bouyantPimpleFluid::meshMotionCorrector()
{
    // Make the fluxes absolute
    fvc::makeAbsolute(phi_(), U_());

    mesh().update();

#       include "pimpleFluidVolContinuity.H"

    if (mesh().changing() && correctPhi_)
    {
#       include "pimpleFluidCorrectPhi.H"
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi_(), U_());

    if (mesh().moving() && checkMeshCourantNo_)
    {
#           include "pimpleFluidMeshCourantNo.H"
    }

    if (mesh().changing())
    {
#           include "pimpleFluidCourantNo.H"
    }
}

// ************************************************************************* //
