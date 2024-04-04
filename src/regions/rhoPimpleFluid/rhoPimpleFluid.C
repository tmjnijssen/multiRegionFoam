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
#include "rhoPimpleFluid.H"
#include "fvCFD.H"
#include "basicPsiThermo.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(rhoPimpleFluid, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        rhoPimpleFluid,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::rhoPimpleFluid::rhoPimpleFluid
(
    const Time& runTime,
    const word& regionName
)
:
    regionType(runTime, regionName),

    regionName_(regionName),

    pimple_(mesh()),

    pThermo_(basicPsiThermo::New(mesh())),

    p_(nullptr),
    h_(nullptr),
    psi_(nullptr),
    rho_(nullptr),
    U_(nullptr),
    phi_(nullptr),
    turbulence_(nullptr),
    T_(nullptr),
    sigma_(nullptr),
    kappaEff_(nullptr),

    DpDt_(nullptr),

    pMin_(pimple_.dict().lookup("pMin")),

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
    p_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "p",
        true,
        true
    );

    h_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "h",
        true,
        true
    );

    psi_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "psi",
        true,
        true
    );

    rho_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "rho",
        false,
        true,
        pThermo_().rho()
    );
    

    U_ = lookupOrRead<volVectorField>
    (
        mesh(),
        "U",
        true,
        true
    );

    phi_ = lookupOrRead<surfaceScalarField>
    (
        mesh(),
        "phi",
        false,
        true,
        linearInterpolate(rho_()*U_()) & mesh().Sf()
    );

    turbulence_ = compressible::turbulenceModel::New
    (
        rho_(),
        U_(),
        phi_(),
        pThermo_()
    );

    T_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "T",
        true,
        true
    );

    sigma_ = lookupOrRead<volSymmTensorField>
    (
        mesh(),
        "sigma",
        false,
        true,
        -p_()*symmTensor(1,0,0,1,0,1) - turbulence_().devRhoReff()
    );

    kappaEff_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "kappaEff",
        false,
        true,
        turbulence_().alphaEff() * pThermo_().Cp()
    );

    DpDt_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "DpDt",
        false,
        true,
        fvc::DDt
        (
            surfaceScalarField("phiU", phi_()/fvc::interpolate(rho_())),
            p_()
        )
    );

    // phid_ = lookupOrRead<surfaceScalarField>
    // (
    //     mesh(),
    //     "phid",
    //     false,
    //     true,
    //     fvc::interpolate(psi)
    //    *(
    //         (fvc::interpolate(U_()) & mesh.Sf())
    //       + fvc::ddtPhiCorr(rAU, rho, U, phi)
    //     )
    // );

    mesh().schemesDict().setFluxRequired(p_().name());
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::rhoPimpleFluid::~rhoPimpleFluid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::rhoPimpleFluid::correct()
{
#       include "rhoPimpleFluidCourantNo.H"
}


Foam::scalar Foam::regionTypes::rhoPimpleFluid::getMinDeltaT()
{
    //- TODO: implement deltaT based on CFL criteria
    return GREAT;
}


void Foam::regionTypes::rhoPimpleFluid::setCoupledEqns()
{

}


void Foam::regionTypes::rhoPimpleFluid::postSolve()
{
    if (myTimeIndex_ < mesh().time().timeIndex())
    {
        mrfZones_.translationalMRFs().correctMRF();

        mrfZones_.translationalMRFs().correctBoundaryVelocity(U_(), phi_());

        myTimeIndex_ = mesh().time().timeIndex();
    }
}


void Foam::regionTypes::rhoPimpleFluid::solveRegion()
{
    // Solve the continuity for density.
    solve(fvm::ddt(rho_()) + fvc::div(phi_()));
}

void Foam::regionTypes::rhoPimpleFluid::prePredictor()
{
    Info<< "Pre-predictor for " << this->typeName
        << " in region " << mesh().name()
        << endl;
        rho_().storePrevIter();
        p_().storePrevIter();
}

void Foam::regionTypes::rhoPimpleFluid::momentumPredictor()
{
    Info<< "Momentum predictor for " << this->typeName
        << " in region " << mesh().name()
        << endl;

    // Convection-diffusion matrix
    tUEqn =
        (
            fvm::ddt(rho_(), U_())
          + fvm::div(phi_(), U_())
          + turbulence_().divDevRhoReff()
        );
    fvVectorMatrix& UEqn = tUEqn();

    mrfZones_.translationalMRFs().addFrameAcceleration(UEqn, rho_());

    UEqn.relax
    (
        mesh().solutionDict().equationRelaxationFactor
        (
            U_().select(pimple_.finalIter())
        )
    );

    volScalarField rAU = 1.0/UEqn.A();

    if (pimple_.momentumPredictor())
    {
        solve
        (
            UEqn == -fvc::grad(p_()),
            mesh().solutionDict().solver((U_().select(pimple_.finalIter())))
        );
    }
    else
    {
        U_() = rAU*(UEqn.H() - fvc::grad(p_()));
        U_().correctBoundaryConditions();
    }

    fvScalarMatrix hEqn
    (
        fvm::ddt(rho_(), h_())
      + fvm::div(phi_(), h_())
      - fvm::laplacian(turbulence_().alphaEff(), h_())
     ==
        DpDt_()
    );

    hEqn.relax
    (
        mesh().solutionDict().equationRelaxationFactor
        (
            h_().select(pimple_.finalIter())
        )
    );
    hEqn.solve
    (
        mesh().solutionDict().solver((h_().select(pimple_.finalIter())))
    );

    pThermo_().correct();
}

void Foam::regionTypes::rhoPimpleFluid::pressureCorrector()
{
    Info<< "Pressure corrector for " << this->typeName
        << " in region " << mesh().name()
        << endl;

    // Get cached matricies from momentum predictor
    fvVectorMatrix& UEqn = tUEqn();

    // --- PISO loop
    while (pimple_.correct())
    {
        rho_() = pThermo_().rho();
        volScalarField rAU = 1.0/UEqn.A();

        // Calculate U from convection-diffusion matrix
        U_() = rAU*UEqn.H();

        if (pimple_.transonic())
        {
            surfaceScalarField phid
            (
                "phid",
                fvc::interpolate(psi_())
               *(
                    (fvc::interpolate(U_()) & mesh().Sf())
                  + fvc::ddtPhiCorr(rAU, rho_(), U_(), phi_())
                )
            );

            while (pimple_.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::ddt(psi_(), p_())
                  + fvm::div(phid, p_())
                  - fvm::laplacian(rho_()*rAU, p_())
                );

                pEqn.solve
                (
                    mesh().solutionDict().solver
                    (
                        p_().select(pimple_.finalInnerIter())
                    )
                );

                if (pimple_.finalNonOrthogonalIter())
                {
                    phi_() == pEqn.flux();
                }
            }
        }
        else
        {
            phi_() =
                fvc::interpolate(rho_())*
                (
                    (fvc::interpolate(U_()) & mesh().Sf())
                );

            while (pimple_.correctNonOrthogonal())
            {
                // Pressure corrector
                fvScalarMatrix pEqn
                (
                    fvm::ddt(psi_(), p_())
                  + fvc::div(phi_())
                  - fvm::laplacian(rho_()*rAU, p_())
                );

                pEqn.solve
                (
                    mesh().solutionDict().solver
                    (
                        p_().select(pimple_.finalInnerIter())
                    )
                );

                if (pimple_.finalNonOrthogonalIter())
                {
                    phi_() += pEqn.flux();
                }
            }
        }

        // Solve continuity for density
        solve(fvm::ddt(rho_()) + fvc::div(phi_()));
        #include "rhoPimpleFluidContinuityErrs.H"

        {
            // Explicitly relax pressure for momentum corrector
            p_().relax();

            rho_() = pThermo_().rho();
            rho_().relax();
            Info<< "  " << mesh().name() << ": rho min/mean/max: " 
                << gMin(rho_()) << "/"
                << gAverage(rho_()) << "/"
                << gMax(rho_()) << endl;
        }

        U_() -= rAU*fvc::grad(p_());
        U_().correctBoundaryConditions();

        DpDt_() = fvc::DDt
            (
                surfaceScalarField("phiU", phi_()/fvc::interpolate(rho_())),
                p_()
            );

        bound(p_(), pMin_);

        // Update sigma field
        sigma_() = -p_()*symmTensor(1,0,0,1,0,1) - turbulence_().devRhoReff();
    }

    turbulence_().correct();

    kappaEff_() = turbulence_().alphaEff() * pThermo_().Cp();

    Info<< "  " << mesh().name() << ": pressure min/mean/max: "
        << gMin(p_()) << "/"
        << gAverage(p_()) << "/"
        << gMax(p_()) << nl
        << "  " << mesh().name() << ": velocity min/mean/max: "
        << gMin(U_()) << "/"
        << gAverage(U_()) << "/"
        << gMax(U_()) << nl
        << "  " << mesh().name() << ": volume: "
        << gSum(mesh().V()) << endl;
}

void Foam::regionTypes::rhoPimpleFluid::meshMotionCorrector()
{
//     // Make the fluxes absolute
//     fvc::makeAbsolute(phi_(), U_());

//     mesh().update();

// #       include "rhoPimpleFluidVolContinuity.H"

//     if (mesh().changing() && correctPhi_)
//     {
// #       include "rhoPimpleFluidCorrectPhi.H"
//     }

//     // Make the fluxes relative to the mesh motion
//     fvc::makeRelative(phi_(), U_());

//     if (mesh().moving() && checkMeshCourantNo_)
//     {
// #           include "rhoPimpleFluidMeshCourantNo.H"
//     }

//     if (mesh().changing())
//     {
// #           include "rhoPimpleFluidCourantNo.H"
//     }
}

// ************************************************************************* //
