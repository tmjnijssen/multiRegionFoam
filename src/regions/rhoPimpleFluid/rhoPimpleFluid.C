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

    laminarTransport_(nullptr),
    
	
    turbulence_(nullptr),

    U_(nullptr),
    phi_(nullptr),
    pKin_(nullptr),
    pAbs_(nullptr),
    p_(nullptr),
    rho_(nullptr),
    thermo.rho_(nullptr),
    psi_(nullptr),
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

    pAbs_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "pAbs",
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
    


    laminarTransport_.set(new singlePhaseTransportModel(U_(), phi_()));

    rho_.set(new dimensionedScalar(laminarTransport_().lookup("rho")));

    turbulence_ = compressible::turbulenceModel::New
    (
        U_(), phi_(), laminarTransport_()
    );

    rho_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "rho",
        true,
        true
    );
    
    p_ = lookupOrRead<volScalarField>
    (
        mesh(),
        "p",
        true,
        true
    );

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
    
    phid_ = lookupOrRead<surfaceScalarField>
    (
        mesh(),
        "phid",
        false,
        true,
        fvc::interpolate(psi)
       *(
            (fvc::interpolate(U_()) & mesh.Sf())
          + fvc::ddtPhiCorr(rUA, rho, U, phi)
        )
    );

    setRefCell(p_(), pimple_.dict(), pRefCell_, pRefValue_);
    mesh().schemesDict().setFluxRequired(p_().name());
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::rhoPimpleFluid::~rhoPimpleFluid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::pimpleFluid::correct()
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
     solve(fvm::ddt(rho) + fvc::div(phi));   // Solve the continuity for density.
}

void Foam::regionTypes::rhoPimpleFluid::prePredictor()
{
    Info<< nl << "Pre-predictor for " << this->typeName
        << " in region " << mesh().name()
        << nl << endl;

    if (mesh().changing() && correctPhi_)
    {
#       include "rhoPimpleFluidCorrectPhi.H"
    }
}

void Foam::regionTypes::rhoPimpleFluid::momentumPredictor()
{
    Info<< nl << "Momentum predictor for " << this->typeName
        << " in region " << mesh().name()
        << nl << endl;

    // Time derivative matrix
    tddtUEqn = fvm::ddt(U_()) + fvm::ddt(rho);
    fvVectorMatrix& ddtUEqn = tddtUEqn();

    // Convection-diffusion matrix
    tHUEqn =
        (
            fvm::ddt(rho_(), U_())
          + fvm::div(phi_(), U_())
          + turbulence_().divDevReff()
        );
    fvVectorMatrix& HUEqn = tHUEqn();

    mrfZones_.translationalMRFs().addFrameAcceleration(ddtUEqn);

    if (pimple_.momentumPredictor())
    {
        solve(relax(ddtUEqn + HUEqn) == -fvc::grad(p_()),
        mesh.solutionDict().solver((U.select(pimple.finalIter()))));
    }
}

void Foam::regionTypes::rhoPimpleFluid::pressureCorrector()
{
    Info<< nl << "Pressure corrector for " << this->typeName
        << " in region " << mesh().name()
        << nl << endl;

    // Get cached matricies from momentum predictor
    fvVectorMatrix& ddtUEqn = tddtUEqn();
    fvVectorMatrix& HUEqn = tHUEqn();
    

	
    // --- PISO loop
    while (pimple_.correct())
    {
        // Update pressure BCs
        p_().boundaryField().updateCoeffs();

        // read rho from thermo base
        
        rho_() = thermo.rho(); 
        
        // Prepare clean 1/a_p without time derivative and under-relaxation
        // contribution
        rAU_() = 1.0/HUEqn.A();

        // Calculate U from convection-diffusion matrix
        U_() = rAU_()*HUEqn.H();

        // Consistently calculate flux
        pimple_.calcTransientConsistentFlux(phi_(), U_(), rAU_(), ddtUEqn);

        // Global flux balance
        adjustPhi(phi_(), U_(), p_());

        while (pimple_.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian
                (
                    fvc::interpolate(rho_()*rAU_())/pimple_.aCoeff(U_().name()),
                    p_(),
                    "laplacian(rAU,p)"
                )
             ==
                fvc::div(phi_())
            );

            pEqn.setReference(pRefCell_, pRefValue_);
            pEqn.solve
            (
                mesh().solutionDict()
                .solver(p_().select(pimple_.finalInnerIter()))
            );

            if (pimple_.finalNonOrthogonalIter())
            {
                phi_() -= pEqn.flux();
            }
        }

        //- Pressure relaxation except for last corrector
        if (!pimple_.finalIter())
        {
            p_().relax();
        }

#       include "rhoPimpleFluidMovingMeshContinuityErrs.H"

        // Consistently reconstruct velocity after pressure equation. Note: flux is
        // made relative inside the function
        pimple_.reconstructTransientVelocity(U_(), phi_(), ddtUEqn, rAU_(), pKin_());

        // Update pressure field
        p_() = pRefValue_()+p_().value()+0.5*rho_()*U_()*U_();

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

void Foam::regionTypes::pimpleFluid::meshMotionCorrector()
{
    // Make the fluxes absolute
    fvc::makeAbsolute(phi_(), U_());

    mesh().update();

#       include "rhoPimpleFluidVolContinuity.H"

    if (mesh().changing() && correctPhi_)
    {
#       include "rhoPimpleFluidCorrectPhi.H"
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi_(), U_());

    if (mesh().moving() && checkMeshCourantNo_)
    {
#           include "rhoPimpleFluidMeshCourantNo.H"
    }

    if (mesh().changing())
    {
#           include "rhoPimpleFluidCourantNo.H"
    }
}

// ************************************************************************* //
