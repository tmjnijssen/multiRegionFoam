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
#include "reactingFluid.H"
#include "fvCFD.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(reactingFluid, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        reactingFluid,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::reactingFluid::reactingFluid
(
    const Time& runTime,
    const word& regionName
)
:
    regionType(runTime, regionName),

    regionName_(regionName),

    pimple_(mesh()),

    pChemistry_(psiChemistryModel::New(mesh())),

    inertSpecie_(pChemistry_().thermo().lookup("inertSpecie")),
    rho_(nullptr),
    U_(nullptr),
    p_(nullptr),
    psi_(nullptr),
    hs_(nullptr),
    T_(nullptr),
    phi_(nullptr),
    chemistryKappa_
    (
        IOobject
        (
            "chemistryKappa",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    turbulence_(nullptr),
    sigma_(nullptr),
    DpDt_(nullptr),
    fields_(),
    turbulentReaction_(false),
    Cmix_(0),
    chemistrySh_
    (
        IOobject
        (
            "chemistry::Sh",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("chemistrySh", dimEnergy/dimTime/dimVolume, 0.0)
    ),
    g_
    (
        IOobject
        (
            "g",
            mesh().time().constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
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
#   include "readChemistryProperties.H"
#   include "createFields.H"
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::reactingFluid::~reactingFluid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::reactingFluid::correct()
{
#   include "reactingFluidCourantNo.H"
}


Foam::scalar Foam::regionTypes::reactingFluid::getMinDeltaT()
{
    //- TODO: implement deltaT based on CFL criteria
    return GREAT;
}


void Foam::regionTypes::reactingFluid::setCoupledEqns()
{

}


void Foam::regionTypes::reactingFluid::postSolve()
{
    /*if (myTimeIndex_ < mesh().time().timeIndex())
    {
        mrfZones_.translationalMRFs().correctMRF();

        mrfZones_.translationalMRFs().correctBoundaryVelocity(U_(), phi_());

        myTimeIndex_ = mesh().time().timeIndex();
    }*/
    
    if (mesh().time().write())
    {
        pChemistry_().dQ()().write();
    }
}


void Foam::regionTypes::reactingFluid::solveRegion()
{
#   include "chemistry.H"
#   include "rhoEqn.H"
}

void Foam::regionTypes::reactingFluid::prePredictor()
{
    Info<< "Pre-predictor for " << this->typeName
        << " in region " << mesh().name()
        << endl;

    loopPIMPLE_ = pimple_.loop();

    if (!loopPIMPLE_)
    {
        Info<< "Breaking PIMPLE loop for " << this->typeName
        << " in region " << mesh().name()
        << endl;
        return;
    }
}

void Foam::regionTypes::reactingFluid::momentumPredictor()
{
    Info<< "Momentum predictor for " << this->typeName
        << " in region " << mesh().name()
        << endl;

#   include "UEqn.H"
#   include "YEqn.H"
#   include "hsEqn.H"
}

void Foam::regionTypes::reactingFluid::pressureCorrector()
{
    Info<< "Pressure corrector for " << this->typeName
        << " in region " << mesh().name()
        << endl;

    // Get cached matricies from momentum predictor
    fvVectorMatrix& UEqn = tUEqn();

    // --- PISO loop
    while (pimple_.correct())
    {
#       include "pEqn.H"

        // Update sigma field
        sigma_() = -p_()*symmTensor(1,0,0,1,0,1) - turbulence_().devRhoReff();
    }

    turbulence_().correct();

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

void Foam::regionTypes::reactingFluid::meshMotionCorrector()
{
//     // Make the fluxes absolute
//     fvc::makeAbsolute(phi_(), U_());

//     mesh().update();

// #       include "reactingFluidVolContinuity.H"

//     if (mesh().changing() && correctPhi_)
//     {
// #       include "reactingFluidCorrectPhi.H"
//     }

//     // Make the fluxes relative to the mesh motion
//     fvc::makeRelative(phi_(), U_());

//     if (mesh().moving() && checkMeshCourantNo_)
//     {
// #           include "reactingFluidMeshCourantNo.H"
//     }

//     if (mesh().changing())
//     {
// #           include "reactingFluidCourantNo.H"
//     }
}

// ************************************************************************* //
