/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Fork:     foam-extend
    \\  /    A nd           | Version:  4.1
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::rhoPimpleFluid.C

SourceFiles
    rhoPimpleFluid.C

Author
    Constantin Habes (constantin.habes@tu-darmstadt.de, Affiliation A)
    Holger Marschall (holger.marschall@tu-darmstadt.de, Affiliation A)

Contact
    Holger Marschall (holger.marschall@tu-darmstadt.de)
    main developer and principal investigator
    TU Darmstadt

Affiliations
    Affiliation A)
    Computational Multiphase Flow
    Department of Mathematics
    Technical University of Darmstadt, Germany

Acknowledgement
    Funded by
    Hessian Ministry of Higher Education, Research, Science and the Arts,
    National High Performance Computing Center for Computational Engineering
    Science (NHR4CES)

Description

    This file is part of the multiRegionFoam library.

    You may refer to this software as :
    https://doi.org/10.48550/arXiv.2306.01924

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/

#ifndef rhoPimpleFluid_H
#define rhoPimpleFluid_H

#include "regionType.H"
#include "MovingReferenceFrameZones.H"
#include "dimensionedScalar.H"
#include "multiRegionSystem.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOReferencer.H"
#include "HashPtrTable.H"
#include "rhoPimpleControl.H"
#include "solutionControl.H"
#include "volFieldsFwd.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{

/*---------------------------------------------------------------------------*\
                           Class rhoPimpleFluid Declaration
\*---------------------------------------------------------------------------*/

class rhoPimpleFluid
:
    public regionType
{

        //- name
        word regionName_;

        //- PIMPLE control
        pimpleControl pimple_;

        //- Transport model
        autoPtr<singlePhaseTransportModel> laminarTransport_;

        //- density
        autoPtr<dimensionedScalar> rho_;

        //- Turbulence model
        autoPtr<compressible::turbulenceModel> turbulence_;

        //- Velocity
        autoPtr<volVectorField> U_;

        //- Face flux
        autoPtr<surfaceScalarField> phi_;
        
        //- compressible term
        autoPtr<surfaceScalarField> psi_;

        //- Kinematic pressure
        autoPtr<volScalarField> pKin_;
        
        //- Absolute
        autoPtr<volScalarField> pAbs_;

        //- Pressure
        autoPtr<volScalarField> p_;
        
        //- Density
        autoPtr<volScalarField> rho_;
        
        //- Thermo density
        autoPtr<volScalarField> thermo.rho_;
        
        //- Temperature
        autoPtr<volScalarField> T_;

        //- forces
        autoPtr<volSymmTensorField> sigma_;

        //- Reziprocal diagonal elements of convection-diffusion matrix
        autoPtr<volScalarField> rAU_;

        //- Reference pressure cell
        label pRefCell_;

        //- Reference pressure value
        scalar pRefValue_;

        //- Moving reference frames
        MovingReferenceFrameZones mrfZones_;
        label myTimeIndex_;

        //- For adjustable time-step, this is the maximum Courant number
        scalar maxCo_;

        //- For adjustable time-step, this is the maximum time-step
        scalar maxDeltaT_;

        //- Correct phi?
        Switch correctPhi_;

        //- Check mesh Courant number?
        Switch checkMeshCourantNo_;

        //- Continuity errors
        scalar sumLocalContErr_;
        scalar globalContErr_;
        scalar cumulativeContErr_;

        //- Cached time derivative matrix
        //  shared between the momentum predictor and pressure corrector
        tmp<fvVectorMatrix> tddtUEqn;

        //- Cached convection-diffusion matrix
        //  shared between the momentum predictor and pressure corrector
        tmp<fvVectorMatrix> tHUEqn;

public:

    TypeName("rhoPimpleFluid");

    // Constructors

        rhoPimpleFluid
        (
            const Time& runTime,
            const word& regionName
        );

    //- Destructor
    virtual ~rhoPimpleFluid();


    // Member Functions

        // Access
        const volVectorField& U() const
        {
            return U_();
        }

        const volScalarField& p() const
        {
            return p_();
        }
        
        const volScalarField& rho() const
        {
            return rho_();
        }

        volVectorField& U()
        {
            return U_();
        }

        volScalarField& p()
        {
            return p_();
        }
        
        volScalarField& p_Abs()
        {
            return p_Abs();
        }
        
        volScalarField& psi_()
        {
            return psi_();
        }
        
        volScalarField& rho()
        {
            return rho_();
        }
        
        volScalarField& rho()
        {
            return thermo.rho_();
        }

        virtual word regionTypeName() const
        {
            return rhoPimpleFluid::typeName;
        }

        virtual bool usesPIMPLE() const
        {
            return true;
        }

        virtual void correct();

        virtual scalar getMinDeltaT();

        virtual void setCoupledEqns();

        virtual void solveRegion();

        virtual void prePredictor();

        virtual void momentumPredictor();

        virtual void pressureCorrector();

        virtual void meshMotionCorrector();

        virtual void postSolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
