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
    Foam::NTGK.C

SourceFiles
    NTGK.C

Author
    Holger Marschall (h.marschall@ianus-simulation.de, Affiliation A)
    Lucas Hüfner (l.huefner@ianus-simulation.de, Affiliation A)

Contact
    Holger Marschall (h.marschall@ianus-simulation.de)
    main developer and principal investigator
    IANUS Simulation GmbH

Affiliations
    Affiliation A)
    IANUS Simulation GmbH
    Department Product Devlopment & Innovation
    Darmstadt, Germany

Description
    
    This file is part of the multiRegionFoam library.

    You may refer to this software as :
    https://doi.org/10.48550/arXiv.2306.01924

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

    Semi-empirical class for Li Ion Battery Modeling solving differential
    equations for charge and energy transport. Based on these sources:

    Kwon, Ki Hyun, et al. "A two-dimensional modeling of a lithium-polymer
    battery." Journal of Power Sources 163.1 (2006): 151-157.
    
    Kim, Gi-Heon, Ahmad Pesaran, and Robert Spotnitz. "A three-dimensional
    thermal abuse model for lithium-ion cells."
    Journal of power sources 170.2 (2007): 476-489.


\*---------------------------------------------------------------------------*/


#ifndef NTGK_H
#define NTGK_H

#include "regionType.H"
#include "dimensionedScalar.H"
#include "multiRegionSystem.H"
#include "volFields.H"
#include "IOReferencer.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{

/*---------------------------------------------------------------------------*\
                       Class NTGK Declaration
\*---------------------------------------------------------------------------*/

class NTGK
:
    public regionType
{
        //- name
        word regionName_;

        //- transport properties dictionary
        IOdictionary transportProperties_;
        
        //- electrochemical properties dictionary
        IOdictionary electrochemicalProperties_;

        //- electric conductivity positive electrode material
        dimensionedScalar sigmaPos_;

        //- electric conductivity negative electrode material
        dimensionedScalar sigmaNeg_;

        //- density of jelly roll
        dimensionedScalar rho_;

        //- specific heat capacity of jelly roll
        dimensionedScalar cp_;

        //- thermal conductivity of jelly roll
        dimensionedScalar k_;

        //- reference temperature
        dimensionedScalar TRef_;

        //- specific area of the battery
        dimensionedScalar spArea_;

        //- fitting parameter a0 for fit U0
        dimensionedScalar a0_;

        //- fitting parameter a1 for fit U0
        dimensionedScalar a1_;

        //- fitting parameter a2 for fit U0
        dimensionedScalar a2_;

        //- fitting parameter a3 for fit U0
        dimensionedScalar a3_;

        //- fitting parameter a4 for fit Y0
        dimensionedScalar a4_;

        //- fitting parameter a5 for fit Y0
        dimensionedScalar a5_;

        //- fitting parameter a6 for fit Y0
        dimensionedScalar a6_;

        //- fitting parameter a7 for fit Y0
        dimensionedScalar a7_;

        //- fitting parameter a8 for fit Y0
        dimensionedScalar a8_;

        //- fitting parameter a9 for fit Y0
        dimensionedScalar a9_;

        //- fitting parameter C1 for fit T-dependency of Y
        dimensionedScalar C1_;

        //- fitting parameter C2 for fit T-dependency of U
        dimensionedScalar C2_;

        //- depth of discharge
        autoPtr<volScalarField> DOD_;

        //- T-dependent open circuit potential
        autoPtr<volScalarField> U_;
        
        //- volumetric exchange current density
        autoPtr<volScalarField> j_;

        //- heat source term
        autoPtr<volScalarField> ST_;

        //- Cached positive electrode potential matrix
        tmp<fvScalarMatrix> faiPosEqn;

        //- Cached negative electrode potential matrix
        tmp<fvScalarMatrix> faiNegEqn;

        //- Cached temperature matrix
        tmp<fvScalarMatrix> TEqn;
        
    //private member functions
        //- calculate electrochemical parameters like current density with semi-empirical approach
        void calculateElectrochemicalParameters();
        
        //- calculate thermal behavior variables
        void calculateThermalBehavior();

        // calculate thermal abuse model variables
        void calculateThermalAbuse();

protected:

        //- positive electrode potential
        autoPtr<volScalarField> faiPos_;

        //- negative electrode potential
        autoPtr<volScalarField> faiNeg_;

        //- temperature
        autoPtr<volScalarField> T_;

public:

    TypeName("NTGK");

    // Constructors

        NTGK
        (
            const Time& runTime,
            const word& regionName
        );

    //- Destructor
    virtual ~NTGK();

    // Member Functions

        // Access
        const volScalarField& faiPos() const
        {
            return faiPos_();
        }

        volScalarField& faiPos()
        {
            return faiPos_();
        }

        const volScalarField& faiNeg() const
        {
            return faiNeg_();
        }

        volScalarField& faiNeg()
        {
            return faiNeg_();
        }
        
        const volScalarField& T() const
        {
            return T_();
        }

        volScalarField& T()
        {
            return T_();
        }

        virtual word regionTypeName() const
        {
            return NTGK::typeName;
        }

        virtual bool usesPIMPLE() const
        {
            return false;
        }

        virtual void correct();

        virtual scalar getMinDeltaT();

        virtual void solveRegion();

        virtual void prePredictor();

        virtual void momentumPredictor();

        virtual void pressureCorrector();

        virtual void meshMotionCorrector();

        virtual void setCoupledEqns();

        virtual void postSolve();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
