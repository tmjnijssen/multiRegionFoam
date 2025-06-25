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
    Foam::separator.C

SourceFiles
    separator.C

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

\*---------------------------------------------------------------------------*/


#ifndef separator_H
#define separator_H

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
                       Class separator Declaration
\*---------------------------------------------------------------------------*/

class separator
:
    public regionType
{
        //- name
        word regionName_;

        //- transport properties dictionary
        IOdictionary transportProperties_;
        
        //- electrochemical properties dictionary
        IOdictionary electrochemicalProperties_;
        
        //- volume fraction electrolyte phase
        dimensionedScalar epsE_;
        
        //- Bruggeman coefficient for porous media
        dimensionedScalar brugg_;
        
        //- transference number
        dimensionedScalar tNo_;
        
        //- diffusional electrolyte conductivity field
        dimensionedScalar Gamma_;
        
        //- dimensions of concentration
        dimensionedScalar dimC_;

        //- density
        dimensionedScalar rho_;

        //- heat capacity
        dimensionedScalar cp_;

        //- heat conductivity
        dimensionedScalar kT_;

        //- electrolyte conductivity field
        autoPtr<volScalarField> kappa_;
         
        //- diffusion coefficient ion in electrolyte
        autoPtr<volScalarField> DE_;

        //- heat source term
        autoPtr<volScalarField> ST_;
        
        //- Cached electrolyte potential matrix
        tmp<fvScalarMatrix> faiEEqn;
        
        //- Cached ion concentration matrix
        tmp<fvScalarMatrix> cEEqn;

        //- Cached temperature matrix
        tmp<fvScalarMatrix> TEqn;
        
    //private member functions
        // calculate transport coefficients
        void calculateTransportCoeffs();
         
        //- calculate heat sources and temperature dependency of open circuit potential 
        void calculateHeatSourceTerms();

protected:

        //- electrolyte potential
        autoPtr<volScalarField> faiE_;
        
        //- ion concentration
        autoPtr<volScalarField> cE_;
        
        //- temperature
        autoPtr<volScalarField> T_;

public:

    TypeName("separator");

    // Constructors

        separator
        (
            const Time& runTime,
            const word& regionName
        );

    //- Destructor
    virtual ~separator();

    // Member Functions

        // Access
        const volScalarField& faiE() const
        {
            return faiE_();
        }

        volScalarField& faiE()
        {
            return faiE_();
        }
        
        const volScalarField& cE() const
        {
            return cE_();
        }

        volScalarField& cE()
        {
            return cE_();
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
            return separator::typeName;
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
