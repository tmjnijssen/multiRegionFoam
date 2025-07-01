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
    Foam::intercalationAnode.C

SourceFiles
    intercalationAnode.C

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


#ifndef intercalationAnode_H
#define intercalationAnode_H

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
                       Class intercalationAnode Declaration
\*---------------------------------------------------------------------------*/

class intercalationAnode
:
    public regionType
{
        //- name
        word regionName_;

        //- transport properties dictionary
        IOdictionary transportProperties_;

        //- electrochemical properties dictionary
        IOdictionary electrochemicalProperties_;

        //- electrode material
        word electrodeMaterial_;

        //- volume fraction fluid phase (porosity)
        dimensionedScalar epsE_;

        //- volume fraction filler material
        dimensionedScalar epsF_;

        //- volume fraction solid phase
        dimensionedScalar epsS_;

        //- particle radius
        dimensionedScalar r_;

        //- Bruggeman coefficient for porous media
        dimensionedScalar brugg_;

        //- specific particle surface area
        dimensionedScalar spArea_;

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

        //- solid diffusion field
        autoPtr<volScalarField> DS_;

        //- solid diffusion prefactor
        dimensionedScalar DSinit_;

        //- activation energy solid state diffusion
        dimensionedScalar EADiff_;

        //- applied current density
        dimensionedScalar IApp_;

        //- anodic transfer coefficient
        dimensionedScalar alphaA_;

        //- cathodic transfer coefficient
        dimensionedScalar alphaC_;

        //-  activation energy reaction
        dimensionedScalar EARct_;

        //- reaction rate constant
        dimensionedScalar kRct_;

        //- reference temperature for electrochemistry
        dimensionedScalar TRef_;

        //- maximum lithium concentration at the solid surface
        autoPtr<volScalarField> cSMax_;

        //- minimum lithium concentration at the solid surface
        autoPtr<volScalarField> cSMin_;

        //- electric conductivity field
        autoPtr<volScalarField> sigma_;

        //- state of charge
        autoPtr<volScalarField> SoC_;

        //- reversible voltage
        autoPtr<volScalarField> EEqRev_;

        //- temperature dependency reversible voltage (entropy change)
        autoPtr<volScalarField> dEEqdT_;

        //- electrode overpotential
        autoPtr<volScalarField> eta_;

        //- volumetric exchange current density
        autoPtr<volScalarField> j_;

        //- heat source term
        autoPtr<volScalarField> ST_;

        //- Cached solid electrode potential matrix
        tmp<fvScalarMatrix> faiSEqn;

        //- Cached electrolyte potential matrix
        tmp<fvScalarMatrix> faiEEqn;

        //- Cached solid surface lithium concentration matrix
        tmp<fvScalarMatrix> cSEqn;

        //- Cached ion concentration matrix
        tmp<fvScalarMatrix> cEEqn;

        //- Cached temperature matrix
        tmp<fvScalarMatrix> TEqn;

    //private member functions
        //- calculate open circuit voltage of electrode and reversible voltage
        void calculateOCV();

        //- calculate Butler-Volmer kinetics and vol. current density
        void calculateButlerVolmer();

        // calculate transport coefficients
        void calculateTransportCoeffs();

        //- calculate lithium concentration at interface of solid and electroylte
       void calculateInterfaceConcentration();

       //- calculate heat sources and temperature dependency of open circuit potential
       void calculateHeatSourceTerms();

protected:

        //- solid electrode potential
        autoPtr<volScalarField> faiS_;

        //- electrolyte potential
        autoPtr<volScalarField> faiE_;

        //- solid surface lithium concentration
        autoPtr<volScalarField> cS_;

        //- ion concentration
        autoPtr<volScalarField> cE_;

        //- lithium concentration at interface of solid and electrolyte
        autoPtr<volScalarField> cSE_;

        //- temperature
        autoPtr<volScalarField> T_;

public:

    TypeName("intercalationAnode");

    // Constructors

        intercalationAnode
        (
            const Time& runTime,
            const word& regionName
        );

    //- Destructor
    virtual ~intercalationAnode();

    // Member Functions

        // Access
        const volScalarField& faiS() const
        {
            return faiS_();
        }

        volScalarField& faiS()
        {
            return faiS_();
        }

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

        const volScalarField& cS() const
        {
            return cS_();
        }

        volScalarField& cS()
        {
            return cS_();
        }

        const volScalarField& cSE() const
        {
            return cSE_();
        }

        volScalarField& cSE()
        {
            return cSE_();
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
            return intercalationAnode::typeName;
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
