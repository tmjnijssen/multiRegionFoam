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
    Foam::diffuseAdsorbSpecie.C

SourceFiles
    diffuseAdsorbSpecie.C

Author
    Holger Marschall (holger.marschall@tu-darmstadt.de, Affiliation A)
    Constantin Habes (constantin.habes@tu-darmstadt.de, Affiliation A)
    Heba Alkafri (heba.alkafri@tu-darmstadt.de, Affiliation A)
    Steffen Hess (s.hess@fz-juelich.de, Affiliation B)

Contact
    Holger Marschall (holger.marschall@tu-darmstadt.de)
    main developer and principal investigator
    TU Darmstadt

Affiliations
    Affiliation A)
    Computational Multiphase Flow
    Department of Mathematics
    Technical University of Darmstadt, Germany

    Affiliation B) 
    Institute of Energy and Climate Research, 
    Forschungszentrum Jülich GmbH, Germany

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


#ifndef diffuseAdsorbSpecie_H
#define diffuseAdsorbSpecie_H

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
                           Class diffuseAdsorbSpecie Declaration
\*---------------------------------------------------------------------------*/

class diffuseAdsorbSpecie
:
    public regionType
{
        //- name
        word regionName_;

        //- transport properties dictionary
        IOdictionary transportProperties_;

        //- sorbent properties dictionary
        IOdictionary sorbentProperties_;        

        //- porosity
        volScalarField eps_;

        //- 1: carbamate reaction
        //- 2: carbonate reaction

        //- Coefficients adsorption equations
        dimensionedScalar A_;   // Fitting parameter Piscina 2024
        dimensionedScalar B_;   // Fitting parameter Piscina 2024 
        dimensionedScalar C_;   // Fitting parameter Piscina 2024  set to zero in sorbentProperties
        dimensionedScalar D_;  // Fitting parameter Piscina 2024
        dimensionedScalar F_;  // Fitting parameter Piscina 2024
        dimensionedScalar G_;  // Fitting parameter Piscina 2024
        dimensionedScalar H_;  // Fitting parameter Piscina 2024
        dimensionedScalar J_;  // Fitting parameter Piscina 2024  set to zero in sorbentProperties
        dimensionedScalar chi_;  // Fitting parameter Piscina 2024
        dimensionedScalar alpha_;  // Fitting parameter Piscina 2024
        dimensionedScalar bT_;     // Temperature dependent affinity Toth parameter Piscina 2024
        dimensionedScalar b0_;     // affinity parameter Toth at reference temperature
        dimensionedScalar tauT_; // Temperature dependent surface heterogenity Toth parameter Piscina 2024
        dimensionedScalar tau0_;  // Reference surface heterogenity at reference temperature 
        dimensionedScalar qCO2inf_;  // Saturation loading Piscina 2024
        dimensionedScalar qH2Oinf_;  // Saturation loading water Young 2021
        dimensionedScalar qCO2inf0_; // Maximum CO2 loading at reference temperature 
        dimensionedScalar dH2O_;        // fitting parameter GAB isotherm Young 2021
        dimensionedScalar fH2O_;     // fitting parameter GAB isotherm Young 2021
        dimensionedScalar gH2O_;        // fitting parameter GAB isotherm Young 2021
        dimensionedScalar qmH2O_;    // Monolayer H2O loading Young 2021 for GAB isotherm

        // Kinetic constants
        dimensionedScalar kCO2_;      // Kinetic Constant CO2 Driessen
        dimensionedScalar kH2O_;      // Kinetic Constant H2O Driessen 2020
        

        // Activation energy
        dimensionedScalar Ea1_;
        dimensionedScalar cH2O_;
        dimensionedScalar dH2O_;

        // Enthalpy of reaction
        dimensionedScalar dH1_;
        dimensionedScalar dH2_;

        // Reference temperature
        dimensionedScalar T0H2O_;
        dimensionedScalar T0CO2_;

        //- pore diffusion coefficients
        volScalarField DporeCO2_;
        volScalarField DporeH2O_;

        //- Relative Humidity
        volScalarField RH_;
        volScalarField pCO2_;
        volScalarField pCO2sat_;
        volScalarField Ctot_;

        //- Adsorption rates
        volScalarField dqdtCO2_;
        volScalarField dqdtCO2ex_;
        volScalarField dqdtCO2im_;
        volScalarField dqdtH2O_;
        volScalarField dqdtH2Oex_;
        volScalarField dqdtH2Oim_;

        //- Concentrations in solid
        volScalarField  qH2Oads_;
        volScalarField  qCO2ads_;

        //- Temperature
        autoPtr<volScalarField> T_;

        //- Heat source
        autoPtr<volScalarField> heatSource_;
         
        //- Cached concentration matrices
        tmp<fvScalarMatrix> CO2Eqn;
        tmp<fvScalarMatrix> H2OEqn;

protected:

        //- Gas phase concentrations
         autoPtr<volScalarField> CO2_;
         autoPtr<volScalarField> H2O_;

public:

    TypeName("diffuseAdsorbSpecie");

    // Constructors

        diffuseAdsorbSpecie
        (
            const Time& runTime,
            const word& regionName
        );

    //- Destructor
    virtual ~diffuseAdsorbSpecie();

    // Member Functions

        // Access
        const volScalarField& CO2() const
        {
            return CO2_();
        }

        volScalarField& CO2()
        {
            return CO2_();
        }

        const volScalarField& H2O() const
        {
            return H2O_();
        }

        volScalarField& H2O()
        {
            return H2O_();
        }

        virtual word regionTypeName() const
        {
            return diffuseAdsorbSpecie::typeName;
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

        virtual void setCoupledEqns(word fieldName = word::null);

        virtual void postSolve();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
