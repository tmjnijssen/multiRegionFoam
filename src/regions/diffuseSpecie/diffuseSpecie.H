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
    Foam::diffuseSpecie.C

SourceFiles
    diffuseSpecie.C

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


#ifndef diffuseSpecie_H
#define diffuseSpecie_H

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
                           Class diffuseSpecie Declaration
\*---------------------------------------------------------------------------*/

class diffuseSpecie
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

        //- carbamate reaction rate
        volScalarField k1_;

        //- carbamate equilibrium constant
        volScalarField K1_;

        //- bicarbonate reaction rate
        volScalarField k2_;

        //- bicarbonate equilibrium constant
        volScalarField K2_;

        //- pore diffusion coefficients
        volScalarField DporeCO2_;
        volScalarField DporeH2O_;

        //- Henry coefficients
        volScalarField HrCO2_;
        volScalarField HrH2O_;

        //- Adsorption rates
        volScalarField dqdtCO2_;
        volScalarField dqdtH2O_;

        //- Concentrations
        volScalarField R2NH_;
        volScalarField R2NH2p_;
        volScalarField R2NCO2m_;
        volScalarField HCO3m_;
         
        //- Cached concentration matrices
        tmp<fvScalarMatrix> CO2Eqn;
        tmp<fvScalarMatrix> H2OEqn;

protected:

        //- Gas phase concentrations
         autoPtr<volScalarField> CO2_;
         autoPtr<volScalarField> H2O_;


public:

    TypeName("diffuseSpecie");

    // Constructors

        diffuseSpecie
        (
            const Time& runTime,
            const word& regionName
        );

    //- Destructor
    virtual ~diffuseSpecie();

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
            return diffuseSpecie::typeName;
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
