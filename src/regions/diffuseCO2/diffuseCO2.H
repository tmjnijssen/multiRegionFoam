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
    Foam::diffuseCO2.C

SourceFiles
    diffuseCO2.C

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


#ifndef diffuseCO2_H
#define diffuseCO2_H

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
                           Class diffuseCO2 Declaration
\*---------------------------------------------------------------------------*/

class diffuseCO2
:
    public regionType
{
        //- name
        word regionName_;

        //- transport properties dictionary
        IOdictionary sorbentProperties_;

        //- specific heat
        // dimensionedScalar t_;
        volScalarField t_;

        //- density
        // dimensionedScalar eps_;
        volScalarField eps_;

        // dimensionedScalar Dpore_;

        volScalarField Dpore_;

        //- Pore diffusion field and adsorption rate
        autoPtr<volScalarField> Dp_;

         autoPtr<volScalarField> q_;

        //- Cached temperature matrix
        tmp<fvScalarMatrix> CO2Eqn;

protected:

        //- temperature
        autoPtr<volScalarField> CO2_;

public:

    TypeName("diffuseCO2");

    // Constructors

        diffuseCO2
        (
            const Time& runTime,
            const word& regionName
        );

    //- Destructor
    virtual ~diffuseCO2();

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

        virtual word regionTypeName() const
        {
            return diffuseCO2::typeName;
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
