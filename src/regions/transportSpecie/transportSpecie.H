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
    Foam::transportSpecie.C

SourceFiles
    transportSpecie.C

Author
    Holger Marschall (holger.marschall@tu-darmstadt.de, Affiliation A)
    Constantin Habes (constantin.habes@tu-darmstadt.de, Affiliation A)
    Heba Alkafri (heba.alkafri@tu-darmstadt.de, Affiliation A)

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


#ifndef transportSpecie_H
#define transportSpecie_H

#include "regionType.H"
#include "dimensionedScalar.H"
#include "multiRegionSystem.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOReferencer.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{

/*---------------------------------------------------------------------------*\
                           Class transportSpecie Declaration
\*---------------------------------------------------------------------------*/

class transportSpecie
:
    public regionType
{
        //- name
        word regionName_;

        //- transport
        //volVectorField U_;

        //surfaceScalarField phi_;

        //- transport properties dictionary
        IOdictionary transportProperties_;

        //- thermal conductivity
        dimensionedScalar D_;

        //- tortuosity
        dimensionedScalar s_;

        //- tortuosity
        dimensionedScalar eps_;

        //- transport
        autoPtr<volVectorField> U_;

        //- transport density
        autoPtr<volVectorField> rho_;

        //- thermal conductivity field
        autoPtr<volScalarField> De_;

        //- flux
        autoPtr<surfaceScalarField> phi_;

        //- Cached temperature matrix
        tmp<fvScalarMatrix> DEqn;

protected:

        //- CO2 field
        autoPtr<volScalarField> CO2_;

public:

    TypeName("transportSpecie");

    // Constructors

        transportSpecie
        (
            const Time& runTime,
            const word& regionName
        );

    //- Destructor
    virtual ~transportSpecie();

    // Member Functions

        // Access
        const volScalarField& T() const
        {
            return CO2_();
        }

        volScalarField& T()
        {
            return CO2_();
        }

        virtual word regionTypeName() const
        {
            return transportSpecie::typeName;
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
