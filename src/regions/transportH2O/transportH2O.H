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
    Foam::transportH2O.C

SourceFiles
    transportH2O.C

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


#ifndef transportH2O_H
#define transportH2O_H

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
                           Class transportH2O Declaration
\*---------------------------------------------------------------------------*/

class transportH2O
:
    public regionType
{
        //- name
        word regionName_;

        //- transport
        //volVectorField U_;

        //surfaceScalarField phi_;

        //- transport
        autoPtr<volVectorField> U_;

        //- Effective diffusivity field
        autoPtr<volScalarField> De_;
        
        //- Effective diffusivity field
        autoPtr<volScalarField> rho_;


        //- flux
        autoPtr<surfaceScalarField> phi_;

        //- Cached temperature matrix
        tmp<fvScalarMatrix> H2OEqn;

protected:

        //- temperature
        autoPtr<volScalarField> H2O_;

public:

    TypeName("transportH2O");

    // Constructors

        transportH2O
        (
            const Time& runTime,
            const word& regionName
        );

    //- Destructor
    virtual ~transportH2O();

    // Member Functions

        // Access
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
            return transportH2O::typeName;
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

