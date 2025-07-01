/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    regionCoupledMassTransferVelocityFlux

Author
    Holger Marschall (marschall@csi.tu-darmstadt.de)

Description
    Neumann boundary condition for velocity within the ALE interface tracking.
    Jump in velocity normal flux is implemented using the normal velocity
    gradient jump from tangential momentum jump condition.

    Example usage:

	    interfacePatchShadow
        {
            type                    regionCoupledMassTransferVelocityFlux;
            k                       "mu";
            neighbourRegionName     gas;
            neighbourPatchName      interfacePatch;
            nonOrthCorr             yes;
            secondOrder             yes;
            value                   uniform (0 0 0);
        }


SourceFiles
    regionCoupledMassTransferVelocityFlux.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledMassTransferVelocityFlux_H
#define regionCoupledMassTransferVelocityFlux_H

#include "genericRegionCoupledFluxFvPatchField.H"
#include "capillaryInterface.H"
#include "massTransferInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class regionCoupledMassTransferVelocityFlux Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledMassTransferVelocityFlux
:
    public genericRegionCoupledFluxFvPatchField<vector>
{

public:

    //- Runtime type information
    TypeName("regionCoupledMassTransferVelocityFlux");


    // Constructors

        //- Construct from patch and internal field
        regionCoupledMassTransferVelocityFlux
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        regionCoupledMassTransferVelocityFlux
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given regionCoupledMassTransferVelocityFlux
        //  onto a new patch
        regionCoupledMassTransferVelocityFlux
        (
            const regionCoupledMassTransferVelocityFlux&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
            	new regionCoupledMassTransferVelocityFlux(*this)
            );
        }

        //- Construct as copy setting internal field reference
        regionCoupledMassTransferVelocityFlux
        (
            const regionCoupledMassTransferVelocityFlux&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new regionCoupledMassTransferVelocityFlux(*this, iF)
            );
        }

    // Member functions

        //- Calculate pressure flux jump
        virtual tmp<vectorField> fluxJump() const;
        virtual word kName() const
        {
            word name = "mu";
            return name;
        };
        const regionInterfaces::capillaryInterface& capInterface() const;
        const regionInterfaces::massTransferInterface& massTrInterface() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
