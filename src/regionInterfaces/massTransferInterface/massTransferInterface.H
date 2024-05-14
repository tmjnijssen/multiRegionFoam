/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::massTransferInterface.C

SourceFiles
    massTransferInterface.C

\*---------------------------------------------------------------------------*/

#ifndef massTransferInterface_H
#define massTransferInterface_H

#include "regionInterfaceType.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFieldsFwd.H"
#include "IOReferencer.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionInterfaces
{

/*---------------------------------------------------------------------------*\
                           Class massTransferInterface Declaration
\*---------------------------------------------------------------------------*/

class massTransferInterface
:
    public regionInterfaceType
{
    // Private data

        // interface dictionary
        dictionary dict_;

        // Surface tension
        dimensionedScalar hlv0_;

        areaScalarField hlv_;

        // Surface tension
        dimensionedScalar TSat0_;
        // Surface tension
        dimensionedScalar mDots0_;

        areaScalarField mDotInterface;


protected:

    // Private data

        //- surface flux field
        mutable autoPtr<areaScalarField> mDotsPtr_;

private:

    // Private Member Functions

        //- Make surface flux
        void makeMDotS() const;

protected:

    // Protected Member Functions

        //- Clear all demand-driven data
        virtual void clearOut() const;

public:

    TypeName("massTransferInterface");

    // Constructors

        massTransferInterface
        (
            const word& type,
            const dictionary& dict,
            const Time& runTime,
            const fvPatch& patchA,
            const fvPatch& patchB
        );

    // Destructor
    virtual ~massTransferInterface()
    {}

    // Public Member Functions

        // Access

        // Correct interface physics - here surface tension
        void correct();

        // Get min time step based on stability criteria
        scalar getMinDeltaT();

        // Edit

        // Update surface flux
        void updateMDotS();

        // Update surface flux
        scalarField saturatedFlux();

        //-- return surface flux
        areaScalarField& mDotS()
        {
            if (mDotsPtr_.empty())
            {
                makeMDotS();
            }

            return mDotsPtr_();
        }

        const areaScalarField& mDotS() const
        {
            if (mDotsPtr_.empty())
            {
               makeMDotS();
            }

            return mDotsPtr_();
        }

        const areaScalarField& hlv() const
        {
            return hlv_;
        }

        // Log interface specific information
        void info() const
        {};
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionInterfaces
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
