/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016-2018 OpenFOAM Foundation
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
    Foam::wallBoilingModels::nucleationSiteModels::LemmertChawlaCustom

Description
    Lemmert & Chawla function for nucleation site density,
    correlation by Egorov & Menter.

    References:
    \verbatim
        Lemmert, M., & Chawla, J. M. (1977).
        Influence of flow velocity on surface boiling heat transfer coefficient.
        Heat Transfer in Boiling, 237, 247.

        Egorov, Y., & Menter, F. (2004).
        Experimental implementation of the RPI wall boiling model in CFX-5.6.
        Staudenfeldweg, 12, 83624.
        Technical Report ANSYS/TR-04-10, ANSYS Gmbh.
    \endverbatim

SourceFiles
    LemmertChawlaCustom.C

\*---------------------------------------------------------------------------*/

#ifndef LemmertChawlaCustom_H
#define LemmertChawlaCustom_H

#include "nucleationSiteModelCustom.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace nucleationSiteModels
{

/*---------------------------------------------------------------------------*\
                         Class LemmertChawlaCustom Declaration
\*---------------------------------------------------------------------------*/

class LemmertChawlaCustom
:
    public nucleationSiteModelCustom
{

    // Private data:

        //- Coefficient for nucleation site density
        scalar Cn_;

public:

    //- Runtime type information
    TypeName("LemmertChawla");

    // Constructors

        //- Construct from a dictionary
        LemmertChawlaCustom(const dictionary& dict);


    //- Destructor
    virtual ~LemmertChawlaCustom();


    // Member Functions

        //- Calculate and return the nucleation-site density
        virtual tmp<scalarField> N
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L,
            const scalarField& dDep,
            const scalarField& fDep
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nucleationSiteModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
