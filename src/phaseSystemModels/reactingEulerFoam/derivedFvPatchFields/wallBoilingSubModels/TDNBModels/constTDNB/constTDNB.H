/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd
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
    Foam::wallBoilingModels:TDNBModels:::constTDNB

Description
    Departure from nulceate boiling constant value.

SourceFiles
    constTDNB.C

\*---------------------------------------------------------------------------*/

#ifndef constTDNB_H
#define constTDNB_H

#include "TDNBModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace TDNBModels
{

/*---------------------------------------------------------------------------*\
                         Class constTDNB Declaration
\*---------------------------------------------------------------------------*/

class constTDNB
:
    public TDNBModel
{

    // Private data:

        scalar TDNB_;

public:

    //- Runtime type information
    TypeName("const");

    // Constructors

        //- Construct from a dictionary
        constTDNB(const dictionary& dict);


    //- Destructor
    virtual ~constTDNB();


    // Member Functions

        //- Calculate and return the nucleation-site density
        virtual tmp<scalarField> TDNB
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L
        ) const;


        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TDNBModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
