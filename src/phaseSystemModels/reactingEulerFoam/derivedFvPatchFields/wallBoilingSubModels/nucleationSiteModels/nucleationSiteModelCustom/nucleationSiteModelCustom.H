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
    Foam::wallBoilingModels::nucleationSiteModelCustom

Description
    Base class for nucleation site density models

SourceFiles
    nucleationSiteModelCustom.C
    newnucleationSiteModelCustom.C

\*---------------------------------------------------------------------------*/

#ifndef nucleationSiteModelCustom_H
#define nucleationSiteModelCustom_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{

/*---------------------------------------------------------------------------*\
                         Class nucleationSiteModelCustom Declaration
\*---------------------------------------------------------------------------*/

class nucleationSiteModelCustom
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        nucleationSiteModelCustom(const nucleationSiteModelCustom&);

        //- Disallow default bitwise assignment
        void operator=(const nucleationSiteModelCustom&);


public:

    //- Runtime type information
    TypeName("nucleationSiteModelCustom");


    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        nucleationSiteModelCustom,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct null
        nucleationSiteModelCustom();


    // Selectors

        //- Select null constructed
        static autoPtr<nucleationSiteModelCustom> New(const dictionary& dict);


    //- Destructor
    virtual ~nucleationSiteModelCustom();


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
        ) const = 0;

        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
