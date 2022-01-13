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

\*---------------------------------------------------------------------------*/

#include "constantTLeid.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace LeidenfrostModels
{
    defineTypeNameAndDebug(constantTLeid, 0);
    addToRunTimeSelectionTable
    (
        LeidenfrostModel,
        constantTLeid,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::LeidenfrostModels::constantTLeid::constantTLeid
(
    const dictionary& dict
)
:
    LeidenfrostModel(),
    TLeid_(dict.lookupOrDefault<scalar>("TLeid", 315.56))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallBoilingModels::LeidenfrostModels::constantTLeid::~constantTLeid()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::LeidenfrostModels::constantTLeid::TLeid
(
    const phaseModel& liquid,
    const phaseModel& vapor,
    const label patchi,
    const scalarField& Tl,
    const scalarField& Tsatw,
    const scalarField& L
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            liquid.thermo().p().boundaryField()[patchi].size(),
            TLeid_
        )
    );
}


void Foam::wallBoilingModels::LeidenfrostModels::constantTLeid::write
(
    Ostream& os
) const
{
    LeidenfrostModel::write(os);
    os.writeKeyword("TLeid") << TLeid_ << token::END_STATEMENT << nl;
}

// ************************************************************************* //
