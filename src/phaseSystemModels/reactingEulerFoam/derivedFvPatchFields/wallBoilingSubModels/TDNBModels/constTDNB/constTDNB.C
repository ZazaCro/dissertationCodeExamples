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

#include "constTDNB.H"
#include "addToRunTimeSelectionTable.H"
#include "physicoChemicalConstants.H"

using Foam::constant::physicoChemical::R;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace TDNBModels
{
    defineTypeNameAndDebug(constTDNB, 0);
    addToRunTimeSelectionTable
    (
        TDNBModel,
        constTDNB,
        dictionary
    );
}
}
}

using Foam::constant::physicoChemical::R;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::TDNBModels::constTDNB::constTDNB
(
    const dictionary& dict
)
:
    TDNBModel(),
    TDNB_(dict.lookupOrDefault<scalar>("TDNB", 300))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallBoilingModels::TDNBModels::constTDNB::~constTDNB()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::TDNBModels::constTDNB::TDNB
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
            TDNB_
        )
    );
}


void Foam::wallBoilingModels::TDNBModels::constTDNB::write
(
    Ostream& os
) const
{
    TDNBModel::write(os);
    os.writeKeyword("TDNB") << TDNB_ << token::END_STATEMENT << nl;
}

// ************************************************************************* //


