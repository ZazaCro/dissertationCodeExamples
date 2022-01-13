/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "yPlusPhase.H"
#include "volFields.H"
#include "turbulenceModel.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(yPlusPhase, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        yPlusPhase,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::yPlusPhase::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "y+ ()");

    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "average");
    os << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::yPlusPhase::yPlusPhase
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    phase_(dict.lookupOrDefault<word>("phase", "None"))
{
    read(dict);

    writeFileHeader(file());
   
    if (phase_ != "None")
    {
        fieldName_ = IOobject::groupName(typeName, phase_);
    }
    else
    {
        fieldName_ = typeName;
    }

    volScalarField* yPlusPhasePtr
    (
        new volScalarField
        (
            IOobject
            (
                fieldName_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        )
    );

    mesh_.objectRegistry::store(yPlusPhasePtr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::yPlusPhase::~yPlusPhase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::yPlusPhase::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    return true;
}


bool Foam::functionObjects::yPlusPhase::execute()
{
    volScalarField& yPlusPhase =
        lookupObjectRef<volScalarField>(fieldName_);

    word turbulenceProperties;
    if (phase_ != "None")
    {
        turbulenceProperties =
            IOobject::groupName(turbulenceModel::propertiesName, phase_);
    }
    else
    {
        turbulenceProperties = turbulenceModel::propertiesName;
    }
        
    if (foundObject<turbulenceModel>(turbulenceProperties))
    {
        volScalarField::Boundary& yPlusPhaseBf = yPlusPhase.boundaryFieldRef();

        const turbulenceModel& model =
            lookupObject<turbulenceModel>
            (
                turbulenceProperties
            );

        const nearWallDist nwd(mesh_);
        const volScalarField::Boundary& d = nwd.y();

        // nut needed for wall function patches
        tmp<volScalarField> tnut = model.nut();
        const volScalarField::Boundary& nutBf = tnut().boundaryField();

        // U needed for plain wall patches
        const volVectorField::Boundary& UBf = model.U().boundaryField();

        const fvPatchList& patches = mesh_.boundary();

        forAll(patches, patchi)
        {
            const fvPatch& patch = patches[patchi];

            if (isA<nutWallFunctionFvPatchScalarField>(nutBf[patchi]))
            {
                const nutWallFunctionFvPatchScalarField& nutPf =
                    dynamic_cast<const nutWallFunctionFvPatchScalarField&>
                    (
                        nutBf[patchi]
                    );

                yPlusPhaseBf[patchi] = nutPf.yPlus();
            }
            else if (isA<wallFvPatch>(patch))
            {
                yPlusPhaseBf[patchi] =
                    d[patchi]
                   *sqrt
                    (
                        model.nuEff(patchi)
                       *mag(UBf[patchi].snGrad())
                    )/model.nu(patchi);
            }
        }
    }
    else
    {
        WarningInFunction
            << "Unable to find turbulence model in the "
            << "database: yPlusPhase will not be calculated" << endl;
        return false;
    }

    return true;
}


bool Foam::functionObjects::yPlusPhase::write()
{
    const volScalarField& yPlusPhase =
        obr_.lookupObject<volScalarField>(fieldName_);

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << yPlusPhase.name() << endl;

    yPlusPhase.write();

    const volScalarField::Boundary& yPlusPhaseBf = yPlusPhase.boundaryField();
    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];

        if (isA<wallFvPatch>(patch))
        {
            const scalarField& yPlusPhasep = yPlusPhaseBf[patchi];

            const scalar minYplus = gMin(yPlusPhasep);
            const scalar maxYplus = gMax(yPlusPhasep);
            const scalar avgYplus = gAverage(yPlusPhasep);

            if (Pstream::master())
            {
                Log << "    patch " << patch.name()
                    << " y+ : min = " << minYplus << ", max = " << maxYplus
                    << ", average = " << avgYplus << nl;

                writeTime(file());
                file()
                    << token::TAB << patch.name()
                    << token::TAB << minYplus
                    << token::TAB << maxYplus
                    << token::TAB << avgYplus
                    << endl;
            }
        }
    }

    return true;
}


// ************************************************************************* //
