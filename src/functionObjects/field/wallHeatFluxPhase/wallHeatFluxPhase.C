/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010, 2016-2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016-2017 OpenFOAM Foundation
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

#include "wallHeatFluxPhase.H"
#include "turbulentFluidThermoModel.H"
#include "solidThermo.H"
#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallHeatFluxPhase, 0);
    addToRunTimeSelectionTable(functionObject, wallHeatFluxPhase, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::wallHeatFluxPhase::writeFileHeader(Ostream& os) const
{
    // Add headers to output data
    writeHeader(os, "Wall heat-flux");
    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "integral");
    os  << endl;
}


void Foam::functionObjects::wallHeatFluxPhase::calcHeatFlux
(
    const volScalarField& alpha,
    const volScalarField& he,
    volScalarField& wallHeatFluxPhase
)
{
    volScalarField::Boundary& wallHeatFluxPhaseBf = wallHeatFluxPhase.boundaryFieldRef();

    const volScalarField::Boundary& heBf = he.boundaryField();

    const volScalarField::Boundary& alphaBf = alpha.boundaryField();

    forAll(wallHeatFluxPhaseBf, patchi)
    {
        if (!wallHeatFluxPhaseBf[patchi].coupled())
        {
            wallHeatFluxPhaseBf[patchi] = alphaBf[patchi]*heBf[patchi].snGrad();
        }
    }

    if (foundObject<volScalarField>(qrName_))
    {
        const volScalarField& qr = lookupObject<volScalarField>(qrName_);

        const volScalarField::Boundary& radHeatFluxBf = qr.boundaryField();

        forAll(wallHeatFluxPhaseBf, patchi)
        {
            if (!wallHeatFluxPhaseBf[patchi].coupled())
            {
                wallHeatFluxPhaseBf[patchi] -= radHeatFluxBf[patchi];
            }
        }
    }
}

void Foam::functionObjects::wallHeatFluxPhase::calcHeatFluxPhase
(
    const volScalarField& alpha,
    const volScalarField& he,
    const volScalarField& phaseFraction,
    volScalarField& wallHeatFluxPhase
)
{
    volScalarField::Boundary& wallHeatFluxPhaseBf = wallHeatFluxPhase.boundaryFieldRef();

    const volScalarField::Boundary& heBf = he.boundaryField();

    const volScalarField::Boundary& alphaBf = alpha.boundaryField();

    const volScalarField::Boundary& phaseFractionBf = phaseFraction.boundaryField();
    
    forAll(wallHeatFluxPhaseBf, patchi)
    {
        if (!wallHeatFluxPhaseBf[patchi].coupled())
        {
            wallHeatFluxPhaseBf[patchi] =
                phaseFractionBf[patchi]*alphaBf[patchi]*heBf[patchi].snGrad();
        }
    }

    if (foundObject<volScalarField>(qrName_))
    {
        const volScalarField& qr = lookupObject<volScalarField>(qrName_);

        const volScalarField::Boundary& radHeatFluxBf = qr.boundaryField();

        forAll(wallHeatFluxPhaseBf, patchi)
        {
            if (!wallHeatFluxPhaseBf[patchi].coupled())
            {
                wallHeatFluxPhaseBf[patchi] -= radHeatFluxBf[patchi];
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallHeatFluxPhase::wallHeatFluxPhase
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    patchSet_(),
    qrName_("qr"),
    phase_(dict.lookupOrDefault<word>("phase", "None"))
{
    if (phase_ != "None")
    {
        fieldName_ = IOobject::groupName(typeName, phase_);
    }
    else
    {
        fieldName_ = typeName;
    }

    volScalarField* wallHeatFluxPhasePtr
    (
        new volScalarField
        (
            IOobject
            (
                fieldName_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimMass/pow3(dimTime), Zero)
        )
    );

    mesh_.objectRegistry::store(wallHeatFluxPhasePtr);

    read(dict);

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::wallHeatFluxPhase::~wallHeatFluxPhase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallHeatFluxPhase::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    dict.readIfPresent("qr", qrName_);

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        for (const label patchi : patchSet_)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall heat-flux on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }

    return true;
}


bool Foam::functionObjects::wallHeatFluxPhase::execute()
{
    volScalarField& wallHeatFluxPhase =
        lookupObjectRef<volScalarField>(fieldName_);

    if 
    (
        foundObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        )
    )
    {
        const compressible::turbulenceModel& turbModel =
            lookupObject<compressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );
        
        calcHeatFlux
        (
            turbModel.alphaEff()(),
            turbModel.transport().he(),
            wallHeatFluxPhase
        );
    }
    else if (phase_ != "None")
    {
        const phaseSystem& fluid =
            lookupObject<phaseSystem>
            (
                "phaseProperties"
            );
    
        const phaseModel& phaseFraction(fluid.phases()[phase_]);

        calcHeatFluxPhase
        (
            phaseFraction.alphaEff(),
            phaseFraction.thermo().he(),
            phaseFraction,
            wallHeatFluxPhase
        );
        
    }
    else if (foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const fluidThermo& thermo =
            lookupObject<fluidThermo>(fluidThermo::dictName);

        calcHeatFlux
        (
            thermo.alpha(),
            thermo.he(),
            wallHeatFluxPhase
        );
    }
    else if (foundObject<solidThermo>(solidThermo::dictName))
    {
        const solidThermo& thermo =
            lookupObject<solidThermo>(solidThermo::dictName);

        calcHeatFlux(thermo.alpha(), thermo.he(), wallHeatFluxPhase);
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find turbulence model in the "
            << "database" << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::wallHeatFluxPhase::write()
{
    const volScalarField& wallHeatFluxPhase =
        lookupObject<volScalarField>(fieldName_);

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << wallHeatFluxPhase.name() << endl;

    wallHeatFluxPhase.write();

    const fvPatchList& patches = mesh_.boundary();

    const surfaceScalarField::Boundary& magSf =
        mesh_.magSf().boundaryField();

    for (const label patchi : patchSet_)
    {
        const fvPatch& pp = patches[patchi];

        const scalarField& hfp = wallHeatFluxPhase.boundaryField()[patchi];

        const scalar minHfp = gMin(hfp);
        const scalar maxHfp = gMax(hfp);
        const scalar integralHfp = gSum(magSf[patchi]*hfp);

        if (Pstream::master())
        {
            writeTime(file());

            file()
                << token::TAB << pp.name()
                << token::TAB << minHfp
                << token::TAB << maxHfp
                << token::TAB << integralHfp
                << endl;
        }

        Log << "    min/max/integ(" << pp.name() << ") = "
            << minHfp << ", " << maxHfp << ", " << integralHfp << endl;
    }

    return true;
}


// ************************************************************************* //
