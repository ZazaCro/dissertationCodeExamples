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

#include "TFluid.H"
#include "volFields.H"
#include "turbulenceModel.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

#include "mappedPatchBase.H"
#include "solidThermo.H"
#include "phaseSystem.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(TFluid, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        TFluid,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::TFluid::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "TFluid ()");

    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "average");
    os << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::TFluid::TFluid
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    phaseName_(dict.lookupOrDefault<word>("phase", "None")),
    otherPhaseName_(dict.lookupOrDefault<word>("otherPhase", "None")),
    TnbrName_(dict.lookupOrDefault<word>("Tnbr", "T"))
{
    read(dict);

    writeFileHeader(file());
   
    if (phaseName_ != "None")
    {
        Tliquid_ = IOobject::groupName("T", phaseName_);
    }
    else
    {
        Tliquid_ = typeName;
    }

    if (otherPhaseName_ != "None")
    {
        Tgas_ = IOobject::groupName("T", otherPhaseName_);
    }
    else
    {
        Tgas_ = typeName;
    }

    volScalarField* TFluidPtr
    (
        new volScalarField
        (
            IOobject
            (
                typeName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimTemperature, Zero)
        )
    );

    mesh_.objectRegistry::store(TFluidPtr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::TFluid::~TFluid()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::TFluid::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    return true;
}


bool Foam::functionObjects::TFluid::execute()
{

    volScalarField& TFluid =
        lookupObjectRef<volScalarField>(typeName);

    volScalarField& TL =
        lookupObjectRef<volScalarField>(Tliquid_);

    /* volScalarField& TV = */
    /*     lookupObjectRef<volScalarField>(Tgas_); */

    volScalarField::Boundary& TFluidBf = TFluid.boundaryFieldRef();
    
    volScalarField::Boundary& TLBf = TL.boundaryFieldRef();
        
    /* volScalarField::Boundary& TVBf = TV.boundaryFieldRef(); */

    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];

        // Check whether the particular patch has neighbouring patch
        if (isA<mappedPatchBase>(patch.patch()))
        {
            const mappedPatchBase& mpp =
                refCast<const mappedPatchBase>(patch.patch());
            const polyMesh& nbrMesh = mpp.sampleMesh();
            const label samplePatchi = mpp.samplePolyPatch().index();
            const fvPatch& nbrPatch =
                refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];


            scalarField KDeltaNbr;
            const solidThermo& thermo = nbrMesh.lookupObject<solidThermo>(basicThermo::dictName);
            KDeltaNbr = thermo.kappa(samplePatchi)*nbrPatch.deltaCoeffs();
            mpp.distribute(KDeltaNbr);
            
            // Fluid region values
            const phaseSystem& fluid = 
            (
                mesh_.lookupObject<phaseSystem>("phaseProperties")
            ); 
           
            const phaseModel& liquid
            (
                fluid.phases()[phaseName_]
            );

           const phaseModel& vapor(fluid.phases()[otherPhaseName_]);

           const fvPatchScalarField& alphav = vapor.boundaryField()[patchi];
           const fvPatchScalarField& alphal = liquid.boundaryField()[patchi];
           
           const scalarField KdeltaVap
           (
            alphav*(vapor.kappaEff(patchi))*patch.deltaCoeffs()
           );  
           const scalarField KdeltaLiq
           (
            alphal*(liquid.kappaEff(patchi))*patch.deltaCoeffs()
           );
         
           const tmp<scalarField> tTLIf = TLBf[patchi].patchInternalField();
           const scalarField& TLIf = tTLIf();
       
           // The temperature field retrieve should be consistent, but currently
           // left due to lack of time. 
           const fvPatchScalarField& TVBf = 
               vapor.thermo().T().boundaryField()[patchi];
           const tmp<scalarField> tTVIf = TVBf.patchInternalField();
           const scalarField& TVIf = tTVIf();

           const scalarField BiV(KdeltaVap/KDeltaNbr);
           const scalarField BiL(KdeltaLiq/KDeltaNbr);  
           Info << gMax(TLIf) <<" " << gMin(TLIf) << endl;
           Info << gMax(TVIf) <<" " << gMin(TVIf) << endl;
           
           TFluidBf[patchi] = (BiL * TLIf + BiV * TVIf) / (BiL + BiV);
        }
    }
    
    return true;
}


bool Foam::functionObjects::TFluid::write()
{
    const volScalarField& TFluid =
        obr_.lookupObject<volScalarField>(typeName);

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << TFluid.name() << endl;

    TFluid.write();

    const volScalarField::Boundary& TFluidBf = TFluid.boundaryField();
    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];

        if (isA<mappedPatchBase>(patch.patch()))
        {
            const scalarField& TFluidp = TFluidBf[patchi];

            const scalar minTFluid = gMin(TFluidp);
            const scalar maxTFluid = gMax(TFluidp);
            const scalar avgTFluid = gAverage(TFluidp);

            if (Pstream::master())
            {
                Log << "    patch " << patch.name()
                    << " TFluid : min = " << minTFluid << ", max = " << maxTFluid
                    << ", average = " << avgTFluid << nl;

                writeTime(file());
                file()
                    << token::TAB << patch.name()
                    << token::TAB << minTFluid
                    << token::TAB << maxTFluid
                    << token::TAB << avgTFluid
                    << endl;
            }
        }
    }

    return true;
}


// ************************************************************************* //
