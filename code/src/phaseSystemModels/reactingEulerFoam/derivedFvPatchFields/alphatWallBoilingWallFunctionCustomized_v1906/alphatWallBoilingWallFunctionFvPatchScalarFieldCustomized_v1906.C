/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"

#include "phaseSystem.H"
#include "compressibleTurbulenceModel.H"
#include "ThermalDiffusivity.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "saturationModel.H"
#include "wallFvPatch.H"
#include "uniformDimensionedFields.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::compressible::
    alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906::phaseType
>
Foam::compressible::
alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906::phaseTypeNames_
{
    { phaseType::vaporPhase, "vapor" },
    { phaseType::liquidPhase, "liquid" },
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906::
alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(p, iF),
    otherPhaseName_("vapor"),
    phaseType_(liquidPhase),
    relax_(),
    AbyV_(p.size(), 0),
    regimeTypes_(p.size(), -1),
    alphatConv_(p.size(), 0),
    V_(p.size(), 0),
    dmdtV_(p.size(), 0),
    Tl_(p.size(), 0),
    dDep_(p.size(), 1e-5),
    fDep_(p.size(), 0),
    N_(p.size(), 0),
    qq_(p.size(), 0),
    K_(4),
    partitioningModel_(nullptr),
    nucleationSiteModel_(nullptr),
    departureDiamModel_(nullptr),
    departureFreqModel_(nullptr),
    filmBoilingModel_(nullptr),
    LeidenfrostModel_(nullptr),
    CHFModel_(nullptr),
    CHFSoobModel_(nullptr),
    MHFModel_(nullptr),
    TDNBModel_(nullptr),
    wp_(1)
{
    AbyV_ = this->patch().magSf();
    
    forAll(AbyV_, facei)
    {
        const label faceCelli = this->patch().faceCells()[facei];
        AbyV_[facei] /= iF.mesh().V()[faceCelli];
    }
    
    // Generate volumes of cells next to patch
    forAll(V_, facei)
    {
        const label faceCelli = this->patch().faceCells()[facei];
        V_[facei] = iF.mesh().V()[faceCelli];
    }
}


alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906::
alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(p, iF, dict),
    otherPhaseName_(dict.lookup("otherPhase")),
    phaseType_(phaseTypeNames_.read(dict.lookup("phaseType"))),
    relax_(Function1<scalar>::New("relax", dict)),
    AbyV_(p.size(), 0),
    regimeTypes_(p.size(), -1),
    alphatConv_(p.size(), 0),
    V_(p.size(), 0),
    dmdtV_(p.size(), 0),
    Tl_(p.size(), 0),
    dDep_(p.size(), 1e-5),
    fDep_(p.size(), 0),
    N_(p.size(), 0),
    qq_(p.size(), 0),
    K_(4),
    partitioningModel_(nullptr),
    nucleationSiteModel_(nullptr),
    departureDiamModel_(nullptr),
    departureFreqModel_(nullptr),
    filmBoilingModel_(nullptr),
    LeidenfrostModel_(nullptr),
    CHFModel_(nullptr),
    CHFSoobModel_(nullptr),
    MHFModel_(nullptr),
    TDNBModel_(nullptr),
    wp_(1)
{

    // Check that otherPhaseName != this phase
    if (internalField().group() == otherPhaseName_)
    {
        FatalErrorInFunction
            << "otherPhase should be the name of the vapor phase that "
            << "corresponds to the liquid base of vice versa" << nl
            << "This phase: " << internalField().group() << nl
            << "otherPhase: " << otherPhaseName_
            << abort(FatalError);
    }

    switch (phaseType_)
    {
        case vaporPhase:
        {
            partitioningModel_ =
                wallBoilingModels::partitioningModel::New
                (
                    dict.subDict("partitioningModel")
                );

            const dictionary* LeidenfrostDict =
                dict.findDict("LeidenfrostModel");

            if (LeidenfrostDict)
            {
                LeidenfrostModel_ =
                    wallBoilingModels::LeidenfrostModel::New(*LeidenfrostDict);
            }

            const dictionary* filmDict = dict.findDict("filmBoilingModel");

            if (filmDict)
            {
                filmBoilingModel_ =
                    wallBoilingModels::filmBoilingModel::New(*filmDict);
            }

            dmdt_ = 0;

            break;
        }
        case liquidPhase:
        {
            partitioningModel_ =
                wallBoilingModels::partitioningModel::New
                (
                    dict.subDict("partitioningModel")
                );

            nucleationSiteModel_ =
                wallBoilingModels::nucleationSiteModelCustom::New
                (
                    dict.subDict("nucleationSiteModel")
                );

            departureDiamModel_ =
                wallBoilingModels::departureDiameterModel::New
                (
                    dict.subDict("departureDiamModel")
                );

            departureFreqModel_ =
                wallBoilingModels::departureFrequencyModel::New
                (
                    dict.subDict("departureFreqModel")
                );

            const dictionary* LeidenfrostDict =
                dict.findDict("LeidenfrostModel");

            if (LeidenfrostDict)
            {
                LeidenfrostModel_ =
                    wallBoilingModels::LeidenfrostModel::New(*LeidenfrostDict);
            }

            const dictionary* CHFDict = dict.findDict("CHFModel");

            if (CHFDict)
            {
                CHFModel_ =
                    wallBoilingModels::CHFModel::New(*CHFDict);
            }

            const dictionary* HFSubCoolDict = dict.findDict("CHFSubCoolModel");

            if (HFSubCoolDict)
            {
                CHFSoobModel_ =
                    wallBoilingModels::CHFSubCoolModel::New(*HFSubCoolDict);
            }

            const dictionary* MHFDict = dict.findDict("MHFModel");

            if (MHFDict)
            {
                MHFModel_ =
                    wallBoilingModels::MHFModel::New(*MHFDict);
            }

            const dictionary* TDNBDict = dict.findDict("TDNBModel");

            if (TDNBDict)
            {
                TDNBModel_ =
                    wallBoilingModels::TDNBModel::New(*TDNBDict);
            }

            const dictionary* filmDict = dict.findDict("filmBoilingModel");

            if (filmDict)
            {
                filmBoilingModel_ =
                    wallBoilingModels::filmBoilingModel::New(*filmDict);
            }

            if (dict.found("dmdtV"))
            {
                dmdtV_ = scalarField("dmdtV", dict, p.size());
            }
            
            if (dict.found("Tl"))
            {
                Tl_ = scalarField("Tl", dict, p.size());
            }
            
            if (dict.found("dDep"))
            {
                dDep_ = scalarField("dDep", dict, p.size());
            }

            if (dict.found("fDep"))
            {
                fDep_ = scalarField("fDep", dict, p.size());
            }
            
            if (dict.found("N"))
            {
                N_ = scalarField("N", dict, p.size());
            }
            
            if (dict.found("K"))
            {
                dict.lookup("K") >> K_;
            }

            if (dict.found("wp"))
            {
                dict.lookup("wp") >> wp_;
            }

            if (dict.found("qQuenching"))
            {
                qq_ = scalarField("qQuenching", dict, p.size());
            }

            break;
        }
    }

    if (dict.found("alphatConv"))
    {
        alphatConv_ = scalarField("alphatConv", dict, p.size());
    }

    AbyV_ = this->patch().magSf();
    forAll(AbyV_, facei)
    {
        const label faceCelli = this->patch().faceCells()[facei];
        AbyV_[facei] /= iF.mesh().V()[faceCelli];
    }
    
    // Generate volumes of cells next to patch
    forAll(V_, facei)
    {
        const label faceCelli = this->patch().faceCells()[facei];
        V_[facei] = iF.mesh().V()[faceCelli];
    }
}


alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906::
alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906
(
    const alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField
    (
        psf,
        p,
        iF,
        mapper
    ),
    otherPhaseName_(psf.otherPhaseName_),
    phaseType_(psf.phaseType_),
    relax_(psf.relax_.clone()),
    AbyV_(psf.AbyV_),
    regimeTypes_(psf.regimeTypes_),
    alphatConv_(psf.alphatConv_, mapper),
    V_(psf.V_, mapper),
    dmdtV_(psf.dmdtV_, mapper),
    Tl_(psf.Tl_, mapper),
    dDep_(psf.dDep_, mapper),
    fDep_(psf.fDep_, mapper),
    N_(psf.N_, mapper),
    qq_(psf.qq_, mapper),
    K_(psf.K_),
    partitioningModel_(psf.partitioningModel_),
    nucleationSiteModel_(psf.nucleationSiteModel_),
    departureDiamModel_(psf.departureDiamModel_),
    filmBoilingModel_(psf.filmBoilingModel_),
    LeidenfrostModel_(psf.LeidenfrostModel_),
    CHFModel_(psf.CHFModel_),
    CHFSoobModel_(psf.CHFSoobModel_),
    MHFModel_(psf.MHFModel_),
    TDNBModel_(psf.TDNBModel_),
    wp_(psf.wp_)
{}


alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906::
alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906
(
    const alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906& psf
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(psf),
    otherPhaseName_(psf.otherPhaseName_),
    phaseType_(psf.phaseType_),
    relax_(psf.relax_.clone()),
    AbyV_(psf.AbyV_),
    regimeTypes_(psf.regimeTypes_),
    alphatConv_(psf.alphatConv_),
    V_(psf.V_),
    dmdtV_(psf.dmdtV_),
    Tl_(psf.Tl_),
    dDep_(psf.dDep_),
    fDep_(psf.fDep_),
    N_(psf.N_),
    qq_(psf.qq_),
    K_(psf.K_),
    partitioningModel_(psf.partitioningModel_),
    nucleationSiteModel_(psf.nucleationSiteModel_),
    departureDiamModel_(psf.departureDiamModel_),
    filmBoilingModel_(psf.filmBoilingModel_),
    LeidenfrostModel_(psf.LeidenfrostModel_),
    CHFModel_(psf.CHFModel_),
    CHFSoobModel_(psf.CHFSoobModel_),
    MHFModel_(psf.MHFModel_),
    TDNBModel_(psf.TDNBModel_),
    wp_(psf.wp_)
{}


alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906::
alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906
(
    const alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField(psf, iF),
    otherPhaseName_(psf.otherPhaseName_),
    phaseType_(psf.phaseType_),
    relax_(psf.relax_.clone()),
    AbyV_(psf.AbyV_),
    regimeTypes_(psf.regimeTypes_),
    alphatConv_(psf.alphatConv_),
    V_(psf.V_),
    dmdtV_(psf.dmdtV_),
    Tl_(psf.Tl_),
    dDep_(psf.dDep_),
    fDep_(psf.fDep_),
    N_(psf.N_),
    qq_(psf.qq_),
    K_(psf.K_),
    partitioningModel_(psf.partitioningModel_),
    nucleationSiteModel_(psf.nucleationSiteModel_),
    departureDiamModel_(psf.departureDiamModel_),
    filmBoilingModel_(psf.filmBoilingModel_),
    LeidenfrostModel_(psf.LeidenfrostModel_),
    CHFModel_(psf.CHFModel_),
    CHFSoobModel_(psf.CHFSoobModel_),
    MHFModel_(psf.MHFModel_),
    TDNBModel_(psf.TDNBModel_),
    wp_(psf.wp_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906::
activePhasePair(const phasePairKey& phasePair) const
{
    if (phasePair == phasePairKey(otherPhaseName_, internalField().group()))
    {
        return true;
    }
    else
    {
        return false;
    }
}

const scalarField& alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906::
dmdt(const phasePairKey& phasePair) const
{
    if (activePhasePair(phasePair))
    {
        return dmdt_;
    }
    else
    {
        FatalErrorInFunction
            << " dmdt requested for invalid phasePair!"
            << abort(FatalError);

        return dmdt_;
    }
}

const scalarField& alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906::
mDotL(const phasePairKey& phasePair) const
{
    if (activePhasePair(phasePair))
    {
        return mDotL_;
    }
    else
    {
        FatalErrorInFunction
            << " mDotL requested for invalid phasePair!"
            << abort(FatalError);

        return mDotL_;
    }
}

void alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906::updateCoeffs()
{

    if (updated())
    {
        return;
    }

    // Check that partitioningModel has been constructed
    if (!partitioningModel_.valid())
    {
        FatalErrorInFunction
            << "partitioningModel has not been constructed!"
            << abort(FatalError);
    }

    // Lookup the fluid model
    const phaseSystem& fluid =
        refCast<const phaseSystem>
        (
            db().lookupObject<phaseSystem>("phaseProperties")
        );

    const saturationModel& satModel =
        db().lookupObject<saturationModel>("saturationModel");

    const label patchi = patch().index();

    const scalar t = this->db().time().timeOutputValue();
    scalar relax = relax_->value(t);

    switch (phaseType_)
    {
        case vaporPhase:
        {
            const phaseModel& vapor
            (
                fluid.phases()[internalField().group()]
            );

            const fvPatchScalarField& hewv =
                vapor.thermo().he().boundaryField()[patchi];

            const phaseModel& liquid(fluid.phases()[otherPhaseName_]);

            const phaseCompressibleTurbulenceModel& turbModel =
                db().lookupObject<phaseCompressibleTurbulenceModel>
                (
                    IOobject::groupName
                    (
                        turbulenceModel::propertiesName,
                        liquid.name()
                    )
                );

            const phaseCompressibleTurbulenceModel& vaporTurbModel =
                db().lookupObject<phaseCompressibleTurbulenceModel>
                (
                    IOobject::groupName
                    (
                        turbulenceModel::propertiesName,
                        vapor.name()
                    )
                );

            const fvPatchScalarField& rhoVaporw =
                vaporTurbModel.rho().boundaryField()[patchi];

            // Vapor Liquid phase fraction at the wall
            const scalarField vaporw(vapor.boundaryField()[patchi]);

            const fvPatchScalarField& Tw =
                liquid.thermo().T().boundaryField()[patchi];
            const scalarField Tc(Tw.patchInternalField());

             // Saturation temperature
            const tmp<volScalarField> tTsat =
                satModel.Tsat(liquid.thermo().p());
            const volScalarField& Tsat = tTsat();
            const fvPatchScalarField& Tsatw(Tsat.boundaryField()[patchi]);
            const scalarField Tsatc(Tsatw.patchInternalField());

            const fvPatchScalarField& hewl =
                liquid.thermo().he().boundaryField()[patchi];

            const fvPatchScalarField& pw =
                liquid.thermo().p().boundaryField()[patchi];

            const fvPatchScalarField& rhow =
                turbModel.rho().boundaryField()[patchi];

            const scalarField hw
            (
                liquid.thermo().he().member() == "e"
              ? hewl.patchInternalField() + pw/rhow.patchInternalField()
              : hewl.patchInternalField()
            );

            const scalarField L
            (
                vapor.thermo().he().member() == "e"
              ? vapor.thermo().he(pw, Tsatc, patchi) + pw/rhoVaporw - hw
              : vapor.thermo().he(pw, Tsatc, patchi) - hw
            );

            // Film boiling models

            scalarField htcFilmBoiling(this->size(), 0);
            scalarField TLeiden(this->size(), GREAT);

            if (filmBoilingModel_.valid() && LeidenfrostModel_.valid())
            {
                // htc for film boiling
                htcFilmBoiling =
                    filmBoilingModel_->htcFilmBoil
                    (
                        liquid,
                        vapor,
                        patchi,
                        Tc,
                        Tsatw,
                        L
                    );

                // Leidenfrost Temperature
                TLeiden =
                    LeidenfrostModel_->TLeid
                    (
                        liquid,
                        vapor,
                        patchi,
                        Tc,
                        Tsatw,
                        L
                    );
            }

            const scalarField qFilm(htcFilmBoiling*max(Tw - Tsatw, scalar(0)));

            // NOTE! Assumes 1-thisPhase for liquid fraction in
            // multiphase simulations
            const scalarField fLiquid
            (
                partitioningModel_->fLiquid(1-vaporw)
            );

            const tmp<scalarField> talphaw = vapor.thermo().alpha(patchi);
            const scalarField& alphaw = talphaw();

            const scalarField heSnGrad(max(hewv.snGrad(), scalar(1e-16)));

            // Convective thermal diffusivity for single phase
            const scalarField alphatv(calcAlphat(*this));

            forAll (*this, i)
            {
                if (Tw[i] > TLeiden[i])
                {
                    this->operator[](i) =
                    (
                        max
                        (
                            (1 - fLiquid[i])
                           *(
                               (qFilm[i]/heSnGrad[i])
                              /max(vaporw[i], scalar(1e-8))
                             - alphaw[i]
                            ),
                            -alphaw[i]
                        )
                    );
                }
                else
                {
                    this->operator[](i) =
                    (
                        (1 - fLiquid[i])*(alphatv[i])
                       /max(vaporw[i], scalar(1e-8))
                    );
                }
            }
        
            /* const phaseModel& vapor */
            /* ( */
            /*     fluid.phases()[internalField().group()] */
            /* ); */
  
            /* const fvPatchScalarField& hewv = */
            /*     vapor.thermo().he().boundaryField()[patchi]; */
               
            /* // Vapor Liquid phase fraction at the wall */
            /* const scalarField vaporw(vapor.boundaryField()[patchi]); */
                          
            /* // NOTE! Assumes 1-thisPhase for liquid fraction in */
            /* // multiphase simulations */
            /* const scalarField fLiquid */
            /* ( */
            /*     partitioningModel_->fLiquid(1-vaporw) */
            /* ); */
                          
            /* const tmp<scalarField> talphaw = vapor.thermo().alpha(patchi); */
            /* const scalarField& alphaw = talphaw(); */
                          
            /* const scalarField heSnGrad(max(hewv.snGrad(), scalar(1e-16))); */
                          
            /* // Convective thermal diffusivity for single phase */
            /* const scalarField alphatv(calcAlphat(*this)); */

            /* forAll(*this, i) */
            /* { */
            /*     this->operator[](i) = */
            /*     ( */
            /*         (1 - fLiquid[i])*(alphatv[i]) */
            /*             /max(vaporw[i], scalar(1e-8)) */
            /*     ); */
            /* } */
            
            if (debug)
            {
                Info<< "alphat for vapour : " << nl << endl;

                Info<< "  alphatEffv: " << gMin(vaporw*(*this + alphaw))
                    << " - " << gMax(vaporw*(*this + alphaw)) << endl;

                const scalarField qEff(vaporw*(*this + alphaw)*hewv.snGrad());

                scalar Qeff = gSum(qEff*patch().magSf());
                Info<< " Effective heat transfer rate to vapor:" << Qeff
                    << nl << endl;
            }
            break;
        }
        case liquidPhase:
        {
            // Check that nucleationSiteModel has been constructed
            if (!nucleationSiteModel_.valid())
            {
                FatalErrorInFunction
                    << "nucleationSiteModel has not been constructed!"
                    << abort(FatalError);
            }

            // Check that departureDiameterModel has been constructed
            if (!departureDiamModel_.valid())
            {
                FatalErrorInFunction
                    << "departureDiameterModel has not been constructed!"
                    << abort(FatalError);
            }

            // Check that nucleationSiteModel has been constructed
            if (!departureFreqModel_.valid())
            {
                FatalErrorInFunction
                    << "departureFrequencyModel has not been constructed!"
                    << abort(FatalError);
            }

            const phaseModel& liquid
            (
                fluid.phases()[internalField().group()]
            );

            const phaseModel& vapor(fluid.phases()[otherPhaseName_]);

            // Retrieve turbulence properties from models
            const phaseCompressibleTurbulenceModel& turbModel =
                db().lookupObject<phaseCompressibleTurbulenceModel>
                (
                    IOobject::groupName
                    (
                        turbulenceModel::propertiesName,
                        liquid.name()
                    )
                );
            const phaseCompressibleTurbulenceModel& vaporTurbModel =
                db().lookupObject<phaseCompressibleTurbulenceModel>
                (
                    IOobject::groupName
                    (
                        turbulenceModel::propertiesName,
                        vapor.name()
                    )
                );

            const tmp<scalarField> tnutw = turbModel.nut(patchi);

            const scalar Cmu25(pow025(Cmu_));

            const scalarField& y = turbModel.y()[patchi];

            const tmp<scalarField> tmuw = turbModel.mu(patchi);
            const scalarField& muw = tmuw();

            const tmp<scalarField> talphaw = liquid.thermo().alphahe(patchi);
            const scalarField& alphaw = talphaw();

            const tmp<volScalarField> tk = turbModel.k();
            const volScalarField& k = tk();
            const fvPatchScalarField& kw = k.boundaryField()[patchi];

            const fvPatchVectorField& Uw =
                turbModel.U().boundaryField()[patchi];
            const scalarField magUp(mag(Uw.patchInternalField() - Uw));
            const scalarField magGradUw(mag(Uw.snGrad()));

            const fvPatchScalarField& rhow =
                turbModel.rho().boundaryField()[patchi];


            const fvPatchScalarField& Tw =
                liquid.thermo().T().boundaryField()[patchi];
            const scalarField Tc(Tw.patchInternalField());

            const scalarField uTau(Cmu25*sqrt(kw));

            const scalarField yPlus(uTau*y/(muw/rhow));

            const scalarField Pr(muw/alphaw);

            // Molecular-to-turbulent Prandtl number ratio
            const scalarField Prat(Pr/Prt_);

            // Thermal sublayer thickness
            const scalarField P(this->Psmooth(Prat));

            const scalarField yPlusTherm(this->yPlusTherm(P, Prat));

            const fvPatchScalarField& rhoVaporw =
                vaporTurbModel.rho().boundaryField()[patchi];

            tmp<volScalarField> tCp = liquid.thermo().Cp();
            const volScalarField& Cp = tCp();
            const fvPatchScalarField& Cpw = Cp.boundaryField()[patchi];

            // Saturation temperature
            const tmp<volScalarField> tTsat =
                satModel.Tsat(liquid.thermo().p());

            const volScalarField& Tsat = tTsat();
            const fvPatchScalarField& Tsatw(Tsat.boundaryField()[patchi]);
            const scalarField Tsatc(Tsatw.patchInternalField());

            const fvPatchScalarField& pw =
                liquid.thermo().p().boundaryField()[patchi];

            const fvPatchScalarField& hew =
                liquid.thermo().he().boundaryField()[patchi];

            const scalarField hw
            (
                liquid.thermo().he().member() == "e"
              ? hew.patchInternalField() + pw/rhow.patchInternalField()
              : hew.patchInternalField()
            );

            const scalarField L
            (
                vapor.thermo().he().member() == "e"
              ? vapor.thermo().he(pw, Tsatc, patchi) + pw/rhoVaporw - hw
              : vapor.thermo().he(pw, Tsatc, patchi) - hw
            );

            // Liquid phase fraction at the wall
            const scalarField liquidw(liquid.boundaryField()[patchi]);

            const scalarField fLiquid(partitioningModel_->fLiquid(liquidw));

            // Liquid temperature at y+=250 is estimated from logarithmic
            // thermal wall function (Koncar, Krepper & Egorov, 2005)
            const scalarField Tplus_y250(Prt_*(log(E_*250)/kappa_ + P));
            const scalarField Tplus(Prt_*(log(E_*yPlus)/kappa_ + P));
            Tl_ = Tw - (Tplus_y250/Tplus)*(Tw - Tc);
            Tl_ = max(Tc - 40, Tl_);

            // Film, transient boiling regimes
            scalarField Qtb(this->size(), 0);
            scalarField tDNB(this->size(), GREAT);
            scalarField TLeiden(this->size(), GREAT);
            scalarField htcFilmBoiling(this->size(), 0);

            if
            (
                CHFModel_.valid()
             && CHFSoobModel_.valid()
             && TDNBModel_.valid()
             && MHFModel_.valid()
             && LeidenfrostModel_.valid()
             && filmBoilingModel_.valid()
            )
            {

                const scalarField CHF
                (
                    CHFModel_->CHF
                    (
                        liquid,
                        vapor,
                        patchi,
                        Tl_,
                        Tsatw,
                        L
                    )
                );

                // Effect of sub-cooling to the CHF in saturated conditions
                const scalarField CHFSubCool
                (
                    CHFSoobModel_->CHFSubCool
                    (
                        liquid,
                        vapor,
                        patchi,
                        Tl_,
                        Tsatw,
                        L
                    )
                );

                const scalarField CHFtotal(CHF*CHFSubCool);

                tDNB =
                    TDNBModel_->TDNB
                    (
                        liquid,
                        vapor,
                        patchi,
                        Tl_,
                        Tsatw,
                        L
                    );

                const scalarField MHF
                (
                    MHFModel_->MHF
                    (
                        liquid,
                        vapor,
                        patchi,
                        Tl_,
                        Tsatw,
                        L
                    )
                );

                TLeiden =
                    LeidenfrostModel_->TLeid
                    (
                        liquid,
                        vapor,
                        patchi,
                        Tl_,
                        Tsatw,
                        L
                    );

                // htc for film boiling
                htcFilmBoiling =
                    filmBoilingModel_->htcFilmBoil
                    (
                        liquid,
                        vapor,
                        patchi,
                        Tl_,
                        Tsatw,
                        L
                    );

                // htc for film transition boiling

                // Indicator between CHF (phi = 0) and MHF (phi = 1)
                const scalarField phi
                (
                    min
                    (
                        max
                        (
                            wp_*(Tw - tDNB)/(TLeiden - tDNB),
                            scalar(0)
                        )
                        , scalar(1)
                    )
                );

                    Qtb = CHFtotal*(1 - phi) + phi*MHF;
            }

            // Bubble departure diameter:
            dDep_ = departureDiamModel_->dDeparture
            (
                liquid,
                vapor,
                patchi,
                Tl_,
                Tsatw,
                L
            );

            // Bubble departure frequency:
            fDep_ = departureFreqModel_->fDeparture
            (
                liquid,
                vapor,
                patchi,
                dDep_
            );

            // Sub-cool boiling Nucleation
            N_ = nucleationSiteModel_->N
            (
                liquid,
                vapor,
                patchi,
                Tl_,
                Tsatw,
                L,
                dDep_,
                fDep_
            );

            // Convective thermal diffusivity for single phase
            alphatConv_ = calcAlphat(alphatConv_);

            // Convective heat transfer area for Sub-cool boiling
            scalarField A1(this->size(), 0);

            const scalarField hewSn(hew.snGrad());

            scalarField alphaFilm(this->size(), 0);

            // Use to identify regimes per face
            /* labelField regimeTypes(A1.size(), -1); */

            forAll (*this, i)
            {
                if (Tw[i] > Tsatw[i])
                {
                    // Sub-cool boiling
                    if (Tw[i] < tDNB[i])
                    {
                        // Sub-cool boiling
                        regimeTypes_[i] = regimeType::subcool;

                        Tl_ = (Tw - (Tplus_y250/Tplus)*(Tw - Tc));
                        Tl_ = max(Tc - 40, Tl_);

                        // Area fractions:

                        /*
                        // Del Valle & Kenning (1985)
                        const scalarField Ja
                        (
                            rhoLiquidw*Cpw*(Tsatw - Tl)/(rhoVaporw*L)
                        );

                        const scalarField Al
                        (
                            fLiquid*4.8*exp(min(-Ja/80, log(VGREAT)))
                        );
                        */

                        // More simple method to calculate area affected by
                        // bubbles
                        const scalar A2
                        (
                            min
                            (
                                fLiquid[i]*pi*sqr(dDep_[i])*N_[i]*K_/4,
                                scalar(1)
                            )
                        );

                        A1[i] = max(1 - A2, 0.0);

                        // Following Bowring(1962)
                        const scalar A2E
                        (
                            min
                            (
                                fLiquid[i]*pi*sqr(dDep_[i])*N_[i],
                                scalar(5)
                            )
                        );

                        // Volumetric mass source in the near wall cell due
                        // to the wall boiling
                        dmdt_[i] =
                            (
                                (1 - relax)*dmdt_[i]
                              + relax*(1.0/6.0)*A2E*dDep_[i]*rhoVaporw[i]
                              * fDep_[i]*AbyV_[i]
                            );

                        // Calculate the mass source multiplied by next
                        // to patch cell volume
                        dmdtV_[i] = dmdt_[i]*V_[i];

                        // Volumetric source in the near wall cell due to
                        // the wall boiling
                        mDotL_[i] = dmdt_[i]*L[i];

                        // Quenching heat transfer coefficient
                        const scalar hQ
                        (
                            2*(alphaw[i]*Cpw[i])*fDep_[i]
                            *sqrt
                            (
                                (0.8/fDep_[i])/(pi*alphaw[i]/rhow[i])
                            )
                        );

                        // Quenching heat flux in Sub-cool boiling
                        qq_[i] =
                            (
                                (1 - relax)*qq_[i]
                              + relax*A2*hQ*max(Tw[i] - Tl_[i], scalar(0))
                            );

                        this->operator[](i) =
                        (
                            max
                            (
                                A1[i]*alphatConv_[i]
                              + (
                                    (qq_[i] + mDotL_[i]/AbyV_[i])
                                  / max(hewSn[i], scalar(1e-16))
                                )
                                /max(liquidw[i], scalar(1e-8)),
                                1e-8
                            )
                        );
                    }
                    else if (Tw[i] > tDNB[i] && Tw[i] < TLeiden[i])
                    {
                        // transient boiling
                        regimeTypes_[i] = regimeType::transient;

                        // No convective heat tranfer
                        alphatConv_[i] = 0.0;

                        // transient boiling
                        dmdt_[i] =
                            fLiquid[i]
                           *(
                                relax*Qtb[i]*AbyV_[i]/L[i]
                              + (1 - relax)*dmdt_[i]
                            );

                        // Calculate the mass source multiplied by next
                        // to patch cell volume
                        dmdtV_[i] = dmdt_[i]*V_[i];
                        
                        mDotL_[i] = dmdt_[i]*L[i];

                        // No quenching flux
                        qq_[i] = 0.0;

                        this->operator[](i) =
                        (
                            max
                            (
                                (
                                    mDotL_[i]/AbyV_[i]
                                    /max(hewSn[i], scalar(1e-16))
                                )/max(liquidw[i], scalar(1e-8)),
                                1e-8
                            )
                        );
                    }
                    else if (Tw[i] > TLeiden[i])
                    {
                        regimeTypes_[i] = regimeType::film; // film boiling

                        // No convective heat tranfer
                        alphatConv_[i] = 0.0;

                        // Film boiling
                        dmdt_[i] =
                            fLiquid[i]
                           *(
                                relax*htcFilmBoiling[i]
                               *max(Tw[i] - Tsatw[i], 0)*AbyV_[i]/L[i]
                             + (1 - relax)*dmdt_[i]
                            );

                        // Calculate the mass source multiplied by next
                        // to patch cell volume
                        dmdtV_[i] = dmdt_[i]*V_[i];

                        mDotL_[i] = dmdt_[i]*L[i];

                        // No quenching flux
                        qq_[i] = 0.0;

                        alphaFilm[i] =
                        (
                            mDotL_[i]/AbyV_[i]/max(hewSn[i], scalar(1e-16))
                        );

                        // alphat is added alphal and multiplied by phase
                        // alphaFilm is lower than alphal. Then we substract
                        // alpha and divide by phase to get alphaFilm
                        this->operator[](i) =
                        (
                            (alphaFilm[i]/max(liquidw[i], scalar(1e-8))
                          - alphaw[i])
                        );
                    }
                }
                else
                {
                    // Tw below Tsat. No boiling single phase convection
                    // Single phase
                    regimeTypes_[i] = regimeType::nonBoiling;
                    A1[i] = 1.0;
                    qq_[i] = 0.0;
                    mDotL_[i] = 0.0;
                    dmdt_[i] = 0.0;

                    // Turbulente thermal diffusivity for single phase.
                    this->operator[](i) =
                    (
                        max
                        (
                            fLiquid[i]*A1[i]*(alphatConv_[i])
                           /max(liquidw[i], scalar(1e-8)),
                            1e-16
                        )
                    );
                }
            }

            const scalarField TsupPrev(Tw - Tsatw);
            Info << "TsupPrev min, max: " << gMin(TsupPrev) <<
                " - " << gMax(TsupPrev) << endl;
                
            Info << "sum(dmdt): " << gSum(dmdt_) << endl;
            Info << "sum(dmdtV): " << gSum(dmdtV_) << endl;
            
            if (debug)
            {
                const scalarField qEff
                (
                    liquidw*(*this + alphaw)*hew.snGrad()
                );

                Info<< "alphat for liquid:  " <<  nl << endl;

                Info<< "  alphatl: " << gMin((*this)) << " - "
                    << gMax((*this)) << endl;

                Info<< "  dmdt: " << gMin((dmdt_)) << " - "
                    << gMax((dmdt_)) << endl;

                Info<< "  alphatlEff: " << gMin(liquidw*(*this + alphaw))
                    << " - " << gMax(liquidw*(*this + alphaw)) << endl;

                scalar Qeff = gSum(qEff*patch().magSf());
                Info<< " Effective heat transfer rate to liquid: " << Qeff
                    << endl << nl;

                if (debug & 2)
                {

                    scalarField nSubCools(this->size(), 0);
                    scalarField nTransients(this->size(), 0);
                    scalarField nFilms(this->size(), 0);
                    scalarField nNonBoilings(this->size(), 0);

                    forAll (*this, i)
                    {
                        switch (regimeTypes_[i])
                        {
                            case regimeType::subcool:
                                nSubCools[i] = 1;
                            break;

                            case regimeType::transient:
                                nTransients[i] = 1;
                            break;

                            case regimeType::film:
                                nFilms[i] = 1;
                            break;

                            case regimeType::nonBoiling:
                                nNonBoilings[i] = 1;
                            break;
                        }
                    }
                    scalar nSubCool(gSum(nSubCools));
                    scalar nTransient(gSum(nTransients));
                    scalar nFilm(gSum(nFilms));
                    scalar nNonBoiling(gSum(nNonBoilings));

                    Info<< "Faces regime :  " <<  nl << endl;

                    Info<< "    sub Cool faces : " << nSubCool << endl;
                    Info<< "    transient faces : " << nTransient << endl;
                    Info<< "    film faces : " << nFilm << endl;
                    Info<< "    non-Boiling faces : " << nNonBoiling << endl;
                    Info<< "    total faces : " << nSubCool + nTransient +
                            nFilm + nNonBoiling << endl << nl;

                    const scalarField qc
                    (
                        nNonBoilings*fLiquid*A1*(alphatConv_ + alphaw)
                       *hew.snGrad()
                    );

                    scalar Qc = gSum(qc*patch().magSf());
                    Info<< " Convective heat transfer: " << Qc << endl;

                    const scalarField qFilm
                    (
                        relax*fLiquid*nFilms*htcFilmBoiling*(Tw - Tsatw)
                    );

                    scalar QFilm = gSum(qFilm*patch().magSf());
                    Info<< " Film boiling heat transfer: " << QFilm << endl;

                    Info<< " Htc Film Boiling coeff: "
                        << gMin(nFilms*htcFilmBoiling)
                        << " - "
                        << gMax(nFilms*htcFilmBoiling) << endl;

                    scalar Qtbtot =
                        gSum(fLiquid*nTransients*Qtb*patch().magSf());
                    Info<< " Transient boiling heat transfer:" << Qtbtot
                        << endl;
                    
                    Info<< " TLeiden: " << gMin(TLeiden) << " - " << gMax(TLeiden)
                        << endl;

                    Info<< " TDNB: " << gMin(tDNB) << " - " << gMax(tDNB)
                        << endl;

                    const scalarField qSubCool
                    (
                        fLiquid*nSubCools*
                        (
                            A1*alphatConv_*hew.snGrad()
                          + qe() + qq()
                        )
                    );

                    scalar QsubCool = gSum(qSubCool*patch().magSf());

                    Info<< " Sub Cool boiling heat transfer: " << QsubCool
                        << endl;

                    Info<< "  N: " << gMin(nSubCools*N_) << " - "
                        << gMax(nSubCools*N_) << endl;
                    Info<< "  dDep: " << gMin(nSubCools*dDep_) << " - "
                        << gMax(nSubCools*dDep_) << endl;
                    Info<< "  fDep: " << gMin(nSubCools*fDep_) << " - "
                        << gMax(nSubCools*fDep_) << endl;
                    Info<< "  A1: " << gMin(nSubCools*A1) << " - "
                        << gMax(nSubCools*A1) << endl;

                    Info<< nl;
                }

            }
        
        }
        break;
        default:
        {
            FatalErrorInFunction
                << "Unknown phase type. Valid types are: "
                << phaseTypeNames_ << nl << exit(FatalError);
        }
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);

    os.writeKeyword("phaseType") << phaseTypeNames_[phaseType_]
        << token::END_STATEMENT << nl;

    relax_->writeData(os);

    switch (phaseType_)
    {
        case vaporPhase:
        {
            os.writeKeyword("partitioningModel") << nl;
            os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
            partitioningModel_->write(os);
            os << decrIndent << indent << token::END_BLOCK << nl;

            if (filmBoilingModel_.valid())
            {
                os.writeKeyword("filmBoilingModel") << nl;
                os << indent << token::BEGIN_BLOCK << incrIndent << nl;
                filmBoilingModel_->write(os);
                os << decrIndent << indent << token::END_BLOCK << nl;
            }

            if (LeidenfrostModel_.valid())
            {
                os.writeKeyword("LeidenfrostModel") << nl;
                os << indent << token::BEGIN_BLOCK << incrIndent << nl;
                LeidenfrostModel_->write(os);
                os << decrIndent << indent << token::END_BLOCK << nl;
            }

            break;
        }
        case liquidPhase:
        {
            os.writeKeyword("partitioningModel") << nl;
            os << indent << token::BEGIN_BLOCK << incrIndent << nl;
            partitioningModel_->write(os);
            os << decrIndent << indent << token::END_BLOCK << nl;

            os.writeKeyword("nucleationSiteModel") << nl;
            os << indent << token::BEGIN_BLOCK << incrIndent << nl;
            nucleationSiteModel_->write(os);
            os << decrIndent << indent << token::END_BLOCK << nl;

            os.writeKeyword("departureDiamModel") << nl;
            os << indent << token::BEGIN_BLOCK << incrIndent << nl;
            departureDiamModel_->write(os);
            os << decrIndent << indent << token::END_BLOCK << nl;

            os.writeKeyword("departureFreqModel") << nl;
            os << indent << token::BEGIN_BLOCK << incrIndent << nl;
            departureFreqModel_->write(os);
            os << decrIndent << indent << token::END_BLOCK << nl;

            if (filmBoilingModel_.valid())
            {
                os.writeKeyword("filmBoilingModel") << nl;
                os << indent << token::BEGIN_BLOCK << incrIndent << nl;
                filmBoilingModel_->write(os);
                os << decrIndent << indent << token::END_BLOCK << nl;
            }

            if (LeidenfrostModel_.valid())
            {
                os.writeKeyword("LeidenfrostModel") << nl;
                os << indent << token::BEGIN_BLOCK << incrIndent << nl;
                LeidenfrostModel_->write(os);
                os << decrIndent << indent << token::END_BLOCK << nl;
            }

            if (CHFModel_.valid())
            {
                os.writeKeyword("CHFModel") << nl;
                os << indent << token::BEGIN_BLOCK << incrIndent << nl;
                CHFModel_->write(os);
                os << decrIndent << indent << token::END_BLOCK << nl;
            }

            if (CHFSoobModel_.valid())
            {
                os.writeKeyword("CHFSubCoolModel") << nl;
                os << indent << token::BEGIN_BLOCK << incrIndent << nl;
                CHFSoobModel_->write(os);
                os << decrIndent << indent << token::END_BLOCK << nl;
            }

            if (MHFModel_.valid())
            {
                os.writeKeyword("MHFModel") << nl;
                os << indent << token::BEGIN_BLOCK << incrIndent << nl;
                MHFModel_->write(os);
                os << decrIndent << indent << token::END_BLOCK << nl;
            }

            if (TDNBModel_.valid())
            {
                os.writeKeyword("TDNBModel") << nl;
                os << indent << token::BEGIN_BLOCK << incrIndent << nl;
                TDNBModel_->write(os);
                os << decrIndent << indent << token::END_BLOCK << nl;
            }

            os.writeKeyword("K") << K_ << token::END_STATEMENT << nl;
            os.writeKeyword("wp") << wp_ << token::END_STATEMENT << nl;
            break;
        }
    }

    os.writeKeyword("otherPhase") << otherPhaseName_
        << token::END_STATEMENT << nl;

    Tl_.writeEntry("Tl", os);
    dmdt_.writeEntry("dmdt", os);
    dmdtV_.writeEntry("dmdtV", os);
    dDep_.writeEntry("dDep", os);
    fDep_.writeEntry("fDep", os);
    N_.writeEntry("N", os);
    qq_.writeEntry("qQuenching", os);
    alphatConv_.writeEntry("alphatConv", os);
    writeEntry("value", os);

    scalarField regimes(this->size(), 0);
    forAll (*this, i)
    {
        switch (regimeTypes_[i])
        {
            case regimeType::nonBoiling:
                regimes[i] = 1;
            break;
            
            case regimeType::subcool:
                regimes[i] = 2;
            break;

            case regimeType::transient:
                regimes[i] = 3;
            break;

            case regimeType::film:
                regimes[i] = 4;
            break;
        }
    }
    regimes.writeEntry("regimeType",os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatWallBoilingWallFunctionFvPatchScalarFieldCustomized_v1906
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
