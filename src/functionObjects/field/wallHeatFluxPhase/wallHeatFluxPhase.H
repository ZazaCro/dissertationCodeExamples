/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

Class
    Foam::functionObjects::wallHeatFluxPhase

Group
    grpForcesFunctionObjects

Description
    Calculates and write the heat-flux at wall patches as the
    volScalarField field 'wallHeatFluxPhase'.

    All wall patches are included by default; to restrict the calculation to
    certain patches, use the optional 'patches' entry.

Usage
    Example of function object specification:
    \verbatim
    wallHeatFluxPhase1
    {
        type        wallHeatFluxPhase;
        libs        ("libfieldFunctionObjectsCustomized");
        ...
        region      fluid;
        patches     (".*Wall");
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property | Description                | Required   | Default value
        type     | type name: wallHeatFluxPhase    | yes        |
        phase    | Name of the phase          | no         |
        patches  | list of patches to process | no         | all wall patches
        qr       | name of radiative heat flux field | no  | qr
        region   | region to be evaluated     | no         | default region
    \endtable

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    wallHeatFluxPhase.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_wallHeatFluxPhase_H
#define functionObjects_wallHeatFluxPhase_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class wallHeatFluxPhase Declaration
\*---------------------------------------------------------------------------*/

class wallHeatFluxPhase
:
    public fvMeshFunctionObject,
    public writeFile
{

protected:

    // Protected data

        //- Optional list of patches to process
        labelHashSet patchSet_;

        //- Name of radiative heat flux name, default = qr
        word qrName_;

        //- Name of the phase 
        const word phase_;

        word fieldName_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Calculate the heat-flux
        void calcHeatFlux
        (
            const volScalarField& alpha,
            const volScalarField& he,
            volScalarField& wallHeatFluxPhase
        );

        //- Calculate phase intensive heat-flux
        void calcHeatFluxPhase
        (
            const volScalarField& alpha,
            const volScalarField& he,
            const volScalarField& phaseFraction,
            volScalarField& wallHeatFluxPhase
        );

private:

    // Private member functions

        //- No copy construct
        wallHeatFluxPhase(const wallHeatFluxPhase&) = delete;

        //- No copy assignment
        void operator=(const wallHeatFluxPhase&) = delete;


public:

    //- Runtime type information
    TypeName("wallHeatFluxPhase");


    // Constructors

        //- Construct from Time and dictionary
        wallHeatFluxPhase
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~wallHeatFluxPhase();


    // Member Functions

        //- Read the wallHeatFluxPhase data
        virtual bool read(const dictionary&);

        //- Calculate the wall heat-flux
        virtual bool execute();

        //- Write the wall heat-flux
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
