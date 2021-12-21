/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
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

Class
    Foam::functionObjects::TF

Group
    grpFieldFunctionObjects

Description
    Evaluates and outputs Temperature of fluid in boundary cell
    base on conjugate heat transfer into two phases. The code is tested
    for conjugate boiling using multiphase Eulerian-Eulerian solvers.
    Although, it has been developed for multiphase solvers,
    it can be used for single phase as well, just omit the 
    phase option.

    The filed is written to the particular boundary mesh instead
    of the boundary cells.

Usage
    Example of function object specification to calculate the TF:
    \verbatim
    Tfluid
    {
        type        TFluid;
        libs        ("libfieldFunctionObjectsCustomized.so");
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: TFluid   | yes         |
        region       | fluid
        phase        | Name of the phase       | no          | None
        otherPhase   | Name of the other phase | no          | None
        Tnbr         | Name of other region T field | no     | T
        result       | Name of TF field        | no          | \<function name\>
        log          | Log to standard output  | no          | yes
    \endtable


See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    TFluid.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_TF_H
#define functionObjects_TF_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                            Class TF Declaration
\*---------------------------------------------------------------------------*/

class TFluid
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- No copy construct
        TFluid(const TFluid&) = delete;

        //- No copy assignment
        void operator=(const TFluid&) = delete;

    // Private variables

        const word phaseName_;

        const word otherPhaseName_;

        const word TnbrName_;

        word Tliquid_;

        word Tgas_;
public:

    //- Runtime type information
    TypeName("TFluid");


    // Constructors

        //- Construct from Time and dictionary
        TFluid
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~TFluid();


    // Member Functions

        //- Read the TFluid data
        virtual bool read(const dictionary&);

        //- Calculate the TFluid field
        virtual bool execute();

        //- Write the TFluid field
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
