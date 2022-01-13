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
    Foam::functionObjects::yPlusPhase

Group
    grpFieldFunctionObjects

Description
    Evaluates and outputs turbulence y+ for turbulence models
    of multiphase Eulerian-Eulerian solvers.
    Although, it has been developed for multiphase solvers,
    it can be used for single phase as well, just omit the 
    phase option.

    The field is stored on the mesh database so that it can be
    retrieved and used for other applications.

Usage
    Example of function object specification to calculate the y+ (LES):
    \verbatim
    yPlusPhase1
    {
        type        yPlusPhase;
        libs        ("libfieldFunctionObjectsCustomized.so");
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: yPlusPhase   | yes         |
        phase        | Name of the phase       | no          |
        result       | Name of y+ field        | no          | \<function name\>
        log          | Log to standard output  | no          | yes
    \endtable


See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    yPlusPhase.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_yPlusPhase_H
#define functionObjects_yPlusPhase_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                            Class yPlusPhase Declaration
\*---------------------------------------------------------------------------*/

class yPlusPhase
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- No copy construct
        yPlusPhase(const yPlusPhase&) = delete;

        //- No copy assignment
        void operator=(const yPlusPhase&) = delete;

    // Private variables

        const word phase_;

        word fieldName_;
public:

    //- Runtime type information
    TypeName("yPlusPhase");


    // Constructors

        //- Construct from Time and dictionary
        yPlusPhase
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~yPlusPhase();


    // Member Functions

        //- Read the yPlusPhase data
        virtual bool read(const dictionary&);

        //- Calculate the yPlusPhase field
        virtual bool execute();

        //- Write the yPlusPhase field
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
