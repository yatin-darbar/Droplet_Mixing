/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2020 OpenFOAM Foundation
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

Class
    Foam::functionObjects::myMixingMetrics

Description
    Evaluates and writes the turbulence intensity field 'I'.

    The turbulence intensity field 'I' is the root-mean-square of the turbulent
    velocity fluctuations normalised by the local velocity magnitude:
    \f[
        I \equiv \frac{\sqrt{\frac{2}{3}\, k}}{U}
    \f]
    To avoid spurious extrema and division by 0 I is limited to 1 where the
    velocity magnitude is less than the turbulent velocity fluctuations.

    Example of function object specification:
    \verbatim
    functions
    {
        .
        .
        .
        myMixingMetrics
        {
            type        myMixingMetrics;
            libs        ("libfieldFunctionObjects.so");
        }
        .
        .
        .
    }
    \endverbatim

    or using the standard configuration file:
    \verbatim
    functions
    {
        .
        .
        .
        #includeFunc myMixingMetrics
        .
        .
        .
    }
    \endverbatim

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::writeLocalObjects
    Foam::functionObjects::timeControl

SourceFiles
    myMixingMetrics.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_myMixingMetrics_H
#define functionObjects_myMixingMetrics_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "writeLocalObjects.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class myMixingMetrics Declaration
\*---------------------------------------------------------------------------*/

class myMixingMetrics
:
    public fvMeshFunctionObject,
    public logFiles,
    public writeLocalObjects
{
protected:

    // Protected Data

	//-- Name of phase
	word phaseName_;

	//-- Full phase string
	word alphaName_;

	//-- Name of diffuse scalar
	word scalarName_;

        //-- MI - for diffusive scalar
        scalar MI1_;   

        //-- Pe_D - dynamic Peclet number based on max droplet velocity
        scalar PeD_;

	//-- L - Length Scale of the problem (droplet radius)
	scalar L_;

	//-- D - diffusion coefficient
	scalar D_;

	//-- No. of interface cells
	scalar interfaceCells_;

	
    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);


public:

    //- Runtime type information
    TypeName("myMixingMetrics");


    // Constructors

        //- Construct from Time and dictionary
        myMixingMetrics
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        myMixingMetrics(const myMixingMetrics&) = delete;


    //- Destructor
    virtual ~myMixingMetrics();


    // Member Functions

        //- Read the myMixingMetrics data
        virtual bool read(const dictionary&);

        //- Calculate the myMixingMetrics field
        virtual bool execute();

        //- Write the myMixingMetrics field
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const myMixingMetrics&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
