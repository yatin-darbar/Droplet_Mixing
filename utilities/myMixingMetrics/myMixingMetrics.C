/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 OpenFOAM Foundation
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

#include "myMixingMetrics.H"
#include "momentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(myMixingMetrics, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        myMixingMetrics,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::myMixingMetrics::writeFileHeader(const label i)
{
    writeHeader(file(), "Mixing Metrics");
    writeCommented(file(), "Time");
    writeTabbed(file(), "MIndex");
    writeTabbed(file(), "Pe_D");
    file() << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::myMixingMetrics::myMixingMetrics
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::myMixingMetrics::~myMixingMetrics()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::myMixingMetrics::read(const dictionary& dict)
{

    fvMeshFunctionObject::read(dict);
    //writeLocalObjects::read(dict);
    dict.readIfPresent("phase", phaseName_);
    word alphastr("alpha.");
    alphaName_ = alphastr + phaseName_;
    
    dict.readIfPresent("scalarName", scalarName_);
    scalarName_ = scalarName_ + "." + phaseName_;

    dict.readIfPresent("L", L_);
    dict.readIfPresent("D", D_);
    
    resetName(typeName);

    return true;
}


bool Foam::functionObjects::myMixingMetrics::execute()
{
    return true;
}


bool Foam::functionObjects::myMixingMetrics::write()
{
    
    //- Initial Stuff to create wrtie files
    Log << type() << " " << name() << " write:" << nl;

    //writeLocalObjects::write();

    logFiles::write();

    //- Calculation of Mixing quantities

    const volScalarField& alpha = mesh_.lookupObject<volScalarField>(alphaName_);
//    const volScalarField& beta = mesh_.lookupObject<volScalarField>("beta");
    const volScalarField& gamma = mesh_.lookupObject<volScalarField>(scalarName_);
    const volVectorField&     U = mesh_.lookupObject<volVectorField>("U");

    scalar MI1_ = 0;   // Standard Dev from mean of gamma
    scalar max_U = 0;

    scalar MIchecker = 0;

    scalar interfaceCells_ = 0;
    scalar alphaCells = 0;
//    scalar betaCells = 0;
    scalar gammaCells = 0;
    int cell = 999;

    scalar gammaBar = 0;

    forAll(mesh_.C(), celli)
    {
	//-- Interface Loop
        if (alpha[celli] > 0.1 && alpha[celli] < 0.9) 	       
        {
	    interfaceCells_ += 1;
	}

	//-- Inside Droplet Loop
        if (alpha[celli] > 0.99999)
        {
	    //-- No. cells contain droplet phase
	    alphaCells += 1;
            
//            if (beta[celli] > 0.5)
//            {
//                betaCells += 1;
//            }
            if (gamma[celli] > 0.5)
            {
                gammaCells += 1;
            }

	    if (gamma[celli] > 0.471 && gamma[celli] < 0.529)
	    {
		MIchecker += 1;
	    }

	    //-- Mixing Metric Calculation
	    gammaBar += gamma[celli];


	    if (mag(U[celli]) > max_U)
	        {
		    max_U = mag(U[celli]);
	   	    cell = celli;
                }

        }
	
    }

Log << "got here 2 " << nl;
    //--Ensure operations done over all processors    
    reduce(interfaceCells_, sumOp<scalar>());
    reduce(alphaCells, sumOp<scalar>());
//    reduce(betaCells, sumOp<scalar>());
    reduce(gammaCells, sumOp<scalar>());

    reduce(MIchecker, sumOp<scalar>());

    reduce(MI1_, sumOp<scalar>());
    reduce(max_U, maxOp<scalar>());
    
    reduce(gammaBar, sumOp<scalar>());
 
Log << "got here 3 " << nl;
 
   gammaBar = gammaBar / alphaCells;

    forAll(mesh_.C(), celli)
        {
        //-- Inside Droplet Loop
            if (alpha[celli] > 0.99999)
            {
	        MI1_ += pow((gamma[celli] - gammaBar),2);
            }
	}
	
    reduce(MI1_, sumOp<scalar>());
    MI1_ = pow(MI1_/alphaCells,0.5);

Log << "got here 4 " << nl;

//Log << "got here MI2 " << nl;
    PeD_ = (max_U*L_)/D_;
//Log << "got here Dym Pe " << nl;

    // Print out to I/O stream the values
    Log << "No. alpha Cells = " << alphaCells << nl;
//    Log << "No. beta Cells = " << betaCells << nl;
    Log << "No. gamma Cells = " << gammaCells << nl;

    Log << "No. MIchecker Cells = " << MIchecker << nl;

    Log << "MI1 = " << MI1_ << nl;	
    Log << "PeD = " << PeD_ << nl;
    Log << "U_max = " << max_U << nl;
    Log << "Cell = " << cell << nl;

    Log << "% of cells in mixing region " << MIchecker / alphaCells << nl;
    if (Pstream::master())
    {
	/*
        Log << "    I : min = " << " 5 "
            << ", max = " << " hello " 
            << ", average = " << " bye " << nl;
	*/
        writeTime(file());
        file()
            << tab << MI1_
            << tab << PeD_
            << endl;

    }

    Log << endl;

    return true;
}
// ************************************************************************* //
