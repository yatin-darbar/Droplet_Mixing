/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "kistlerAlphaContactAngleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "fvPatchFields.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::kistlerAlphaContactAngleFvPatchScalarField::
    convertToDeg = 180.0/constant::mathematical::pi;

const Foam::scalar Foam::kistlerAlphaContactAngleFvPatchScalarField::
    convertToRad = constant::mathematical::pi/180.0;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kistlerAlphaContactAngleFvPatchScalarField::
kistlerAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleFvPatchScalarField(p, iF),
    thetaA_(0.0),
    thetaR_(0.0),
    theta0_(0.0),
    muName_("undefined"),
    sigmaName_("undefined")
{}

Foam::kistlerAlphaContactAngleFvPatchScalarField::
kistlerAlphaContactAngleFvPatchScalarField
(
    const kistlerAlphaContactAngleFvPatchScalarField& acpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleFvPatchScalarField(acpsf, p, iF, mapper),
    thetaA_(acpsf.thetaA_),
    thetaR_(acpsf.thetaR_),
    theta0_(acpsf.theta0_),
    muName_(acpsf.muName_),
    sigmaName_(acpsf.sigmaName_)
{}


Foam::kistlerAlphaContactAngleFvPatchScalarField::
kistlerAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleFvPatchScalarField(p, iF),
    thetaA_(readScalar(dict.lookup("thetaA"))),
    thetaR_(readScalar(dict.lookup("thetaR"))),
    theta0_(readScalar(dict.lookup("theta0"))),
    muName_(dict.lookup("muKistler")),
    sigmaName_(dict.lookup("sigmaKistler"))
{
    evaluate();
}

/*
Foam::kistlerAlphaContactAngleFvPatchScalarField::
kistlerAlphaContactAngleFvPatchScalarField
(
    const kistlerAlphaContactAngleFvPatchScalarField& acpsf
)
:
    alphaContactAngleFvPatchScalarField(acpsf),
    thetaA_(acpsf.thetaA_),
    thetaR_(acpsf.thetaR_),
    theta0_(acpsf.theta0_),
    muName_(acpsf.muName_),
    sigmaName_(acpsf.sigmaName_)
{}
*/

Foam::kistlerAlphaContactAngleFvPatchScalarField::
kistlerAlphaContactAngleFvPatchScalarField
(
    const kistlerAlphaContactAngleFvPatchScalarField& acpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleFvPatchScalarField(acpsf, iF),
    thetaA_(acpsf.thetaA_),
    thetaR_(acpsf.thetaR_),
    theta0_(acpsf.theta0_),
    muName_(acpsf.muName_),
    sigmaName_(acpsf.sigmaName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Function is read into alphaContactAngle class to calculate the body force
Foam::tmp<Foam::scalarField>
Foam::kistlerAlphaContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField& Up,
    const fvsPatchVectorField& nHat
) const
{
    // Check patchFields for viscosity and surface tension are ok (add)
    if((muName_ != "muKistler") || (sigmaName_ != "sigmaKistler"))
    {
        FatalErrorIn
        (
            "kistlerAlphaContactAngleFvPatchScalarField"
        )   << " muKistler or sigma set inconsitently, muKistler = "
            << muName_ << ", sigmaKistler = " << sigmaName_ << '.' << nl
            << "    Set both muKistler and sigmaKistler according to the "
            << "definition of kistlerAlphaContactAngle"
            << exit(FatalError);
    }

    // Assign mup as the mixture viscosity, taking the value from muKistler
    // which is defined as muName_ in a constructor above (add)
    const fvPatchField<scalar>& mup =
        patch().lookupPatchField<volScalarField, scalar>(muName_);

    // Assign sigmap as the surface tension, taking the value from sigmaKistler
    // which is defined as sigmaName_ in a constructor above (add)
    const fvPatchField<scalar>& sigmap =
        patch().lookupPatchField<volScalarField, scalar>(sigmaName_);

    // (orig)
    vectorField nf = patch().nf();

    // Calculate the component of the velocity parallel to the wall (orig)
    vectorField Uwall(Up.patchInternalField() - Up);
    Uwall -= (nf & Uwall)*nf;

    // Find the direction of the interface parallel to the wall (orig)
    vectorField nWall(nHat - (nf & nHat)*nf);

    // Normalise nWall (orig)
    nWall /= (mag(nWall) + SMALL);

    // Calculate Uwall resolved normal to the interface parallel to
    // the wall (orig)
    scalarField uwall(nWall & Uwall);

    // Calculate local Capillary number (add)
    scalarField Ca(mup*mag(uwall)/sigmap);

    // Define InverseHoffmanFunction function object (f_iH) (add)
    // Advancing case
    kistlerAlphaContactAngleFvPatchScalarField::InverseHoffmanFunction
    InvHoffFuncThetaA
    (
        convertToRad*thetaA_
    );
    // Receding case
    kistlerAlphaContactAngleFvPatchScalarField::InverseHoffmanFunction
    InvHoffFuncThetaR
    (
        convertToRad*thetaR_
    );

    // Calculate InverseHoffmanFunction values using RiddersRoot (add)
    // Advancing case
    RiddersRoot RRInvHoffFuncThetaA(InvHoffFuncThetaA, 1.e-10);
    scalar InvHoffFuncThetaAroot = RRInvHoffFuncThetaA.root(0,65);
    // Receding case
    RiddersRoot RRInvHoffFuncThetaR(InvHoffFuncThetaR, 1.e-10);
    scalar InvHoffFuncThetaRroot = RRInvHoffFuncThetaR.root(0,65);

    // Calculate and return the value of contact angle on patch faces (add)
    //     Approach: the product of Uwall and nWall is negative for advancing
    //     positive for receding motion. thetaDp is initalised according to
    //     theta0.
    scalarField thetaDp(patch().size(), convertToRad*theta0_);
    forAll(uwall, pfacei)
    {
        if(uwall[pfacei] < 0.0)
        {
            thetaDp[pfacei] = HoffmanFunction(   Ca[pfacei]
                                               + InvHoffFuncThetaAroot);
        }
        else if (uwall[pfacei] > 0.0)
        {
            thetaDp[pfacei] = HoffmanFunction(   Ca[pfacei]
                                               + InvHoffFuncThetaRroot);
        }
    }
    
    return convertToDeg*thetaDp; // return the dynamic CA 
}

// Function returns the value of the Hoffman function for input x
Foam::scalar
Foam::kistlerAlphaContactAngleFvPatchScalarField::HoffmanFunction
(
    const scalar& x
) const
{
    return acos(1 - 2*tanh(5.16*pow(x/(1+1.31*pow(x,0.99)),0.706)));
}

// Function writes out the appropriate values to the OS
void
Foam::kistlerAlphaContactAngleFvPatchScalarField::write(Ostream& os) const
{
    alphaContactAngleFvPatchScalarField::write(os);
    writeEntry(os, "thetaA", thetaA_);
    writeEntry(os, "thetaR", thetaR_);
    writeEntry(os, "theta0", theta0_);
    writeEntry(os, "muKistler", muName_);
    writeEntry(os, "sigmaKistler", sigmaName_);
    writeEntry(os, "value", *this);
}    

/*{
    fvPatchScalarField::write(os);
    os.writeKeyword("thetaA") << thetaA_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaR") << thetaR_ << token::END_STATEMENT << nl;
    os.writeKeyword("theta0") << theta0_ << token::END_STATEMENT << nl;
    os.writeKeyword("muKistler") << muName_ << token::END_STATEMENT << nl;
    os.writeKeyword("sigmaKistler") << sigmaName_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        kistlerAlphaContactAngleFvPatchScalarField
    );
}


// ************************************************************************* //
