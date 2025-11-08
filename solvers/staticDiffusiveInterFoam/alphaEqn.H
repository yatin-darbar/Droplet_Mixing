{
    #include "alphaScheme.H"
    // Check subcycling alphaEqn is not being attempted
    if (nAlphaSubCycles > 1)
    {
	FatalErrorInFunction
	    << "Sub-cycling is not supported with ScalarInterFoam. Either use"
	       "an original solver or set nAlphaSubCycles=1 in"
               "system/fvSolution"
	       << exit(FatalError);

    }

    // Check semi-implict MULES has not been requested
    if (MULESCorr)
    {
        FatalErrorInFunction
            << "Semi-implicit variant of MULES is not supported with "
               "ScalarInterFoam. Either set MULESCorr=FALSE or remove the "
               "reference to MULESCorr from system/fvSolution"
                << exit(FatalError);
    }

    // Set the off-centering coefficient according to ddt scheme
    scalar ocCoeff = 0;
    {
        tmp<fv::ddtScheme<scalar>> tddtAlpha
        (
            fv::ddtScheme<scalar>::New
            (
                mesh,
                mesh.ddtScheme("ddt(alpha)")
            )
        );
        const fv::ddtScheme<scalar>& ddtAlpha = tddtAlpha();

        if
        (
            isType<fv::EulerDdtScheme<scalar>>(ddtAlpha)         
        )
        {
            ocCoeff = 0;
        }
        else if (isType<fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha))
        {
            if
            (
                alphaRestart
             || mesh.time().timeIndex() > mesh.time().startTimeIndex() + 1
            )
            {
                ocCoeff =
                    refCast<const fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha)
                   .ocCoeff();
            }
        }
        else
        {
            FatalErrorInFunction
                << "Only Euler and CrankNicolson ddt schemes are supported"
                << exit(FatalError);
        }
    }

    // Set the time blending factor, 1 for Euler
    scalar cnCoeff = 1.0/(1.0 + ocCoeff);

    tmp<surfaceScalarField> phiCN(phi);

    // Calculate the Crank-Nicolson off-centred volumetric flux
    if (ocCoeff > 0)
    {
        phiCN = surfaceScalarField::New
        (
            "phiCN",
            cnCoeff*phi + (1.0 - cnCoeff)*phi.oldTime()
        );
    }

    #include "alphaSuSp.H"  

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        // Split operator
        tmp<surfaceScalarField> talphaPhi1Un
        (
            fvc::flux
            (
                phiCN(),
                (cnCoeff*alpha1 + (1.0 - cnCoeff)*alpha1.oldTime())(),
                compressionScheme.rewind()
            )
        );
   
        alphaPhi10 = talphaPhi1Un;

        if (divU.valid())
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phiCN,
                alphaPhi10,
                Sp(),
                (Su() + divU()*min(alpha1(), scalar(1)))(),
                oneField(),
                zeroField()
            );
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phiCN,
                alphaPhi10,
                oneField(),
                zeroField()
            );
        }

        alpha2 = 1.0 - alpha1;

        mixture.correct();
    }

    talphaPhi1Corr0.clear();
 
    #include "rhofs.H"
    // Calculate end-of-time-step fluxes for use in UEqn.H
    if
    (
        word(mesh.ddtScheme("ddt(rho,U)"))
     == fv::EulerDdtScheme<vector>::typeName
    )
    {
        rhoPhi = alphaPhi10*(rho1f - rho2f) + phiCN*rho2f;
    }
    else
    {
        if (ocCoeff > 0)
        {
            // Calculate the end-of-time-step alpha flux
            alphaPhi10 =
                (alphaPhi10 - (1.0 - cnCoeff)*alphaPhi10.oldTime())/cnCoeff;
        }

        // Calculate the end-of-time-step mass flux
        rhoPhi = alphaPhi10*(rho1f - rho2f) + phi*rho2f;
    }

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;

    rho == alpha1*rho1 + alpha2*rho2;  // Originally in alphaEqnSubCycle.H
}
