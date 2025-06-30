{
    pcorr_() *= 0;
    // [MP] Phi formulation taken from interTrackFluid
    if (mesh().changing())
    {
        forAll(U_().boundaryField(), patchi)
        {
            if (U_().boundaryField()[patchi].fixesValue())
            {
                U_().boundaryField()[patchi].initEvaluate();
            }
        }

        forAll(U_().boundaryField(), patchi)
        {
            if (U_().boundaryField()[patchi].fixesValue())
            {
                // Fluxes updated in the BC
                U_().boundaryField()[patchi].evaluate();
                phi_().boundaryField()[patchi] =
                    U_().boundaryField()[patchi] 
                    & mesh().Sf().boundaryField()[patchi];
            }
        }
    }

    Info << "U BC Pimple PrePredictor ";
    forAll(U_().boundaryField(), patchi)
    {
        Info << " " << mesh().boundaryMesh()[patchi].name() << " " 
            << gSum(mesh().time().deltaT().value()*((U_().boundaryField()[patchi] & mesh().boundary()[patchi].Sf()))) << ", ";
    }
    Info << endl;

    // Initialise flux with interpolated velocity
    // [MP] Old Phi formulation taken 

    //phi_() = fvc::interpolate(U_()) & mesh().Sf();

    if (closedVolume_ && pKin_().needReference())
    {
        if (mesh().moving() && mesh().time().timeIndex() > 10)
        {
            fvc::makeRelative(phi_(), U_());
            forAll (phi_().boundaryField(), patchi)
            {
                const fvPatchVectorField& Up = U_().boundaryField()[patchi];
                if
                (
                    isA<regionCoupledMassTransferVelocityValue>(Up)
                )
                {
                    scalar dV = sum(fvc::meshPhi(U_())().boundaryField()[patchi]);

                    scalar correction = -dV/sum(phi_().boundaryField()[patchi]);

                    phi_().boundaryField()[patchi] *= correction;
                }
            }
            fvc::makeAbsolute(phi_(), U_());
        }
    }
    else
    {
        //adjustPhi(phi_(), U_(), pcorr_());
    }

    mesh().schemesDict().setFluxRequired(pcorr_().name());

    while (pimple_.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAU_(), pcorr_()) == fvc::div(phi_())
        );

        pcorrEqn.setReference(pRefCell_, pRefValue_);
        pcorrEqn.solve();

        if (pimple_.finalNonOrthogonalIter())
        {
            phi_() -= pcorrEqn.flux();
        }

        // Fluxes are corrected to absolute velocity and further corrected
        // later.  HJ, 6/Feb/2009
    }

#   include "pimpleFluidContinuityErrs.H"
}
