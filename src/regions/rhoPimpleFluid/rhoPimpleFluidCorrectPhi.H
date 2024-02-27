{
    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("pcorr", p.dimensions(), 0.0),
        pcorrTypes
    );
    
    pcorr == p;

    // Initialise flux with interpolated velocity
    phi_() = fvc::interpolate(rho_())
        *((fvc::interpolate(U_()) & mesh.Sf()) - fvc::meshPhi(rho_(), U));

    adjustPhi(phi_(), psi_(), U_(), rho_(), pcorr);

    mesh().schemesDict().setFluxRequired(pcorr.name());

    while (pimple_.correctNonOrthogonal())
    {        
        fvScalarMatrix pcorrEqn
        (
            fvm::ddt(psi_(), pcorr_()) + fvc::div(phi_()) - fvm::laplacian(rho_()*rAU_(), pcorr) == fvc::div(phi_())
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

#   include "rhoPimpleFluidContinuityErrs.H"
}
