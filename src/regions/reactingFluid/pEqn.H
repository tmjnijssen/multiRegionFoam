rho_() = pChemistry_().thermo().rho();

volScalarField rUA = 1.0/UEqn.A();
U_() = rUA*UEqn.H();

if (pimple_.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi_())
       *(
            (fvc::interpolate(U_()) & mesh().Sf())
          + fvc::ddtPhiCorr(rUA, rho_(), U_(), phi_())
        )
    );

    while (pimple_.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi_(), p_())
          + fvm::div(phid, p_())
          - fvm::laplacian(rho_()*rUA, p_())
        );

        pEqn.solve();

        if (pimple_.finalNonOrthogonalIter())
        {
            phi_() == pEqn.flux();
        }
    }
}
else
{
    phi_() =
        fvc::interpolate(rho_())
       *(
            (fvc::interpolate(U_()) & mesh().Sf())
          + fvc::ddtPhiCorr(rUA, rho_(), U_(), phi_())
        );

    while (pimple_.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi_(), p_())
          + fvc::div(phi_())
          - fvm::laplacian(rho_()*rUA, p_())
        );

        pEqn.solve();

        if (pimple_.finalNonOrthogonalIter())
        {
            phi_() += pEqn.flux();
        }
    }
}

#include "rhoEqn.H"
#include "reactingFluidContinuityErrs.H"

U_() -= rUA*fvc::grad(p_());
U_().correctBoundaryConditions();

DpDt_() = fvc::DDt(surfaceScalarField("phiU", phi_()/fvc::interpolate(rho_())), p_());
