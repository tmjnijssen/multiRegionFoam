PtrList<volScalarField>& Y = pChemistry_().thermo().composition().Y();

tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh(),
        fields_,
        phi_(),
        mesh().schemesDict().divScheme("div(phi,Yi_h)")
    )
);

{
    label inertIndex = -1;
    volScalarField Yt("Yt", 0.0*Y[0]);

    for (label i = 0; i < Y.size(); i++)
    {
        if (Y[i].name() != inertSpecie_)
        {
            volScalarField& Yi = Y[i];

            solve
            (
                fvm::ddt(rho_(), Yi)
              + mvConvection->fvmDiv(phi_(), Yi)
              - fvm::laplacian(turbulence_().muEff(), Yi)
             ==
                chemistryKappa_*pChemistry_().RR(i),
                mesh().solutionDict().solver("Yi")
            );

            Yi.max(0.0);
            Yt += Yi;
        }
        else
        {
            inertIndex = i;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
