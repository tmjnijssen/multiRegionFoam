{
    Info<< "Solving chemistry" << endl;

    pChemistry_().solve
    (
        mesh().time().value() - mesh().time().deltaT().value(),
        mesh().time().deltaT().value()
    );

    // turbulent time scale
    if (turbulentReaction_)
    {
        volScalarField tk =
                Cmix_*sqrt(turbulence_().muEff()/rho_()/turbulence_().epsilon());
        volScalarField tc = pChemistry_().tc();

        // Chalmers PaSR model
        chemistryKappa_ = (mesh().time().deltaT() + tc)/(mesh().time().deltaT() + tc + tk);
    }
    else
    {
        chemistryKappa_ = 1.0;
    }

    chemistrySh_ = chemistryKappa_*pChemistry_().Sh()();
}
