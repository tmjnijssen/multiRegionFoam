{
    fvScalarMatrix hsEqn
    (
        fvm::ddt(rho_(), hs_())
      + mvConvection->fvmDiv(phi_(), hs_())
      - fvm::laplacian(turbulence_().alphaEff(), hs_())
//      - fvm::laplacian(turbulence_().muEff(), hs_())  // unit lewis no.
     ==
        DpDt_()
      + chemistrySh_
    );

    hsEqn.relax();
    hsEqn.solve();

    pChemistry_().thermo().correct();

    Info<< "  " << mesh().name() << ": T min/mean/max: " 
      << gMin(T_()) << "/"
      << gAverage(T_()) << "/"
      << gMax(T_()) << endl;
}
