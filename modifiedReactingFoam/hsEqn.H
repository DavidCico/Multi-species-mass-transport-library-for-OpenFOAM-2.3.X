
{
    fvScalarMatrix hsEqn
            (
                fvm::ddt(rho, hs)
                + mvConvection->fvmDiv(phi, hs)
                - fvm::laplacian(turbulence->alphaEff(), hs)
                + mstm().multiSpeciesHeatSource()
                ==
                dpdt
                // ------------------------------------------------------------------------- //
                + chemistrySh
                // ------------------------------------------------------------------------- //
                );

    hsEqn.relax();
    hsEqn.solve();

    thermo.correct();

    // A.Alexiou 2015 - Debugging
    //Info << "Check T" << thermo.T();

    Info<< "T gas min/max   = " << min(T).value() << ", "
        << max(T).value() << endl;
}
