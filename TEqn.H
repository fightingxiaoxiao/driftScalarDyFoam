
    fvScalarMatrix TEqn
    (
//      fvm::ddt(T)
        fvm::div(phi, T)        // 被动输运  // passive transport
    +   fvm::div(phiWf, T)      // 以速度w_f下落 // fall down with velocity w_f
    -   fvm::laplacian(nut, T)  // 湍流扩散  // turbulent diffusion
    ==
        fvOptions(T)
    );

    TEqn.relax();
    fvOptions.constrain(TEqn);
    TEqn.solve();
    fvOptions.correct(T);
