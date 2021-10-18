
    fvScalarMatrix TEqn
    (
      fvm::ddt(T)
    + fvm::div(phi, T)              // 被运         // passive transport
    + fvm::div(phiWf, T)            // 以速度w下落    // fall down with velocity w_f
    - fvm::laplacian(nut/S_ct, T)   // 湍散         // turbulent diffusion
    ==
      fvOptions(T)
    );
    TEqn.relax();
    fvOptions.constrain(TEqn);
    TEqn.solve();
    fvOptions.correct(T);
