
    fvScalarMatrix TEqn
    (
      fvm::div(phi, T)              // 被运         // passive transport
      + fvm::div(phiWf, T)            // 以速度w下落    // fall down with velocity w_f
      - fvm::laplacian(turbulence->nut()/S_ct, T)   // 湍散         // turbulent diffusion
      ==
        fvOptions(T)
    );

    TEqn.relax();
    fvOptions.constrain(TEqn);
    solverPerformance TInfo = TEqn.solve();
    TResidual = TInfo.initialResidual();
    fvOptions.correct(T);
/*
    forAll(T, i)
    {
      if(T[i] < 0)
      {
        T[i] = 0;
      }
    }
*/
