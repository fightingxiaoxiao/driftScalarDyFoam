/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | https://github.com/fightingxiaoxiao
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 Chenxiaoxiao
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    driftScalarUnsteadyDyFoam

Group
    extendSolvers

Description
    Passive transport solver for snow drifting.

Solver details
    The equation is given by:

    Where:
        T       | Passive scalar

    Required fields
        T       | Passive scalar
        U       | Velocity [m/s]
        nut     | Turbulent viscosity [m^2/s]


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Passive scalar transport equation solver."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
        
    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    #include "findSnowPatch.H"

    const scalar rhoAir = readScalar(erosionDepositionProperties.lookup("rhoAir"));
    const scalar rhoSnow = readScalar(erosionDepositionProperties.lookup("rhoSnow"));
    const scalar ca = readScalar(erosionDepositionProperties.lookup("ca"));
    const scalar Uthreshold = readScalar(erosionDepositionProperties.lookup("Uthreshold"));

    volSymmTensorField Reff = turbulence->devReff();
    

    // 修正质量交换率M
    vector zNormal = vector(0., 0., -1.);
    for (label patchi : snowPatchList)
    {
        // 获取雪面的剪切应力
        // get shear stress on snow surface
        const vectorField &Sfp = mesh.Sf().boundaryField()[patchi]; // 面积矢量
        const scalarField &magSfp = mesh.magSf().boundaryField()[patchi];   // 面积矢量模长
        const scalarField &Tp = T.boundaryField()[patchi]; // 边界处的雪漂浓度
        const symmTensorField& Reffp = Reff.boundaryField()[patchi];

        vectorField &ssp = wallShearStress.boundaryFieldRef()[patchi]; // 剪切应力
        

        ssp = (-Sfp/magSfp) & Reffp;    // 剪切应力

        scalarField UShear = sqrt(mag(ssp) / rhoAir); // 剪切风速模量
        scalarField &Mp = M.boundaryFieldRef()[patchi];
        scalarField &deltaHp = deltaH.boundaryFieldRef()[patchi];

        // 更新质量交换率（侵蚀/沉积）
        // update mass exchange rate on snow surface (erosion & deposition)
        //Info << "Update mass exchange rate." << endl;
        forAll(Mp, i)
        {
            const scalar zArea = Sfp[i] & zNormal;
            if (UShear[i] > Uthreshold) // 侵蚀
            {
                //Mp[i] = -ca * rhoSnow * UShear[i] * (1.-sqr(Uthreshold)/sqr(UShear[i])) *zArea;
                Mp[i] = -ca * rhoSnow * (sqr(UShear[i]) - sqr(Uthreshold)) * zArea;
            }
            else // 沉积
            {
                Mp[i] = Tp[i] * mag(wf.value()) * (1 - sqr(UShear[i]) / sqr(Uthreshold)) * zArea;
            }
            //tmpDeltaH[i] = Mp[i] / rhoSnow / zArea * runTime.deltaTValue();
            deltaHp[i] = Mp[i] / rhoSnow / zArea;
        }
    }

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T)
                + fvm::div(phi, T)                              // 被动输运         // passive transport
                + fvm::div(phiWf, T)                            // 以速度wf下落      // fall down with velocity w_f
                - fvm::laplacian(turbulence->nut()/S_ct, T)     // 湍散             // turbulent diffusion
                ==
                fvOptions(T)
            );

            //TEqn.relax();
            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }
        
        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
