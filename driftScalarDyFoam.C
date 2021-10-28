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
    driftScalarDyFoam

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
#include "dynamicFvMesh.H"
#include "SolverPerformance.H"

#include "primitivePatchInterpolation.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "passive transport solver for snow drifting."
    );
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    //Foam::Info<< "Create runStage\n" << Foam::endl;

    //Foam::Time runStage(Foam::Time::controlDictName, args);
    //#include "createMesh.H"
    #include "createNamedDynamicFvMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "findSnowPatch.H"

    const scalar rhoAir = readScalar(erosionDepositionProperties.lookup("rhoAir"));
    const scalar rhoSnow = readScalar(erosionDepositionProperties.lookup("rhoSnow"));
    const scalar ca = readScalar(erosionDepositionProperties.lookup("ca"));
    const scalar Uthreshold = readScalar(erosionDepositionProperties.lookup("Uthreshold"));

    const scalar residualT = readScalar(erosionDepositionProperties.lookup("residualT"));
    const label nSubCycles = readLabel(erosionDepositionProperties.lookup("nSubCycles"));

    label nStage = 0;

    while (runTime.run())
    {
        Info << nl << "-----------------------" << endl;
        Info << "Stage = " << runTime.timeIndex() << endl;
        Info << "Physical Time = " << runTime.timeName() << endl;
        Info << "-----------------------" << endl;

        TimeState subCycleTimeState = runTime.subCycle(nSubCycles);

        // 启动单个阶段内的子循环迭代
        for (label cycleI = 0; cycleI < nSubCycles; cycleI++)
        {
            Info<< "\n\nsubCycles = " << runTime.timeName() << nl << endl;
            runTime++;
            // --- Pressure-velocity SIMPLE corrector
            {
                #include "UEqn.H"
                #include "pEqn.H"
            }

            // RANS
            laminarTransport.correct();
            turbulence->correct();

            // 标量输运
            // scalar transport
            solverPerformance Tres;
            while (simple.correctNonOrthogonal())
            {
                #include "TEqn.H"
            }
            
            // 由于扩散需要获取湍流粘度nut,故检查湍流模型
            if (!turbulence)
            {
                FatalErrorInFunction
                << "Unable to find turbulence model in the "
                << "database" << exit(FatalError);
                return 1;
            }

            volSymmTensorField Reff = turbulence->devReff();

            // 修正质量交换率M
            vector zNormal = vector(0.,0.,-1.);
            for (label patchi : snowPatchList)
            {
                // 获取雪面的剪切应力
                // get shear stress on snow surface
                const vectorField &Sfp = mesh.Sf().boundaryField()[patchi];         // 面积矢量
                const scalarField &magSfp = mesh.magSf().boundaryField()[patchi];   // 面积矢量模长
                const scalarField &Tp = T.boundaryField()[patchi];                  // 边界处的雪漂浓度

                const symmTensorField& Reffp = Reff.boundaryField()[patchi];        

                vectorField& ssp = wallShearStress.boundaryFieldRef()[patchi];

                ssp = (-Sfp/magSfp) & Reffp;    // 剪切应力 

                
                scalarField UShear = sqrt(mag(ssp)/rhoAir);         // 剪切风速模量

                scalarField& Mp = M.boundaryFieldRef()[patchi];

                // 更新质量交换率（侵蚀/沉积）
                // update mass exchange rate on snow surface (erosion & deposition)
                Info << "Update mass exchange rate." << endl; 
                forAll(Mp, i)
                {
                    const scalar zArea = Sfp[i] & zNormal;
                    if (UShear[i] > Uthreshold) // 侵蚀
                    {
                        //Mp[i] = -ca * rhoSnow * UShear[i] * (1.-sqr(Uthreshold)/sqr(UShear[i])) * zArea;
                        Mp[i] = -ca * rhoSnow * (sqr(UShear[i]) - sqr(Uthreshold)) * zArea;
                    }
                    else // 沉积
                    {
                        Mp[i] = Tp[i] * mag(wf.value()) * zArea;
                    }

                    //tmpDeltaH[i] = Mp[i] / rhoSnow / zArea * runTime.deltaTValue();
                }
            }

            runTime.printExecutionTime(Info);   // 打印子循环的运行时

            // 检查是否达到收敛
            if (Tres.initialResidual() < residualT && nStage > 0)
            {   
                Info << "Subcycle converged." << endl;
                break; // 跳出子循环
            }
        }
        runTime.endSubCycle(); //结束子循环

        // 计算根据M计算deltaH
        vector zNormal = vector(0.,0.,-1.);
        for (label patchi : snowPatchList)
        {
            scalarField &deltaHp = deltaH.boundaryFieldRef()[patchi];

            const scalarField &Mp = M.boundaryField()[patchi];
            const vectorField &Sfp = mesh.Sf().boundaryField()[patchi]; 

            forAll(Mp, i)
            {
                const scalar zArea = Sfp[i] & zNormal;
                deltaHp[i] = Mp[i] / rhoSnow / zArea;
            }
        }
        
        runTime.write();
        ++runTime;

        if (!runTime.run())
        {
            break;
        }
        
        Info << "\nUpdate Mesh" << nl << endl;
        mesh.update();

        ++nStage;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
