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

    const label nSubCycles = readLabel(erosionDepositionProperties.lookup("nSubCycles"));

    label nStage = 0;

    while (runTime.run())
    {
        Info << nl << "-----------------------" << endl;
        Info << "Stage = " << runTime.timeIndex() << endl;
        Info << "Physical Time = " << runTime.timeName() << endl;
        Info << "-----------------------" << endl;

        Info << "\nUpdate Mesh" << nl << endl;
        mesh.update();
        
        TimeState subCycleTimeState = runTime.subCycle(nSubCycles);

        for (label cycleI = 0; cycleI < nSubCycles; cycleI++)
        {
            Info<< "\n\nsubCycles = " << runTime.timeName() << nl << endl;
            runTime++;
            // --- Pressure-velocity SIMPLE corrector
            {
                #include "UEqn.H"
                #include "pEqn.H"
            }
            laminarTransport.correct();
            turbulence->correct();

            solverPerformance Tres;
            while (simple.correctNonOrthogonal())
            {
                #include "TEqn.H"
            }
            
            if (Tres.initialResidual() < 1e-6 && nStage > 0)
            {
                runTime.printExecutionTime(Info);
                Info << "nStage = " <<nStage << endl;
                Info << "Subcycle converged." << endl;
                break;
            }
            runTime.printExecutionTime(Info);
        }
        runTime.endSubCycle();

        if (!turbulence)
        {
            FatalErrorInFunction
            << "Unable to find turbulence model in the "
            << "database" << exit(FatalError);
            return 1;
        }

        volSymmTensorField Reff = turbulence->devReff();

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

            ssp = (-Sfp/magSfp) & Reffp;

            // 剪切应力 
            scalarField UShear = sqrt(mag(ssp)/rhoAir);         // 剪切风速模量

            scalarField& Mp = M.boundaryFieldRef()[patchi];
            scalarField &deltaHp = deltaH.boundaryFieldRef()[patchi];

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

                deltaHp[i] += Mp[i] / rhoSnow / zArea * runTime.deltaTValue();
            }
        }
        runTime.write();
        runTime.printExecutionTime(Info);
        ++runTime;
        ++nStage;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
