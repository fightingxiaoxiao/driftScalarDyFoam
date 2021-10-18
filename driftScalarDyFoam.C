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
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    // 获取雪面的Patch名
    // get snow surface patch id
    std::vector<label> snowPatchList;
    const fvPatchList& patches = mesh.boundary();
    forAll(patches,i)
    {
        const std::string name = static_cast<std::string>(patches[i].name()); // 强制转换std::string
        std::smatch match;
        std::regex e("(.snow)");   // 匹配包含".snow"字符串的边界
        if(std::regex_search(name, match, e))
        {
            Info << "Recognized snow surface: \"" << patches[i].name() << "\"." << endl;
            snowPatchList.push_back(patches[i].index());
            phiWf.boundaryFieldRef()[patches[i].index()] = Zero; // 雪面处下落速度为0
        }
    }

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "pEqn.H"
        }

        laminarTransport.correct();
        turbulence->correct();

        while (simple.correctNonOrthogonal())
        {
            //必要时更新下落速度的面通量场
            // update face flux of w_f
            // phiWf = fvc::flux(Wf);

            fvScalarMatrix TEqn
            (
              fvm::ddt(T)
            + fvm::div(phi, T)              // 被动输运         // passive transport
            + fvm::div(phiWf, T)            // 以速度w_f下落    // fall down with velocity w_f
            - fvm::laplacian(nut/S_ct, T)   // 湍流扩散         // turbulent diffusion
            ==
              fvOptions(T)
            );

            TEqn.relax();
            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }

        if (!turbulence)
        {
            FatalErrorInFunction
            << "Unable to find turbulence model in the "
            << "database" << exit(FatalError);
            return 1;            
        }

        Info << "Find turbulence model." <<endl;
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

            const scalar rhoAir = 1.225;
            const scalar rhoSnow = 150;

            vectorField& ssp = wallShearStress.boundaryFieldRef()[patchi];

            ssp = (-Sfp/magSfp) & Reffp;

            // 剪切应力 
            scalarField UShear = sqrt(mag(ssp)/rhoAir);         // 剪切风速模量

            scalarField& Mp = M.boundaryFieldRef()[patchi];   

            // 更新质量交换率（侵蚀/沉积）
            // update mass exchange rate on snow surface (erosion & deposition)
            Info << "Update mass exchange rate." << endl; 
            const scalar UThreshold = 0.2;
            forAll(Mp, i)
            {
                const scalar zArea = Sfp[i] & zNormal;
                if (UShear[i] > UThreshold) // 侵蚀
                {
                    Mp[i] = -5e-4 * rhoSnow * UShear[i] * (1.-sqr(UThreshold)/sqr(UShear[i])) * zArea;
                }
                else // 沉积
                {
                    Mp[i] = Tp[i] * mag(wf.value()) * zArea;
                }
            }
        }

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
