#include "fvMesh.H"
#include "Time.H"
#include "argList.H" // TODO: add some cmd arguments
#include "findRefCell.H"
#include "adjustPhi.H"

#include "fvcGrad.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"

#include "solidcloud.h"

int main(int argc, char* argv[])
{
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.h"
#include "initContinuityErrs.H"

    std::string dictfile;

    // if start-time > 0, read from start-time-folder for solidDict, otherwise read from case root
    if (runTime.time().value() > 0)
    {
        if (!Foam::Pstream::parRun())
            dictfile = mesh.time().name() + "/solidDict";
        else
            dictfile = "processor0/" + mesh.time().name() + "/solidDict";
    }
    else
    {
        dictfile = "solidDict";
    }

    sdfibm::SolidCloud solidcloud(runTime.globalPath() + "/" + dictfile, U, runTime.value());
    solidcloud.saveState(); // write the initial condition

    while (runTime.loop())
    {
        Foam::Info << "Time = " << runTime.name() << Foam::endl;

#include "CourantNo.H"
        Foam::dimensionedScalar dt = runTime.deltaT();

        if (solidcloud.isOnFluid())
        {
            Foam::fvVectorMatrix UEqn(fvm::ddt(U) + 1.5 * fvc::div(phi, U) -
                                          0.5 * fvc::div(phi.oldTime(), U.oldTime()) ==
                                      0.5 * fvm::laplacian(nu, U) + 0.5 * fvc::laplacian(nu, U));
            UEqn.solve();

            phi = linearInterpolate(U) & mesh.Sf();
            Foam::fvScalarMatrix pEqn(fvm::laplacian(p) == fvc::div(phi) / dt - fvc::div(Fs));
            pEqn.solve();

            U = U - dt * fvc::grad(p);
            phi = phi - dt * fvc::snGrad(p) * mesh.magSf();

            Foam::fvScalarMatrix TEqn(fvm::ddt(T) + fvm::div(phi, T) == fvm::laplacian(alpha, T));
            TEqn.solve();
        }

        solidcloud.interact(runTime.value(), dt.value());

        if (solidcloud.isOnFluid())
        {
            U = U - Fs * dt;
            phi = phi - dt * (linearInterpolate(Fs) & mesh.Sf());

            U.correctBoundaryConditions();
            adjustPhi(phi, U, p);

            T = (1.0 - As) * T + Ts;
            T.correctBoundaryConditions();

#include "continuityErrs.H"
        }

        solidcloud.evolve(runTime.value(), dt.value());
        solidcloud.saveState();

        if (solidcloud.isOnFluid())
        {
            solidcloud.fixInternal(dt.value());
        }

        if (runTime.writeTime())
        {
            runTime.write();

            if (Foam::Pstream::master())
            {
                std::string file_name;
                if (Foam::Pstream::parRun())
                    file_name = "./processor0/" + runTime.name() + "/solidDict";
                else
                    file_name = "./" + runTime.name() + "/solidDict";
                solidcloud.saveRestart(file_name);
            }
        }
    }

    Foam::Info << "DONE\n" << Foam::endl;
    return 0;
}
