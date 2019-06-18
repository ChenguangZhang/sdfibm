#include "fvCFD.H"
#include "solidcloud.h"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.h"
    #include "initContinuityErrs.H"

    SolidCloud solidcloud;

    // link solid cloud with the fluid
    solidcloud.linkMesh(mesh);
    if(solidcloud.isOnFluid())
    {
        solidcloud.linkFluid(U);
    }

    while (runTime.loop())
    {
        Foam::Info << "Time = " << runTime.timeName() << Foam::endl;

        #include "CourantNo.H"
        Foam::dimensionedScalar dt = runTime.deltaT();

        if(solidcloud.isOnFluid())
        {
            Foam::fvVectorMatrix UEqn(
                fvm::ddt(U)
              + 1.5*fvc::div(phi, U) - 0.5*fvc::div(phi.oldTime(), U.oldTime())
              ==0.5*fvm::laplacian(nu, U) + 0.5*fvc::laplacian(nu, U));
            UEqn.solve();

            Foam::surfaceScalarField phiI("phiI", linearInterpolate(U) & mesh.Sf());
            Foam::fvScalarMatrix pEqn(fvm::laplacian(p) == fvc::div(phiI)/dt);
            pEqn.setReference(0, 0); // method of fvMatrix
            pEqn.solve();

            U   = U    - dt*fvc::grad(p);
            phi = phiI - dt*fvc::snGrad(p)*mesh.magSf();

            // solve tracer after U
            Foam::fvScalarMatrix TEqn(
                fvm::ddt(T)
              + fvm::div(phi, T)
              ==fvm::laplacian(alpha, T));
            TEqn.solve();

            // above is fluid-only, now account for solid effect
            solidcloud.interact(runTime.value(), dt.value());
            U = U - Fs*dt;

            U.correctBoundaryConditions();
            adjustPhi(phi, U, p);
            
            T = (1.0 - As)*T + Ts;
            T.correctBoundaryConditions();
            Tphi = phi*fvc::interpolate(T);

            #include "continuityErrs.H"
        }

        solidcloud.evolve(runTime.value(), dt.value());

        if(solidcloud.isOnFluid())
            solidcloud.fixInternal();

        if(runTime.outputTime())
        {
            // write flow fields
            runTime.write();

            if(Foam::Pstream::master())
            {
                std::string file_name;
                if(Foam::Pstream::parRun())
                    file_name = "./processor0/" + runTime.timeName() + "/restart";
                else
                    file_name = "./" + runTime.timeName() + "/restart";
                solidcloud.saveRestart(file_name);
            }
        }
    }

    Info << "DONE\n" << endl;
    return 0;
}
