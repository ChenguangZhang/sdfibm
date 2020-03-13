#include "fvCFD.H"
#include "solidcloud.h"

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addOption("name", "word", "Name of field to be initialized");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if (!args.checkRootCase())
        Foam::FatalError.exit();

    Foam::word field_name;

    if(args.options().empty())
    {
        field_name = "alpha.water";
        Foam::Info << "* No field name provided, default to " << field_name << "\n";
    }
    else
        args.optionReadIfPresent("name", field_name);

    volScalarField      As(IOobject(field_name, runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
    volScalarField      Ct(IOobject(      "Ct", runTime.timeName(), mesh, IOobject::NO_READ,   IOobject::NO_WRITE  ), mesh, dimensionedScalar("Ct",  dimless, 0.0), zeroGradientFvPatchScalarField::typeName);
    surfaceScalarField phi(IOobject(     "phi", runTime.timeName(), mesh, IOobject::NO_READ,   IOobject::NO_WRITE  ), mesh, dimensionedScalar("tmp", dimless, 0.0), "calculated");

    sdfibm::SolidCloud solidcloud;
    solidcloud.setFieldName(field_name);

    solidcloud.writeVOF(mesh);

    return 0;
}
