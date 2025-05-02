using namespace Foam;

volScalarField p(IOobject("p", runTime.name(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
volVectorField U(IOobject("U", runTime.name(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
surfaceScalarField phi(IOobject("phi", runTime.name(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                       linearInterpolate(U) & mesh.Sf());
label pRefCell = 0;
scalar pRefValue = 0;
setRefCell(p, mesh.solution().dict().subDict("PISO"), pRefCell, pRefValue);

IOdictionary transportProperties(
    IOobject("transportProperties", runTime.constant(), mesh, IOobject::MUST_READ_IF_MODIFIED, IOobject::NO_WRITE));
dimensionedScalar nu(transportProperties.lookup("nu"));

// sdfibm related
dimensionedScalar rho(transportProperties.lookup("rho"));
dimensionedScalar alpha(transportProperties.lookup("alpha"));

volScalarField As(IOobject("As", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
                  mesh,
                  dimensionedScalar("As", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0),
                  "zeroGradient");

volScalarField Ct(IOobject("Ct", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
                  mesh,
                  dimensionedScalar("Ct", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0),
                  "zeroGradient");
volVectorField Fs(IOobject("Fs", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                  mesh,
                  dimensionedVector("Fs", dimAcceleration, vector::zero),
                  "fixedValue");
volScalarField T(IOobject("T", runTime.name(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
volScalarField Ts(IOobject("Ts", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                  mesh,
                  dimensionedScalar("NULL", dimensionSet(0, 0, 0, 1, 0, 0, 0), 0.0));
