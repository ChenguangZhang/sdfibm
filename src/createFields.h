volScalarField p ( IOobject ( "p", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE ), mesh );
volVectorField U ( IOobject ( "U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE ), mesh );
surfaceScalarField phi ( IOobject ( "phi", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE ),
                         linearInterpolate(U) & mesh.Sf() );
label  pRefCell  = 0;
scalar pRefValue = 0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

IOdictionary transportProperties
(
    IOobject(
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
dimensionedScalar nu   (transportProperties.lookup("nu") );

// sdfibm related
dimensionedScalar rho  (transportProperties.lookup("rho"));
dimensionedScalar alpha(transportProperties.lookup("alpha"));

volScalarField As (IOobject ("As", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE ), mesh,
                    dimensionedScalar("As", dimless, 0.0), zeroGradientFvPatchScalarField::typeName);
volScalarField Ct (IOobject ("Ct", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE ), mesh,
                    dimensionedScalar("Ct", dimless, 0.0), zeroGradientFvPatchScalarField::typeName);
volVectorField Fs (IOobject ("Fs", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ), mesh,
                    dimensionedVector("Fs", dimAcceleration, vector::zero), zeroGradientFvPatchVectorField::typeName );
volScalarField T  (IOobject ( "T", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE ), mesh );
volScalarField Ts (IOobject ("Ts", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ), mesh,
                   dimensionedScalar("NULL", dimensionSet(0, 0, 0, 1, 0, 0, 0), 0.0) );
