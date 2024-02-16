
    Info<< "Reading physicalProperties\n" << endl;

    IOdictionary physicalProperties
    (
        IOobject
        (
            "physicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading field phiE\n" << endl;
    volScalarField phiE
    (
        IOobject
        (
            "phiE",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field JE\n" << endl;
    volVectorField JE
    (
        IOobject
        (
            "JE",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

Info<< "Reading field UxB\n" << endl;
volVectorField UxB
(
    IOobject
    (
        "UxB",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

dimensionedScalar corrUnits("corrUnits", dimensionSet(-2,-2,5,1,0,2,0), scalar(1.0));

#include "createFvOptions.H"
