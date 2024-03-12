Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar rho
(
    "rho",
    dimDensity,
    transportProperties
);

dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties
);

dimensionedScalar cp
(
    "cp",
    dimensionSet(0,2,-2,-1,0,0,0),
    transportProperties
);

dimensionedScalar mu
(
    "mu",
    dimensionSet(1, -1, -1, 0, 0, 0, 0),
    transportProperties
);

dimensionedScalar sigma
(
    "sigma",
    dimensionSet(-1, -3, 3, 0, 0, 2, 0),
    transportProperties
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Reading field B\n" << endl;
volVectorField B
(
    IOobject
    (
        "B",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
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

Info<< "Reading field phiEold\n" << endl;
volScalarField phiEold
(
    IOobject
    (
        "phiEold",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
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

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
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
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Lorentz\n" << endl;
volVectorField Lorentz
(
    IOobject
    (
        "Lorentz",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//dimensionedScalar corrUnits("corrUnits", dimensionSet(0,4,-1,1,0,-2,0), scalar(1.0));

Info<< "Reading field JouleHeating\n" << endl;
volScalarField JouleHeating
(
    IOobject
    (
        "JouleHeating",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

/*
Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
Info<< "Reading diffusivity DT\n" << endl;
dimensionedScalar DT("DT", dimViscosity, transportProperties);
#include "createFvOptions.H"
*/
#include "createFvOptions.H"
dimensionedScalar DT("DT", dimViscosity, transportProperties);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, piso.dict(), pRefCell, pRefValue);

mesh.setFluxRequired(p.name());
