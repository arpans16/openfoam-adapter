
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
    dimensionSet(1,-3,0,0,0,0,0),
    transportProperties
);

dimensionedScalar cp
(
    "cp",
    dimensionSet(0,2,-2,-1,0,0,0),
    transportProperties
);

dimensionedScalar sigma
(
    "sigma",
    dimensionSet(-1, -3, 3, 0, 0, 2, 0),
    transportProperties
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

Info<< "Reading diffusivity DT\n" << endl;

volScalarField DT
(
    IOobject
    (
        "DT",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimViscosity, Zero)
);

//dimensionedScalar corrUnits("corrUnits", dimensionSet(-2,-2,5,1,0,2,0), scalar(1.0));

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

if (!DT.headerOk())
{
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
    DT = dimensionedScalar("DT", dimViscosity, transportProperties);
}

#include "createFvOptions.H"
