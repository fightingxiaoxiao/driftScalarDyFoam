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

Info<< "Reading field nut\n" << endl;

volScalarField nut
(
    IOobject
    (
        "nut",
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

dimensionedVector wf("wf", dimVelocity, transportProperties);

// (Tominaga, 2011)
// (Kang, 2018)
dimensionedScalar S_ct("S_ct", 
                        dimensionSet(0, 0, 0, 0, 0, 0 ,0), 
                        transportProperties); 

volVectorField Wf
(
    IOobject
    (
        "Wf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    wf
);

surfaceScalarField phiWf
(
    IOobject
    (
        "phiWf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(Wf)
);

dimensionedScalar M_("M_", 
                    dimensionSet(0, 0, 0, 0, 0, 0 ,0), 
                    Zero); 

surfaceScalarField M
(
    IOobject
    (
        "M",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    M_
);


Info<< "Reading erosionDepositionProperties\n" << endl;

IOdictionary erosionDepositionProperties
(
    IOobject
    (
        "erosionDepositionProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


#include "createPhi.H"

#include "createFvOptions.H"
