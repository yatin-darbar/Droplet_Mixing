Info<< "Creating dummy fields for Kistler CA model\n" << endl;

// Dummy viscosity field (adds to the object registry as volScalarField)
volScalarField muKistler
(
    IOobject
    (
	"muKistler",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("muKistler",dimensionSet(1,-1,-1,0,0,0,0),10000.0)
);

// Dummy surface tension field (adds to the object registry as volScalarField)
volScalarField sigma
(
    IOobject
    (
        "sigmaKistler",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(transportProperties.lookup("sigma"))
);

