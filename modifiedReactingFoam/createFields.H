

//Info<< nl << "Reading thermophysicalProperties" << endl;
Info<< nl << "Creating combustion model" << endl;

// A.Alexiou 2014
// autoPtr<combustionModels:: psiChemistryCombustionModel> combustion // OF 2.1
autoPtr<combustionModels::psiCombustionModel> combustion // OF 2.3
(
 // A.Alexiou 2014
 //combustionModels::psiChemistryCombustionModel::New(mesh) // OF 2.1
 combustionModels::psiCombustionModel::New(mesh) // OF 2.3
);

// M. Lindner, A.Alexiou 2014
// psiChemistryModel& chemistry = combustion->pChemistry(); // OF 2.1
autoPtr<psiChemistryModel> pChemistry // OF 2.3
(
    psiChemistryModel::New(mesh) // OF 2.3
);

// M. Lindner 2014
psiChemistryModel& chemistry = pChemistry(); // OF 2.3

// A.Alexiou 2014
// hsCombustionThermo& thermo = chemistry.thermo(); // OF 2.1
psiReactionThermo& thermo = chemistry.thermo(); // OF 2.3

basicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();



word inertSpecie(thermo.lookup("inertSpecie"));

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
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

volScalarField& p = thermo.p();
const volScalarField& psi = thermo.psi();

// A.Alexiou 2014
// volScalarField& hs = thermo.hs();
volScalarField& hs = thermo.he(); // .hs();

// A.Alexiou 2014 - runtime.write() does not work for field T for reason unknown
// Instance of an object T didn't help change that
//    volScalarField T
//    (
//        IOobject
//        (
//            "T",
//            runTime.timeName(),
//            mesh,
//            IOobject::MUST_READ,
//            IOobject::AUTO_WRITE
//        ),
//        thermo.T()
//    );

const volScalarField& T = thermo.T();

#include "compressibleCreatePhi.H"

volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

// Set the turbulence into the combustion model
combustion->setTurbulence(turbulence());

/*Info<< "Creating field DpDt\n" << endl;
volScalarField DpDt =
fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);*/

Info<< "Creating field dpdt\n" << endl; //Mohsen, may have to use DpDt!
volScalarField dpdt("dpdt", fvc::ddt(p));

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(hs);

DimensionedField<scalar, volMesh> chemistrySh
(
    IOobject
    (
        "chemistry::Sh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("chemistrySh", dimEnergy/dimTime/dimVolume, 0.0)
);

volScalarField dQ
(
    IOobject
    (
        "dQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dQ", dimEnergy/dimTime, 0.0)
);


// ------------------------------------------------------------------------- //
autoPtr<multiSpeciesTransportModel> mstm
(
    multiSpeciesTransportModel::New
    (
        thermo,
        turbulence()
    )
);
// ------------------------------------------------------------------------- //
