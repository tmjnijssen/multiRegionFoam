Info<< nl << "Reading thermophysicalProperties" << endl;

hsCombustionThermo& thermo = pChemistry_().thermo();
PtrList<volScalarField>& Y = thermo.composition().Y();

rho_ = lookupOrRead<volScalarField>
(
    mesh(),
    "rho",
    false,
    true,
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
U_ = lookupOrRead<volVectorField>
(
    mesh(),
    "U",
    true,
    true
);

p_ = lookupOrRead<volScalarField>
(
    mesh(),
    "p",
    true,
    true,
    thermo.p()
);

psi_ = lookupOrRead<volScalarField>
(
    mesh(),
    "psi",
    true,
    true,
    thermo.psi()
);

hs_ = lookupOrRead<volScalarField>
(
    mesh(),
    "hs",
    true,
    true,
    thermo.hs()
);

T_ = lookupOrRead<volScalarField>
(
    mesh(),
    "T",
    true,
    true,
    thermo.T()
);

phi_ = lookupOrRead<surfaceScalarField>
(
    mesh(),
    "phi",
    false,
    true,
    linearInterpolate(rho_()*U_()) & mesh().Sf()
);

Info << "Creating turbulence model.\n" << nl;
turbulence_ = compressible::turbulenceModel::New
(
    rho_(),
    U_(),
    phi_(),
    thermo
);

sigma_ = lookupOrRead<volSymmTensorField>
(
    mesh(),
    "sigma",
    false,
    true,
    -p_()*symmTensor(1,0,0,1,0,1) - turbulence_().devRhoReff()
); 

Info<< "Creating field DpDt\n" << endl;
DpDt_ = lookupOrRead<volScalarField>
(
    mesh(),
    "DpDt",
    false,
    true,
    fvc::DDt
    (
        surfaceScalarField("phiU", phi_()/fvc::interpolate(rho_())),
        p_()
    )
);

forAll(Y, i)
{
    fields_.add(Y[i]);
}
fields_.add(hs_());

mesh().schemesDict().setFluxRequired(p_().name());
