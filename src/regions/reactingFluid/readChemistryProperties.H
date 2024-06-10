Info<< "Reading chemistry properties\n" << endl;

IOdictionary chemistryProperties
(
    IOobject
    (
        "chemistryProperties",
        mesh().time().constant(),
        mesh(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false
    )
);

turbulentReaction_ = Switch(chemistryProperties.lookup("turbulentReaction"));

if (turbulentReaction_)
{
    Cmix_ = readScalar(chemistryProperties.lookup("Cmix"));
}
