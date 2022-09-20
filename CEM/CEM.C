#include "CEM.H"

#include "Utilities.H"

using namespace Foam;

preciceAdapter::CEM::ConjugateElectroMagnetics::ConjugateElectroMagnetics(
    const Foam::fvMesh& mesh)
: mesh_(mesh)
{
}

bool preciceAdapter::CEM::ConjugateElectroMagnetics::configure(const IOdictionary& adapterConfig)
{
    DEBUG(adapterInfo("Configuring the CEM module..."));

    // Read the CEM-specific options from the adapter's configuration file
    if (!readConfig(adapterConfig))
    {
        return false;
    }

    DEBUG(adapterInfo("Nothing to configure, only basic incompressible implemented."));

    return true;
}

bool preciceAdapter::CEM::ConjugateElectroMagnetics::readConfig(const IOdictionary& adapterConfig)
{
    const dictionary& CEMdict = adapterConfig.subOrEmptyDict("CEM");

    // Read the name of the potential field (if different)
    namePhiE_ = CEMdict.lookupOrDefault<word>("namePhiE", "phiE");
    DEBUG(adapterInfo("    potential field name : " + namePhiE_));

    return true;
}

bool preciceAdapter::CEM::ConjugateElectroMagnetics::addWriters(std::string dataName, Interface* interface)
{
    bool found = true; // Set to false later, if needed.

    if (dataName.find("Potential") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new Potential(mesh_, namePhiE_));
        DEBUG(adapterInfo("Added writer: Potential."));
    }
    else if (dataName.find("ElectricFlux") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new ElectricFlux(mesh_, namePhiE_));
        DEBUG(adapterInfo("Added writer: Electric Flux."));
    }
    else
    {
        found = false;
    }

    // NOTE: If you want to couple another variable, you need
    // to add your new coupling data user as a coupling data
    // writer here (and as a reader below).
    // The argument of the dataName.compare() needs to match
    // the one provided in the adapter's configuration file.

    return found;
}

bool preciceAdapter::CEM::ConjugateElectroMagnetics::addReaders(std::string dataName, Interface* interface)
{
    bool found = true; // Set to false later, if needed.

    if (dataName.find("Potential") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new Potential(mesh_, namePhiE_));
        DEBUG(adapterInfo("Added reader: Potential."));
    }
    else if (dataName.find("ElectricFlux") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new ElectricFlux(mesh_, namePhiE_));
        DEBUG(adapterInfo("Added reader: Electric Flux."));
    }
    else
    {
        found = false;
    }

    // NOTE: If you want to couple another variable, you need
    // to add your new coupling data user as a coupling data
    // reader here (and as a writer above).
    // The argument of the dataName.compare() needs to match
    // the one provided in the adapter's configuration file.

    return found;
}
