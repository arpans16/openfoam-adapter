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
    namePhiEold_ = CEMdict.lookupOrDefault<word>("namePhiEold", "phiEold");
    DEBUG(adapterInfo("    potential field name : " + namePhiEold_));

    // Read the name of the current field (if different)
    nameJE_ = CEMdict.lookupOrDefault<word>("nameJE", "JE");
    DEBUG(adapterInfo("    Current (J) field name : " + nameJE_));

    // Read the name of the u cross b field (if different)
    nameuxb_ = CEMdict.lookupOrDefault<word>("nameuxb", "UxB");
    DEBUG(adapterInfo("    u cross b field name : " + nameuxb_));

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
    else if (dataName.find("Current") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new Current(mesh_, nameJE_, namePhiE_, nameuxb_)); //Arpan - remove namePhiE and nameuxb later
        DEBUG(adapterInfo("Added writer: Current."));
    }
    else if (dataName.find("CurrentRobin") == 0)
    {
        interface->addCouplingDataWriter(
            dataName,
            new CurrentRobin(mesh_, nameJE_, namePhiE_, namePhiEold_, nameuxb_));
        DEBUG(adapterInfo("Added writer: CurrentRobin."));
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
    else if (dataName.find("Current") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new Current(mesh_, nameJE_, namePhiE_, nameuxb_)); //Arpan - remove nameJE later
        DEBUG(adapterInfo("Added reader: Current."));
    }
    else if (dataName.find("CurrentRobin") == 0)
    {
        interface->addCouplingDataReader(
            dataName,
            new CurrentRobin(mesh_, nameJE_, namePhiE_, namePhiEold_, nameuxb_));
        DEBUG(adapterInfo("Added reader: CurrentRobin."));
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
