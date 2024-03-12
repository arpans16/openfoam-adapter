#ifndef CEM_H
#define CEM_H

#include "Interface.H"

#include "CEM/Potential.H"
#include "CEM/Current.H"
#include "CEM/CurrentRobin.H"

#include "fvCFD.H"

namespace preciceAdapter
{
namespace CEM
{

class ConjugateElectroMagnetics
{

protected:
    //- OpenFOAM fvMesh object
    const Foam::fvMesh& mesh_;

    // TODO: Create a struct for all the parameter names

    //- Name of the electrical conductivity for a basic solver
    std::string nameSigma_ = "sigma";

    //- Name of the potential field
    std::string namePhiE_ = "phiE";
    std::string namePhiEold_ = "phiEold";

    //- Name of the current field
    std::string nameJE_ = "JE";

    //- Name of the u cross b field
    std::string nameuxb_ = "UxB";

protected:
    //- Read the CEM-related options from the adapter's configuration file
    bool readConfig(const IOdictionary& adapterConfig);

public:
    //- Constructor
    ConjugateElectroMagnetics(const Foam::fvMesh& mesh);

    //- Configure
    bool configure(const IOdictionary& adapterConfig);

    //- Add coupling data writers
    bool addWriters(std::string dataName, Interface* interface);

    //- Add coupling data readers
    bool addReaders(std::string dataName, Interface* interface);
};

}
}

#endif
