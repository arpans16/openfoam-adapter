#ifndef CEM_GET_SIGMA_H
#define CEM_GET_SIGMA_H

#include "fvCFD.H"
#include "fvOptions.H"
#include "pisoControl.H"
//#include "addCheckCaseOptions.H"
//#include "setRootCaseLists.H"
//#include "createTime.H"
//#include "createMesh.H"
//#include "createControl.H"

#include "turbulentTransportModel.H"

namespace preciceAdapter
{
namespace CEM
{

//- Class that extracts the electric conductivity.
class getSigma
{

protected:
    //- OpenFOAM fvMesh object
    const Foam::fvMesh& mesh_;
    //- Electrical conductivity
    Foam::dimensionedScalar sigma_;
    //- OpenFOAM transportProperties
    // const Foam::IOdictionary& transProps_;
    //- Name of the user-provided conductivity (in the preciceDict)
    const std::string nameSigma_;

public:
    //- Constructor
    getSigma(const Foam::fvMesh& mesh, const std::string nameKappa);

    //- Extract the electrical conductivity and store it.
    void extract();

    //- Get the value of sigma
    Foam::scalar getValue();
};

}
}

#endif
