#ifndef CEM_POTENTIAL_H
#define CEM_POTENTIAL_H

#include "CouplingDataUser.H"

#include "fvCFD.H"

namespace preciceAdapter
{
namespace CEM
{

//- Class that writes and reads potential
class Potential : public CouplingDataUser
{

private:
    //- Potential field
    Foam::volScalarField* phiE_;
    const Foam::fvMesh& mesh_;

public:
    //- Constructor
    Potential(
        const Foam::fvMesh& mesh,
        const std::string namePhiE);

    //- Write the potential values into the buffer
    void write(double* buffer, bool meshConnectivity, const unsigned int dim);

    //- Read the potential values from the buffer
    void read(double* buffer, const unsigned int dim);

    bool isLocationTypeSupported(const bool meshConnectivity) const override;

    //- Get the name of the current data field
    std::string getDataName() const override;
};

}
}

#endif
