#ifndef CEM_CURRENT_H
#define CEM_CURRENT_H

#include "CouplingDataUser.H"

#include "fvCFD.H"

namespace preciceAdapter
{
namespace CEM
{

//- Abstract class that writes and reads current (electric flux).
class Current : public CouplingDataUser
{

protected:
    //- Potential field
    Foam::volScalarField* phiE_;
    Foam::volVectorField* uxb_;
    const Foam::fvMesh& mesh_;

public:
    //- Constructor
    Current(const Foam::fvMesh& mesh, const std::string namePhiE, const std::string nameuxb);

    //- Compute current values from the potential field
    //  and write them into the buffer
    virtual void write(double* buffer, bool meshConnectivity, const unsigned int dim);

    //- Read current values from the buffer and assign them to
    //  the gradient of the potential field
    virtual void read(double* buffer, const unsigned int dim);

    bool isLocationTypeSupported(const bool meshConnectivity) const override;

    //- Get the name of the current data field
    std::string getDataName() const override;

    //- Destructor
    virtual ~Current() {};
};

}
}


#endif
