#ifndef CEM_ELECTRICFLUX_H
#define CEM_ELECTRICFLUX_H

#include "CouplingDataUser.H"

#include "fvCFD.H"

namespace preciceAdapter
{
namespace CEM
{

//- Abstract class that writes and reads electric flux.
class ElectricFlux : public CouplingDataUser
{

protected:
    //- Potential field
    Foam::volScalarField* phiE_;
    const Foam::fvMesh& mesh_;

public:
    //- Constructor
    ElectricFlux(const Foam::fvMesh& mesh, const std::string namePhiE);

    //- Compute electric flux values from the potential field
    //  and write them into the buffer
    virtual void write(double* buffer, bool meshConnectivity, const unsigned int dim);

    //- Read electric flux values from the buffer and assign them to
    //  the gradient of the potential field
    virtual void read(double* buffer, const unsigned int dim);

    bool isLocationTypeSupported(const bool meshConnectivity) const override;

    //- Get the name of the current data field
    std::string getDataName() const override;

    //- Destructor
    virtual ~ElectricFlux() {};
};

}
}


#endif
