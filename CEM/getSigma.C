#include "getSigma.H"
#include "primitivePatchInterpolation.H"
#include "Utilities.H"

using namespace Foam;

//----- preciceAdapter::CEM::getSigma ------------------

preciceAdapter::CEM::getSigma::getSigma(
  const Foam::fvMesh& mesh,
  const std::string nameSigma)
: mesh_(mesh),
  nameSigma_(nameSigma)
{
    DEBUG(adapterInfo("Constructed getSigma."));

    // Get the preciceDict/CHT dictionary
    const dictionary& CEMDict = mesh_.lookupObject<IOdictionary>("preciceDict").subOrEmptyDict("CEM");

    // Read the conductivity
    if (!CEMDict.readIfPresent<dimensionedScalar>(nameSigma_, sigma_))
    {
        adapterInfo("Cannot find the conductivity in preciceDict/CEM using the name " + nameSigma_, "error");
    }
    else
    {
        DEBUG(adapterInfo("sigma = " + std::to_string(sigma_.value())));
    }    
}

void preciceAdapter::CEM::getSigma::extract()
{
    // Already extracted in the constructor
}

scalar preciceAdapter::CEM::getSigma::getValue()
{
    return sigma_.value();
}
