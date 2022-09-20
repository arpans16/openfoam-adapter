#include "Potential.H"
#include "primitivePatchInterpolation.H"


using namespace Foam;

preciceAdapter::CEM::Potential::Potential(
    const Foam::fvMesh& mesh,
    const std::string namePhiE)
: phiE_(
    const_cast<volScalarField*>(
        &mesh.lookupObject<volScalarField>(namePhiE))),
  mesh_(mesh)
{
    dataType_ = scalar;
}

void preciceAdapter::CEM::Potential::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        const scalarField& phiEPatch(
            phiE_->boundaryField()[patchID]);

        //If we use the mesh connectivity, we interpolate from the centres to the nodes
        if (meshConnectivity)
        {
            //Create an Interpolation object at the boundary Field
            primitivePatchInterpolation patchInterpolator(mesh_.boundaryMesh()[patchID]);

            //Interpolate from centers to nodes
            scalarField phiEPoints(patchInterpolator.faceToPointInterpolate(phiEPatch));

            forAll(phiEPoints, i)
            {
                // Copy the potential into the buffer
                buffer[bufferIndex++] = phiEPoints[i];
            }
        }
        else
        {
            forAll(phiEPatch, i)
            {
                // Copy the potential into the buffer
                buffer[bufferIndex++] = phiEPatch[i];
            }
        }
    }
}

void preciceAdapter::CEM::Potential::read(double* buffer, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(phiE_->boundaryField()[patchID], i)
        {
            // Set the potential as the buffer value
            phiE_->boundaryFieldRef()[patchID][i] = buffer[bufferIndex++];
        }
    }
}

bool preciceAdapter::CEM::Potential::isLocationTypeSupported(const bool meshConnectivity) const
{
    // For cases with mesh connectivity, we support:
    // - face nodes, only for writing
    // - face centers, only for reading
    // However, since we do not distinguish between reading and writing in the code, we
    // always return true and offload the handling to the user.
    if (meshConnectivity)
    {
        return true;
    }
    else
    {
        return (this->locationType_ == LocationType::faceCenters);
    }
}

std::string preciceAdapter::CEM::Potential::getDataName() const
{
    return "Potential";
}
