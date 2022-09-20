#include "ElectricFlux.H"
#include "primitivePatchInterpolation.H"

//#include "fvCFD.H"

using namespace Foam;

//----- preciceAdapter::CEM::ElectricFlux -----------------------------------------

preciceAdapter::CEM::ElectricFlux::ElectricFlux(
    const Foam::fvMesh& mesh,
    const std::string namePhiE)
: phiE_(
    const_cast<volScalarField*>(
        &mesh.lookupObject<volScalarField>(namePhiE))),
  mesh_(mesh)
{
    dataType_ = scalar;
}

void preciceAdapter::CEM::ElectricFlux::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        const scalarField gradientPatch(
            (phiE_->boundaryField()[patchID])
                .snGrad());

        // If we use the mesh connectivity, we interpolate from the centres to the nodes
        if (meshConnectivity)
        {
            //Setup Interpolation object
            primitivePatchInterpolation patchInterpolator(mesh_.boundaryMesh()[patchID]);

            scalarField gradientPoints;

            //Interpolate
            gradientPoints = patchInterpolator.faceToPointInterpolate(gradientPatch);

            // For every cell of the patch
            forAll(gradientPoints, i)
            {
                // Copy the heat flux into the buffer
                // Q = - k * gradient(T)
                //TODO: Interpolate kappa in case of a turbulent calculation
                buffer[bufferIndex++] = - gradientPoints[i];
            }
        }
        else
        {
            // For every cell of the patch
            forAll(gradientPatch, i)
            {
                // Copy the heat flux into the buffer
                // Q = - k * gradient(T)
                //TODO: Interpolate kappa in case of a turbulent calculation
                buffer[bufferIndex++] = - gradientPatch[i];
            }
        }
    }
}

void preciceAdapter::CEM::ElectricFlux::read(double* buffer, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // Get the potential gradient boundary patch
        scalarField& gradientPatch(
            refCast<fixedGradientFvPatchScalarField>(
                phiE_->boundaryFieldRef()[patchID])
                .gradient());

        // For every cell of the patch
        forAll(gradientPatch, i)
        {
            // Compute and assign the gradient from the buffer.
            // The sign of the heat flux needs to be inversed,
            // as the buffer contains the flux that enters the boundary:
            // gradient(T) = -Q / -k
            gradientPatch[i] = buffer[bufferIndex++];
        }
    }
}

bool preciceAdapter::CEM::ElectricFlux::isLocationTypeSupported(const bool meshConnectivity) const
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

std::string preciceAdapter::CEM::ElectricFlux::getDataName() const
{
    return "ElectricFlux";
}
