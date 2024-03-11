#include "CurrentRobin.H"
#include "primitivePatchInterpolation.H"

#include "fvCFD.H"

using namespace Foam;

//----- preciceAdapter::CEM::Current -----------------------------------------

preciceAdapter::CEM::CurrentRobin::CurrentRobin(
    const Foam::fvMesh& mesh,
    const std::string nameJE,
    const std::string namePhiE,
    const std::string namePhiEold,
    const std::string nameuxb)
: JE_(const_cast<volVectorField*>(&mesh.lookupObject<volVectorField>(nameJE))), 
  phiE_(const_cast<volScalarField*>(&mesh.lookupObject<volScalarField>(namePhiE))),
  phiEold_(const_cast<volScalarField*>(&mesh.lookupObject<volScalarField>(namePhiEold))),
  uxb_(const_cast<volVectorField*>(&mesh.lookupObject<volVectorField>(nameuxb))),
  mesh_(mesh)
{
    dataType_ = scalar;
}

void preciceAdapter::CEM::CurrentRobin::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

	const vectorField norm = mesh_.Sf().boundaryField()[patchID] / mesh_.magSf().boundaryField()[patchID];
	const scalarField JE_wn(JE_->boundaryField()[patchID] & norm);

        // If we use the mesh connectivity, we interpolate from the centres to the nodes
        if (meshConnectivity)
        {
            //Setup Interpolation object
            primitivePatchInterpolation patchInterpolator(mesh_.boundaryMesh()[patchID]);

	    //Interpolate on patches
	    scalarField JEPoints = patchInterpolator.faceToPointInterpolate(JE_wn);

            // For every cell of the patch
            forAll(JEPoints, i)
            {
                buffer[bufferIndex++] = JEPoints[i];
            }
        }
        else
        {
            // For every cell of the patch
            forAll(JE_wn, i)
            {
                buffer[bufferIndex++] = JE_wn[i];
            }
        }
    }
}

void preciceAdapter::CEM::CurrentRobin::read(double* buffer, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
	double sigma = 400.0; //Arpan - read in this later
        int patchID = patchIDs_.at(j);

        // Get the potential gradient boundary patch
        scalarField& gradientPatch(refCast<fixedGradientFvPatchScalarField>(phiE_->boundaryFieldRef()[patchID]).gradient());
        const scalarField& delCoef(refCast<fixedGradientFvPatchScalarField>(phiE_->boundaryFieldRef()[patchID]).patch().deltaCoeffs());

        vectorField norm = mesh_.Sf().boundaryField()[patchID] / mesh_.magSf().boundaryField()[patchID];
        scalarField uxb_scalar(uxb_->boundaryField()[patchID] & norm);
	//fvPatchField del = phiE_->boundaryField()[patchID].delta();

        // For every cell of the patch
        forAll(gradientPatch, i)
        {
	    //Arpan changes March 08 2024
	    //new Robin - attempt
	    double fluxin 	= uxb_scalar[i] - buffer[bufferIndex++]/sigma;
	    double alpha 	= 0.2*(delCoef[i]*2); //phiE_->boundaryField()[patchID].deltaCoeffs()[i];
	    double phi_new	= phiE_->boundaryFieldRef()[patchID][i];
	    double phi_old	= phiEold_->boundaryFieldRef()[patchID][i];
	    gradientPatch[i] 	= alpha * (phi_new - phi_old) + fluxin;
        }
    }
}

bool preciceAdapter::CEM::CurrentRobin::isLocationTypeSupported(const bool meshConnectivity) const
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

std::string preciceAdapter::CEM::CurrentRobin::getDataName() const
{
    return "CurrentRobin";
}
