#include "readTimeControls.H"

//correctPhi = pimple.dict().getOrDefault("correctPhi", false);
correctPhi = pimple.dict().lookupOrDefault<bool>("correctPhi", false);

//checkMeshCourantNo = pimple.dict().getOrDefault("checkMeshCourantNo", false);
checkMeshCourantNo = pimple.dict().lookupOrDefault<bool>("checkMeshCourantNo", false);

//massFluxInterpolation = pimple.dict().getOrDefault("massFluxInterpolation", false);
massFluxInterpolation = pimple.dict().lookupOrDefault<bool>("massFluxInterpolation", false);

//ddtCorr = pimple.dict().getOrDefault("ddtCorr", true);
ddtCorr = pimple.dict().lookupOrDefault<bool>("ddtCorr", true);
