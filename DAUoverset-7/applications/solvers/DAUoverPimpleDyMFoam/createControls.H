#include "createTimeControls.H"

bool correctPhi
(
    //pimple.dict().getOrDefault("correctPhi", false)
    pimple.dict().lookupOrDefault<bool>("correctPhi", false)
);

bool checkMeshCourantNo
(
    //pimple.dict().getOrDefault("checkMeshCourantNo", false)
    pimple.dict().lookupOrDefault<bool>("checkMeshCourantNo", false)
);

bool massFluxInterpolation
(
    //pimple.dict().getOrDefault("massFluxInterpolation", false)
    pimple.dict().lookupOrDefault<bool>("massFluxInterpolation", false)
);

bool adjustFringe
(
    //pimple.dict().getOrDefault("oversetAdjustPhi", false)
    pimple.dict().lookupOrDefault<bool>("oversetAdjustPhi", false)
);

bool ddtCorr
(
    //pimple.dict().getOrDefault("ddtCorr", true)
    pimple.dict().lookupOrDefault<bool>("ddtCorr", true)
);
