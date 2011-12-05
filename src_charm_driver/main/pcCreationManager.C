#include "pcCreationManager.h"
#include "paircalc/pcBuilder.h"
#include "orthog_ctrl/orthoBuilder.h"
#include "load_balance/IntMap.h"
#include "ortho.decl.h"

namespace cp {
    namespace startup {

PCCreationManager::PCCreationManager(const paircalc::pcConfig &_symmCfg, const paircalc::pcConfig &_asymmCfg, const ortho::orthoConfig &_orthoCfg):
    symmCfg(_symmCfg), asymmCfg(_asymmCfg), orthoCfg(_orthoCfg)
{
    if (symmCfg.orthoGrainSize != asymmCfg.orthoGrainSize)
        CkAbort("Ortho grain size mismatch in supllied configs\n");
    if (symmCfg.instanceIndex != asymmCfg.instanceIndex)
        CkAbort("Cannot wire together two PCs from different instances!!");
}




void PCCreationManager::build(CkCallback cb, const PCMapConfig mapCfg)
{
    /// Create the message that will hold the handles to created chares
    pcSetupMsg *msg = new pcSetupMsg();
    msg->symmCfg    = symmCfg;
    msg->asymmCfg   = asymmCfg;

    // Create the symmetric (psi) and asymmetric (lambda) paircalc instances
    CkPrintf("\n\nCreating the symmetric (psi) and asymmetric (lambda) paircalculators\n");
    cp::paircalc::Builder symmBuilder(symmCfg), asymmBuilder(asymmCfg);
    msg->symmIDs  = symmBuilder.build (mapCfg.boxSize, mapCfg.getPeList, mapCfg.gSpaceMap);
    msg->asymmIDs = asymmBuilder.build(mapCfg.boxSize, mapCfg.getPeList, mapCfg.gSpaceMap);

    // Spawn the ortho array and its world of chares/classes (CLA_Matrix, OrthoHelper etc.)
    CkPrintf("Creating the ortho array\n");
    cp::ortho::Builder orthoBuilder(orthoCfg);
    msg->orthoAID = orthoBuilder.build(msg->asymmIDs, mapCfg.getPeList);

    // Ask ortho to setup its communication sections of paircalcs
    CkPrintf("Setting up communication between gspace <--> paircalc <--> ortho\n");
    CProxy_Ortho ortho = CProxy_Ortho(msg->orthoAID);
    ortho.makeSections(symmCfg, asymmCfg, msg->symmIDs.pcAID, msg->asymmIDs.pcAID);

    /// Send handles to the created chares via a callback
    cb.send(msg);
}

    } // end namespace startup
} // end namespace cp

