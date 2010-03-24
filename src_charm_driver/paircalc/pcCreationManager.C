#include "pcCreationManager.h"
#include "pcBuilder.h"
#include "orthog_ctrl/orthoBuilder.h"
#include "orthog_ctrl/orthoConfig.h"
#include "load_balance/IntMap.h"
#include "ortho.decl.h"

namespace cp {
    namespace paircalc {

CreationManager::CreationManager(const pcConfig &_symmCfg, const pcConfig &_asymmCfg, const ortho::orthoConfig &_orthoCfg):
    symmCfg(_symmCfg), asymmCfg(_asymmCfg), orthoCfg(_orthoCfg)
{
    if (symmCfg.orthoGrainSize != asymmCfg.orthoGrainSize)
        CkAbort("Ortho grain size mismatch in supllied configs\n");
}




void CreationManager::build(CkCallback cb, const int boxSize, PeListFactory getPeList, MapType2 *gSpaceMap)
{
    /// Create the message that will hold the handles to created chares
    pcSetupMsg *msg = new pcSetupMsg();
    msg->symmCfg    = symmCfg;
    msg->asymmCfg   = asymmCfg;

    // Create the symmetric (psi) and asymmetric (lambda) paircalc instances
    CkPrintf("\n\nCreating the symmetric (psi) and asymmetric (lambda) paircalculators\n");
    cp::paircalc::Builder symmBuilder(symmCfg), asymmBuilder(asymmCfg);
    msg->symmIDs  = symmBuilder.build (boxSize, getPeList, gSpaceMap);
    msg->asymmIDs = asymmBuilder.build(boxSize, getPeList, gSpaceMap);

    // Spawn the ortho array and its world of chares/classes (CLA_Matrix, OrthoHelper etc.)
    CkPrintf("Creating the ortho array\n");
    cp::ortho::Builder orthoBuilder(orthoCfg);
    msg->orthoAID = orthoBuilder.build(msg->asymmIDs, getPeList);

    // Ask ortho to setup its communication sections of paircalcs
    CkPrintf("Setting up communication between gspace <--> paircalc <--> ortho\n");
    CProxy_Ortho ortho = CProxy_Ortho(msg->orthoAID);
    ortho.makeSections(symmCfg, asymmCfg, msg->symmIDs.pcAID, msg->asymmIDs.pcAID);

    /// Send handles to the created chares via a callback
    cb.send(msg);
}

    } // end namespace paircalc
} // end namespace cp

