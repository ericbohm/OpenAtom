#include "startupMessages.decl.h"

#include "paircalc/pcConfig.h"
#include "paircalc/pcInstanceIDs.h"
#include "charm++.h"

#ifndef STARTUP_MESSAGES_H
#define STARTUP_MESSAGES_H

/// paircalc::CreationManager returns relevant chare array handles via this msg
class pcSetupMsg: public CMessage_pcSetupMsg
{
  public:
    CkArrayID gspAID;
    cp::paircalc::pcConfig symmCfg, asymmCfg;
    cp::paircalc::InstanceIDs symmIDs, asymmIDs;
    CkArrayID orthoAID;
};

#endif // STARTUP_MESSAGES_H

