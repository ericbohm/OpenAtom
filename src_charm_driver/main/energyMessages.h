#ifndef ENERGYMESSAGES_H
#define ENERGYMESSAGES_H

#include "energyMessages.decl.h"


class ECookieMsg : public CkMcastBaseMsg, public CMessage_ECookieMsg {
  public:
    int junk;
};

class EnergyStructMsg : public CkMcastBaseMsg, public CMessage_EnergyStructMsg {
  public:
  EnergyStruct es;
  UberCollection sender;
};

#endif
