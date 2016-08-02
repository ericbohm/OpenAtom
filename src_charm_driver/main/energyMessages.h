#ifndef ENERGYMESSAGES_H
#define ENERGYMESSAGES_H

#include "energyMessages.decl.h"


struct ECookieMsg : public CkMcastBaseMsg, public CMessage_ECookieMsg {
  public:
    int junk;
   void pup(PUP::er &p){
	  CMessage_ECookieMsg::pup(p);
	  p|junk;
   }

};

class EnergyStructMsg : public CkMcastBaseMsg, public CMessage_EnergyStructMsg {
  public:
  EnergyStruct es;
  UberCollection sender;
};

#endif
