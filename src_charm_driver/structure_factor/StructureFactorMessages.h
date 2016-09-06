#include "structureFactor.decl.h"
#include "structureFactorCache.decl.h"

#ifndef STRUCTURE_FACTOR_MESSAGES_H
#define STRUCTURE_FACTOR_MESSAGES_H


class SFDummyMsg: public CkMcastBaseMsg, public CMessage_SFDummyMsg  {
  public:
    int iteration_src;
};



class PPDummyMsg: public CMessage_PPDummyMsg {
  public:
    int atmGrp;
    int sfindex;
};



class StructFactorMsg: public CkMcastBaseMsg, public CMessage_StructFactorMsg {
  public:
    int datalen;
    int atmGrpIndex;
    int gsSize;
    int planeIndex; 
    complex *structFactor; 
    complex *structFactor_fx, *structFactor_fy, *structFactor_fz;
};

#endif // STRUCTURE_FACTOR_MESSAGES_H

