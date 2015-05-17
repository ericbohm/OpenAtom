#ifndef GW_BSE_H
#define GW_BSE_H

#include "gw_bse.decl.h"

class GWDriver : public CBase_GWDriver {
  public:
    GWDriver(CkArgMsg* msg);

    void readConfig();
    void readState();
};

class RowMessage : public CMessage_RowMessage {
  public:
    RowMessage(unsigned i, unsigned s) : index(i), size(s) {}
    unsigned index, size;
    double* row;
};

class PMatrix : public CBase_PMatrix {
  public:
    PMatrix() {}
    PMatrix(CkMigrateMessage* msg) {}

    void receiveRowContribution(RowMessage* msg) { delete msg; }
};

/*readonly*/extern GWConfig config;
/*readonly*/extern CkGroupID mcast_ID;
/*readonly*/extern CProxy_Psi kpsi;
/*readonly*/extern CProxy_Psi qpsi;
/*readonly*/extern CProxy_FCalculator fcalc;
/*readonly*/extern CProxy_PMatrix pmatrix;

#endif
