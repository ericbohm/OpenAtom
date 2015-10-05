#ifndef GW_BSE_H
#define GW_BSE_H

#include "gw_bse.decl.h"

class GWDriver : public CBase_GWDriver {
  public:
    GWDriver(CkArgMsg* msg);

    void readConfig();
    void readState();
};

/*readonly*/extern GWConfig config;
/*readonly*/extern CProxy_PsiCache psicache;
/*readonly*/extern CProxy_Psi psi;
/*readonly*/extern CProxy_PMatrix pmatrix;

#endif
