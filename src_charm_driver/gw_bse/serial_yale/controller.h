#ifndef __CONTROLLER_H__
#define __CONTROLLER_H__

#include "ckcomplex.h"

#include "controller.decl.h"

class Controller : public CBase_Controller {
  public:
    Controller() { CkPrintf("HI! %d\n", CkMyPe()); }
};

class PsiCache : public CBase_PsiCache {
  public:
    PsiCache();

    void receivePsi(PsiMessage*);
    complex* getPsi(unsigned) const;
  private:
    unsigned psi_count, psi_size, received_psis, pipeline_stages;
    complex** psis;
};

extern /* readonly */ CProxy_Controller controller_proxy;
extern /* readonly */ CProxy_PsiCache psi_cache_proxy;

#endif
