#ifndef __CONTROLLER_H__
#define __CONTROLLER_H__

#include "ckcomplex.h"

#include "controller.decl.h"

class Controller : public CBase_Controller {
  public:
    Controller();
    void cachesFilled();
    void psiComplete();

  private:
    unsigned L, M, pipeline_stages, next_state;
};

class PsiCache : public CBase_PsiCache {
  public:
    PsiCache();

    void receivePsi(PsiMessage*);
    complex* getPsi(unsigned) const;
  private:
    unsigned psi_count, psi_size, received_psis;
    complex** psis;
};

extern /* readonly */ CProxy_Controller controller_proxy;
extern /* readonly */ CProxy_PsiCache psi_cache_proxy;

#endif
