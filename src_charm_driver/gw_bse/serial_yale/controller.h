#ifndef __CONTROLLER_H__
#define __CONTROLLER_H__

#include "ckcomplex.h"

#include "controller.decl.h"

class Controller : public CBase_Controller {
  Controller_SDAG_CODE
  public:
    Controller();
    void pComplete(); // Just for debugging

  private:
    unsigned K, L, M, pipeline_stages;
    unsigned next_K, next_state, total_sent, total_complete;
};

class PsiCache : public CBase_PsiCache {
  public:
    PsiCache();

    void receivePsi(PsiMessage*);
    complex* getPsi(unsigned, unsigned, unsigned, unsigned) const;
  private:
    unsigned psi_count, psi_size, received_psis;
    complex** psis;
};

extern /* readonly */ CProxy_Controller controller_proxy;
extern /* readonly */ CProxy_PsiCache psi_cache_proxy;

#endif
