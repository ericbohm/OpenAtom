#ifndef __CONTROLLER_H__
#define __CONTROLLER_H__

#include <cstdlib>
#include "ckcomplex.h"

#include "controller.decl.h"

class Controller : public CBase_Controller {
  Controller_SDAG_CODE
  public:
    Controller();

  private:
    bool do_output;
    double start, end;
    unsigned K, L, M, pipeline_stages;
    unsigned next_K, next_state, total_sent, total_complete;
};

struct FComputePacket {
  unsigned int ikq;
  unsigned int size;
  complex* unocc_psi;
  complex* umklapp_factor;
  complex** fs;
  complex*** occ_psis;
};

class PsiCache : public CBase_PsiCache {
  public:
    PsiCache();

    void receivePsi(PsiMessage*);
    void computeFs(PsiMessage*);
    complex* getPsi(unsigned, unsigned, unsigned) const;
    complex* getF(unsigned) const;
  private:
    void kqIndex(unsigned, unsigned&, int*);
    void computeUmklappFactor(int*);

    // Used for CkLoop parameters
    FComputePacket f_packet;

    unsigned K, L, psi_size, received_psis, qindex;
    complex*** psis;
    complex** fs;
    complex* umklapp_factor;
};

extern /* readonly */ CProxy_Controller controller_proxy;
extern /* readonly */ CProxy_PsiCache psi_cache_proxy;

#endif
