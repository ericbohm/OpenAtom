#ifndef PSI_H
#define PSI_H

#include "gw_bse.decl.h"

// Message sent from psi used to compute f. Some of these are cached on each
// node, and others are streamed in to the PMatrix as needed.
class PsiMessage : public CMessage_PsiMessage {
  public:
    PsiMessage(unsigned s, double* p) : size(s) {
      std::copy(p, p+size, psi);
    }
    unsigned k_index, state_index, size;
    double* psi;
};

class Psi : public CBase_Psi {
  public:
    Psi(unsigned);
    Psi(CkMigrateMessage* msg) {}

    void sendToP();

  private:
    void calculatePsiR();
    void sendToCache();

    // Information about the kind of psi we are
    int size;
    double* psi;
};

class PsiCache: public CBase_PsiCache {
  public:
    PsiCache(unsigned, unsigned);

    void receivePsi(PsiMessage*);
    double* getPsi(unsigned) const;
  private:
    unsigned psi_count, psi_size, received_psis;
    double** psis;
};

#endif
