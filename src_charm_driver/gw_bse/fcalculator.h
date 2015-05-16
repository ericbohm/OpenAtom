#ifndef FCALCULATOR_H
#define FCALCULATOR_H

#include "gw_bse.decl.h"
#include "ckmulticast.h"

// Message used in setting up section reductions
class PlaneMsg : public CkMcastBaseMsg, public CMessage_PlaneMsg {
  public:
    PlaneMsg(unsigned idx) : plane_index(idx) {}
    unsigned plane_index;
};

// Message sent from psi used to compute f
class PsiMsg : public CkMcastBaseMsg, public CMessage_PsiMsg {
  public:
    PsiMsg(unsigned s, double* p) : size(s) {
      std::copy(p, p+size, psi);
    }
    unsigned size;
    double* psi;
};

class FCalculator : public CBase_FCalculator {
  FCalculator_SDAG_CODE
  public:
    FCalculator();
    FCalculator(CkMigrateMessage* msg) {}

  private:
    // Setup an l-plane section and broadcast with it to get section info
    void setupSections() const;

    // Transform a row index using chare index to stagger sends to PMatrix
    unsigned rowOffset(unsigned) const;

    // Computations for f and p
    void computeF(PsiMsg*, PsiMsg*);
    void computeAndSendRow(unsigned row) const;

    // Size of f and array pointer for f
    unsigned size;
    double* f;

    // Section reduction variables
    CkMulticastMgr* mcast_ptr;
    CkSectionInfo sec_info;
};

#endif
