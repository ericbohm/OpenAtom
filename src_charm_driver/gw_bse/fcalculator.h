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

    void setupSections(); // Setup plane sections for the section reductions

    void computeF(PsiMsg*, PsiMsg*);  // Compute f from the two psis
    void sendPContribution();         // Compute parts of P and send them

  private:
    double* f;

    CkMulticastMgr* mcast_ptr;
    CkSectionInfo sec_info;
};
