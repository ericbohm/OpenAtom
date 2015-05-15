#include "gw_bse.decl.h"

#include "ckmulticast.h"

class GWDriver : public CBase_GWDriver {
  public:
    GWDriver(CkArgMsg* msg);

    void readConfig();
    void readState();
};

class Psi : public CBase_Psi {
  Psi_SDAG_CODE
  public:
    Psi(bool occupied);
    Psi(CkMigrateMessage* msg) {}

    void setupSections(); // Create sections for broadcasting psi to fcalc
    void calculatePsiR(); // Do an IFFT to take PsiG to PsiR

  private:
    // Information about the kind of psi we are
    bool occupied;
    int size;
    double* psi;

    // Indices used for sdag loops
    unsigned l_index, section_index;

    // Variables for section management
    CkMulticastMgr* mcast_ptr;
    std::vector<std::pair<unsigned, CProxySection_FCalculator> > sections;
};

class PlaneMsg : public CkMcastBaseMsg, public CMessage_PlaneMsg {
  public:
    PlaneMsg(unsigned idx) : plane_index(idx) {}
    unsigned plane_index;
};

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

class PMatrix : public CBase_PMatrix {
  public:
    PMatrix() {}
    PMatrix(CkMigrateMessage* msg) {}
};
