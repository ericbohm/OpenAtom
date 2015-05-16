#include "gw_bse.decl.h"
#include "ckmulticast.h"

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
