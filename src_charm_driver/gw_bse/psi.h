#ifndef PSI_H
#define PSI_H

#include "gw_bse.decl.h"
#include "ckmulticast.h"

class Psi : public CBase_Psi {
  Psi_SDAG_CODE
  public:
    Psi(bool occupied);
    Psi(CkMigrateMessage* msg) {}

  private:
    void setupSections(unsigned); // Create bcast sections for a given q iter
    void calculatePsiR();         // Do an IFFT to take PsiG to PsiR

    // Information about the kind of psi we are
    bool occupied;
    int size;
    double* psi;

    // Indices used for sdag loops
    unsigned q_index, l_index, section_index;

    // Variables for section management
    CkMulticastMgr* mcast_ptr;
    std::vector<std::pair<unsigned, CProxySection_FCalculator> > sections;
};

#endif
