#include "gw_bse.decl.h"

#include "ckmulticast.h"

class GWBSEDriver : public CBase_GWBSEDriver {
  public:
    GWBSEDriver(CkArgMsg* msg);

    void readState();
};

class Psi : public CBase_Psi {
  Psi_SDAG_CODE
  public:
    Psi(bool occupied);
    Psi(CkMigrateMessage* msg) {}
    void createSections();
    void calculatePsiR();

  private:
    bool occupied;
    unsigned l, section_index;
    double* psi;

    std::vector<std::pair<unsigned, CProxySection_FCalculator> > sections;
};

class PlaneMsg : public CkMcastBaseMsg, public CMessage_PlaneMsg {
  public:
    unsigned plane_index;
};

class PsiMsg : public CkMcastBaseMsg, public CMessage_PsiMsg {
  public:
    double* psi;
};

class FCalculator : public CBase_FCalculator {
  FCalculator_SDAG_CODE
  public:
    FCalculator();
    FCalculator(CkMigrateMessage* msg) {}

    void createSections();
    void computeF(double*, double*);

  private:
    unsigned l;
    double* f;

    CkMulticastMgr* mcast_ptr;
    CkSectionInfo sec_info;
    CProxySection_Psi bcast_section;
};

class PMatrix : public CBase_PMatrix {
  public:
    PMatrix() {}
    PMatrix(CkMigrateMessage* msg) {}
};
