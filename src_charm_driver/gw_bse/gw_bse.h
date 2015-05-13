#include "gw_bse.decl.h"

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

class FCalculator : public CBase_FCalculator {
  FCalculator_SDAG_CODE
  public:
    FCalculator();
    FCalculator(CkMigrateMessage* msg) {}

    void computeF(double*, double*);

  private:
    unsigned l;
    double* f;
};

class PMatrix : public CBase_PMatrix {
  public:
    PMatrix() {}
    PMatrix(CkMigrateMessage* msg) {}
};
