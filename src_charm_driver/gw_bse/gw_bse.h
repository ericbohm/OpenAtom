#include "gw_bse.decl.h"

class GWBSEDriver : public CBase_GWBSEDriver {
  public:
    GWBSEDriver(CkArgMsg* msg);

    void readState();
};

class KPsi : public CBase_KPsi {
  KPsi_SDAG_CODE
  public:
    KPsi();
    KPsi(CkMigrateMessage* msg) {}
    void createSections();
    void calculatePsiR();

  private:
    unsigned l;
    double* psi;

    CProxySection_FCalculator section;
};

class QPsi : public CBase_QPsi {
  QPsi_SDAG_CODE
  public:
    QPsi();
    QPsi(CkMigrateMessage* msg) {}

    void createSections();
    void calculatePsiR();

  private:
    unsigned l;
    double* psi;

    std::vector<CProxySection_FCalculator> sections;
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
