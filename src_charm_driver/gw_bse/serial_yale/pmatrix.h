#ifndef PMATRIX_H
#define PMATRIX_H

#include "pmatrix.decl.h"

class PMatrix : public CBase_PMatrix {
  PMatrix_SDAG_CODE
  public:
    PMatrix();
    PMatrix(CkMigrateMessage* msg) {}

    void receivePsi(PsiMessage*);
    void printRowAndExit(int);

  private:
    // TODO: These will be moved to the parallel controller
    unsigned pipeline_stages, L, M;
    unsigned num_rows, num_cols, start_row, start_col, done_count;
    complex** data;
};

extern /* readonly */ CProxy_PMatrix pmatrix_proxy;

#endif
