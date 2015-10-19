#ifndef PMATRIX_H
#define PMATRIX_H

#include "states.decl.h"

class PMatrix : public CBase_PMatrix {
  PMatrix_SDAG_CODE
  public:
    PMatrix();
    PMatrix(CkMigrateMessage* msg) {}

    void receivePsi(PsiMessage*);

  private:
    unsigned num_rows, num_cols, start_row, start_col, done_count;
    double** data;
};

#endif
