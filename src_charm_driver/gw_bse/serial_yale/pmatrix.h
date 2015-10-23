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
    unsigned L; // Number of occupied psis
    unsigned num_rows, num_cols, start_row, start_col; // The shape of our data
    complex** data;

    unsigned qindex;

    void kqIndex(unsigned, unsigned&, int*);
};

extern /* readonly */ CProxy_PMatrix pmatrix_proxy;

#endif
