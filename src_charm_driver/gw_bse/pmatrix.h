#ifndef PMATRIX_H
#define PMATRIX_H

#include "gw_bse.decl.h"

class RowMessage : public CMessage_RowMessage {
  public:
    RowMessage(unsigned i, unsigned s) : index(i), size(s) {}
    unsigned index, size;
    double* row;
};

class PMatrix : public CBase_PMatrix {
  public:
    PMatrix();
    PMatrix(CkMigrateMessage* msg) {}

    void receiveRowContribution(RowMessage* msg);

  private:
    unsigned num_rows, row_size, start_row;
    double** rows;
};

#endif
