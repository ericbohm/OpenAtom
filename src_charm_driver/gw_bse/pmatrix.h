#ifndef PMATRIX_H
#define PMATRIX_H

#include "pmatrix.decl.h"

#ifdef FORTRANUNDERSCORE
#define ZGERC zgerc_
#else
#define ZGERC zgerc
#endif

extern "C" {
  void ZGERC (int*, int*, complex*, complex*, int*, complex*, int*, complex*, int*);
}

class FFTController;

class PMatrix : public CBase_PMatrix {
  PMatrix_SDAG_CODE
  public:
    PMatrix();
    PMatrix(CkMigrateMessage* msg) {}

    void receivePsi(PsiMessage*);
    void applyFs(int, int, int, int);
    void fftRows(int);
    void printRows(int, const char*);

  private:
    unsigned L; // Number of occupied psis
    unsigned num_rows, num_cols, start_row, start_col; // The shape of our data
    unsigned trans_count, num_chares; // SDAG variables
    complex* data;
    int* nfft; // number of fft grids in each direction
    unsigned qindex;
    FFTController* fft_controller;

    void kqIndex(unsigned, unsigned&, int*);
    complex* umklapp_factor;
    void getUmklappFactor(complex*, int[3]);
};

extern /* readonly */ CProxy_PMatrix pmatrix_proxy;

#endif
