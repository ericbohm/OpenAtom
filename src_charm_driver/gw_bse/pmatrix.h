#ifndef PMATRIX_H
#define PMATRIX_H

#include "pmatrix.decl.h"

#ifdef FORTRANUNDERSCORE
#define ZGEMM zgemm_
#define ZGERC zgerc_
#else
#define ZGEMM zgemm
#define ZGERC zgerc
#endif

extern "C" {
  void ZGERC (int*, int*, complex*, complex*, int*, complex*, int*, complex*, int*);
  void ZGEMM (char*, char*, int*, int*, int*, complex*, complex*, int*, complex*, int*, complex*, complex*, int*);
}

class FFTController;

class PMatrix : public CBase_PMatrix {
  PMatrix_SDAG_CODE
  public:
    PMatrix();
    PMatrix(CkMigrateMessage* msg) {}

    void applyFs();
    void fftRows(int);
    void printRows(int, const char*);
    void reportPTime();

  private:
    unsigned L; // Number of occupied psis
    unsigned matrix_dimension; // Size of the entire matrix
    unsigned num_rows, num_cols, start_row, start_col; // The shape of our data
    unsigned trans_count, num_chares; // SDAG variables
    complex* data;
    int* nfft; // number of fft grids in each direction
    unsigned qindex;
    FFTController* fft_controller;
    unsigned local_mtx_size_1d_y;
    int receive_counter;

    void kqIndex(unsigned, unsigned&, int*);
    complex* umklapp_factor;
    void getUmklappFactor(complex*, int[3]);

    double total_time;
};

extern /* readonly */ CProxy_PMatrix pmatrix_proxy;

#endif
