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

class PMatrix2D : public CBase_PMatrix2D {
  PMatrix2D_SDAG_CODE
  public:
    PMatrix2D();
    PMatrix2D(CkMigrateMessage* msg) {}

    void applyFs();
    void sendTo1D();
    void receiveChunk(Phase2Message*);
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

class PMatrix1D : public CBase_PMatrix1D {
  PMatrix1D_SDAG_CODE
  public:
    PMatrix1D(int, int);
    PMatrix1D(CkMigrateMessage* msg) {}

    void fftRows(int);
    void sendTo2D();
    void receiveRow(Phase2Message*);

  private:
    unsigned L; // Number of occupied psis
    unsigned matrix_dimension; // Size of the entire matrix
    unsigned num_rows, num_cols, start_row, start_col; // The shape of our data
    unsigned trans_count, num_chares; // SDAG variables
    complex* data;
    int* nfft; // number of fft grids in each direction
    unsigned qindex;
    FFTController* fft_controller;

    unsigned local_mtx_size_1d_x;
    unsigned local_mtx_size_1d_y;

    unsigned local_mtx_size_2d_x;
    unsigned local_mtx_size_2d_y;

    unsigned number_of_chares_1d;
    unsigned number_of_chares_2d_x;
    unsigned number_of_chares_2d_y;

    unsigned start_index, end_index;
    unsigned chare_chunk;
    unsigned send_count, recv_count;
    unsigned iteration;
    unsigned max_iterations;

    int arrival_counter;
};

extern /* readonly */ CProxy_PMatrix2D pmatrix2D_proxy;
extern /* readonly */ CProxy_PMatrix1D pmatrix1D_proxy;

#endif
