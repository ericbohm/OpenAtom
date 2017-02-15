#ifndef PMATRIX_H
#define PMATRIX_H

#include "pmatrix.decl.h"

#include "mylapack.h"
#include "CLA_Matrix.h"
#include "ckcomplex.h"

class FFTController;

class PMatrix2D : public CBase_PMatrix2D {
  PMatrix2D_SDAG_CODE
  public:
    PMatrix2D();
    PMatrix2D(CkMigrateMessage* msg) {}

    void checkReady();
    void applyFs();
    void sendTo1D();
    void sendTo1D_tmp();
    void receiveChunk(Phase2Message*);
    void calc_vcoulb();
    void calc_Eps(Phase3Message* msg);

    void reportPTime();
    void generateEpsilon(std::vector<double> vcoulb, std::vector<int> accept, int inew, int jnew, int size, int max_inew, int max_jnew);   
 
    virtual void pup(PUP::er &p){
      CBase_PMatrix2D::pup(p);
      p | num_rows;
      p | num_cols;
      if(p.isUnpacking())
        data = new complex[num_rows * num_cols];
      PUParray(p, data, num_cols * num_rows);
    }
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
    complex total[144];
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
    void receiveRow_tmp(Phase2Message*);
    void generateEpsilon(std::vector<int> accept);

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

extern /* readonly */ CProxy_EpsMatrix2D eps_matrix2D_proxy;
#endif