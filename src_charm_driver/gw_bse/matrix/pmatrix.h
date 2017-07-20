#ifndef PMATRIX_H
#define PMATRIX_H

#include "matrix.h"
#include "pmatrix.decl.h"

#include "mylapack.h"
#include "CLA_Matrix.h"
#include "ckcomplex.h"

class FFTController;

class PMatrix : public CBase_PMatrix {
  PMatrix_SDAG_CODE
  public:
    PMatrix(MatrixConfig config);
    PMatrix(CkMigrateMessage* msg) {}

    void fftRows(int);
    void generateEpsilon(CProxy_EpsMatrix proxy, std::vector<int> accept);
    void applyFs();
    void calc_vcoulb();
    void calc_Eps(Phase3Message* msg);

    void reportPTime();
    void generateEpsilon(std::vector<double> vcoulb, std::vector<int> accept, int inew, int jnew, int size, int max_inew, int max_jnew);   
    void registerTileSections();
 
  private:
    unsigned L; // Number of occupied psis
    unsigned trans_count, num_chares; // SDAG variables
    unsigned completed_chunks;
    int* nfft; // number of fft grids in each direction
    unsigned qindex;
    FFTController* fft_controller;

    unsigned start_index, end_index;
    unsigned chare_chunk;
    unsigned send_count, recv_count;
    unsigned iteration;
    unsigned max_iterations;

    int arrival_counter;
    int receive_counter;

    complex total[144];
    void kqIndex(unsigned, unsigned&, int*);
    complex* umklapp_factor;
    void getUmklappFactor(complex*, int[3]);

    double total_time;
};

extern /* readonly */ CProxy_PMatrix pmatrix2D_proxy;
extern /* readonly */ CProxy_PMatrix pmatrix1D_proxy;
#endif
