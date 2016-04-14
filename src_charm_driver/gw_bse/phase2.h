#ifndef PMATRIX_H_PHASE2__
#define PMATRIX_H_PHASE2__

#include <vector>
#include "ckcomplex.h"

#define IDX(r,c) ((r)*num_cols + (c))

class Phase2 : public CBase_Phase2 {
    Phase2_SDAG_CODE
public:
    Phase2();
    Phase2(CkMigrateMessage *msg);
    Phase2(int local_x, int local_y);

    void fftRows(int direction);
private:
    unsigned L; // Number of occupied psis
    unsigned matrix_dimension; // Size of the entire matrix
    unsigned num_rows, num_cols, start_row, start_col; // The shape of our data
    unsigned trans_count, num_chares; // SDAG variables
    std::vector<complex> data;
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

    int arrival_counter;
};

extern /* readonly */ CProxy_Phase2 phase2_proxy;
extern /* readonly */ CProxy_PMatrix pmatrix_proxy;
#endif //PMATRIX_H_PHASE2__
