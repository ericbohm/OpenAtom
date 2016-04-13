#include "phase2.decl.h"

#include "standard_include.h"
#include "allclass_gwbse.h"
#include "pmatrix.h"
#include "messages.h"
#include "controller.h"
#include "states.h"
#include "fft_routines.h"
#include "fft_controller.h"

#include "phase2.h"


Phase2::Phase2() : Phase2(10,10)
{
    CkPrintf("Default constructor of Phase 2 called!\n");
    CkExit();
}

Phase2::Phase2(int local_size_x, int local_size_y) : data(local_size_x*local_size_y)
{
    GWBSE* gwbse = GWBSE::get();

    // Set some constants
    L = gwbse->gw_parallel.L;
    nfft = gwbse->gw_parallel.fft_nelems;
    qindex = Q_IDX; // Eventually the controller will set this

    // Grab a local pointer to the fft controller for fft-ing our rows
    // TODO: Is this guaranteed to be safe (is the local branch created for sure)?
    fft_controller = fft_controller_proxy.ckLocalBranch();

    // Figure out what part of the matrix we have and allocate space
    matrix_dimension = gwbse->gw_parallel.n_elems;
    num_rows = gwbse->gw_parallel.rows_per_chare;
    num_cols = gwbse->gw_parallel.cols_per_chare;
    num_chares = (matrix_dimension / num_rows) * (matrix_dimension / num_cols);

    // Do the same for 1D:
    local_mtx_size_1d_x = local_size_x;
    local_mtx_size_1d_y = local_size_y;

    local_mtx_size_2d_x = num_cols;
    local_mtx_size_2d_y = num_rows;

    number_of_chares_1d = matrix_dimension / local_mtx_size_1d_y;
    number_of_chares_2d_x = matrix_dimension / local_mtx_size_2d_x;
    number_of_chares_2d_y = matrix_dimension / local_mtx_size_2d_y;
}

Phase2::Phase2(CkMigrateMessage* msg) {}

void Phase2::fftRows(int direction) {
  // FFT each row stored in this chare
  for (int i=0; i < local_mtx_size_1d_y; i++){
    // First set up the data structures in the FFTController
    fft_controller->setup_fftw_3d(nfft, direction);
    fftw_complex* in_pointer = fft_controller->get_in_pointer();
    fftw_complex* out_pointer = fft_controller->get_out_pointer();

    // Pack our data, do the fft, then get the output
    put_into_fftbox(nfft, &data[IDX(i,0)], in_pointer);
    fft_controller->do_fftw();
    fftbox_to_array(num_cols, out_pointer, &data[IDX(i,0)], 1);
  }

  // Let the controller know we have completed the fft
  contribute(CkCallback(CkReductionTarget(Controller, fftComplete), controller_proxy));
}

#include "phase2.def.h"