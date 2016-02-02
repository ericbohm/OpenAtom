#include "fft_controller.h"
#include "fft_routines.h"

FFTController::FFTController() {
  first_time = true;
  in_pointer = out_pointer = NULL;
}

void FFTController::do_fftw() {
  //setup_fftw_3d(nfft, direction);
  //put_into_fftbox(nfft, in_data, in_pointer);
  fftw_execute(plan);
  //fftbox_to_array(out_size, out_pointer, out_data, 1);
}

void FFTController::destroy_fftw_stuff() {
  fftw_destroy_plan(plan);
  fftw_free(in_pointer);
  fftw_free(out_pointer);
  in_pointer = out_pointer = NULL;
}

void FFTController::setup_fftw_3d(int nfft[3], int direction) {
  // check for some dumb input values
  if (nfft[0]<=0 || nfft[1] <=0 || nfft[2] <=0) {
    CkPrintf("setup_fftw_3d routine received illegal value for nfft. \
              nfft should be positive number.");
    CkExit();
  }
  if (!(direction==-1 || direction==1)) {
    CkPrintf("setup_fftw_3d routine received illegal value for direction. \
              FFTW direction must either 1 or -1");
    CkExit();
  }

  // if not the first time and any parameters (size, direction) mismatch
  // we need to destroy old plans to set up new ones
  if (!first_time && (direction != old_direction ||
      nfft[0] != old_nfft[0] ||
      nfft[1] != old_nfft[1] ||
      nfft[2] != old_nfft[2])) {
    destroy_fftw_stuff();
    first_time = true;
  }

  // if first time, we need to set up
  if(first_time) {
    const int ndata = nfft[0]*nfft[1]*nfft[2];
    in_pointer = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ndata);
    out_pointer = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ndata);
    plan = fftw_plan_dft_3d(nfft[0], nfft[1], nfft[2],
        in_pointer, out_pointer, direction, FFTW_ESTIMATE);
  }

  // now the old value changes to new value
  first_time = false;
  old_nfft[0] = nfft[0];
  old_nfft[1] = nfft[1];
  old_nfft[2] = nfft[2];
  old_direction = direction;
}
