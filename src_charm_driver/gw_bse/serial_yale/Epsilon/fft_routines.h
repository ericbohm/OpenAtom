#include <cstdlib>
#include "include/ckcomplex.h"
#include "/sw/include/fftw3.h"

void gidx_to_fftidx(int, int**, int [3], int**);
void put_into_fftbox(int, complex*, int**, int [3], fftw_complex*);
void fftbox_to_array(int, fftw_complex*, complex*, double);
void fftidx_to_gidx(int*, int*, int*, int[3]);

