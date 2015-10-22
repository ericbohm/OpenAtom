#include "/sw/include/fftw3.h"
#include "util.h"

void setup_fftw_3d(int [3], int);
void destroy_fftw_stuff();
void do_fftw();

// global variables
extern fftw_complex *in_pointer;
extern fftw_complex *out_pointer;
extern fftw_plan my_plan;
