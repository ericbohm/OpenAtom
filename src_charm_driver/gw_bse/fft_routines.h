#include <cstdlib>
#include "ckcomplex.h"
#include "fftw3.h"

void set_radix(int nrad_in,int *nrad_ret, int *krad);
void set_fftsize(int, int, int*, int*, int*, int*);

void gidx_to_fftidx(int, int**, int [3], int**);
void put_into_fftbox(int, complex*, int**, int [3], fftw_complex*, bool);
void put_into_fftbox(int [3], complex*, fftw_complex*);
void fftbox_to_array(int, fftw_complex*, complex*, double);
void fftidx_to_gidx(int*, int*, int*, int[3]);

void setup_fftw_3d(int [3], int);
void destroy_fftw_stuff();
void do_fftw();

// global variables
extern fftw_complex *in_pointer;
extern fftw_complex *out_pointer;
extern fftw_plan my_plan;
