// lapack interface
#include <cstdlib>
#include "ckcomplex.h"

#ifdef USE_LAPACK

#ifdef USE_FORTRAN_UNDERSCORE

#define ZGEMM zgemm_
#define ZGERC zgerc_
#define ZHEEV zheev_

#else

#define ZGEMM zgemm
#define ZGERC zgerc
#define ZHEEV zheev

#endif

#else

#error "Currently need LAPACK for CLA_Matrix"

#endif

// ZGEMM: complex matrix-matrix multiplication
extern "C" {
  void ZGEMM (char *, char *, int *, int *, int *,complex *,complex *, int *, complex *, int *, complex *, complex *, int * );
  void ZGERC (int*, int*, complex*, complex*, int*, complex*, int*, complex*, int*);
}


void myGEMM (char *opA, char *opB, int *m, int *n, int *k, double *alpha, complex *A, int *lda, complex *B, int *ldb, double *beta, complex *C, int *ldc);


