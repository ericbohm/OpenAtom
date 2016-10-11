// lapack interface
#include <cstdlib>
#include "ckcomplex.h"

#define ZGEMM zgemm_
#define ZHEEV zheev_

// ZGEMM: complex matrix-matrix multiplication
extern "C" {
  void ZGEMM (char *, char *, int *, int *, int *,complex *,complex *, int *, complex *, int *, complex *, complex *, int * );
}


void myGEMM (char *opA, char *opB, int *m, int *n, int *k, double *alpha, complex *A, int *lda, complex *B, int *ldb, double *beta, complex *C, int *ldc);


