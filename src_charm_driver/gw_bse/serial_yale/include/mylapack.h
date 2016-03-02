// lapack interface
#include <cstdlib>
#include "ckcomplex.h"

#ifdef FORTRANUNDERSCORE
#define ZGEMM zgemm_
#define ZHEEV zheev_
#else
#define ZGEMM zgemm
#define ZHEEV zheev
#endif

// ZGEMM: complex matrix-matrix multiplication
extern "C" {
  void ZGEMM (char *, char *, int *, int *, int *,complex *,complex *, int *, complex *, int *, complex *, complex *, int * );
}


void myGEMM (char *opA, char *opB, int *m, int *n, int *k, double *alpha, complex *A, int *lda, complex *B, int *ldb, double *beta, complex *C, int *ldc);




// ZHEEV: computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices
extern "C" {
  void ZHEEV (char *, char *, int *, complex *, int *, double *, complex *, int *, double *, int *);
}

void myHEEV(char *jobz, char *uplo, int *n, complex *A, int *lda, double *W, complex *WORK, int *lwork, double *RWORK, int *info);
