#include "include/mylapack.h"

// from ckPairCalculator.C source code
// GEMM: matrix matrix multiplication
void myGEMM(char *opA, char *opB, int *m, int *n, int *k, double *alpha, complex *A, int *lda, complex *B, int *ldb, double *beta, complex *C, int *ldc){
  complex cAlpha(*alpha,0.), cBeta(*beta,0.);
  ZGEMM(opA, opB, m, n, k, &cAlpha, A, lda, B, ldb, &cBeta, C, ldc);
}


// HEEV: eigenvalues and eigenvectors
void myHEEV(char *jobz, char *uplo, int *n, complex *A, int *lda, double *W, complex *WORK, int *lwork, double *RWORK, int *info){
  ZHEEV(jobz, uplo, n, A, lda, W, WORK, lwork, RWORK, info);
}
