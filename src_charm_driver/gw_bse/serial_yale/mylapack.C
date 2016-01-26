#include "include/mylapack.h"

// from ckPairCalculator.C source code

void myGEMM(char *opA, char *opB, int *m, int *n, int *k, double *alpha, complex *A, int *lda, complex *B, int *ldb, double *beta, complex *C, int *ldc){
    complex cAlpha(*alpha,0.), cBeta(*beta,0.);
    ZGEMM(opA, opB, m, n, k, &cAlpha, A, lda, B, ldb, &cBeta, C, ldc);
}
