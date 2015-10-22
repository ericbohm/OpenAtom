// lapack interface
#include <cstdlib>
#include "ckcomplex.h"

#ifdef FORTRANUNDERSCORE
#define ZGEMM zgemm_
#else
#define ZGEMM zgemm
#endif



extern "C" {
    void ZGEMM (char *, char *, int *, int *, int *,complex *,complex *, int *, complex *, int *, complex *, complex *, int * );
}



void myGEMM(char *opA, char *opB, int *m, int *n, int *k, double *alpha, complex *A, int *lda, complex *B, int *ldb, double *beta, complex *C, int *ldc);
