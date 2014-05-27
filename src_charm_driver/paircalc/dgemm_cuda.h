void dgemm_cuda(char* opA, char* opB, int *m, int *n, int *k, 
        double* alpha, double* A, int* lda, double* B, int* ldb,
        double* beta, double* C, int* ldc);
