#include <cuda_runtime.h>
#include <cublas_v2.h>

#include <stdlib.h>
#include <iostream>

void dgemm_cuda(char* opA, char* opB, int *m, int *n, int *k,
        double* alpha, double* A, int* lda, double* B, int* ldb,
        double* beta, double* C, int* ldc) {
    double *d_A = 0;
    double *d_B = 0;
    double *d_C = 0;

    cublasStatus_t status;
    cublasHandle_t handle;

    cublasOperation_t d_opA;
    if (*opA == 'N' || *opA == 'n') {
        d_opA = CUBLAS_OP_N;
    } else {
        d_opA = CUBLAS_OP_T;
    }
    cublasOperation_t d_opB;
    if (*opB == 'N' || *opB == 'n') {
        d_opB = CUBLAS_OP_N;
    } else {
        d_opB = CUBLAS_OP_T;
    }

    status = cublasCreate(&handle);
    if (status != CUBLAS_STATUS_SUCCESS) {
        //ckout << "Failed to create handle\n";
        return;
    }

	if (cudaMalloc((void**)&d_A, *m * *k * sizeof(d_A[0])) != cudaSuccess) {
        //ckout << "Failed to malloc d_A\n";
        return;
    }
    if (cudaMalloc((void**)&d_B, *k * *n * sizeof(d_B[0])) != cudaSuccess) {
        //ckout << "Failed to malloc d_B\n";
        return;
    }
    if (cudaMalloc((void**)&d_C, *m * *n * sizeof(d_C[0])) != cudaSuccess) {
        //ckout << "Failed to malloc d_C\n";
        return;
    }

    status = cublasSetVector(*m * *k, sizeof(A[0]), A, 1, d_A, 1);
    if (status != CUBLAS_STATUS_SUCCESS) {
        //ckout << "Failed to set d_A\n";
        return;
    }
    status = cublasSetVector(*k * *n, sizeof(B[0]), B, 1, d_B, 1);
    if (status != CUBLAS_STATUS_SUCCESS) {
        //ckout << "Failed to set d_B\n";
        return;
    }
    status = cublasSetVector(*m * *n, sizeof(C[0]), C, 1, d_C, 1);
    if (status != CUBLAS_STATUS_SUCCESS) {
        //ckout << "Failed to set d_C\n";
        return;
    }

    status = cublasDgemm(handle, d_opA, d_opB, *m, *n, *k, alpha, d_A, *lda, d_B, *ldb, beta, d_C, *ldc);
    if (status != CUBLAS_STATUS_SUCCESS) {
        std::cout << "Failed to do DGEMM\n";
        return;
    }

    status = cublasGetVector(*m * *n, sizeof(C[0]), d_C, 1, C, 1);
    if (status != CUBLAS_STATUS_SUCCESS) {
        //ckout << "Failed to get result\n";
        return;
    }

    if (cudaFree(d_A) != cudaSuccess) {
        //ckout << "Failed to free d_A\n";
        return;
    }
    if (cudaFree(d_B) != cudaSuccess) {
        //ckout << "Failed to free d_B\n";
        return;
    }
    if (cudaFree(d_C) != cudaSuccess) {
        //ckout << "Failed to free d_C\n";
        return;
    }

    status = cublasDestroy(handle);
    if (status != CUBLAS_STATUS_SUCCESS) {
        //ckout << "Failed to destroy handle\n";
        return;
    }
}
