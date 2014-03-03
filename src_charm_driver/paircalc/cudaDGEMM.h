#include <cuda_runtime.h>
#include <cublas_v2.h>

void testDGEMM() {
	float* h_A;
	float* h_B;
	float* h_C;
	float* d_A = 0;
	float* d_B = 0;
	float* d_C = 0;
	float alpha = 1.0;
	float beta = 0.0;
	const int N = 100;
	const int n2 = N * N;

	int devCount;
	cudaGetDeviceCount(&devCount);
	if (devCount == 0) {
		std::cout << "Error finding device\n";
		return;
	}
	cudaSetDevice(0);

	cublasStatus_t status;
	cublasHandle_t handle;

	cublasCreate(&handle);
	if (status != CUBLAS_STATUS_SUCCESS) {
		std::cout << "Couldn't create handle\n";
	}

	h_A = (float*)malloc(n2 * sizeof(float));
	h_B = (float*)malloc(n2 * sizeof(float));
	h_C = (float*)malloc(n2 * sizeof(float));
	for (int i = 0; i < n2; i++) {
		h_A[i] = rand() / (float)RAND_MAX;
		h_B[i] = rand() / (float)RAND_MAX;
		h_C[i] = 0.0;
	}
	cudaMalloc((void **)&d_A, n2 * sizeof(float));
	cudaMalloc((void **)&d_B, n2 * sizeof(float));
	cudaMalloc((void **)&d_C, n2 * sizeof(float));
	cublasSetVector(n2, sizeof(h_A[0]), h_A, 1, d_A, 1);
	cublasSetVector(n2, sizeof(h_B[0]), h_B, 1, d_B, 1);
	cublasSetVector(n2, sizeof(h_C[0]), h_C, 1, d_C, 1);

	status = cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, N, N, &alpha, d_A, N, d_B, N, &beta, d_C, N);
	if (status == CUBLAS_STATUS_NOT_INITIALIZED) {
		std::cout << "Not initialized!\n";
	}
	if (status == CUBLAS_STATUS_ALLOC_FAILED) {
		std::cout << "Alloc failed!\n";
	}
	if (status == CUBLAS_STATUS_INVALID_VALUE) {
		std::cout << "Invalid value!\n";
	}
	if (status == CUBLAS_STATUS_ARCH_MISMATCH) {
		std::cout << "Arch mismatch!\n";
	}
	if (status == CUBLAS_STATUS_EXECUTION_FAILED) {
		std::cout << "Execution failed!\n";
	}
	if (status == CUBLAS_STATUS_INTERNAL_ERROR) {
		std::cout << "Internal error!\n";
	}
	if (status != CUBLAS_STATUS_SUCCESS) {
		std::cout << "Error with dummy DGEMM\n";
	}

	cublasGetVector(n2, sizeof(h_C[0]), d_C, 1, h_C, 1);
	free(h_A);
	free(h_B);
	free(h_C);
	cudaFree(d_A);
	cudaFree(d_B);
	cudaFree(d_C);
	cublasDestroy(handle);
	std::cout << "Lovely! We've finished!\n";
}

/*void cudaDGEMM(char* opA, char* opB, int *m, int *n, int *k,
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
		if (*lda != *m) {
			ckout << "SHIT!!!!!\n";
		}
		//ckout << "SHIT!!!!!\n";
	} else {
		d_opA = CUBLAS_OP_T;
		if (*lda != *k) {
			ckout << "SHIT!!!!!\n";
		}
		//ckout << "SHIT!!!!!\n";
	}
	cublasOperation_t d_opB;
	if (*opB == 'N' || *opB == 'n') {
		d_opB = CUBLAS_OP_N;
		if (*ldb != *k) {
			ckout << "SHIT!!!!!\n";
		}
	} else {
		d_opB = CUBLAS_OP_T;
		if (*ldb != *n) {
			ckout << "SHIT!!!!!\n";
		}
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
		//ckout << "Failed to do DGEMM\n";
		return;
	}

	status = cublasGetVector(*ldc, sizeof(C[0]), d_C, 1, C, 1);
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
}*/
