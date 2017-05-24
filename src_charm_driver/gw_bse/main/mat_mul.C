#include <iostream>
#include <sstream>
#include <math.h>

#include "pmatrix.h"
#include "eps_matrix.h"
#include "fft_routines.h"
#include "controller.h"
#include "mat_mul.h"
#include "CLA_Matrix.h"
#include "mat_mul.decl.h"
#include "ckmulticast.h"
#include "ckarray.h"
#include "ckcomplex.h"


/*readonly*/ int M_stride;
/*readonly*/ int K_stride;
/*readonly*/ int N_stride;

#define NUM_ITERATIONS      4//5
#define MIGRATE_FREQ        10

class MatMul : public CBase_MatMul{
private:
  int msg_received;

public:
   MatMul() {}

   void setup(CLA_Matrix_interface matA,
              CLA_Matrix_interface matB,
              CLA_Matrix_interface matC,
              CProxy_EpsMatrix aproxy,
              CProxy_EpsMatrix bproxy,
              CProxy_EpsMatrix cproxy,
              int dimension, int rows) {
    int M, N, K;
    int m, n, k;
    int M_chunks, N_chunks, K_chunks;
    int algorithm = MM_ALG_2D;
    M = K = N = dimension;//10;
    m = k = n = rows;//2;
    M_chunks = (M + m - 1) / m;
    K_chunks = (K + k - 1) / k;
    N_chunks = (N + n - 1) / n;
    msg_received = 0;

    M_stride = 1; //2;
    K_stride = 1; //3;
    N_stride = 1; //5;

    CkGroupID gid = CProxy_CkMulticastMgr::ckNew();

    CkCallback cb(CkIndex_MatMul::chunk_inited(), thisProxy);

    /* create matrix multiplication objects */
    make_multiplier(&matA, &matB, &matC, aproxy, bproxy, cproxy,
        M, K, N, m, k, n, M_stride, K_stride, N_stride,
        cb, cb, cb, gid, algorithm);

    /* A Matrix */
    for(int i = 0; i < M_stride * M_chunks; i += M_stride){
      for(int j = 0; j < N_stride * N_chunks; j += N_stride){
        aproxy(i, j).setI(matA, false);
      }
    }

    /* B Matrix */
    for(int i = 0; i < K_stride * K_chunks; i += K_stride){
      for(int j = 0; j < N_stride * N_chunks; j += N_stride){
        bproxy(i, j).setI(matB, false);
      }
    }

    /* C Matrix */
    for(int i = 0; i < M_stride * M_chunks; i += M_stride){
      for(int j = 0; j < N_stride * N_chunks; j += N_stride){
        cproxy(i, j).setI(matC, true); //since output, it's values need to be reset
      }
    }
  }
    

  void quiet(){
    CkPrintf("quiecense detected!!!!!!!!\n");
    CkExit();
  }

  void do_multiply(CProxy_EpsMatrix aproxy, CProxy_EpsMatrix bproxy, CProxy_EpsMatrix cproxy, double alpha) {
    msg_received = 0;

    aproxy.multiply(1.0, 0);
    bproxy.multiply(1.0, 0);
    cproxy.multiply(alpha, 0);
  }

  void chunk_inited() {
    msg_received++;
    if (msg_received==3) {
      msg_received = 0;
      controller_proxy.matrixReady();
    }
  }
};

#include "mat_mul.def.h"
