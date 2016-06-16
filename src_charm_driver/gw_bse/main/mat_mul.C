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


/* readonly */ CProxy_Controller controller_proxy;

/*readonly*/ int M_stride;
/*readonly*/ int K_stride;
/*readonly*/ int N_stride;

#define NUM_ITERATIONS      4//5
#define MIGRATE_FREQ        10

class MatMul : public CBase_MatMul{
  public:
   MatMul(){
      msg_received = 0;
   }

   void setup(CLA_Matrix_interface matA, CLA_Matrix_interface matB, CLA_Matrix_interface matC, CProxy_EpsMatrix2D aproxy, CProxy_EpsMatrix2D bproxy, CProxy_EpsMatrix2D cproxy,
        int dimension, int rows,
        bool insertA, bool insertB, bool insertC){
//      CLA_Matrix_interface matA, matB, matC;
      //matProxy = thishandle;
      int algorithm = MM_ALG_2D;
        M = K = N = dimension;//10;
        m = k = n = rows;//2;
      M_chunks = (M + m - 1) / m;
      K_chunks = (K + k - 1) / k;
      N_chunks = (N + n - 1) / n;
      iteration = 0;
      msg_received = 0;

      M_stride = 1; //2;
      K_stride = 1; //3;
      N_stride = 1; //5;

//      CkPrintf("\nNumber of chares = %d\n",K_chunks);

      
      CkGroupID gid = CProxy_CkMulticastMgr::ckNew();

      /* create matrix multiplication objects */
      CkCallback *cb = new CkCallback(CkIndex_MatMul::chunk_inited(), thisProxy);

      make_multiplier(&matA, &matB, &matC,
      aproxy, bproxy, cproxy,
       M, K, N, m, k, n, M_stride, K_stride, N_stride, *cb, *cb, *cb, gid,
       algorithm);

  /* A Matrix */
      if(insertA){
        for(int i = 0; i < M_stride * M_chunks; i += M_stride){
          for(int j = 0; j < N_stride * N_chunks; j += N_stride){
            aproxy(i, j).insert(matA);
          }
        }
        aproxy.doneInserting();
      } else{
        for(int i = 0; i < M_stride * M_chunks; i += M_stride){
          for(int j = 0; j < N_stride * N_chunks; j += N_stride){
            aproxy(i, j).setI(matA, false);
          }
        }
      }

  /* B Matrix */
      if(insertB){
        for(int i = 0; i < K_stride * K_chunks; i += K_stride){
          for(int j = 0; j < N_stride * N_chunks; j += N_stride){
            bproxy(i, j).insert(matB);
          }
        }
        bproxy.doneInserting();
      } else{
        for(int i = 0; i < K_stride * K_chunks; i += K_stride){
          for(int j = 0; j < N_stride * N_chunks; j += N_stride){
            bproxy(i, j).setI(matB, false);
          }
        }
      }

   /* C Matrix */ 
      if(insertC){
        for(int i = 0; i < M_stride * M_chunks; i += M_stride){
          for(int j = 0; j < N_stride * N_chunks; j += N_stride){
            cproxy(i, j).insert(matC);
          }
        }
        cproxy.doneInserting();
      }else{
        for(int i = 0; i < M_stride * M_chunks; i += M_stride){
          for(int j = 0; j < N_stride * N_chunks; j += N_stride){
            cproxy(i, j).setI(matC, true); //since output, it's values need to be reset
          }
        }
      }
   }
    

   void quiet(){
      CkPrintf("quiecense detected!!!!!!!!\n");
      CkExit();
    }

    void do_multiply(CProxy_EpsMatrix2D aproxy, CProxy_EpsMatrix2D bproxy, CProxy_EpsMatrix2D cproxy){
/*      if(++msg_received != 3)
        return;
*/      msg_received = 0;

      aproxy.multiply(1, 0);
      bproxy.multiply(1, 0);
      cproxy.multiply(1, 0);

    }


/*
    void migration_done(CkReductionMsg *msg){
      delete msg;
      if(++msg_received != 3)
        return;
//      CkPrintf("Migration finished\n");
      msg_received = 0;
      pmatrix2D_proxy.multiply(1, 0);
      pmatrix2D_bproxy.multiply(1, 0);
      pmatrix2D_cproxy.multiply(1, 0);
    }
*/


    void chunk_inited(){
    
    msg_received++;
    if(msg_received==3){
      msg_received = 0;
      controller_proxy.matrixReady();
    }

    }
  private:
    int runs;
    int M, N, K;
    int m, n, k;
    int M_chunks, N_chunks, K_chunks;
    CkCallback *cb_ready, *cb_done;
    int iteration;
    int msg_received;
    double start_t;
};



#include "mat_mul.def.h"
