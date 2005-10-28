#include "CLA_Matrix.h"

#if 0
#include <sstream>
using std::ostringstream;
using std::endl;
#endif

#ifdef FORTRANUNDERSCORE
#define DGEMM dgemm_
#else
#define DGEMM dgemm
#endif


extern "C" {void  DGEMM(char *, char *, int *, int *, int *, double *,
			   double *, int *, double *, int *, double *,
			   double *, int *);}

#define MULTARG_A	0
#define MULTARG_B	1
#define MULTARG_C	2

/******************************************************************************/
/* helper functions */

/* Should be called by user to create matrices. Documented in header file. */
int make_multiplier(CLA_Matrix_interface *A, CLA_Matrix_interface *B,
 CLA_Matrix_interface *C, CProxy_ArrayElement bindA,
 CProxy_ArrayElement bindB, CProxy_ArrayElement bindC,
 int M, int K, int N, int m, int k, int n, CkCallback cbA, CkCallback cbB,
 CkCallback cbC, int algorithm){
  /* validate arguments */
  if(algorithm < MM_ALG_MIN || MM_ALG_MAX < algorithm)
    return ERR_INVALID_ALG;

  if(m > M || k > K || n > N)
    return ERR_INVALID_DIM;

  /* create arrays */
  CkArrayOptions optsA(0);
  CkArrayOptions optsB(0);
  CkArrayOptions optsC(0);

  optsA.bindTo(bindA);
  optsB.bindTo(bindB);
  optsC.bindTo(bindC);
  CProxy_CLA_Matrix pa = CProxy_CLA_Matrix::ckNew(optsA);
  CProxy_CLA_Matrix pb = CProxy_CLA_Matrix::ckNew(optsB);
  CProxy_CLA_Matrix pc = CProxy_CLA_Matrix::ckNew(optsC);
  A->setProxy(pa);
  B->setProxy(pb);
  C->setProxy(pc);

  /* populate arrays */
  int M_chunks = (M + m - 1) / m; // ceil(1.0 * M / m)
  int K_chunks = (K + k - 1) / k; // ceil(1.0 * K / k)
  int N_chunks = (N + n - 1) / n; // ceil(1.0 * N / n)

  if(algorithm == MM_ALG_2D){
    for(int i = 0; i < M_chunks; i++)
      for(int j = 0; j < K_chunks; j++)
        (A->p(i, j)).insert(M, K, N, m, k, n, MULTARG_A, B->p, C->p, bindA,cbA);
    A->p.doneInserting();

    for(int i = 0; i < K_chunks; i++)
      for(int j = 0; j < N_chunks; j++)
        (B->p(i, j)).insert(M, K, N, m, k, n, MULTARG_B, A->p, C->p, bindB,cbB);
    B->p.doneInserting();

    for(int i = 0; i < M_chunks; i++)
      for(int j = 0; j < N_chunks; j++)
        (C->p(i, j)).insert(M, K, N, m, k, n, MULTARG_C, A->p, B->p, bindC,cbC);
    C->p.doneInserting();
  }

  return 0;
}

/* Transpose data, which has dimension m x n */
void transpose(double *data, int m, int n){
  if(m == n){
    /* transpose square matrix in place */
    for(int i = 0; i < m; i++)
      for(int j = i + 1; j < n; j++){
        double tmp = data[i * n + j];
        data[i * n + j] = data[j * m + i];
        data[j * m + i] = tmp;
      }
  }
  else {
    double *tmp = new double[m * n];
    memcpy(tmp, data, m * n * sizeof(double));
    for(int i = 0; i < m; i++)
      for(int j = 0; j < n; j++)
        data[j * m + i] = tmp[i * n + j];
    delete [] tmp;
  }
}

/******************************************************************************/
/* CLA_Matrix */

CLA_Matrix::CLA_Matrix(int M, int K, int N, int m, int k, int n,
 int part, CProxy_CLA_Matrix other1, CProxy_CLA_Matrix other2,
 CProxy_ArrayElement bound, CkCallback ready){
    /* initialize simple members */
    this->M = M; this->K = K; this->N = N;
    this->um = m; this->uk = k; this->un = n;
    this->part = part;
    this->algorithm = MM_ALG_2D;
    this->other1 = other1; this->other2 = other2;
    this->bound = bound;
    M_chunks = (M + m - 1) / m;
    K_chunks = (K + k - 1) / k;
    N_chunks = (N + n - 1) / n;

    /* figure out size of our sections */
    if(part == MULTARG_A){
      if(thisIndex.x == M_chunks - 1){
        this->m = M % m;
        if(this->m == 0) this->m = m;
      }
      else this->m = m;
      if(thisIndex.y == K_chunks - 1){
        this->k = K % k;
        if(this->k == 0) this->k = k;
      }
      else this->k = k;
      this->n = n;
    } else if(part == MULTARG_B) {
      if(thisIndex.x == K_chunks - 1){
        this->k = K % k;
        if(this->k == 0) this->k = k;
      }
      else this->k = k;
      if(thisIndex.y == N_chunks - 1){
        this->n = N % n;
        if(this->n == 0) this->n = n;
      }
      else this->n = n;
      this->m = m;
    } else if(part == MULTARG_C) {
      if(thisIndex.x == M_chunks - 1){
        this->m = M % m;
        if(this->m == 0) this->m = m;
      }
      else this->m = m;
      if(thisIndex.y == N_chunks - 1){
        this->n = N % n;
        if(this->n == 0) this->n = n;
      }
      else this->n = n;
      this->k = k;
      got_start = false;
      row_count = col_count = 0;
    }

    /* make communication group for A, B, destination arrays for C */
    if(part == MULTARG_A){
      commGroup = CProxySection_CLA_Matrix::ckNew(other2, thisIndex.x,
       thisIndex.x, 1, 0, N_chunks - 1, 1);
      tmpA = tmpB = NULL;
    } else if(part == MULTARG_B) {
      commGroup = CProxySection_CLA_Matrix::ckNew(other2, 0, M_chunks - 1, 1,
       thisIndex.y, thisIndex.y, 1);
      tmpA = tmpB = NULL;
    } else if(part == MULTARG_C) {
      tmpA = new double[this->m * K];
      tmpB = new double[K * this->n];
    }

    /* let the caller know we are ready */
    contribute(0, NULL, CkReduction::sum_int, ready);
}

void CLA_Matrix::multiply(double alpha, double beta, double *data,
 void (*fptr) (void*), void *usr_data){
  /* A and B send out their chunks, ignoring alpha, beta, ftpr, and usr_data */
  if(part == MULTARG_A){
    CLA_Matrix_msg *msg = new (m * k) CLA_Matrix_msg(data, m, k, thisIndex.x,
     thisIndex.y);
    commGroup.receiveA(msg);
  } else if(part == MULTARG_B){
    CLA_Matrix_msg *msg = new (k * n) CLA_Matrix_msg(data, k, n, thisIndex.x,
     thisIndex.y);
    commGroup.receiveB(msg);
  }
  /* C stores the paramters for the multiplication */
  else if(part == MULTARG_C){
    fcb = fptr;
    user_data = usr_data;
    dest = data;
    this->alpha = alpha;
    this->beta = beta;
    got_start = true;
    /* Check if we were slow to arrive */
    if(row_count == K_chunks && col_count == K_chunks)
      multiply();
  }
  else
    CmiAbort("CLA_Matrix internal error");
}

void CLA_Matrix::receiveA(CLA_Matrix_msg *msg){
  /* store current part */
  row_count++;
  for(int i = 0; i < m; i++)
    memcpy(&tmpA[K * i + uk * msg->fromY], &msg->data[i * msg->d2],
     msg->d2 * sizeof(double));
  delete msg;

  /* If we have all the parts, multiply */
  if(row_count == K_chunks && col_count == K_chunks && got_start)
    multiply();
}

void CLA_Matrix::receiveB(CLA_Matrix_msg *msg){
  /* store current part */
  col_count++;
  memcpy(&tmpB[n * uk * msg->fromX], msg->data,
   msg->d1 * msg->d2 * sizeof(double));
  delete msg;

  /* If we have all the parts, multiply */
  if(row_count == K_chunks && col_count == K_chunks && got_start)
    multiply();
}

void CLA_Matrix::multiply(){
  /* reset counters */
  row_count = col_count = 0;
  got_start = false;

  /* transpose result matrix (if beta != 0) */
  if(beta != 0)
    transpose(dest, m, n);

  /* multiply */
  char trans = 't';
  DGEMM(&trans, &trans, &m, &n, &K, &alpha, tmpA, &K, tmpB, &n, &beta,
   dest, &m);
  
  /* transpose result */
  transpose(dest, n, m);

  /* tell caller we are done */
  fcb(user_data);
}

/******************************************************************************/
/* CLA_Matrix_msg */
CLA_Matrix_msg::CLA_Matrix_msg(double *data, int d1, int d2, int fromX,
 int fromY){
  memcpy(this->data, data, d1 * d2 * sizeof(double));
  this->d1 = d1; this->d2 = d2;
  this->fromX = fromX; this->fromY = fromY;
}

#include "CLA_Matrix.def.h"
