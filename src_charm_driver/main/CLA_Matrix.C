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
//#define PRINT_DGEMM_PARAMS

extern "C" {void  DGEMM(char *, char *, int *, int *, int *, double *,
			   double *, int *, double *, int *, double *,
			   double *, int *);}

#define MULTARG_A	0
#define MULTARG_B	1
#define MULTARG_C	2
extern CkReduction::reducerType sumFastDoubleType;
/******************************************************************************/
/* helper functions */

/* Should be called by user to create matrices. Documented in header file. */
int make_multiplier(CLA_Matrix_interface *A, CLA_Matrix_interface *B,
 CLA_Matrix_interface *C, CProxy_ArrayElement bindA,
 CProxy_ArrayElement bindB, CProxy_ArrayElement bindC,
 int M, int K, int N, int m, int k, int n, int strideM, int strideK,
 int strideN, CkCallback cbA, CkCallback cbB,
 CkCallback cbC, CkGroupID gid, int algorithm, int gemmSplitOrtho){
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
  int M_chunks = (M + m - 1) / m; // same as ceil(1.0 * M / m)
  int K_chunks = (K + k - 1) / k; // same as ceil(1.0 * K / k)
  int N_chunks = (N + n - 1) / n; // same as ceil(1.0 * N / n)

  if(algorithm == MM_ALG_2D){
    for(int i = 0; i < M_chunks; i++)
      for(int j = 0; j < K_chunks; j++)
        (A->p(i * strideM, j * strideK)).insert(M, K, N, m, k, n, strideM,
         strideK, strideN, MULTARG_A, B->p, C->p, cbA, gemmSplitOrtho);
    A->p.doneInserting();

    for(int i = 0; i < K_chunks; i++)
      for(int j = 0; j < N_chunks; j++)
        (B->p(i * strideK, j * strideN)).insert(M, K, N, m, k, n, strideM,
         strideK, strideN, MULTARG_B, A->p, C->p, cbB, gemmSplitOrtho);
    B->p.doneInserting();

    for(int i = 0; i < M_chunks; i++)
      for(int j = 0; j < N_chunks; j++)
        (C->p(i * strideM, j * strideN)).insert(M, K, N, m, k, n, strideM,
         strideK, strideN, MULTARG_C, A->p, B->p, cbC, gemmSplitOrtho);
    C->p.doneInserting();
  }
  else if(algorithm == MM_ALG_3D){
    CProxy_CLA_MM3D_multiplier mult = CProxy_CLA_MM3D_multiplier::ckNew();
    int curpe = 0;
    int totpe = CkNumPes();
    for(int i = 0; i < M_chunks; i++){
      int mm = m;
      if(i == M_chunks - 1) {
        mm = M % m;
        if(mm == 0)
          mm = m;
      }
      for(int j = 0; j < N_chunks; j++){
        int nn = n;
        if(j == N_chunks - 1) {
          nn = N % n;
          if(nn == 0)
            nn = n;
        }
        for(int l = 0; l < K_chunks; l++){
          int kk = k;
          if(l == K_chunks - 1) {
            kk = K % k;
            if(kk == 0)
              kk = k;
          }
          mult(i, j, l).insert(mm, kk, nn, curpe);
          curpe = (curpe + 1) % totpe;
        }
      }
    }
    mult.doneInserting();

    for(int i = 0; i < M_chunks; i++)
      for(int j = 0; j < K_chunks; j++)
        (A->p(i * strideM, j * strideK)).insert(mult, M, K, N, m, k, n,
         strideM, strideK, strideN, MULTARG_A, cbA, gid, gemmSplitOrtho);
    A->p.doneInserting();

    for(int i = 0; i < K_chunks; i++)
      for(int j = 0; j < N_chunks; j++)
        (B->p(i * strideK, j * strideN)).insert(mult, M, K, N, m, k, n,
         strideM, strideK, strideN, MULTARG_B, cbB, gid, gemmSplitOrtho);
    B->p.doneInserting();

    for(int i = 0; i < M_chunks; i++)
      for(int j = 0; j < N_chunks; j++)
        (C->p(i * strideM, j * strideN)).insert(mult, M, K, N, m, k, n,
         strideM, strideK, strideN, MULTARG_C, cbC, gid, gemmSplitOrtho);
    C->p.doneInserting();
  }

  return SUCCESS;
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

/* constructor for 2D algorithm */
CLA_Matrix::CLA_Matrix(int M, int K, int N, int m, int k, int n,
 int strideM, int strideK, int strideN, int part,
 CProxy_CLA_Matrix other1, CProxy_CLA_Matrix other2, CkCallback ready, int _gemmSplitOrtho){
  /* initialize simple members */
  this->M = M; this->K = K; this->N = N;
  this->um = m; this->uk = k; this->un = n;
  this->part = part;
  this->algorithm = MM_ALG_2D;
  this->other1 = other1; this->other2 = other2;
  this->M_stride = strideM;
  this->K_stride = strideK;
  this->N_stride = strideN;
  gemmSplitOrtho=_gemmSplitOrtho;
  M_chunks = (M + m - 1) / m;
  K_chunks = (K + k - 1) / k;
  N_chunks = (N + n - 1) / n;
  algorithm = MM_ALG_2D;
  usesAtSync = CmiFalse;
  setMigratable(false);
  /* figure out size of our sections */
  if(part == MULTARG_A){
    if(thisIndex.x == (M_chunks - 1) * strideM){
      this->m = M % m;
      if(this->m == 0) this->m = m;
    }
    else this->m = m;
    if(thisIndex.y == (K_chunks - 1) * strideK){
      this->k = K % k;
      if(this->k == 0) this->k = k;
    }
    else this->k = k;
    this->n = n;
  } else if(part == MULTARG_B) {
    if(thisIndex.x == (K_chunks - 1) * strideK){
      this->k = K % k;
      if(this->k == 0) this->k = k;
    }
    else this->k = k;
    if(thisIndex.y == (N_chunks - 1) * strideN){
      this->n = N % n;
      if(this->n == 0) this->n = n;
    }
    else this->n = n;
    this->m = m;
  } else if(part == MULTARG_C) {
    if(thisIndex.x == (M_chunks - 1) * strideM){
      this->m = M % m;
      if(this->m == 0) this->m = m;
    }
    else this->m = m;
    if(thisIndex.y == (N_chunks - 1) * strideN){
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
    commGroup2D = CProxySection_CLA_Matrix::ckNew(other2, thisIndex.x,
     thisIndex.x, 1, 0, (N_chunks - 1) * strideN, strideN);
    tmpA = tmpB = NULL;
  } else if(part == MULTARG_B) {
    commGroup2D = CProxySection_CLA_Matrix::ckNew(other2, 0,
     (M_chunks - 1) * strideM, strideM, thisIndex.y, thisIndex.y, 1);
    tmpA = tmpB = NULL;
  } else if(part == MULTARG_C) {
    tmpA = new double[this->m * K];
    tmpB = new double[K * this->n];
  }

  /* let the caller know we are ready */
  contribute(0, NULL, CkReduction::sum_int, ready);
}

/* constructor for 3D algorithm */
CLA_Matrix::CLA_Matrix(CProxy_CLA_MM3D_multiplier p, int M, int K, int N,
 int m, int k, int n, int strideM, int strideK, int strideN, int part,
 CkCallback cb, CkGroupID gid, int _gemmSplitOrtho){
  /* set up easy variable */
  this->M = M; this->K = K; this->N = N;
  this->um = m; this->uk = k; this->un = n;
  this->part = part;
  this->algorithm = MM_ALG_2D;
  this->other1 = other1; this->other2 = other2;
  this->M_stride = strideM;
  this->K_stride = strideK;
  this->N_stride = strideN;
  gemmSplitOrtho=_gemmSplitOrtho;
  M_chunks = (M + m - 1) / m;
  K_chunks = (K + k - 1) / k;
  N_chunks = (N + n - 1) / n;
  got_data = got_start = false;
  res_msg = NULL;
  algorithm = MM_ALG_3D;
  usesAtSync = CmiFalse;
  setMigratable(false);
  /* figure out size of our sections */
  if(part == MULTARG_A){
    if(thisIndex.x == (M_chunks - 1) * strideM){
      this->m = M % m;
      if(this->m == 0) this->m = m;
    }
    else this->m = m;
    if(thisIndex.y == (K_chunks - 1) * strideK){
      this->k = K % k;
      if(this->k == 0) this->k = k;
    }
    else this->k = k;
    this->n = n;
  } else if(part == MULTARG_B) {
    if(thisIndex.x == (K_chunks - 1) * strideK){
      this->k = K % k;
      if(this->k == 0) this->k = k;
    }
    else this->k = k;
    if(thisIndex.y == (N_chunks - 1) * strideN){
      this->n = N % n;
      if(this->n == 0) this->n = n;
    }
    else this->n = n;
    this->m = m;
  } else if(part == MULTARG_C) {
    if(thisIndex.x == (M_chunks - 1) * strideM){
      this->m = M % m;
      if(this->m == 0) this->m = m;
    }
    else this->m = m;
    if(thisIndex.y == (N_chunks - 1) * strideN){
      this->n = N % n;
      if(this->n == 0) this->n = n;
    }
    else this->n = n;
    this->k = k;
  }

  /* make communication groups, C also has to initialize reduction sections */
  if(part == MULTARG_A){
    int x = thisIndex.x / strideM;
    int y = thisIndex.y / strideK;
    commGroup3D = CProxySection_CLA_MM3D_multiplier::ckNew(p, x, x, 1, 0,
     N_chunks - 1, 1, y, y, 1);
    contribute(0, NULL, CkReduction::sum_int, cb);
  } else if(part == MULTARG_B) {
    int x = thisIndex.x / strideK;
    int y = thisIndex.y / strideN;
    commGroup3D = CProxySection_CLA_MM3D_multiplier::ckNew(p, 0,
     M_chunks - 1, 1, y, y, 1, x, x, 1);
    contribute(0, NULL, CkReduction::sum_int, cb);
  } else if(part == MULTARG_C) {
    init_cb = cb;
    int x = thisIndex.x / strideM;
    int y = thisIndex.y / strideN;
    commGroup3D = CProxySection_CLA_MM3D_multiplier::ckNew(p, x, x, 1, y, y, 1,
     0, K_chunks - 1, 1);
    commGroup3D.ckSectionDelegate(CProxy_CkMulticastMgr(gid).ckLocalBranch());
    CLA_MM3D_mult_init_msg *m = new CLA_MM3D_mult_init_msg(gid,
     CkCallback(CkIndex_CLA_Matrix::readyC(NULL),
     thisProxy(thisIndex.x, thisIndex.y)), CkCallback(
     CkIndex_CLA_Matrix::mult_done(NULL), thisProxy(thisIndex.x,
     thisIndex.y)));
    commGroup3D.initialize_reduction(m);
  }
}

CLA_Matrix::~CLA_Matrix(){
  if(algorithm == MM_ALG_2D){
    delete [] tmpA;
    delete [] tmpB;
  }
  else if(algorithm == MM_ALG_3D){
    if(res_msg != NULL)
      delete res_msg;
  }
}

void CLA_Matrix::pup(PUP::er &p){
  /* make sure we are not using the 3D algorithm */
  if(algorithm == MM_ALG_3D){
    CmiAbort("3D algorithm does not currently support migration.\n");
  }

  /* pup super class */
  CBase_CLA_Matrix::pup(p);

  /* pup shared vars */
  p | M; p | K; p | N; p | m; p | k; p | n;  p | um; p | uk; p | un;
  p | M_chunks; p | K_chunks; p | N_chunks;
  p | M_stride; p | K_stride; p | N_stride;
  p | part; p | algorithm;
  p | alpha; p | beta;
  p | gemmSplitOrtho;
  /* pup vars used by each algorithm */
  if(algorithm == MM_ALG_2D){
    p | row_count; p | col_count;
    p | other1; p | other2;
    if(part == MULTARG_C){
      if(p.isUnpacking()){
        tmpA = new double[m * K];
        tmpB = new double[K * n];
      }
      PUParray(p, tmpA, m * K);
      PUParray(p, tmpB, K * n);
    }
  }
  else if(algorithm == MM_ALG_3D){
    p | init_cb;
    p | got_start; p | got_data;
    p | commGroup3D;
  }
}

void CLA_Matrix::ResumeFromSync(void){
  /* recreate the section proxies */
  if(algorithm == MM_ALG_2D){
    if(part == MULTARG_A){
      commGroup2D = CProxySection_CLA_Matrix::ckNew(other2, thisIndex.x,
       thisIndex.x, 1, 0, (N_chunks - 1) * N_stride, N_stride);
      tmpA = tmpB = NULL;
    } else if(part == MULTARG_B) {
      commGroup2D = CProxySection_CLA_Matrix::ckNew(other2, 0,
       (M_chunks - 1) * M_stride, M_stride, thisIndex.y, thisIndex.y, 1);
      tmpA = tmpB = NULL;
    }
  } else if(algorithm == MM_ALG_3D){
#if 0
    if(part == MULTARG_A){
      int x = thisIndex.x / M_stride;
      int y = thisIndex.y / K_stride;
      commGroup3D = CProxySection_CLA_MM3D_multiplier::ckNew(p, x, x, 1, 0,
       N_chunks - 1, 1, y, y, 1);
      contribute(0, NULL, CkReduction::sum_int, cb);
    } else if(part == MULTARG_B) {
      int x = thisIndex.x / K_stride;
      int y = thisIndex.y / N_stride;
      commGroup3D = CProxySection_CLA_MM3D_multiplier::ckNew(p, 0,
       M_chunks - 1, 1, y, y, 1, x, x, 1);
      contribute(0, NULL, CkReduction::sum_int, cb);
    } else if(part == MULTARG_C) {
      init_cb = cb;
      int x = thisIndex.x / M_stride;
      int y = thisIndex.y / N_stride;
      commGroup3D = CProxySection_CLA_MM3D_multiplier::ckNew(p, x, x, 1, y, y,
       1, 0, K_chunks - 1, 1);
/*
      commGroup3D.ckSectionDelegate(CProxy_CkMulticastMgr(gid).ckLocalBranch());
      CLA_MM3D_mult_init_msg *m = new CLA_MM3D_mult_init_msg(gid,
       CkCallback(CkIndex_CLA_Matrix::readyC(NULL),
       thisProxy(thisIndex.x, thisIndex.y)), CkCallback(
       CkIndex_CLA_Matrix::mult_done(NULL), thisProxy(thisIndex.x,
       thisIndex.y)));
      commGroup3D.initialize_reduction(m);
*/
    }
#endif
  }
}

void CLA_Matrix::multiply(double alpha, double beta, double *data,
 void (*fptr) (void*), void *usr_data){
  if(algorithm == MM_ALG_2D){
    // A and B send out their chunks, ignoring alpha, beta, ftpr, and usr_data
    if(part == MULTARG_A){
      CLA_Matrix_msg *msg = new (m * k) CLA_Matrix_msg(data, m, k, thisIndex.x,
       thisIndex.y);
      commGroup2D.receiveA(msg);
    } else if(part == MULTARG_B){
      CLA_Matrix_msg *msg = new (k * n) CLA_Matrix_msg(data, k, n, thisIndex.x,
       thisIndex.y);
      commGroup2D.receiveB(msg);
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
  } else if(algorithm == MM_ALG_3D){
    if(part == MULTARG_A){
      CLA_Matrix_msg *msg = new (m * k) CLA_Matrix_msg(data, m, k, thisIndex.x,
       thisIndex.y);
      commGroup3D.receiveA(msg);
    } else if(part == MULTARG_B){
      CLA_Matrix_msg *msg = new (k * n) CLA_Matrix_msg(data, k, n, thisIndex.x,
       thisIndex.y);
      commGroup3D.receiveB(msg);
    } else if(part == MULTARG_C){
      fcb = fptr;
      user_data = usr_data;
      dest = data;
      this->alpha = alpha;
      this->beta = beta;
      got_start = true;
      if(got_data){
        got_start = got_data = false;
        /* transpose reduction msg and do the alpha and beta multiplications */
        double *data = (double*) res_msg->getData();
        transpose(data, n, m);
        for(int i = 0; i < m; i++)
          for(int j = 0; j < n; j++)
            dest[i * n + j] = beta * dest[i * n + j] + alpha * data[i * n + j];
        delete res_msg;
        res_msg = NULL;
        (*fcb)(user_data);
      }
    }
  }
}

void CLA_Matrix::receiveA(CLA_Matrix_msg *msg){
  /* store current part */
  row_count++;
  for(int i = 0; i < m; i++)
    memcpy(&tmpA[K * i + uk * (msg->fromY / K_stride)], &msg->data[i * msg->d2],
     msg->d2 * sizeof(double));
  delete msg;

  /* If we have all the parts, multiply */
  if(row_count == K_chunks && col_count == K_chunks && got_start)
    multiply();
}

void CLA_Matrix::receiveB(CLA_Matrix_msg *msg){
  /* store current part */
  col_count++;
  memcpy(&tmpB[n * uk * (msg->fromX / K_stride)], msg->data,
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
#define ORTHO_DGEMM_SPLIT 

#define BUNDLE_USER_EVENTS
#ifdef ORTHO_DGEMM_SPLIT 
  double betap = 1.0;
  int Ksplit_m=gemmSplitOrtho;  
  int Ksplit   = (K > Ksplit_m) ? Ksplit_m : K;
  int Krem     = (K % Ksplit);
  int Kloop    = K/Ksplit-1;

#ifndef CMK_OPTIMIZE
  double StartTime=CmiWallTimer();
#endif

#ifdef TEST_ALIGN
  CkAssert((unsigned int) tmpA %16==0);
  CkAssert((unsigned int) tmpB %16==0);
  CkAssert((unsigned int) dest %16==0);
#endif

#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", trans, trans, m, n, Ksplit, alpha, beta, K, n, m);
#endif
  DGEMM(&trans, &trans, &m, &n, &Ksplit, &alpha, tmpA, &K, tmpB, &n, &beta,
        dest,&m);

#ifndef BUNDLE_USER_EVENTS
#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(401, StartTime, CmiWallTimer());
#endif
#endif

  CmiNetworkProgress();
  
  for(int i=1;i<=Kloop;i++){
    int aoff = Ksplit*i;
    int boff = n*i*Ksplit;
    if(i==Kloop){Ksplit+=Krem;}
#ifndef BUNDLE_USER_EVENTS
#ifndef CMK_OPTIMIZE
    StartTime=CmiWallTimer();
#endif
#endif

#ifdef TEST_ALIGN
    CkAssert((unsigned int) &(tmpA[aoff]) %16==0);
    CkAssert((unsigned int) &(tmpB[boff]) %16==0);
    CkAssert((unsigned int) dest %16==0);
#endif

#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", trans, trans, m, n, Ksplit, alpha, betap, K, n, m);
#endif
    DGEMM(&trans, &trans, &m, &n, &Ksplit, &alpha, &tmpA[aoff], &K, 
          &tmpB[boff], &n, &betap,dest, &m);
#ifndef BUNDLE_USER_EVENTS
#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(401, StartTime, CmiWallTimer());
#endif
#endif
    CmiNetworkProgress();
  }//endfor

#ifdef BUNDLE_USER_EVENTS
#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(401, StartTime, CmiWallTimer());
#endif
#endif

#else
  /* old unsplit version */
#ifndef CMK_OPTIMIZE
    double StartTime=CmiWallTimer();
#endif

#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", trans, trans, m, n, K, alpha, beta, K, n, m);
#endif
     DGEMM(&trans, &trans, &m, &n, &K, &alpha, tmpA, &K, tmpB, &n, &beta,
   dest, &m);
#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(401, StartTime, CmiWallTimer());
#endif

#endif
  /* transpose result */
  transpose(dest, n, m);

  /* tell caller we are done */
  fcb(user_data);
}

void CLA_Matrix::readyC(CkReductionMsg *msg){
  CkCallback cb(CkIndex_CLA_Matrix::ready(NULL), thisProxy(0, 0));
  contribute(0, NULL, CkReduction::sum_int, cb);
  delete msg;
}

void CLA_Matrix::ready(CkReductionMsg *msg){
  init_cb.send();
  delete msg;
}

void CLA_Matrix::mult_done(CkReductionMsg *msg){
  if(got_start){
    got_start = got_data = false;
    /* transpose reduction msg and do the alpha and beta multiplications */
    double *data = (double*) msg->getData();
    transpose(data, n, m);
    for(int i = 0; i < m; i++)
      for(int j = 0; j < n; j++)
        dest[i * n + j] = beta * dest[i * n + j] + alpha * data[i * n + j];
    delete msg;
    msg = NULL;
    (*fcb)(user_data);
  }
  else{
    got_data = true;
    res_msg = msg;
  }
}

/******************************************************************************/
/* CLA_Matrix_msg */
CLA_Matrix_msg::CLA_Matrix_msg(double *data, int d1, int d2, int fromX,
 int fromY){
  memcpy(this->data, data, d1 * d2 * sizeof(double));
  this->d1 = d1; this->d2 = d2;
  this->fromX = fromX; this->fromY = fromY;
}

/******************************************************************************/
/* CLA_MM3D_multiplier */
CLA_MM3D_multiplier::CLA_MM3D_multiplier(int m, int k, int n){
  this->m = m; this->k = k; this->n = n;
  data_msg = NULL;
  gotA = gotB = false;
}

void CLA_MM3D_multiplier::initialize_reduction(CLA_MM3D_mult_init_msg *m){
  reduce_CB = m->reduce;
  CkGetSectionInfo(sectionCookie, m);
  redGrp = CProxy_CkMulticastMgr(m->gid).ckLocalBranch();
  redGrp->contribute(0, NULL, CkReduction::sum_int, sectionCookie, m->ready);
  delete m;
}

void CLA_MM3D_multiplier::receiveA(CLA_Matrix_msg *msg){
  gotA = true;
  if(gotB){
    multiply(msg->data, data_msg->data);
    delete msg;
    delete data_msg;
  }
  else
    data_msg = msg;
}

void CLA_MM3D_multiplier::receiveB(CLA_Matrix_msg *msg){
  gotB = true;
  if(gotA){
    multiply(data_msg->data, msg->data);
    delete msg;
    delete data_msg;
  }
  else
    data_msg = msg;
}

void CLA_MM3D_multiplier::multiply(double *A, double *B){
  double alpha = 1, beta = 0;
  gotA = gotB = false;
  char trans = 't';
  double *C = new double[m * n];
#ifndef CMK_OPTIMIZE
  double  StartTime=CmiWallTimer();
#endif
#ifdef TEST_ALIGN
  CkAssert((unsigned int) A %16==0);
  CkAssert((unsigned int) B %16==0);
  CkAssert((unsigned int) C %16==0);
#endif

#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", trans, trans, m, n, k, alpha, beta, k, n, m);
#endif
  DGEMM(&trans, &trans, &m, &n, &k, &alpha, A, &k, B, &n, &beta, C, &m);
#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(402, StartTime, CmiWallTimer());
#endif
  CmiNetworkProgress();
  //    redGrp->contribute(m * n * sizeof(double), C, CkReduction::sum_double,
  redGrp->contribute(m * n * sizeof(double), C, sumFastDoubleType,
   sectionCookie, reduce_CB);
  delete [] C;
}

#include "CLA_Matrix.def.h"
