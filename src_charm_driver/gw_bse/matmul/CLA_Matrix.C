#include "CLA_Matrix.h"
#include "ckcomplex.h"

#include "mylapack.h"
#if 0
#include <sstream>
using std::ostringstream;
using std::endl;
#endif

#ifdef FORTRANUNDERSCORE
#define DGEMM dgemm_
#else
#define DGEMM dgemm_
#endif


extern "C" {void  DGEMM(char *, char *, int *, int *, int *, double *,
               double *, int *, double *, int *, double *,
               double *, int *);}

#define MULTARG_A   0
#define MULTARG_B   1
#define MULTARG_C   2

/* readonly */ CkGroupID CLA_GID;
/* readonly */ CProxy_CLA_MM3D_multiplier CLA_3D_multiplier_proxy;

/******************************************************************************/
/* helper functions */

/* Should be called by user to create matrices. Documented in header file. */
int make_multiplier(CLA_Matrix_interface *A, CLA_Matrix_interface *B,
 CLA_Matrix_interface *C, CProxy_ArrayElement bindA,
 CProxy_ArrayElement bindB, CProxy_ArrayElement bindC,
 int M, int K, int N, int m, int k, int n, int strideM, int strideK,
 int strideN, CkCallback cbA, CkCallback cbB,
 CkCallback cbC, CkGroupID gid, int algorithm){

//  CkPrintf("\ncalling multiplier\n");

  /* validate arguments */
  if(algorithm < MM_ALG_MIN || MM_ALG_MAX < algorithm)
    return ERR_INVALID_ALG;

  if(m > M || k > K || n > N)
    return ERR_INVALID_DIM;

  /* initialize readonly vars */
  CLA_GID = gid;

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
  A->init(pa);
  B->init(pb);
  C->init(pc);

  /* populate arrays */
  int M_chunks = (M + m - 1) / m; // same as ceil(1.0 * M / m)
  int K_chunks = (K + k - 1) / k; // same as ceil(1.0 * K / k)
  int N_chunks = (N + n - 1) / n; // same as ceil(1.0 * N / n)

  if(algorithm == MM_ALG_2D){
//    CkPrintf("\nPerforming 2D MM\n");

    for(int i = 0; i < M_chunks; i++)
      for(int j = 0; j < K_chunks; j++)
        (A->p(i * strideM, j * strideK)).insert(M, K, N, m, k, n, strideM,
         strideK, strideN, MULTARG_A, B->p, C->p, cbA);
    A->p.doneInserting();

    for(int i = 0; i < K_chunks; i++)
      for(int j = 0; j < N_chunks; j++)
        (B->p(i * strideK, j * strideN)).insert(M, K, N, m, k, n, strideM,
         strideK, strideN, MULTARG_B, A->p, C->p, cbB);
    B->p.doneInserting();

    for(int i = 0; i < M_chunks; i++)
      for(int j = 0; j < N_chunks; j++)
        (C->p(i * strideM, j * strideN)).insert(M, K, N, m, k, n, strideM,
         strideK, strideN, MULTARG_C, A->p, B->p, cbC);
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
         strideM, strideK, strideN, MULTARG_A, cbA);//, gid);
    A->p.doneInserting();

    for(int i = 0; i < K_chunks; i++)
      for(int j = 0; j < N_chunks; j++)
        (B->p(i * strideK, j * strideN)).insert(mult, M, K, N, m, k, n,
         strideM, strideK, strideN, MULTARG_B, cbB);//, gid);
    B->p.doneInserting();

    for(int i = 0; i < M_chunks; i++)
      for(int j = 0; j < N_chunks; j++)
        (C->p(i * strideM, j * strideN)).insert(mult, M, K, N, m, k, n,
         strideM, strideK, strideN, MULTARG_C, cbC);//, gid);
    C->p.doneInserting();
  }

  return SUCCESS;
}

/* Transpose data, which has dimension m x n */
void transpose(complex *data, int m, int n){
  if(m == n){
    /* transpose square matrix in place */
    for(int i = 0; i < m; i++)
      for(int j = i + 1; j < n; j++){
        complex tmp = data[i * n + j];
        data[i * n + j] = data[j * m + i];
        data[j * m + i] = tmp;
      }
  }
  else {
    complex *tmp = new complex[m * n];
    memcpy(tmp, data, m * n * sizeof(complex));
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
 CProxy_CLA_Matrix other1, CProxy_CLA_Matrix other2, CkCallback ready){
  /* initialize simple members */
  this->M = M; this->K = K; this->N = N;
  this->um = m; this->uk = k; this->un = n;
  this->part = part;
  this->algorithm = MM_ALG_2D;
  this->other1 = other1; this->other2 = other2;
  this->M_stride = strideM;
  this->K_stride = strideK;
  this->N_stride = strideN;
  M_chunks = (M + m - 1) / m;
  K_chunks = (K + k - 1) / k;
  N_chunks = (N + n - 1) / n;
  algorithm = MM_ALG_2D;
//  usesAtSync = CmiTrue;

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
    tmpA = new complex[this->m * K];
    tmpB = new complex[K * this->n];
  }

  /* let the caller know we are ready */
  contribute(0, NULL, CkReduction::sum_int, ready);
}

/* constructor for 3D algorithm */
CLA_Matrix::CLA_Matrix(CProxy_CLA_MM3D_multiplier p, int M, int K, int N,
 int m, int k, int n, int strideM, int strideK, int strideN, int part,
 CkCallback cb){//, CkGroupID gid){
  /* set up easy variable */
  this->M = M; this->K = K; this->N = N;
  this->um = m; this->uk = k; this->un = n;
  this->part = part;
  this->algorithm = MM_ALG_2D;
  this->other1 = other1; this->other2 = other2;
  this->M_stride = strideM;
  this->K_stride = strideK;
  this->N_stride = strideN;
  M_chunks = (M + m - 1) / m;
  K_chunks = (K + k - 1) / k;
  N_chunks = (N + n - 1) / n;
  got_data = got_start = false;
  res_msg = NULL;
  algorithm = MM_ALG_3D;
//  usesAtSync = CmiTrue;

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
    commGroup3D.ckSectionDelegate(
     CProxy_CkMulticastMgr(CLA_GID).ckLocalBranch());
//    CLA_MM3D_mult_init_msg *m = new CLA_MM3D_mult_init_msg(gid,
    CLA_MM3D_mult_init_msg *m = new CLA_MM3D_mult_init_msg(
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
//  if(thisIndex.x == 0 && thisIndex.y == 0){
//    CkPrintf("[0, 0] got pup command %d %d %d\n", p.isUnpacking()?1:0,
//             p.isPacking()?1:0, p.isSizing()?1:0);
//  }
  /* make sure we are not using the 3D algorithm */
//  if(algorithm == MM_ALG_3D){
//    CmiAbort("3D algorithm does not currently support migration.\n");
//  }

  /* pup super class */
  CBase_CLA_Matrix::pup(p);

  /* pup shared vars */
  p | M; p | K; p | N; p | m; p | k; p | n;  p | um; p | uk; p | un;
  p | M_chunks; p | K_chunks; p | N_chunks;
  p | M_stride; p | K_stride; p | N_stride;
  p | part; p | algorithm;
  p | alpha; p | beta;

  /* pup vars used by each algorithm */
  if(algorithm == MM_ALG_2D){
    p | row_count; p | col_count;
    p | other1; p | other2;
    if(part == MULTARG_C){
      if(p.isUnpacking()){
        tmpA = new complex[m * K];
        tmpB = new complex[K * n];
      }
      PUParray(p, tmpA, m * K);
      PUParray(p, tmpB, K * n);
    }
    p | commGroup2D;
  }
  else if(algorithm == MM_ALG_3D){
    p | init_cb;
    p | got_start; p | got_data;
    p | commGroup3D;
    if(p.isUnpacking()){
      res_msg = NULL;
      if(part == MULTARG_C){
        CkMulticastMgr *mgr = CProxy_CkMulticastMgr(CLA_GID).ckLocalBranch();
        mgr->resetSection(commGroup3D);
      }
    }
  }
}

void CLA_Matrix::ResumeFromSync(void){
  /* recreate the section proxies */
//  CkPrintf("multiplier [%d %d] of %d resumed\n", thisIndex.x, thisIndex.y, part);
  if(algorithm == MM_ALG_2D){
    if(part == MULTARG_A){
      tmpA = tmpB = NULL;
    } else if(part == MULTARG_B) {
      tmpA = tmpB = NULL;
    }
  } else if(algorithm == MM_ALG_3D){
    if(part == MULTARG_C){
      CLA_MM3D_mult_reinit_msg *msg = new CLA_MM3D_mult_reinit_msg(
       CkCallback(CkCallback::ignore));
      commGroup3D.reinitialize_reduction(msg);
    }
  }
}

void CLA_Matrix::multiply(double alpha, double beta, complex *data,
 void (*fptr) (void*), void *usr_data){
  if(algorithm == MM_ALG_2D){
//      CkPrintf("\nCalling multiply\n");
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
        complex *data = (complex*) res_msg->getData();
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
     msg->d2 * sizeof(complex));
  delete msg;

  /* If we have all the parts, multiply */
  if(row_count == K_chunks && col_count == K_chunks && got_start)
    multiply();
}

void CLA_Matrix::receiveB(CLA_Matrix_msg *msg){
  /* store current part */
  col_count++;
  memcpy(&tmpB[n * uk * (msg->fromX / K_stride)], msg->data,
   msg->d1 * msg->d2 * sizeof(complex));
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

  myGEMM(&trans, &trans, &m, &n, &K, &alpha, tmpA, &K, tmpB, &n, &beta, dest, &m);
 
//Assuming m,n,k are same
#if 0
  for (int c = 0; c < m; c++) {
      for (int d = 0; d < n; d++) {
        for (int k = 0; k < K; k++) {
          dest[c*n+d] += tmpA[c*K+k]* tmpB[k*n+d];
        }
      }
    }
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
    complex *data = (complex*) msg->getData();
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

void CLA_Matrix::synchronize(){
//  CkPrintf("[%d %d] of %d got sync command\n", thisIndex.x, thisIndex.y, part);
  if(algorithm == MM_ALG_2D)
;//    AtSync();
  else if(algorithm == MM_ALG_3D){
    if(part == MULTARG_C)
      commGroup3D.synchronize(new CLA_MM3D_migrate_msg());
//    AtSync();
  }
  else
    CkAbort("This matrix multiplication algorithm does not support migration");
}

/******************************************************************************/
/* CLA_Matrix_msg */
CLA_Matrix_msg::CLA_Matrix_msg(complex *data, int d1, int d2, int fromX,
 int fromY){
  memcpy(this->data, data, d1 * d2 * sizeof(complex));
  this->d1 = d1; this->d2 = d2;
  this->fromX = fromX; this->fromY = fromY;
}

/******************************************************************************/
/* CLA_MM3D_multiplier */
CLA_MM3D_multiplier::CLA_MM3D_multiplier(int m, int k, int n){
  this->m = m; this->k = k; this->n = n;
  data_msg = NULL;
  gotA = gotB = false;
  //usesAtSync = CmiTrue;
}

void CLA_MM3D_multiplier::initialize_reduction(CLA_MM3D_mult_init_msg *m){
  reduce_CB = m->reduce;
  CkGetSectionInfo(sectionCookie, m);
  redGrp = CProxy_CkMulticastMgr(CLA_GID).ckLocalBranch();
//  redGrp = CProxy_CkMulticastMgr(m->gid).ckLocalBranch();
  redGrp->contribute(0, NULL, CkReduction::sum_int, sectionCookie, m->ready);
  delete m;
}

void CLA_MM3D_multiplier::reinitialize_reduction(CLA_MM3D_mult_reinit_msg *m){
  CkGetSectionInfo(sectionCookie, m);
  redGrp = CProxy_CkMulticastMgr(CLA_GID).ckLocalBranch();
//  redGrp = CProxy_CkMulticastMgr(m->gid).ckLocalBranch();
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

void CLA_MM3D_multiplier::multiply(complex *A, complex *B){
//  CkPrintf("(%d %d %d) lives on %d\n", thisIndex.x, thisIndex.y, thisIndex.z, CkMyPe());
  double alpha = 1, beta = 0;
  gotA = gotB = false;
  char trans = 't';
  complex *C = new complex[m * n];
  myGEMM(&trans, &trans, &m, &n, &k, &alpha, A, &k, B, &n, &beta, C, &m);
  redGrp->contribute(m * n * sizeof(double), C, CkReduction::sum_double,
   sectionCookie, reduce_CB);
  delete [] C;
}

void CLA_MM3D_multiplier::synchronize(CLA_MM3D_migrate_msg *msg){
  delete msg;
//  CkPrintf("(%d %d %d) got sync\n", thisIndex.x, thisIndex.y, thisIndex.z);
  AtSync();
}

void CLA_MM3D_multiplier::pup(PUP::er &p){
  CBase_CLA_MM3D_multiplier::pup(p);
  p | m; p | k; p | n;
  p | gotA; p | gotB;
  p | sectionCookie;
  p | reduce_CB;
  if(p.isUnpacking())
    redGrp = CProxy_CkMulticastMgr(CLA_GID).ckLocalBranch();
}

void CLA_MM3D_multiplier::ResumeFromSync(){
//  CkPrintf("(%d %d %d) has reduction group %p\n", thisIndex.x, thisIndex.y,
//    thisIndex.z, redGrp);
}

#include "CLA_Matrix.def.h"
