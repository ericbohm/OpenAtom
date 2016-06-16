#include "standard_include.h"
#include "allclass_gwbse.h"
#include "eps_matrix.h"
#include "messages.h"
#include "pmatrix.h"
#include "controller.h"
#include "states.h"
#include "fft_routines.h"
#include "fft_controller.h"

#define eps_rows 20
#define eps_cols 20
#define IDX_eps(r,c) ((r)*eps_cols + (c))

EpsMatrix2D::EpsMatrix2D(){
//CkPrintf("\nCalling default construction for EpsMatrix2D(%d, %d)\n", thisIndex.x, thisIndex.y);
  GWBSE* gwbse = GWBSE::get();

  // Set some constants
  L = gwbse->gw_parallel.L;
  nfft = gwbse->gw_parallel.fft_nelems;
  qindex = Q_IDX; // Eventually the controller will set this

  // Grab a local pointer to the fft controller for fft-ing our rows
  // TODO: Is this guaranteed to be safe (is the local branch created for sure)?

  // Figure out what part of the matrix we have and allocate space
  num_rows = eps_rows;//120;//gwbse->gw_parallel.rows_per_chare;
  num_cols = eps_cols;//120;//gwbse->gw_parallel.cols_per_chare;
  num_chares = (matrix_dimension / num_rows) * (matrix_dimension / num_cols);
  start_row = thisIndex.x * num_rows;
  start_col = thisIndex.y * num_cols;
  local_mtx_size_1d_y = 1;//1728; // TODO: where from?
  receive_counter = 0;

  data = new complex[num_rows * num_cols];

  total_time = 0.0;
}

void EpsMatrix2D::setSize(int matrix_dimension_in){
  matrix_dimension = matrix_dimension_in;
}

EpsMatrix2D::EpsMatrix2D(CLA_Matrix_interface mat){
//  CkPrintf("\nCalling CLA_Matrix_interface construction for EpsMatrix2D(%d, %d)\n", thisIndex.x, thisIndex.y);   
  GWBSE* gwbse = GWBSE::get();

  // Set some constants
  L = gwbse->gw_parallel.L;
  nfft = gwbse->gw_parallel.fft_nelems;
  qindex = Q_IDX; // Eventually the controller will set this

  // Grab a local pointer to the fft controller for fft-ing our rows
  // TODO: Is this guaranteed to be safe (is the local branch created for sure)?

  // Figure out what part of the matrix we have and allocate space
  num_rows = eps_rows;//gwbse->gw_parallel.rows_per_chare;
  num_cols = eps_cols;//gwbse->gw_parallel.cols_per_chare;
  num_chares = (matrix_dimension / num_rows) * (matrix_dimension / num_cols);
  start_row = thisIndex.x * num_rows;
  start_col = thisIndex.y * num_cols;
  local_mtx_size_1d_y = 1;//1728; // TODO: where from?
  receive_counter = 0;

   
  data = new complex[num_rows * num_cols];

  for(int i=0;i<num_rows * num_cols;i++){
    data[i].re = 0;
    data[i].im = 0;
  }
  total_time = 0.0;
  
  this->matrix = mat;
}

void EpsMatrix2D::setI(CLA_Matrix_interface mat, bool clean){
  this->matrix = mat;
  if(clean)
    data = new complex[num_rows * num_cols];
}



void EpsMatrix2D::receiveFs(complex eps_data[140*140]){

  int start_i = (matrix_dimension/num_cols)*thisIndex.x;
  int start_j = (matrix_dimension/num_rows)*thisIndex.y;

  for(int i=0;i<num_rows;i++)
    for(int j=0;j<num_cols;j++)
        data[IDX_eps(i,j)] = eps_data[start_row*num_cols+ start_col];
  //print_res();

    int i = 0;//count;
    CkCallback *cb = new CkCallback(CkReductionTarget(Controller, epsilon_created), controller_proxy);
    contribute(sizeof(int), &i, CkReduction::sum_int, *cb);
  
}


 
void EpsMatrix2D::multiply(double alpha, double beta){
//   CkPrintf("\n accessing %d,%d\n",thisIndex.x, thisIndex.y);
   matrix.multiply(alpha, beta, data, EpsMatrix2D::done_cb, (void*) this,
       thisIndex.x, thisIndex.y);

 }


void EpsMatrix2D::complement_two(){
//     for(int i=0; i<N; i++)
    //      M1[i*N+i] += compl_two;

    int i = 0;
    complex compl_two(2.0, 0.0);
    for(int i=0;i<17;i++)//num_rows;i++)
      data[IDX_eps(i,i)] += compl_two;

  CkCallback *cb = new CkCallback(CkReductionTarget(Controller, complement_multiplied), controller_proxy);
      contribute(sizeof(int), &i, CkReduction::sum_int,
        *cb
      );
}
 void inline EpsMatrix2D::round_done(void){
    int i=0;
    //print_res();
    CmiMemoryCheck();
    CkCallback *cb = new CkCallback(CkReductionTarget(Controller, m_multiplied), controller_proxy);
      contribute(sizeof(int), &i, CkReduction::sum_int,
        *cb
      );
    }


long double EpsMatrix2D::max_fn(int size){

    long double max = -1;

    for(int i =0;i<size;i++){
      double temp = abs(total[i]);
      if(temp < 0)
        temp *= -1;
      if(temp>max)
        max = temp;
    }

    return max;
}

void EpsMatrix2D::scalar_multiply(double alpha){
    alpha = 1/alpha;
#if 0
     std::ostringstream sout;
     sout << alpha << "\n";
     CkPrintf("\nSending alpha = %s",sout.str().c_str());
#endif
    for(int i=0;i<num_rows;i++)
        for(int j=0;j<num_cols;j++)
            data[IDX_eps(i,j)] = alpha*data[IDX_eps(i,j)]; 

 //   print_res();
    CkCallback *cb = new CkCallback(CkReductionTarget(Controller, scalar_multiplied), controller_proxy);
    contribute(*cb);
}

void EpsMatrix2D::findAlpha(){


    //print_res();
    for(int i=0;i<num_rows;i++)
        for(int j=0;j<num_cols;j++)
        total[i] += data[IDX_eps(i,j)];

    long double max = 0;
    int chunk_count = matrix_dimension/num_rows;
    if(thisIndex.y==chunk_count-1){
        max = max_fn(num_rows);   //return the result;

#if 0
        std::ostringstream sout;
        sout << max << "\n";
        CkPrintf("\nSending max = %s",sout.str().c_str());
#endif 

    }else
        thisProxy(thisIndex.x, thisIndex.y+1).findAlpha();
 
    max = 5e+105; 
    contribute(sizeof(long double), &max, CkReduction::max_double,
        CkCallback(CkReductionTarget(Controller, found_alpha), controller_proxy)
      );
}


void EpsMatrix2D::convergence_check(CProxy_EpsMatrix2D cproxy){//complex *M1, complex *M0, int N

    
    std::vector<complex> data_out(num_rows*num_cols);
    for(int i=0;i<num_rows*num_cols;i++)
      data_out[i] = data[i];

    cproxy(thisIndex.x, thisIndex.y).receiveConvCheck(data_out);
}

void EpsMatrix2D::receiveConvCheck(std::vector<complex> data_in){

   complex Rmax=0;  // the largest element
   complex tmp;
   for(int i=0; i<num_rows*num_cols; i++){
        tmp = data[i] - data_in[i] ;
        if( tmp.re > Rmax.re ){ Rmax = tmp; }
    }

#if 0
   std::ostringstream sout;
   sout << Rmax << "\n";
   CkPrintf("\nSending max = %s",sout.str().c_str());
#endif

   contribute(sizeof(complex), &Rmax, CkReduction::max_double,
        CkCallback(CkReductionTarget(Controller, converge_results), controller_proxy)
      );
}


void EpsMatrix2D::print_res(){
      std::ostringstream sout;
      complex total = 0;
      for(int i = 0; i < num_rows; i++){
        for(int j = 0; j < num_cols; j++)
          total+= data[i*num_rows+j];
      }
      sout << total << " \n";
      if(thisIndex.x==0&&thisIndex.y==4)
      CkPrintf("(%d, %d) result:\n%s", thisIndex.x, thisIndex.y,
       sout.str().c_str());
    }

void EpsMatrix2D::reportPTime() {
  CkReduction::statisticsElement stats(total_time);
  int tuple_size = 2;
  CkReduction::tupleElement tuple_reduction[] = {
    CkReduction::tupleElement(sizeof(double), &total_time, CkReduction::sum_double),
    CkReduction::tupleElement(sizeof(CkReduction::statisticsElement), &stats, CkReduction::statistics) };

  CkReductionMsg* msg = CkReductionMsg::buildFromTuple(tuple_reduction, tuple_size);
  msg->setCallback(CkCallback(CkIndex_Controller::reportPTime(NULL), controller_proxy));
  contribute(msg);
}

void EpsMatrix2D::applyFs() {
  double start = CmiWallTimer();
//  print_res();
  PsiCache* psi_cache = psi_cache_proxy.ckLocalBranch();
//CkPrintf("\napplyFs entered\n");
#ifdef USE_LAPACK
  // Common variables for both ZGERC and ZGEMM
  int M = num_rows, N = num_cols;
  complex alpha = -1.0;
#ifdef USE_ZGEMM
  int K = L; // If using ZGEMM, we compute all outer products with one call
  int LDF = matrix_dimension; // Leading dimension of fs
  complex beta = 1.0;
  char opA = 'N', opB = 'C';
  complex* fs = psi_cache->getF(0);
  ZGEMM(&opA, &opB, &N, &M, &K,
    &alpha, &(fs[start_col]), &LDF,
    &(fs[start_row]), &LDF,
    &beta, data, &N);
#else
  int K = 1; // If using ZGERC, we compute each outer product one at a time
  for (int l = 0; l < L; l++) {
    complex* f = psi_cache->getF(l);
    ZGERC(&N, &M, &alpha, &(f[start_col]), &K, &(f[start_row]), &K, data, &N);
  }
#endif
#else
  for (int l = 0; l < L; l++) {
    complex* f = psi_cache->getF(l);
    for (int r = 0; r < num_rows; r++) {
      for (int c = 0; c < num_cols; c++) {
        data[IDX_eps(r,c)] += f[r+start_row]*f[c+start_col].conj() * -1.0;
//        CkPrintf("\nMultiplying %lf\n", f[r+start_row].re);        
      }
    }
  }
//  CkPrintf("\nApply Fs got fs that were %d", psi_cache->getWrote());
#endif // endif for ifdef USE_LAPACK
//CkPrintf("\napplyFs exited\n");
  //print_res();
  contribute(CkCallback(CkReductionTarget(Controller, psiComplete), controller_proxy));
  total_time += CmiWallTimer() - start;
}

void EpsMatrix2D::checkReady(){
//    contribute(CkCallback(CkReductionTarget(Controller, allReady), controller_proxy));
}

void EpsMatrix2D::receiveChunk(Phase2Message* msg) {
  int local_mtx_size_x = num_cols;
  int local_mtx_size_y = num_rows;
  int local_y = msg->global_y - thisIndex.y * local_mtx_size_y;
  int local_x = msg->global_x - thisIndex.x * local_mtx_size_x;
  for(unsigned i=0; i<msg->size; ++i){
    data[IDX_eps(local_y, local_x) + i] = msg->data[i];
  }

  if(++receive_counter == num_rows){
    contribute(CkCallback(CkReductionTarget(Controller, phase2_complete), controller_proxy));
  }
  delete msg;
}

void EpsMatrix2D::createTranspose(){
    std::vector<complex> incoming;
    unsigned n = 0;
    complex* new_data;
    new_data = new complex[num_rows*num_cols];
    for(int i=0;i<num_rows;i++){
        for(int j=0;j<num_cols;j++){
            new_data[n] =  data[IDX_eps(j,i)];//msg->data[n++] = data[IDX(j,i)];//.conj();    //transposed, need to conjugate as well
            incoming.push_back(new_data[n++]);
        }
    }
#if 0
    std::ostringstream sout;
    sout << msg->data[4] << "\n";
    CkPrintf("\nSending msg->data[4] = %s",sout.str().c_str());
#endif
    pmatrix2D_bproxy(thisIndex.y, thisIndex.x).receiveTranspose(incoming);//new_data);
}

void EpsMatrix2D::receiveTranspose(std::vector<complex> new_data){

#if 0 
   std::ostringstream sout;
    sout << msg->data[4] << "\n";
    CkPrintf("\nReceiving msg->data[4] = %s",sout.str().c_str());
#endif
  unsigned n = 0;
  for(int i=0;i<num_rows;i++)
        for(int j=0;j<num_cols;j++)
            data[IDX_eps(i,j)] = new_data[n++];//msg->data[n++];
    
  contribute(CkCallback(CkReductionTarget(Controller, transpose_complete), controller_proxy));  
}


void EpsMatrix2D::sendTo(CProxy_EpsMatrix2D receiver_proxy){
    std::vector<complex> incoming;
    for(int i=0;i<num_rows;i++){
        for(int j=0;j<num_cols;j++){
            incoming.push_back(data[IDX_eps(i,j)]);
        }
    }
    receiver_proxy(thisIndex.y, thisIndex.x).receiveData(incoming);
}

void EpsMatrix2D::receiveData(std::vector<complex> new_data){

  unsigned n = 0;
  for(int i=0;i<num_rows;i++)
        for(int j=0;j<num_cols;j++)
            data[IDX_eps(i,j)] = new_data[n++];//msg->data[n++];
   
  contribute(CkCallback(CkReductionTarget(Controller, finished_copy), controller_proxy));
}



// TODO: These methods shouldn't be part of PMatrix, and should also just be
// computed once at startup.
void EpsMatrix2D::kqIndex(unsigned ikpt, unsigned& ikq, int* uklapp){
  GWBSE* gwbse = GWBSE::get();

  // temporary space to save k/q/k+q vectors
  double *this_k, *this_q;
  double k_plus_q[3], k_plus_q_orig[3];

  this_k = gwbse->gwbseopts.kvec[ikpt];
  this_q = gwbse->gwbseopts.qvec[qindex];

  for (int i=0; i<3; i++) {
    // calculate k+q vector
    k_plus_q[i] = this_k[i] + this_q[i]; // k+q vector
    k_plus_q_orig[i] = k_plus_q[i]; // save it for Umklapp process
    // if not 0 =< k+q [i] <1, adjust k+q so that k+q[i] is in the Brillouine zone
    if ( k_plus_q[i] >= 1 ) {
      k_plus_q[i] -= 1;
    }
    else if( k_plus_q[i] < 0 ){
      k_plus_q[i] += 1;
    }
  }

  // find k+q vector index
  for (int kk=0; kk < gwbse->gwbseopts.nkpt; kk++) {
    bool match = true;
    this_k = gwbse->gwbseopts.kvec[kk];
    //this_k is now a difference between k and k+q
    for (int i=0; i<3; i++) {
      if (this_k[i] != k_plus_q[i]) {
        match = false;
        break;
      }
    }
    if (match) {
      ikq = kk;
      break;
    }
  }
  // save umklapp scattering information
  for (int i=0; i<3; i++) {
    uklapp[i] = int( k_plus_q_orig[i] - k_plus_q[i] );
  }

}

void EpsMatrix2D::getUmklappFactor(complex* umklapp_factor, int uklpp[3]){

  if (uklpp[0]==0 && uklpp[1]==0 && uklpp[2]==0){
    // do nothing
  }
  else{
    GWBSE *gwbse = GWBSE::get();
    int* nfft;
    nfft = gwbse->gw_parallel.fft_nelems;
    double *a1, *a2, *a3, *b1, *b2, *b3;
    a1 = gwbse->gwbseopts.a1;
    a2 = gwbse->gwbseopts.a2;
    a3 = gwbse->gwbseopts.a3;
    b1 = gwbse->gwbseopts.b1;
    b2 = gwbse->gwbseopts.b2;
    b3 = gwbse->gwbseopts.b3;
    double lattconst = gwbse->gwbseopts.latt;

    double rijk, G0, phase;
    unsigned counter = 0;
    for(int i=0; i<nfft[0]; i++){
      for(int j=0; j<nfft[1]; j++){
        for(int k=0; k<nfft[2]; k++){
          phase = 0;
          for (int l=0; l<3; l++){
            rijk = a1[l]*i/nfft[0] + a2[l]*j/nfft[1] + a3[l]*k/nfft[2];
            G0 = b1[l]*uklpp[0] + b2[l]*uklpp[1] + b3[l]*uklpp[2];
            G0 *= -2*M_PI/lattconst;
            phase += rijk*G0;
          }
          umklapp_factor[counter].re = cos(phase);
          umklapp_factor[counter].im = sin(phase);
          counter += 1;
        }// end k loop
      }// end j loop
    }// end i loop
  }//end if-else statement

}//end function

#include "eps_matrix.def.h"
