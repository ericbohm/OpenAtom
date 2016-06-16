#include "standard_include.h"
#include "allclass_gwbse.h"
#include "pmatrix.h"
#include "messages.h"
#include "eps_matrix.h"
#include "controller.h"
#include "states.h"
#include "fft_routines.h"
#include "fft_controller.h"

#define CHARE_NUM 10
#define IDX(r,c) ((r)*num_cols + (c))
#define IDX_eps(r,c) ((r)*eps_cols + (c))
#define eps_chares_x 10
#define eps_chares_y 10


PMatrix2D::PMatrix2D(){
//CkPrintf("\nCalling default construction for PMatrix2D(%d, %d)\n", thisIndex.x, thisIndex.y);
  GWBSE* gwbse = GWBSE::get();

  // Set some constants
  L = gwbse->gw_parallel.L;
  nfft = gwbse->gw_parallel.fft_nelems;
  qindex = Q_IDX; // Eventually the controller will set this

  // Grab a local pointer to the fft controller for fft-ing our rows
  // TODO: Is this guaranteed to be safe (is the local branch created for sure)?
  fft_controller = fft_controller_proxy.ckLocalBranch();

  // Figure out what part of the matrix we have and allocate space
  matrix_dimension = gwbse->gw_parallel.n_elems;
  num_rows = gwbse->gw_parallel.rows_per_chare;
  num_cols = gwbse->gw_parallel.cols_per_chare;
  num_chares = (matrix_dimension / num_rows) * (matrix_dimension / num_cols);
  start_row = thisIndex.x * num_rows;
  start_col = thisIndex.y * num_cols;
  local_mtx_size_1d_y = 1;//1728; // TODO: where from?
  receive_counter = 0;

  data = new complex[num_rows * num_cols];

  total_time = 0.0;
}

void PMatrix2D::generateEpsilon(std::vector<double> vcoulb, std::vector<int> accept, int inew, int jnew, int size, int max_inew, int max_jnew){
#if 1
  int n = 0; 
//  std::vector<complex> eps_data_tmp(num_rows*num_cols);
  complex *eps_data_tmp;
  eps_data_tmp = new complex[num_rows*num_cols];
  for(int i=0;i<num_rows;i++)
    for(int j=0;j<num_cols;j++)
      eps_data_tmp[n++] = data[IDX(i,j)];

  int *accept_arr;
  double *vcoulb_arr;
  accept_arr = new int[accept.size()];  
  vcoulb_arr = new double[vcoulb.size()];
  for(int i=0;i<accept.size();i++){
    accept_arr[i] = accept[i];
    vcoulb_arr[i] = vcoulb[i];
  }
 
  controller_proxy.inputData(thisIndex.x*12+thisIndex.y, eps_data_tmp, vcoulb_arr, accept_arr, num_rows, num_cols, accept.size());             
 // eps_matrix2D_proxy(thisIndex.x, thisIndex.y).receiveFs(eps_data_tmp, 0, 0, num_rows, num_cols,0);
#endif
  
}

 
void PMatrix2D::reportPTime() {
  CkReduction::statisticsElement stats(total_time);
  int tuple_size = 2;
  CkReduction::tupleElement tuple_reduction[] = {
    CkReduction::tupleElement(sizeof(double), &total_time, CkReduction::sum_double),
    CkReduction::tupleElement(sizeof(CkReduction::statisticsElement), &stats, CkReduction::statistics) };

  CkReductionMsg* msg = CkReductionMsg::buildFromTuple(tuple_reduction, tuple_size);
  msg->setCallback(CkCallback(CkIndex_Controller::reportPTime(NULL), controller_proxy));
  contribute(msg);
}

void PMatrix2D::applyFs() {
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
        data[IDX(r,c)] += f[r+start_row]*f[c+start_col].conj() * -1.0;
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

void PMatrix2D::checkReady(){
    //contribute(CkCallback(CkReductionTarget(Controller, allReady), controller_proxy));
}
void PMatrix2D::sendTo1D() {
//  CkPrintf("\nIn sendTo1D\n");
  int local_mtx_size_x = num_cols;
  int local_mtx_size_y = num_rows;
  int global_x = local_mtx_size_x*thisIndex.x;
  for (unsigned i=0;i<local_mtx_size_y;++i) {
    int global_y = local_mtx_size_y*thisIndex.y + i;
    int dest_chare = (int) ((double)global_y/local_mtx_size_1d_y);
    Phase2Message* msg;
    msg = new (local_mtx_size_x) Phase2Message();
    msg->global_x = global_x;
    msg->global_y = global_y;
    msg->size = local_mtx_size_x;
    for(unsigned c = 0; c < local_mtx_size_x; ++c){
      msg->data[c] = data[IDX(i,c)];
    }
    pmatrix1D_proxy[dest_chare].receiveRow(msg);
  }
}

void PMatrix2D::receiveChunk(Phase2Message* msg) {
  int local_mtx_size_x = num_cols;
  int local_mtx_size_y = num_rows;
  int local_y = msg->global_y - thisIndex.y * local_mtx_size_y;
  int local_x = msg->global_x - thisIndex.x * local_mtx_size_x;
  for(unsigned i=0; i<msg->size; ++i){
    data[IDX(local_y, local_x) + i] = msg->data[i];
  }

  if(++receive_counter == num_rows){
    contribute(CkCallback(CkReductionTarget(Controller, phase2_complete), controller_proxy));
  }
  delete msg;
}


PMatrix1D::PMatrix1D(int local_size_x, int local_size_y) {
  GWBSE* gwbse = GWBSE::get();

  // Set some constants
  L = gwbse->gw_parallel.L;
  nfft = gwbse->gw_parallel.fft_nelems;
  qindex = Q_IDX; // Eventually the controller will set this

  // Grab a local pointer to the fft controller for fft-ing our rows
  // TODO: Is this guaranteed to be safe (is the local branch created for sure)?
  fft_controller = fft_controller_proxy.ckLocalBranch();

  // Figure out what part of the matrix we have and allocate space
  matrix_dimension = gwbse->gw_parallel.n_elems;
  num_rows = gwbse->gw_parallel.rows_per_chare;
  num_cols = gwbse->gw_parallel.cols_per_chare;
  num_chares = (matrix_dimension / num_rows) * (matrix_dimension / num_cols);

  // Do the same for 1D:
  local_mtx_size_1d_x = local_size_x;
  local_mtx_size_1d_y = local_size_y;

  local_mtx_size_2d_x = num_cols;
  local_mtx_size_2d_y = num_rows;

  number_of_chares_1d = matrix_dimension / local_mtx_size_1d_y;
  number_of_chares_2d_x = matrix_dimension / local_mtx_size_2d_x;
  number_of_chares_2d_y = matrix_dimension / local_mtx_size_2d_y;

  arrival_counter = 0;

  data = new complex[local_mtx_size_1d_x * local_mtx_size_1d_y];
}

void PMatrix1D::fftRows(int direction) {
  // FFT each row stored in this chare
  for (int i=0; i < local_mtx_size_1d_y; i++){
    // First set up the data structures in the FFTController
    fft_controller->setup_fftw_3d(nfft, direction);
    fftw_complex* in_pointer = fft_controller->get_in_pointer();
    fftw_complex* out_pointer = fft_controller->get_out_pointer();

    // Pack our data, do the fft, then get the output
    put_into_fftbox(nfft, &data[IDX(i,0)], in_pointer);
    fft_controller->do_fftw();
    fftbox_to_array(num_cols, out_pointer, &data[IDX(i,0)], 1);
  }

  // Let the controller know we have completed the fft
  contribute(CkCallback(CkReductionTarget(Controller, fftComplete), controller_proxy));
}



void PMatrix1D::sendTo2D() {
  int n_col1 = number_of_chares_2d_x;
  int num_rows=local_mtx_size_1d_y;
  int chunksize = local_mtx_size_1d_x / n_col1;
  Phase2Message* msg;
  for(unsigned rows=0;rows<local_mtx_size_1d_y; ++rows){
    for(unsigned dest_x=0;dest_x<n_col1;++dest_x){
      int global_x = dest_x*chunksize;
      int global_y = rows + thisIndex * local_mtx_size_1d_y;
      int dest_y = (int) ((double)global_y / local_mtx_size_2d_y);

      msg = new (chunksize) Phase2Message();
      for(unsigned x=0;x<chunksize;++x){
        msg->data[x] = data[IDX(rows,x+global_x)];
      }
      msg->global_y = global_y;
      msg->global_x = global_x;
      msg->size = chunksize;
      pmatrix2D_proxy(dest_x,dest_y).receiveChunk(msg);
    }
  }
}


void PMatrix1D::receiveRow(Phase2Message* msg) {
  int local_y = msg->global_y - thisIndex * local_mtx_size_1d_y;
  for(unsigned i=0; i< msg->size; ++i){
    data[IDX(local_y,msg->global_x + i)] = msg->data[i];
  }
  if(++arrival_counter == local_mtx_size_1d_y) {
    contribute(CkCallback(CkReductionTarget(Controller, dataSendComplete), controller_proxy));
  }
  delete msg;
}

// TODO: These methods shouldn't be part of PMatrix, and should also just be
// computed once at startup.
void PMatrix2D::kqIndex(unsigned ikpt, unsigned& ikq, int* uklapp){
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

void PMatrix2D::getUmklappFactor(complex* umklapp_factor, int uklpp[3]){

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

#include "pmatrix.def.h"
