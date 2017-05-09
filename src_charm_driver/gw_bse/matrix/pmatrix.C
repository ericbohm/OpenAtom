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
#define IDX(r,c) ((r)*config.tile_cols + (c))
#define IDX_eps(r,c) ((r)*eps_cols + (c))
#define eps_chares_x 10
#define eps_chares_y 10
#define eps_rows 20
#define eps_cols 20

PMatrix::PMatrix(MatrixConfig config) : CBase_PMatrix(config) {
  GWBSE* gwbse = GWBSE::get();

  // Set some constants
  L = gwbse->gw_parallel.L;
  nfft = gwbse->gw_parallel.fft_nelems;
  qindex = Q_IDX; // Eventually the controller will set this

  // Grab a local pointer to the fft controller for fft-ing our rows
  // TODO: Is this guaranteed to be safe (is the local branch created for sure)?
  fft_controller = fft_controller_proxy.ckLocalBranch();

  max_iterations = gwbse->gw_parallel.transpose_stages;
  iteration = 0;
  arrival_counter = 0;
  receive_counter = 0;
  total_time = 0.0;
}

void PMatrix::generateEpsilon(std::vector<int> accept){
  int inew = 0;
  for(int i=0;i<thisIndex.x;i++)
    if(accept[i])
      inew++;

  if(accept[thisIndex.x]){// && thisIndex==0){
    int jnew_local = 0;
    int jnew = 0;
    Phase3Message *msg;
    msg = new(eps_cols)Phase3Message();  
    for(int j=0;j<config.tile_cols;j++){

      if(accept[j]){
        msg->data[jnew_local++] = data[j];
        if ( inew == jnew )
          msg->data[jnew_local-1] += double(1);
        jnew++;
        if(jnew_local == eps_cols){
          int dest_chare_x = inew/eps_rows;
          int dest_chare_y = (jnew/eps_cols)-1;
          msg->start_i = inew%eps_rows;
          msg->start_j = 0;
          msg->end_i = inew%eps_rows;
          msg->end_j = (jnew-1)%eps_cols;
//          CkPrintf("\nSending data to (%d,%d)\n", dest_chare_x, dest_chare_y);
          eps_matrix2D_proxy(dest_chare_x,dest_chare_y).receiveFs(msg);
          jnew_local = 0;
          msg = new(eps_cols)Phase3Message();  
        }
      }
    } 
    if(jnew_local != eps_cols){
      for(int i=jnew_local; i<eps_cols; i++){
        msg->data[i] = 0;
        jnew++;
      }
        
      int dest_chare_x = inew/eps_rows;
      int dest_chare_y = (jnew/eps_cols)-1;
          msg->start_i = inew%eps_rows;
          msg->start_j = 0;
          msg->end_i = inew%eps_rows;
          msg->end_j = (jnew-1)%eps_cols;
          eps_matrix2D_proxy(dest_chare_x,dest_chare_y).receiveFs(msg);
    }
  }
 
//  CkPrintf("\n local_mtx_size_2d_x = %d, config.tile_cols = %d\n", local_mtx_size_2d_x, config.tile_cols); 
//  int num = 1728;
  int counter = 0;
  if(thisIndex.x == config.tile_cols-1){
    
    int padded_send_size = inew + (eps_cols - (inew%eps_cols));
    int remainder = eps_cols - (inew+1)%eps_cols;
//    CkPrintf("\nSending to i=%d,j=%d\n", remainder, padded_send_size);
    for(int i=inew+1;counter<remainder;i++){
      counter++;
      Phase3Message *msg;
      msg = new(eps_cols)Phase3Message();
      int jnew_local = 0;
      int jnew = 0;
      msg->start_i = i%eps_rows;

      for(int j=0;j<padded_send_size;j++){
        msg->data[jnew_local++] = 0;
        jnew++;
        if(jnew_local == eps_cols){
          int dest_chare_x = i/eps_rows;
          int dest_chare_y = (jnew/eps_cols)-1;
          msg->start_i = i%eps_rows;
          msg->start_j = 0;
          msg->end_i = i%eps_rows;
          msg->end_j = (jnew_local-1)%eps_cols;
 //         CkPrintf("\nSending data to (%d,%d)\n", dest_chare_x, dest_chare_y);
          eps_matrix2D_proxy(dest_chare_x,dest_chare_y).receiveFs(msg);
          jnew_local = 0;
          msg = new(eps_cols)Phase3Message();
        }
      }
    }
  }
}
 
void PMatrix::reportPTime() {
  CkReduction::statisticsElement stats(total_time);
  int tuple_size = 2;
  CkReduction::tupleElement tuple_reduction[] = {
    CkReduction::tupleElement(sizeof(double), &total_time, CkReduction::sum_double),
    CkReduction::tupleElement(sizeof(CkReduction::statisticsElement), &stats, CkReduction::statistics) };

  CkReductionMsg* msg = CkReductionMsg::buildFromTuple(tuple_reduction, tuple_size);
  msg->setCallback(CkCallback(CkIndex_Controller::reportPTime(NULL), controller_proxy));
  contribute(msg);
}

void PMatrix::applyFs() {
  double start = CmiWallTimer();

  PsiCache* psi_cache = psi_cache_proxy.ckLocalBranch();

#ifdef USE_LAPACK
  // Common variables for both ZGERC and ZGEMM
  int M = config.tile_rows, N = config.tile_cols;
  complex alpha = -1.0;
#ifdef USE_ZGEMM
  int K = L; // If using ZGEMM, we compute all outer products with one call
  int LDF = config.mat_rows; // Leading dimension of fs
  complex beta = 1.0;
  char opA = 'N', opB = 'C';
  complex* fs = psi_cache->getF(0,completed_chunks);
  ZGEMM(&opA, &opB, &N, &M, &K,
    &alpha, &(fs[start_col]), &LDF,
    &(fs[start_row]), &LDF,
    &beta, data, &N);
#else
  int K = 1; // If using ZGERC, we compute each outer product one at a time
  for (int l = 0; l < L; l++) {
    complex* f = psi_cache->getF(l,completed_chunks);
    ZGERC(&N, &M, &alpha, &(f[start_col]), &K, &(f[start_row]), &K, data, &N);
  }
#endif // endif for ifdef USE_ZGEMM
#else
  for (int l = 0; l < L; l++) {
    complex* f = psi_cache->getF(l,completed_chunks);
    for (int r = 0; r < config.tile_rows; r++) {
      for (int c = 0; c < config.tile_cols; c++) {
        data[IDX(r,c)] += f[r+start_row]*f[c+start_col].conj() * -1.0;
      }
    }
  }
#endif // endif for ifdef USE_LAPACK
  completed_chunks++;
  contribute(CkCallback(CkReductionTarget(Controller, psiComplete), controller_proxy));
  total_time += CmiWallTimer() - start;
}

void PMatrix::fftRows(int direction) {
  if (config.chareCols() != 1) {
    CkAbort("FFT not supported for 2D decompositions\n");
  }
  // FFT each row stored in this chare
  for (int i=0; i < config.tile_rows; i++){
    // First set up the data structures in the FFTController
    fft_controller->setup_fftw_3d(nfft, direction);
    fftw_complex* in_pointer = fft_controller->get_in_pointer();
    fftw_complex* out_pointer = fft_controller->get_out_pointer();

    // Pack our data, do the fft, then get the output
    put_into_fftbox(nfft, &data[IDX(i,0)], in_pointer);
    fft_controller->do_fftw();
    fftbox_to_array(config.tile_cols, out_pointer, &data[IDX(i,0)], 1);
  }

  // Let the controller know we have completed the fft
  contribute(CkCallback(CkReductionTarget(Controller, fftComplete), controller_proxy));
}

// TODO: These methods shouldn't be part of PMatrix, and should also just be
// computed once at startup.
void PMatrix::kqIndex(unsigned ikpt, unsigned& ikq, int* uklapp){
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

void PMatrix::getUmklappFactor(complex* umklapp_factor, int uklpp[3]){

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
