#include "controller.h"

#include "standard_include.h"
#include "allclass_gwbse.h"
#include "messages.h"
#include "eps_matrix.h"
#include "pmatrix.h"
#include "mat_mul.h"
#include "main.h"
#include "states.h"
#include "fft_controller.h"
#include "CkLoopAPI.h"

#define eps_rows 20
#define eps_cols 20

void init_plan_lock();

Controller::Controller() {
  GWBSE *gwbse = GWBSE::get();

  // Set our class variables
  K = gwbse->gw_parallel.K;
  L = gwbse->gw_parallel.L;
  M = gwbse->gw_parallel.M;
  pipeline_stages = gwbse->gw_parallel.pipeline_stages;

  next_K = next_state = total_sent = total_complete = next_report_threshold = 0;

  dimension = gwbse->gw_parallel.n_elems;
  rows = gwbse->gw_parallel.rows_per_chare;

  epsCut = 5;
  alat = 10.261200; 
  vol = 10;//how to read this?
  shift[0] = 0;
  shift[1] = 0;
  shift[2] = 0.001;
  // TODO: Make these config options
  do_output = true;
  max_sends = M*K;  // For debugging this can be changed to a smaller number
  maxiter = 1;
  msg_received = 0;
  global_inew = 0;
  max_local_inew = global_inew;
  global_jnew = 0;
}


void Controller::prep(){

  GWBSE *gwbse = GWBSE::get();

  double *this_q, *b1, *b2, *b3;
  b1 = gwbse->gwbseopts.b1;
  b2 = gwbse->gwbseopts.b2;
  b3 = gwbse->gwbseopts.b3;

  int qindex = Q_IDX;

  this_q = gwbse->gwbseopts.qvec[qindex];


  int* nfft;
  nfft = gwbse->gw_parallel.fft_nelems;

  int ndata = nfft[0]*nfft[1]*nfft[2];
  FFTController* fft_controller = fft_controller_proxy.ckLocalBranch();

//output geps, accept
  fft_controller->get_geps(epsCut, this_q, b1, b2, b3,
                            alat, nfft);
  }

  void Controller::got_geps(std::vector<int> accept, int epsilon_size){
        //CkPrintf("\nGot geps");fflush(stdout);
        accept_result = accept;
        thisProxy.done_geps(epsilon_size);
  }
  
  void Controller::calc_Geps(){

  
  GWBSE *gwbse = GWBSE::get();

  double *this_q, *b1, *b2, *b3;
  b1 = gwbse->gwbseopts.b1;
  b2 = gwbse->gwbseopts.b2;
  b3 = gwbse->gwbseopts.b3;

  int qindex = Q_IDX;

  this_q = gwbse->gwbseopts.qvec[qindex];


  int* nfft;
  nfft = gwbse->gw_parallel.fft_nelems;

  int ndata = nfft[0]*nfft[1]*nfft[2];
  FFTController* fft_controller = fft_controller_proxy.ckLocalBranch();


//output - vcoulb
  fft_controller->calc_vcoulb(this_q, b1, b2, b3, shift, alat, vol, gwbse->gwbseopts.nkpt, qindex);
}

void Controller::got_vcoulb(std::vector<double> vcoulb_in){

  vcoulb = vcoulb_in;
  thisProxy.prepare_epsilon();
}

PsiCache::PsiCache() {
  GWBSE *gwbse = GWBSE::get();
  K = gwbse->gw_parallel.K;
  L = gwbse->gw_parallel.L;
  qindex = Q_IDX;
  psi_size = gwbse->gw_parallel.n_elems;
  pipeline_stages = gwbse->gw_parallel.pipeline_stages;
  received_psis = 0;
  psis = new complex**[K];
  for (int k = 0; k < K; k++) {
    psis[k] = new complex*[L];
    for (int l = 0; l < L; l++) {
      psis[k][l] = new complex[psi_size];
    }
  }
  // shifted k grid psis. Need this for qindex=0
  psis_shifted = new complex**[K];
  for (int k = 0; k < K; k++) {
    psis_shifted[k] = new complex*[L];
    for (int l = 0; l < L; l++) {
      psis_shifted[k][l] = new complex[psi_size];
    }
  }

  fs = new complex[L*psi_size*pipeline_stages];

  umklapp_factor = new complex[psi_size];

  total_time = 0.0;
}

void PsiCache::reportFTime() {
  CkReduction::statisticsElement stats(total_time);
  int tuple_size = 2;
  CkReduction::tupleElement tuple_reduction[] = {
    CkReduction::tupleElement(sizeof(double), &total_time, CkReduction::max_double),
    CkReduction::tupleElement(sizeof(CkReduction::statisticsElement), &stats, CkReduction::statistics) };

  CkReductionMsg* msg = CkReductionMsg::buildFromTuple(tuple_reduction, tuple_size);
  msg->setCallback(CkCallback(CkIndex_Controller::reportFTime(NULL), controller_proxy));
  contribute(msg);
}

void PsiCache::receivePsi(PsiMessage* msg) {
  if (msg->spin_index != 0) {
    CkAbort("Error: We don't support multiple spins yet!\n");
  }
  CkAssert(msg->k_index < K);
  CkAssert(msg->state_index < L);
  CkAssert(msg->size == psi_size);
  if(msg->shifted==false){std::copy(msg->psi, msg->psi+psi_size, psis[msg->k_index][msg->state_index]);}
  if(msg->shifted==true){std::copy(msg->psi, msg->psi+psi_size, psis_shifted[msg->k_index][msg->state_index]);}
  delete msg;

  // Once the cache has received all of it's data start the sliding pipeline
  // sending of psis to P to start the accumulation of fxf'.
  int expected_psis = K*L;
  if(qindex == 0)
    expected_psis += K*L;
  if (++received_psis == expected_psis) {
    //CkPrintf("[%d]: Cache filled\n", CkMyPe());
    contribute(CkCallback(CkReductionTarget(Controller,cachesFilled), controller_proxy));
  }
}

// Called by CkLoop to spread the computation of f vectors across the node
void computeF(int first, int last, void* result, int count, void* params) {
  FComputePacket* f_packet = (FComputePacket*)params;
  unsigned psi_size = f_packet->size;
  complex* psi_unocc = f_packet->unocc_psi;
  complex* umklapp_factor = f_packet->umklapp_factor;
  double* e_occ = f_packet->e_occ;
  double e_unocc = f_packet->e_unocc;
  complex* fs = f_packet->fs;

  for (int l = first; l <= last; l++) {
    complex* f = &(fs[l*psi_size]);
    complex* psi_occ = f_packet->occ_psis[l];
    double scaling_factor = 2/sqrt(e_unocc - e_occ[l]);

    for (int i = 0; i < psi_size; i++) {
      f[i] = psi_occ[i] * psi_unocc[i].conj() * scaling_factor;
      if (umklapp_factor) {
        f[i] *= umklapp_factor[i];
      }
#ifdef USE_LAPACK
      // BLAS calls compute the complex conjugate of P, which is hermitian. This
      // change to f corrects that so we get the correct P.
      f[i] = f[i].conj();
#endif
    }
  }
}

// Receive an unoccupied psi, and split off the computation of all associated f
// vectors across the node using CkLoop.
void PsiCache::computeFs(PsiMessage* msg) {
  double start = CmiWallTimer();

  if (msg->spin_index != 0) {
    CkAbort("Error: We don't support multiple spins yet!\n");
  }
  CkAssert(msg->size == psi_size);

  // Compute ikq index and the associated umklapp factor
  // TODO: This should just be a table lookup
  unsigned ikq;
  int umklapp[3];
  kqIndex(msg->k_index, ikq, umklapp);

  bool uproc = false;
  if (umklapp[0] != 0 || umklapp[1] != 0 || umklapp[2] != 0) {
    uproc = true;
    computeUmklappFactor(umklapp);
  }

  GWBSE* gwbse = GWBSE::get();
  double*** e_occ = gwbse->gw_epsilon.Eocc;
  double*** e_occ_shifted = gwbse->gw_epsilon.Eocc_shifted;
  double*** e_unocc = gwbse->gw_epsilon.Eunocc;

  // Create the FComputePacket for this set of f vectors and start CkLoop
  f_packet.size = psi_size;
  f_packet.unocc_psi = msg->psi;
  if ( qindex == 0 ) { 
    f_packet.occ_psis = psis_shifted[ikq]; 
    f_packet.e_occ = e_occ_shifted[msg->spin_index][ikq];
  }
  else { 
    f_packet.occ_psis = psis[ikq];
    f_packet.e_occ = e_occ[msg->spin_index][ikq]; 
  }
  f_packet.e_unocc = e_unocc[msg->spin_index][msg->k_index][msg->state_index-L];
  f_packet.fs = fs + (L*psi_size*(received_chunks%pipeline_stages));

  if (uproc) { f_packet.umklapp_factor = umklapp_factor; }
  else { f_packet.umklapp_factor = NULL; }

#ifdef USE_CKLOOP
  CkLoop_Parallelize(computeF, 1, &f_packet, L, 0, L - 1);
#else
  for (int l = 0; l < L; l++) {
    computeF(l,l,NULL,1,&f_packet);
  }
#endif
  received_chunks++;

  FVectorCache *fvec_cache = fvector_cache_proxy.ckLocalBranch();
  Phase4Message *fvec_msg;
  fvec_msg = new(L*psi_size) Phase4Message();
  fvec_msg->n = msg->state_index-L;
  fvec_msg->data = fs;
  fvec_cache->putFVec(fvec_msg);
  // Let the matrix chares know that the f vectors are ready
  CkCallback cb(CkReductionTarget(PMatrix2D, applyFs), pmatrix2D_proxy);
  contribute(cb);

  // Cleanup
  delete msg;
  total_time += CmiWallTimer() - start;
}

complex* PsiCache::getPsi(unsigned ispin, unsigned ikpt, unsigned istate) const {
  if (ispin != 0) {
    CkAbort("Error: We don't support multiple spins yet!\n");
  }
  CkAssert(ikpt >= 0 && ikpt < K);
  CkAssert(istate >= 0 && istate < L);
  return psis[ikpt][istate];
}

complex* PsiCache::getF(unsigned idx, unsigned req_no) const {
  CkAssert(idx >= 0 && idx < L);
  CkAssert(req_no < received_chunks && req_no >= received_chunks - pipeline_stages);
  return &(fs[idx*psi_size+(L*psi_size*(req_no%pipeline_stages))]);
}

void PsiCache::kqIndex(unsigned ikpt, unsigned& ikq, int* uklapp){
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


void PsiCache::computeUmklappFactor(int uklpp[3]){

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

FVectorCache::FVectorCache() {
  GWBSE *gwbse = GWBSE::get();
  L = gwbse->gw_parallel.L;
  psi_size = gwbse->gw_parallel.n_elems;
  fcount = 0;
  n_list_size = 10; //TODO: to be greater than 0, read in from user input file
  fs = new complex[L*psi_size*n_list_size];//
}

void FVectorCache::putFVec(Phase4Message *msg) {
  if(fcount==0)
    CkPrintf("\nCaching f[%d*L] vectors\n", msg->n);
//Simple scheme where all Fs are stored in each node
//  if(inList(msg->n))
  { //Store all n*l's for all l in L
    complex* f = &(fs[fcount*L*psi_size]);
    f = msg->data;
  }

  fcount++;
//TODO: Partition Fs and duplicate to reduce access latency
}

#include "psi_cache.def.h"
#include "fvector_cache.def.h"
#include "fft_controller.def.h"
#include "controller.def.h"
