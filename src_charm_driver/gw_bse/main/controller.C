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
#include "fft_routines.h"
#include "CkLoopAPI.h"

#define eps_rows 20
#define eps_cols 20
#define NSIZE 4

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
  contribute(CkCallback(CkReductionTarget(Controller,psiCacheReady), controller_proxy));
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


#ifdef TESTING
{
  FVectorCache *fvec_cache = fvector_cache_proxy.ckLocalBranch();
  fvec_cache->computeFTilde(fs);
//  fvec_cache->applyCutoff(msg->accept_size, msg->accept);
//  fvec_cache->init(140);
//compute ftilde first - similar to ckloop above for all L's
  fvec_cache->putFVec(msg->state_index-L, fs);
}
#endif

  // Let the matrix chares know that the f vectors are ready
  CkCallback cb(CkReductionTarget(PMatrix, applyFs), pmatrix2D_proxy);
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

void FVectorCache::findIndices(){

  int count = 0;
  for(int i=0;i<eps_chares_x;i++){
    for(int j=0;j<eps_chares_y;j++){
      count++;
      if(count == my_chare_start+1){
        eps_start_chare_x = i;
        eps_start_chare_y = j;
      }
      if(count == my_chare_start+my_chare_count){
        eps_end_chare_x = i;
        eps_end_chare_y = j;
        return;
      }
    }
  }

  return;
}

FVectorCache::FVectorCache() {
  eps_chares_x = 7;
  eps_chares_y = 7;
  totalSize = 0;
  GWBSE *gwbse = GWBSE::get();
  L = gwbse->gw_parallel.L;
  int total_eps_chares = eps_chares_x*eps_chares_y;

  my_chare_count = total_eps_chares/CkNumNodes();

  my_chare_start = CkMyNode()*my_chare_count;
  int remaining = total_eps_chares%CkNumNodes();

  if(CkMyNode()>0)
    my_chare_start += remaining;

  if(CkMyNode()==0)
    my_chare_count += remaining;

  my_eps_chare_indices_x = new int[my_chare_count];
  my_eps_chare_indices_y = new int[my_chare_count];

  findIndices();
  int count = 0;
  for(int i=eps_start_chare_x;i<=eps_end_chare_x;i++){
    int j = 0;
    if(i==eps_start_chare_x)
      j = eps_start_chare_y;
    int j_end = eps_chares_y-1;
    if(i==eps_end_chare_x)
      j_end = eps_end_chare_y;
    while(j<=j_end){
      my_eps_chare_indices_x[count] = i;
      my_eps_chare_indices_y[count++] = j;
      j++;
    }
  }

  ndata = gwbse->gw_parallel.n_elems;
  data_size_x = ndata/eps_chares_x;
  if(ndata%eps_chares_x > 0)
    data_size_x += 2;
  data_size_y = ndata/eps_chares_y;
    if(ndata%eps_chares_y > 0)
      data_size_y += 2;
  data_offset_x = new int[my_chare_count];
  data_offset_y = new int[my_chare_count];

  for(int i=0;i<my_chare_count;i++){
    data_offset_x[i] = my_eps_chare_indices_x[i]*data_size_x;
    data_offset_y[i] = my_eps_chare_indices_y[i]*data_size_y;
  }

  int size_x = data_size_x;
  int size_y = data_size_y;
  local_offset =  new int[my_chare_count*2];
  global_offset = new int[my_chare_count*2];
  for(int i=0;i<my_chare_count;i++){
    global_offset[2*i] = data_offset_x[i];//totalSize;
    local_offset[2*i] = totalSize;
    totalSize += size_x;

    global_offset[2*i+1] = data_offset_y[i];//totalSize;
    local_offset[2*i+1] = totalSize;
    totalSize += size_y;
  }

  fs = new complex[NSIZE*L*totalSize];

  contribute(CkCallback(CkReductionTarget(Controller,fCacheReady), controller_proxy));
}

void FVectorCache::init(int size_xy){
}

void FVectorCache::putFVec(int n, complex* fs_input){ //fs_input has all L's corresponding to n
 sum = (0.0,0.0);
 for(int i=0;i<my_chare_count;i++){
    if(my_eps_chare_indices_x[i] == my_eps_chare_indices_y[i])
    for(int l=0;l<L;l++){
      int global_x = global_offset[2*i];
      global_x += l*ndata;
      int local_x = local_offset[2*i];
      local_x += n*L*totalSize + l*totalSize;

      complex *store_x = &fs[local_x];
      complex *load_x = &fs_input[global_x];
      for(int ii=0;ii<data_size_x;ii++)
        store_x[ii] = load_x[ii];

      int global_y = global_offset[2*i+1];
      global_y += l*ndata;
      int local_y = local_offset[2*i+1];
      local_y += n*L*totalSize + l*totalSize;

      complex *store_y = &fs[local_y];
      complex *load_y = &fs_input[global_y];
      for(int ii=0;ii<data_size_x;ii++)
        store_y[ii] = load_y[ii];
    }
  }
}

complex* FVectorCache::getFVec(int n, int l, int chare_start_index, int size){
  for(int i=0;i<my_chare_count;i++){
    if(my_eps_chare_indices_x[i] == my_eps_chare_indices_y[i] && chare_start_index == my_eps_chare_indices_x[i]){
      int local_x = local_offset[2*i];
      local_x += n*L*totalSize + l*totalSize;
      complex *f = &fs[local_x];
      return f;
    }
  }
  return NULL;
}


// Called by CkLoop to spread the computation of f vectors across the node
void fTildeWorkUnit(int first, int last, void* result, int count, void* params) {

  FComputePacket* f_packet = (FComputePacket*)params;
  complex* fs = f_packet->fs;
  int psi_size = f_packet->size;
  GWBSE *gwbse = GWBSE::get();
  int* nfft;
  nfft = gwbse->gw_parallel.fft_nelems;
  int vector_count = 1;
  int direction = -1;
  int L = gwbse->gw_parallel.L;
  FFTController* fft_controller = fft_controller_proxy.ckLocalBranch();


  for (int i=0; i < L; i++){ //for all the L*n_list_size, L are computed in this node
    // First set up the data structures in the FFTController
    fft_controller->setup_fftw_3d(nfft, direction);
    fftw_complex* in_pointer = fft_controller->get_in_pointer();
    fftw_complex* out_pointer = fft_controller->get_out_pointer();

    // Pack our data, do the fft, then get the output
    put_into_fftbox(nfft, &fs[i*psi_size], in_pointer);
    fft_controller->do_fftw();
    fftbox_to_array(psi_size, out_pointer, &fs[i*psi_size], direction); //Now cached on the same partitions
    // replace f_vector to f_tilde_vector
  }
}

//Each node calculates its own ftilde
void FVectorCache::computeFTilde(complex *fs_in){

  // Create the FComputePacket for this set of f vectors and start CkLoop
  f_packet.size = ndata;
  f_packet.fs = fs_in;
  
  
#ifdef USE_CKLOOP
  CkLoop_Parallelize(fTildeWorkUnit, 1, &f_packet, n_list_size, 0, n_list_size - 1);
#else
    fTildeWorkUnit(0,0,NULL,1,&f_packet);
#endif
}

void FVectorCache::applyCutoff(int size, int* accept){

  int inew = 0;

  for(int i=0;i<size;i++){
    if(accept[i]){
      fs[inew] = fs[i];
      inew++;
    }
  }
//  fs.resize(inew);
}

#include "psi_cache.def.h"
#include "fvector_cache.def.h"
#include "fft_controller.def.h"
#include "controller.def.h"
