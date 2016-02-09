#include "controller.h"

#include "standard_include.h"
#include "allclass_gwbse.h"
#include "messages.h"
#include "pmatrix.h"
#include "states.h"
#include "fft_controller.h"

void init_plan_lock();

Controller::Controller() {
  GWBSE *gwbse = GWBSE::get();

  // Set our class variables
  K = gwbse->gw_parallel.K;
  L = gwbse->gw_parallel.L;
  M = gwbse->gw_parallel.M;
  pipeline_stages = gwbse->gw_parallel.pipeline_stages;

  next_K = next_state = total_sent = total_complete = 0;
}

PsiCache::PsiCache() {
  GWBSE *gwbse = GWBSE::get();
  K = gwbse->gw_parallel.K;
  L = gwbse->gw_parallel.L;
  qindex = 0;
  psi_size = gwbse->gw_parallel.n_elems;
  received_psis = 0;
  psis = new complex**[K];
  for (int k = 0; k < K; k++) {
    psis[k] = new complex*[L];
    for (int l = 0; l < L; l++) {
      psis[k][l] = new complex[psi_size];
    }
  }

  fs = new complex*[L];
  for (int l = 0; l < L; l++) {
    fs[l] = new complex[psi_size];
  }
}

void PsiCache::receivePsi(PsiMessage* msg) {
  if (msg->spin_index != 0) {
    CkAbort("Error: We don't support multiple spins yet!\n");
  }
  CkAssert(msg->k_index < K);
  CkAssert(msg->state_index < L);
  CkAssert(msg->size == psi_size);
  std::copy(msg->psi, msg->psi+psi_size, psis[msg->k_index][msg->state_index]);
  delete msg;

  // Once the cache has received all of it's data start the sliding pipeline
  // sending of psis to P to start the accumulation of fxf'.
  if (++received_psis == K*L) {
    CkPrintf("[%d]: Cache filled\n", CkMyPe());
    contribute(CkCallback(CkReductionTarget(Controller,cachesFilled), controller_proxy));
  }
}

void PsiCache::computeFs(PsiMessage* msg) {
  double end, start = CmiWallTimer();
  if (msg->spin_index != 0) {
    CkAbort("Error: We don't support multiple spins yet!\n");
  }
  CkAssert(msg->size == psi_size);

  unsigned ikq;
  int umklapp[3];
  kqIndex(msg->k_index, ikq, umklapp);

  bool uproc = false;
  complex* umklapp_factor = NULL;
  if (umklapp[0] != 0 || umklapp[1] != 0 || umklapp[2] != 0) {
    uproc = true;
    umklapp_factor = new complex[psi_size];
    getUmklappFactor(umklapp_factor, umklapp);
  }

  complex* psi_occ;
  complex* psi_unocc = msg->psi;
  for (int l = 0; l < L; l++) {
    psi_occ = psis[ikq][l];
    for (int i = 0; i < psi_size; i++) {
      fs[l][i] = psi_occ[i] * psi_unocc[i].conj();
      if (uproc) {
        fs[l][i] *= umklapp_factor[i];
      }
    }
  }

  int tmp[4] = { msg->spin_index, msg->k_index, msg->state_index - L, ikq };
  CkCallback cb(CkReductionTarget(PMatrix, applyFs), pmatrix_proxy);
  contribute(4*sizeof(int), tmp, CkReduction::nop, cb);

  delete [] umklapp_factor;
  delete msg;

  end = CmiWallTimer();
  if (CkMyNode() == 0) {
    CkPrintf("[PSICACHE] Computed fs in %fs\n", end - start);
  }
}

complex* PsiCache::getPsi(unsigned ispin, unsigned ikpt, unsigned istate) const {
  if (ispin != 0) {
    CkAbort("Error: We don't support multiple spins yet!\n");
  }
  CkAssert(ikpt >= 0 && ikpt < K);
  CkAssert(istate >= 0 && istate < L);
  return psis[ikpt][istate];
}

complex* PsiCache::getF(unsigned idx) const {
  CkAssert(idx >= 0 && idx < L);
  return fs[idx];
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


void PsiCache::getUmklappFactor(complex* umklapp_factor, int uklpp[3]){

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

#include "controller.def.h"
