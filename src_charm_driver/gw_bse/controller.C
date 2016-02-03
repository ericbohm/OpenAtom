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
  psi_size = gwbse->gw_parallel.n_elems;
  received_psis = 0;
  psis = new complex**[K];
  for (int k = 0; k < K; k++) {
    psis[k] = new complex*[L];
    for (int l = 0; l < L; l++) {
      psis[k][l] = new complex[psi_size];
    }
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

complex* PsiCache::getPsi(unsigned ispin, unsigned ikpt, unsigned istate) const {
  if (ispin != 0) {
    CkAbort("Error: We don't support multiple spins yet!\n");
  }
  CkAssert(ikpt >= 0 && ikpt < K);
  CkAssert(istate >= 0 && istate < L);
  return psis[ikpt][istate];
}

#include "controller.def.h"
