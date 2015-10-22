#include "controller.h"

#include "standard_include_gwbse.h"
#include "allclass_gwbse.h"
#include "messages.h"
#include "states.h"

Controller::Controller() {
  GWBSE *gwbse = GWBSE::get();

  // Set our class variables
  L = gwbse->gw_parallel.L;
  M = gwbse->gw_parallel.M;
  pipeline_stages = gwbse->gw_parallel.pipeline_stages;

  next_state = 0;

  // First fill the cache with the occupied psis
  for (next_state = 0; next_state < L; next_state++) {
    states_proxy(0,0,next_state).sendToCache();
  }
}

void Controller::cachesFilled() {
  for (int i = 0; i < pipeline_stages; i++) {
    states_proxy(0,0,next_state++).sendToP();
  }
}

void Controller::psiComplete() {
  if (next_state < L+M) {
    states_proxy(0, 0, next_state++).sendToP();
  } else {
    CkExit();
  }
}

PsiCache::PsiCache() {
  GWBSE *gwbse = GWBSE::get();
  psi_count = gwbse->gw_parallel.L;
  psi_size = gwbse->gw_parallel.n_elems;
  received_psis = 0;
  psis = new complex*[psi_count];
  for (int i = 0; i < psi_count; i++) {
    psis[i] = new complex[psi_size];
  }
}

void PsiCache::receivePsi(PsiMessage* msg) {

  CkAssert(msg->state_index < psi_count);
  CkAssert(msg->size == psi_size);
  std::copy(msg->psi, msg->psi+psi_size, psis[msg->state_index]);
  delete msg;

  // Once the cache has received all of it's data start the sliding pipeline
  // sending of psis to P to start the accumulation of fxf'.
  if (++received_psis == psi_count) {
    CkPrintf("[%d]: Cache filled\n", CkMyPe());
    contribute(CkCallback(CkReductionTarget(Controller,cachesFilled), controller_proxy));
  }

}

complex* PsiCache::getPsi(unsigned index) const {
  CkAssert(index < psi_count);
  return psis[index];
}

#include "controller.def.h"
