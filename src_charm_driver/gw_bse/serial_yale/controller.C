#include "controller.h"

#include "standard_include_gwbse.h"
#include "allclass_gwbse.h"
#include "messages.h"
#include "pmatrix.h"
#include "states.h"

Controller::Controller() {
  GWBSE *gwbse = GWBSE::get();

  // Set our class variables
  K = gwbse->gw_parallel.K;
  L = gwbse->gw_parallel.L;
  M = gwbse->gw_parallel.M;
  pipeline_stages = gwbse->gw_parallel.pipeline_stages;

  next_K = next_state = total_sent = total_complete = 0;
}

void Controller::pComplete() {
  pmatrix_proxy.printRowAndExit(0);
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

complex* PsiCache::getPsi(unsigned q, unsigned ispin, unsigned ikpt, unsigned istate) const {
  // TODO: Minjung will use q and the three indices to decide which psi to return
  return psis[istate];
}

#include "controller.def.h"
