#include "controller.h"

#include "standard_include_gwbse.h"
#include "allclass_gwbse.h"
#include "messages.h"
#include "states.h"

PsiCache::PsiCache() {
  GWBSE *gwbse = GWBSE::get();
  psi_count = gwbse->gw_parallel.L;
  psi_size = gwbse->gw_parallel.n_elems;
  pipeline_stages = gwbse->gw_parallel.pipeline_stages;
  received_psis = 0;
  psis = new complex*[psi_count];
  for (int i = 0; i < psi_count; i++) {
    psis[i] = new complex[psi_size];
  }
}

void PsiCache::receivePsi(PsiMessage* msg) {

  CkAssert(msg->state_index < psi_count);
  CkAssert(msg->size == psi_size);
  std::copy(psis[msg->state_index], psis[msg->state_index]+psi_size, msg->psi);
  delete msg;

  // Once the cache has received all of it's data start the sliding pipeline
  // sending of psis to P to start the accumulation of fxf'.
  if (++received_psis == psi_count) {
    CkPrintf("[%d]: Cache filled\n", CkMyPe());
    for (int i = 0; i < pipeline_stages; i++) {
      contribute(CkCallback(CkReductionTarget(States,sendToP), states_proxy(0,0,psi_count + i)));
    }
    //contribute(CkCallback(CkCallback::ckExit));
  }

}

complex* PsiCache::getPsi(unsigned index) const {
  CkAssert(index < psi_count);
  return psis[index];
}

#include "controller.def.h"
