#include "psi.h"
#include "gw_bse.h"

Psi::Psi(unsigned s) : size(s) {
  psi = new double[size];
  // TODO: Read in the Psi somehow or get it from the driver?
  // NOTE: For now, using randomly generated psi for performance impact
  for (int i = 0; i < size; i++) {
    psi[i] = ((double) std::rand() / RAND_MAX);
  }

  calculatePsiR();
  if (thisIndex.y < config.L) {
    sendToCache();
  }
}

void Psi::calculatePsiR() {
  // TODO (Yale): Take the psi array in this object and IFFT it to take it to
  // real space
  CkPrintf("[%i,%i]: Calculating real-space psi...\n",thisIndex.x, thisIndex.y);
}

void Psi::sendToCache() {
  CkPrintf("[%i,%i]: Sending psi to node cache...\n",thisIndex.x, thisIndex.y);
  PsiMessage* msg = new (size) PsiMessage(size, psi);
  msg->k_index = thisIndex.x;
  msg->state_index = thisIndex.y;
  psicache.receivePsi(msg);
}

void Psi::sendToP() {
  CkPrintf("[%i,%i]: Sending psi to P matrix...\n",thisIndex.x, thisIndex.y);
  PsiMessage* msg = new (size) PsiMessage(size, psi);
  msg->k_index = thisIndex.x;
  msg->state_index = thisIndex.y;
  pmatrix.receivePsi(msg);
}

PsiCache::PsiCache(unsigned c, unsigned s) : psi_count(c), psi_size(s) {
  received_psis = 0;
  psis = new double*[psi_count];
  for (int i = 0; i < psi_count; i++) {
    psis[i] = new double[psi_size];
  }
}

void PsiCache::receivePsi(PsiMessage* msg) {
  CkPrintf("[%i]: Receiving psi from [%i,%i]\n",CkMyPe(), msg->k_index, msg->state_index);
  CkAssert(msg->state_index < psi_count);
  CkAssert(msg->size == psi_size);
  std::copy(psis[msg->state_index], psis[msg->state_index]+psi_size, msg->psi);

  // Once the cache has received all of it's data start the sliding pipeline
  // sending of psis to P to start the accumulation of fxf'.
  if (++received_psis == psi_count) {
    for (int i = 0; i < config.pipeline_stages; i++) {
      contribute(CkCallback(CkReductionTarget(Psi,sendToP), psi(0,config.L + i)));
    }
  }

  delete msg;
}

double* PsiCache::getPsi(unsigned index) const {
  CkAssert(index < psi_count);
  return psis[index];
}
