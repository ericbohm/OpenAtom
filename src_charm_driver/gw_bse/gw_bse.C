#include "gw_bse.h"

/*readonly*/int K;
/*readonly*/int Q;
/*readonly*/int L;
/*readonly*/int M;

/*readonly*/int psi_size;
/*readonly*/int pipeline_stages;

/*readonly*/CProxy_KPsi kpsi;
/*readonly*/CProxy_QPsi qpsi;
/*readonly*/CProxy_FCalculator fcalc;

GWBSEDriver::GWBSEDriver(CkArgMsg* msg) {
  readState();

  K = 4;
  Q = 4;
  L = 10;
  M = 100;

  psi_size = 16;
  pipeline_stages = 1;

  kpsi = CProxy_KPsi::ckNew(K, L);
  qpsi = CProxy_KPsi::ckNew(Q, M);

  fcalc = CProxy_FCalculator::ckNew();
  for (int k = 0; k < K; k++) {
    for (int q = 0; q < Q; q++) {
      for (int l = 0; l < L; l++) {
        for (int m = 0; m < M; m++) {
          fcalc(k,q,l,m).insert();
        }
      }
    }
  }
  fcalc.doneInserting();

  kpsi.sendPsi();
  //qpsi.sendPsi();
  fcalc.run();
}

// TODO: This function should maybe be distributed to the Psi chares in some way
void GWBSEDriver::readState() {
}

KPsi::KPsi() {
  calculatePsiR();
  CkPrintf("[%d,%d]: Constructed\n", thisIndex.x, thisIndex.y);
}

void KPsi::createSections() {
  CkPrintf("[%d,%d]: Creating sections\n", thisIndex.x, thisIndex.y);
  int k   = thisIndex.x;
  int l   = thisIndex.y;
  section = CProxySection_FCalculator::ckNew( fcalc.ckGetArrayID(),
                                              k, k, 1,
                                              0, Q, 1,
                                              l, l, 1,
                                              0, M, 1 );
}

void KPsi::calculatePsiR() {
 // TODO (Yale): Take the psi array in this object and IFFT it to take it to
 // real space
} 

QPsi::QPsi() {
  calculatePsiR();
}

void QPsi::createSections() {
  sections.resize(L);
  int q   = thisIndex.x;
  int m   = thisIndex.y;
  for (int l = 0; l < L; l++) {
    sections[l] = CProxySection_FCalculator::ckNew( fcalc.ckGetArrayID(),
                                                    0, K, 1,
                                                    q, q, 1,
                                                    l, l, 1,
                                                    m, m, 1 );
  }
};

void QPsi::calculatePsiR() {
 // TODO (Yale): Take the psi array in this object and IFFT it to take it to
 // real space
}

FCalculator::FCalculator() {}

void FCalculator::computeF(double* psi1, double* psi2) {
  // TODO (Yale): Take the two psis and compute f and store it in the field
}

#include "gw_bse.def.h"
