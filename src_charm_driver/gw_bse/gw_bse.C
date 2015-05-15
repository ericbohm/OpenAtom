#include "gw_bse.h"

/*readonly*/int K;
/*readonly*/int Q;
/*readonly*/int L;
/*readonly*/int M;

/*readonly*/int psi_size;
/*readonly*/int pipeline_stages;

/*readonly*/CkGroupID mcast_ID;
/*readonly*/CProxy_Psi kpsi;
/*readonly*/CProxy_Psi qpsi;
/*readonly*/CProxy_FCalculator fcalc;

GWBSEDriver::GWBSEDriver(CkArgMsg* msg) {
  readState();

  K = 4;
  Q = 4;
  L = 8;
  M = 16;

  psi_size = 16;
  pipeline_stages = 2;

  mcast_ID = CProxy_CkMulticastMgr::ckNew();
  kpsi = CProxy_Psi::ckNew(true, K, L);
  qpsi = CProxy_Psi::ckNew(false, Q, M);

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
  qpsi.sendPsi();
  fcalc.run();
}

// TODO: This function should maybe be distributed to the Psi chares in some way
void GWBSEDriver::readState() {
}

Psi::Psi(bool occupied) : occupied(occupied) {
  psi = new double[psi_size];
  calculatePsiR();
}

void Psi::createSections() {
  CkMulticastMgr* mcast_ptr = CProxy_CkMulticastMgr(mcast_ID).ckLocalBranch();
  CProxySection_FCalculator section;
  unsigned k_lower = 0, k_upper = K - 1;
  unsigned q_lower = 0, q_upper = Q - 1;
  unsigned l_lower = 0, l_upper = L - 1;
  unsigned m_lower = 0, m_upper = M - 1;
  if (occupied) {
    k_lower = k_upper = thisIndex.x;
    l_lower = l_upper = thisIndex.y;
  } else {
    q_lower = q_upper = thisIndex.x;
    m_lower = m_upper = thisIndex.y;
  }
  for (section_index = l_lower; section_index <= l_upper; section_index++) {
    section = CProxySection_FCalculator::ckNew( fcalc.ckGetArrayID(),
                                                k_lower, k_upper, 1,
                                                q_lower, q_upper, 1,
                                                section_index, section_index, 1,
                                                m_lower, m_upper, 1 );
    section.ckSectionDelegate(mcast_ptr);
    sections.push_back(std::make_pair(section_index, section));
  }
  section_index = 0;
}

void Psi::calculatePsiR() {
 // TODO (Yale): Take the psi array in this object and IFFT it to take it to
 // real space
} 

FCalculator::FCalculator() {
  mcast_ptr = CProxy_CkMulticastMgr(mcast_ID).ckLocalBranch();
}

void FCalculator::createSections() {
  CProxySection_FCalculator red_section;
  red_section = CProxySection_FCalculator::ckNew( thisProxy.ckGetArrayID(),
                                                  0, K-1, 1,
                                                  0, Q-1, 1,
                                                  thisIndex.y, thisIndex.y, 1,
                                                  0, M-1, 1 );
  red_section.ckSectionDelegate(mcast_ptr);
  PlaneMsg* msg = new PlaneMsg();
  msg->plane_index = thisIndex.y;
  red_section.receiveSectionInfo(msg);
}

void FCalculator::computeF(double* psi1, double* psi2) {
  // TODO (Yale): Take the two psis and compute f and store it in the field
}

#include "gw_bse.def.h"
