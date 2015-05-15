#include "gw_bse.h"

/*readonly*/GWConfig config;

/*readonly*/CkGroupID mcast_ID;
/*readonly*/CProxy_Psi kpsi;
/*readonly*/CProxy_Psi qpsi;
/*readonly*/CProxy_FCalculator fcalc;

GWDriver::GWDriver(CkArgMsg* msg) {
  readConfig();
  readState();

  mcast_ID = CProxy_CkMulticastMgr::ckNew();
  kpsi = CProxy_Psi::ckNew(true, config.K, config.L);
  qpsi = CProxy_Psi::ckNew(false, config.Q, config.M);

  fcalc = CProxy_FCalculator::ckNew();
  for (int k = 0; k < config.K; k++) {
    for (int q = 0; q < config.Q; q++) {
      for (int l = 0; l < config.L; l++) {
        for (int m = 0; m < config.M; m++) {
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

// Read in configuration data from the file system
void GWDriver::readConfig() {
  config.K = 4;
  config.Q = 4;
  config.L = 8;
  config.M = 16;

  config.occupied_size = 32;
  config.unoccupied_size = 16;

  config.pipeline_stages = 2;
}

// Read in state date from the file system
void GWDriver::readState() {
  // TODO: State input should be parallel
}

Psi::Psi(bool occupied) : occupied(occupied) {
  mcast_ptr = CProxy_CkMulticastMgr(mcast_ID).ckLocalBranch();

  if (occupied) {
    size = config.occupied_size;
  } else {
    size = config.unoccupied_size;
  }
  psi = new double[size];
  // TODO: Read in the Psi somehow or get it from the driver?
  calculatePsiR();
}

void Psi::setupSections() {
  CProxySection_FCalculator section;
  unsigned k_lower = 0, k_upper = config.K - 1;
  unsigned q_lower = 0, q_upper = config.Q - 1;
  unsigned l_lower = 0, l_upper = config.L - 1;
  unsigned m_lower = 0, m_upper = config.M - 1;
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
  // Reset section index to 0. It is used in SDAG control flow.
  section_index = 0;
}

void Psi::calculatePsiR() {
 // TODO (Yale): Take the psi array in this object and IFFT it to take it to
 // real space
} 

FCalculator::FCalculator() {
  mcast_ptr = CProxy_CkMulticastMgr(mcast_ID).ckLocalBranch();
}

// Setup the reduction section for our plane and broadcast a message to the
// section, allowing it to get the section info needed for the reduction.
// This method should only be called on one member of each l-plane.
void FCalculator::setupSections() {
  CkAssert(thisIndex.w == 0 && thisIndex.x == 0 && thisIndex.z == 0);

  PlaneMsg* msg = new PlaneMsg(thisIndex.y);
  CProxySection_FCalculator section;
  section = CProxySection_FCalculator::ckNew( thisProxy.ckGetArrayID(),
                                              0, config.K-1, 1,
                                              0, config.Q-1, 1,
                                              thisIndex.y, thisIndex.y, 1,
                                              0, config.M-1, 1 );
  section.ckSectionDelegate(mcast_ptr);
  section.receiveSectionInfo(msg);
}

void FCalculator::computeF(PsiMsg* m1, PsiMsg* m2) {
  // TODO (Yale): Take the two psis and compute f and store it in the field
}

void FCalculator::sendPContribution() {
  // TODO: This will probably be called from some more complicated SDAG flow
}

#include "gw_bse.def.h"
