#include "psi.h"
#include "gw_bse.h"
#include "fcalculator.h"

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
