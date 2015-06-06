#include "psi.h"
#include "gw_bse.h"
#include "fcalculator.h"

Psi::Psi(bool occupied) : occupied(occupied), section_index(0) {
  mcast_ptr = CProxy_CkMulticastMgr(mcast_ID).ckLocalBranch();

  size = config.n_elems;
  psi = new double[size];
  // TODO: Read in the Psi somehow or get it from the driver?
  // NOTE: For now, using randomly generated psi for performance impact
  for (int i = 0; i < size; i++) {
    psi[i] = ((double) std::rand() / RAND_MAX);
  }

  calculatePsiR();
}

void Psi::setupSections(unsigned q) {
  std::vector<unsigned> k_indices, l_indices, m_indices;

  if (occupied) {
    k_indices.push_back(thisIndex.x);
    l_indices.push_back(thisIndex.y);
    for (int m = 0; m < config.M; m++) m_indices.push_back(m);
  } else {
    int k = (thisIndex.x + q) % config.K;
    k_indices.push_back(k);
    for (int l = 0; l < config.L; l++) l_indices.push_back(l);
    m_indices.push_back(thisIndex.y);
  }

  // Build the sections from the indices constructed above
  for (unsigned l : l_indices) {
    CProxySection_FCalculator section;
    std::vector<CkArrayIndex3D> indices;
    for (unsigned k : k_indices) {
      for (unsigned m : m_indices) {
        indices.push_back(CkArrayIndex3D(k,l,m));
      }
    }
    section = CProxySection_FCalculator::ckNew( fcalc.ckGetArrayID(),
                                                indices.data(), indices.size());
    section.ckSectionDelegate(mcast_ptr);
    sections.push_back(std::make_pair(l, section));
  }
}

void Psi::calculatePsiR() {
  // TODO (Yale): Take the psi array in this object and IFFT it to take it to
  // real space
} 
