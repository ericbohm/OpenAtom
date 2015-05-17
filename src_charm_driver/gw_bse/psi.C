#include "psi.h"
#include "gw_bse.h"
#include "fcalculator.h"

Psi::Psi(bool occupied) : occupied(occupied), section_index(0) {
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
  // K,Q indices are not independent of one another in the unoccupied case
  std::vector<std::pair<unsigned, unsigned> > k_q_indices;
  // L,M indices are always independent of one another
  std::vector<unsigned> l_indices, m_indices;

  if (occupied) {
    // Occupied has the Psi for a single k and all q
    for (int q = 0; q < config.Q; q++) {
      k_q_indices.push_back(std::make_pair(thisIndex.x, q));
    }

    // Occupied has the Psi for a single l and all m
    l_indices.push_back(thisIndex.y);
    for (int m = 0; m < config.M; m++) m_indices.push_back(m);
  } else {
    // Unoccupied has the Psi for a single (k+q)
    unsigned k_plus_q = thisIndex.x;
    for (int k = 0; k < config.K; k++) {
      for (int q = 0; q < config.Q; q++) {
        if ((k + q - 1) % config.K == k_plus_q) {
          k_q_indices.push_back(std::make_pair(k,q));
        }
      }
    }

    // Unoccupied has the Psi for a single m and all l
    for (int l = 0; l < config.L; l++) l_indices.push_back(l);
    m_indices.push_back(thisIndex.y);
  }

  // Build the sections from the indices constructed above
  for (unsigned l : l_indices) {
    CProxySection_FCalculator section;
    std::vector<CkArrayIndex4D> indices;
    for (std::pair<unsigned, unsigned> k_q : k_q_indices) {
      for (unsigned m : m_indices) {
        indices.push_back(CkArrayIndex4D(k_q.first, k_q.second,l,m));
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
