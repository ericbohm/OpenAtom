#include "fcalculator.h"
#include "gw_bse.h"

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
