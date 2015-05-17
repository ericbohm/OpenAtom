#include "fcalculator.h"
#include "gw_bse.h"

FCalculator::FCalculator() {
  mcast_ptr = CProxy_CkMulticastMgr(mcast_ID).ckLocalBranch();
  size = config.rows_of_p;
}

// Setup the reduction section for our plane and broadcast a message to the
// section, allowing it to get the section info needed for the reduction.
// This method should only be called on one member of each l-plane.
void FCalculator::setupSections() const {
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

unsigned FCalculator::rowOffset(unsigned r) const {
  // TODO (UIUC): Make it so this staggers the rows by using our index
  return r;
}

void FCalculator::computeF(PsiMsg* m1, PsiMsg* m2) {
  f = new double[size];
  // TODO (Yale): Take the two psis and compute f and store it in the field
}

void FCalculator::computeAndSendRow(unsigned row_index) const {
  unsigned chare_index = row_index / config.rows_per_chare;
  RowMessage* msg = new (size) RowMessage(row_index, size);
  double* row = msg->row;
  // TODO (Yale): Compute the rth row of P using f
  pmatrix[chare_index].receiveRowContribution(msg);
}
