#include "gw_bse.h"
#include "psi.h"
#include "fcalculator.h"
#include "ckmulticast.h"

/*readonly*/GWConfig config;
/*readonly*/CkGroupID mcast_ID;
/*readonly*/CProxy_Psi kpsi;
/*readonly*/CProxy_Psi qpsi;
/*readonly*/CProxy_FCalculator fcalc;
/*readonly*/CProxy_PMatrix pmatrix;

GWDriver::GWDriver(CkArgMsg* msg) {
  readConfig();
  readState();

  mcast_ID = CProxy_CkMulticastMgr::ckNew();
  kpsi = CProxy_Psi::ckNew(true, config.K, config.L);
  qpsi = CProxy_Psi::ckNew(false, config.K, config.M);
  pmatrix = CProxy_PMatrix::ckNew(config.n_elems / config.rows_per_chare);

  fcalc = CProxy_FCalculator::ckNew(config.K, config.L, config.M);

  kpsi.sendPsi();
  qpsi.sendPsi();
  fcalc.run();
}

// Read in configuration data from the file system
void GWDriver::readConfig() {
  config.K = 4;
  config.L = 1;
  config.M = 4;

  config.n_elems = 256;
  config.rows_per_chare = 16;

  config.pipeline_stages = 1;
}

// Read in state date from the file system
void GWDriver::readState() {
  // TODO: State input should be parallel
}

#include "gw_bse.def.h"
