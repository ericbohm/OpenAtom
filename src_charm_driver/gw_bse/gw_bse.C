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
  qpsi = CProxy_Psi::ckNew(false, config.Q, config.M);
  pmatrix = CProxy_PMatrix::ckNew(config.n_elems / config.rows_per_chare);

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

  config.n_elems = 32;
  config.rows_per_chare = 8;

  config.pipeline_stages = 2;
}

// Read in state date from the file system
void GWDriver::readState() {
  // TODO: State input should be parallel
}

#include "gw_bse.def.h"
