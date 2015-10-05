#include "gw_bse.h"
#include "psi.h"
#include "pmatrix.h"

/*readonly*/GWConfig config;
/*readonly*/CProxy_PsiCache psicache;
/*readonly*/CProxy_Psi psi;
/*readonly*/CProxy_PMatrix pmatrix;

GWDriver::GWDriver(CkArgMsg* msg) {
  readConfig();
  readState();

  psicache = CProxy_PsiCache::ckNew(config.L, config.n_elems);
  pmatrix = CProxy_PMatrix::ckNew(config.n_elems / config.rows_per_chare);
  psi = CProxy_Psi::ckNew(config.n_elems, config.K, config.L + config.M);
}

// Read in configuration data from the file system
void GWDriver::readConfig() {
  config.K = 1;
  config.L = 4;
  config.M = 16;

  config.n_elems = 256;
  config.rows_per_chare = 16;
  config.matrix_nchares= config.n_elems / config.rows_per_chare;

  config.pipeline_stages = 4;
}

// Read in state data from the file system
void GWDriver::readState() {
  // TODO: State input should be parallel
}

#include "gw_bse.def.h"
