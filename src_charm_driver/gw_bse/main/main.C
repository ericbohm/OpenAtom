//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// The main module for GWBSE (Up to phase 3)
// There are 2 state 3D chare arrays (spin, kpt, state), occ and unocc
// There is a single NodeCache nodegroup for storing occupied psis on each node
// There is a polarization matrix in realspace, decomposed by rows
// There is a polarization matrix in GSpace, decompsed by rows - requires a transpose
// There is a dielectric matrix in GSpace, decomposed best for iterative inverse
// There is the inverse of the dielectric matrix, created iteratively
// End of Phase 3.
//==========================================================================
#include "main.decl.h"

#define CHARM_ON
#define PUP_ON
#include "standard_include.h"
#include "main.h"
#include "states.h"
#include "pmatrix.h"
#include "eps_matrix.h"
#include "controller.h"
#include "fft_controller.h"
#include "allclass_gwbse.h"
#include "configure_gwbse.h"
#include "states.h"
#include "CkLoopAPI.h"

// =========================================================================
// states_occ_proxy is declared readonly in states.ci file that makes them parallel global.
// This declaration in main where we instantiated makes them c++ global.
// module states is declared external in main.ci inside module main.
// module states is defined in states.ci file, no external.
/* readonly */ GWBSE readonly_gwbse;

// Controller proxies
/* readonly */ CProxy_Main mainProxy;
/* readonly */ CProxy_Controller controller_proxy;
/* readonly */ CProxy_FFTController fft_controller_proxy;
/* readonly */ CProxy_PsiCache psi_cache_proxy;
/* readonly */ CProxy_FVectorCache fvector_cache_proxy;

// States (Phase 1)
/* readonly */ CProxy_States states_proxy;

// PMatrix Proxies (Phase 1 & 2)
/* readonly */ CProxy_PMatrix pmatrix2D_proxy;
/* readonly */ CProxy_PMatrix pmatrix1D_proxy;

// Epsilon matrices (Phase 3 & 4)
/* readonly */ CProxy_EpsMatrix2D eps_matrix2D_proxy;
/* readonly */ CProxy_EpsMatrix1D eps_proxy1D;
/* readonly */ CProxy_EpsMatrix2D eps_matrix2D_bproxy;
/* readonly */ CProxy_EpsMatrix2D eps_matrix2D_cproxy;
/* readonly */ CProxy_EpsMatrix2D eps_matrix2D_bbproxy;
/* readonly */ CProxy_EpsMatrix2D eps_matrix2D_ccproxy;
/* readonly */ CProxy_MatMul mat_mul_proxy;

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Entry point of Charm++ application
Main::Main(CkArgMsg* msg) {
  // -------------------------------------------------------------------
  // Initialize CkLoop, which is used by the PsiCache
#ifdef USE_CKLOOP
  CkLoop_Init(1);
#endif

  // -------------------------------------------------------------------
  // Get the readonly pointer and declare the Config class
  GWBSE *gwbse = GWBSE::get();
  char input_name[10000];
  Config gw_configure;
    
  // -------------------------------------------------------------------
  // Read in all of the configurable parameters (take name from CkArgMsg)
  sprintf(input_name, msg->argv[1]);
  delete msg;
  gw_configure.readConfig(input_name, gwbse);

  // -------------------------------------------------------------------
  // Dump the contents of the readonlies on PE 0
  gwbse->state_class_out();

  // -------------------------------------------------------------------
  // Set the mainProxy readonly to point to a
  //   proxy for the Main chare object (this
  //   chare object).
  mainProxy = thisProxy;

  // -------------------------------------------------------------------
  // Create the controller chares
  controller_proxy = CProxy_Controller::ckNew();
  fft_controller_proxy = CProxy_FFTController::ckNew();
  psi_cache_proxy = CProxy_PsiCache::ckNew();
  fvector_cache_proxy = CProxy_FVectorCache::ckNew();

  // -------------------------------------------------------------------
  // Create the array of state chare objects.
  int nspin = gwbse->gwbseopts.nspin;
  int nkpt = gwbse->gwbseopts.nkpt;
  int nocc = gwbse->gwbseopts.nocc;
  int nunocc = gwbse->gwbseopts.nunocc;
  states_proxy = CProxy_States::ckNew(nspin, nkpt, nocc + nunocc);
  CkPrintf("[MAIN] Creating %i occupied states and %i unoccupied states\n",
      nspin*nkpt*nocc, nspin*nkpt*nunocc);
  
  // -------------------------------------------------------------------
  // Create the arrays of P matrix chare objects.
  int dimension = gwbse->gw_parallel.n_elems;
  MatrixConfig pmatrix2D_config, pmatrix1D_config;
  pmatrix2D_config.mat_rows = pmatrix2D_config.mat_cols = dimension;
  pmatrix2D_config.tile_rows = gwbse->gw_parallel.rows_per_chare;
  pmatrix2D_config.tile_cols = gwbse->gw_parallel.cols_per_chare;
  pmatrix1D_config = convertTo1D(pmatrix2D_config, 1);

  CkPrintf("[MAIN] PMatrix size: %ix%i\n", dimension, dimension);
  CkPrintf("[MAIN] 2D decomposition: %ix%i chares with %ix%i tiles\n",
      pmatrix2D_config.chareRows(), pmatrix2D_config.chareCols(),
      pmatrix2D_config.tile_rows, pmatrix2D_config.tile_cols);
  CkPrintf("[MAIN] 1D decomposition: %ix%i chares with %ix%i tiles\n",
      pmatrix1D_config.chareRows(), pmatrix1D_config.chareCols(),
      pmatrix1D_config.tile_rows, pmatrix1D_config.tile_cols);

  pmatrix2D_proxy = CProxy_PMatrix::ckNew(pmatrix2D_config,
      pmatrix2D_config.chareRows(), pmatrix2D_config.chareCols());
  pmatrix1D_proxy = CProxy_PMatrix::ckNew(pmatrix1D_config,
      pmatrix1D_config.chareRows(), pmatrix1D_config.chareCols());

  // -------------------------------------------------------------------
  // Create the arrays of epsilon matrix chare objects.
  CProxy_EpsMap map = CProxy_EpsMap::ckNew();
  CkArrayOptions opts;
  opts.setMap(map);

  eps_proxy1D = CProxy_EpsMatrix1D::ckNew();
  eps_matrix2D_proxy = CProxy_EpsMatrix2D::ckNew(opts);
  eps_matrix2D_bproxy = CProxy_EpsMatrix2D::ckNew();
  eps_matrix2D_cproxy = CProxy_EpsMatrix2D::ckNew();
  eps_matrix2D_ccproxy = CProxy_EpsMatrix2D::ckNew();
  eps_matrix2D_bbproxy = CProxy_EpsMatrix2D::ckNew();
  mat_mul_proxy = CProxy_MatMul::ckNew();

  // -------------------------------------------------------------------
  // Start the computation
  controller_proxy.run();
}// end routine
// End of the Main constructor
// =========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Constructor needed for chare object migration (ignore for now)
// NOTE: This constructor does not need to appear in the ".ci" file
Main::Main(CkMigrateMessage* msg) { }

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// When called, the "done()" entry method will cause the program
//   to exit.
void Main::done() {
  CkExit();
}


#include "main.def.h"
