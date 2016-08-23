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
#include "standard_include.h"
#include "allclass_gwbse.h"
#include "configure_gwbse.h"
#include "states.h"
#include "CkLoopAPI.h"

// =========================================================================
// states_occ_proxy is declared readonly in states.ci file that makes them parallel global.
// This declaration in main where we instantiated makes them c++ global.
// module states is declared external in main.ci inside module main.
// module states is defined in states.ci file, no external.
/* readonly */ CProxy_Main mainProxy;
/* readonly */ CProxy_MatMul mat_mul_proxy;
/* readonly */ GWBSE readonly_gwbse;
/* readonly */ CProxy_States states_proxy;
/* readonly */ CProxy_PsiCache psi_cache_proxy;

/* readonly */ CProxy_EpsMatrix2D eps_matrix2D_proxy;

/* readonly */ CProxy_PMatrix2D pmatrix2D_proxy;
/* readonly */ CProxy_EpsMatrix2D pmatrix2D_bproxy;
/* readonly */ CProxy_EpsMatrix2D pmatrix2D_cproxy;

/* readonly */ CProxy_EpsMatrix1D eps_proxy1D;

/* readonly */ CProxy_EpsMatrix2D pmatrix2D_aaproxy;
/* readonly */ CProxy_EpsMatrix2D pmatrix2D_bbproxy;
/* readonly */ CProxy_EpsMatrix2D pmatrix2D_ccproxy;

/* readonly */ CProxy_PMatrix1D pmatrix1D_proxy;
///* readonly */ CProxy_Controller controller_proxy;
/* readonly */ CProxy_FFTController fft_controller_proxy;


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Entry point of Charm++ application
Main::Main(CkArgMsg* msg) {

  // -------------------------------------------------------------------
  // Get the readonly pointer and declare the Config class
  GWBSE *gwbse = GWBSE::get();
  char input_name[10000];
  Config gw_configure;
    
  // -------------------------------------------------------------------
  // Read in all of the configurable parameters (take name from CkArgMsg)
  sprintf(input_name, msg->argv[1]);
  //sprintf (input_name,"minjung_test_file");
  delete msg;
  gw_configure.readConfig (input_name, gwbse);

//  gw_configure.read_usrinput();    
  // -------------------------------------------------------------------
  // Dump the contents of the readonlies on PE 0
  gwbse->state_class_out();

  // -------------------------------------------------------------------
  // Display some info about this execution for the user
  CkPrintf("Running GWBSE_V0 using %d processors.\n\n",CkNumPes());

  // -------------------------------------------------------------------
  // Initialize CkLoop, which is used by the PsiCache
#ifdef USE_CKLOOP
  CkLoop_Init(1);
#endif

  // -------------------------------------------------------------------
  // Set the mainProxy readonly to point to a
  //   proxy for the Main chare object (this
  //   chare object).
  mainProxy = thisProxy;

  // -------------------------------------------------------------------
  // Create the controller chares


  mat_mul_proxy = CProxy_MatMul::ckNew();


  controller_proxy = CProxy_Controller::ckNew();
  fft_controller_proxy = CProxy_FFTController::ckNew();
  psi_cache_proxy = CProxy_PsiCache::ckNew();
   

  // -------------------------------------------------------------------
  // Create the array of P matrix chare objects.
  int dimension = gwbse->gw_parallel.n_elems;
  int nrows = gwbse->gw_parallel.rows_per_chare;
  int ncols = gwbse->gw_parallel.cols_per_chare;
  CkPrintf("[MAIN] Creating %ix%i matrix chares with %i rows and %i cols each\n",
      dimension/nrows, dimension/ncols, nrows, ncols);


  eps_matrix2D_proxy = CProxy_EpsMatrix2D::ckNew();
  eps_proxy1D = CProxy_EpsMatrix1D::ckNew(140);

  pmatrix2D_proxy = CProxy_PMatrix2D::ckNew(dimension/nrows, dimension/ncols);
  pmatrix2D_bproxy = CProxy_PMatrix2D::ckNew();
  pmatrix2D_cproxy = CProxy_PMatrix2D::ckNew();
 
  pmatrix2D_ccproxy = CProxy_PMatrix2D::ckNew();
  pmatrix2D_aaproxy = CProxy_PMatrix2D::ckNew();
  pmatrix2D_bbproxy = CProxy_PMatrix2D::ckNew();

  //eps_matrix2D_proxy = CProxy_EpsMatrix2D::ckNew();

  int local_mtx_size_1d_x = dimension;
  int local_mtx_size_1d_y = 1;//1728; 
  int number_of_chares_1d = dimension / local_mtx_size_1d_y;
  CkPrintf("[MAIN] Creating %i array chares with %i rows and %i cols each\n",
      number_of_chares_1d, local_mtx_size_1d_y, local_mtx_size_1d_x);

  pmatrix1D_proxy = CProxy_PMatrix1D::ckNew(local_mtx_size_1d_x, local_mtx_size_1d_y, number_of_chares_1d);



  int nspin = gwbse->gwbseopts.nspin;
  int nkpt = gwbse->gwbseopts.nkpt;
  int nocc = gwbse->gwbseopts.nocc;
  int nunocc = gwbse->gwbseopts.nunocc;
  states_proxy = CProxy_States::ckNew(nspin, nkpt, nocc + nunocc);
  CkPrintf("[MAIN] Creating %i occupied states and %i unoccupied states\n",
      nspin*nkpt*nocc, nspin*nkpt*nunocc);
  // -------------------------------------------------------------------
  // Create the array of state chare objects.
  
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
