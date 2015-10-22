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
#define _DEBUG_PUP_
#include "main.h"
#include "states.h"
#include "pmatrix.h"
#include "controller.h"
#include "standard_include_gwbse.h"
#include "allclass_gwbse.h"
#include "configure_gwbse.h"
#include "states.h"

// =========================================================================
// states_occ_proxy is declared readonly in states.ci file that makes them parallel global.
// This declaration in main where we instantiated makes them c++ global.
// module states is declared external in main.ci inside module main.
// module states is defined in states.ci file, no external.
/* readonly */ CProxy_Main mainProxy;
/* readonly */ GWBSE readonly_gwbse;
/* readonly */ CProxy_States states_proxy;
/* readonly */ CProxy_PsiCache psi_cache_proxy;
/* readonly */ CProxy_PMatrix pmatrix_proxy;
/* readonly */ CProxy_Controller controller_proxy;


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
  sprintf (input_name,"minjung_test_file");
  delete msg;
  gw_configure.readConfig (input_name, gwbse);
    
  // -------------------------------------------------------------------
  // Dump the contents of the readonlies on PE 0
  gwbse->state_class_out();

  // -------------------------------------------------------------------
  // Display some info about this execution for the user
  CkPrintf("Running GWBSE_V0 using %d processors.\n\n",CkNumPes());

  // -------------------------------------------------------------------
  // Set the mainProxy readonly to point to a
  //   proxy for the Main chare object (this
  //   chare object).
  mainProxy = thisProxy;

  // -------------------------------------------------------------------
  // Create the controller chares
  controller_proxy = CProxy_Controller::ckNew();
  psi_cache_proxy = CProxy_PsiCache::ckNew();

  // -------------------------------------------------------------------
  // Create the array of P matrix chare objects.
  int nchares = gwbse->gw_parallel.matrix_nchares;
  pmatrix_proxy = CProxy_PMatrix::ckNew(nchares); 

  // -------------------------------------------------------------------
  // Create the array of state chare objects.
  int nspin = gwbse->gwbseopts.nspin;
  int nkpt = gwbse->gwbseopts.nkpt;
  int nocc = gwbse->gwbseopts.nocc;  
  int nunocc = gwbse->gwbseopts.nunocc;  
  states_proxy = CProxy_States::ckNew(nspin, nkpt, nocc + nunocc);

  // -------------------------------------------------------------------
  // Tell the controller to start the computation
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
