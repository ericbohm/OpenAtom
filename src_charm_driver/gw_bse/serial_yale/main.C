#include "main.decl.h"

#define CHARM_ON
#define PUP_ON
#define _DEBUG_PUP_
#include "main.h"
#include "states.decl.h"
#include "standard_include_gwbse.h"
#include "allclass_gwbse.h"
#include "configure_gwbse.h"
#include "states.h"

// states_occ_proxy is declared readonly in states.ci file that makes them parallel global.
// This declaration in main where we instantiated makes them c++ global.
// module states is declared external in main.ci inside module main.
// module states is defined in states.ci file, no external.
/* readonly */ CProxy_Main mainProxy;
/* readonly */ GWBSE readonly_gwbse;
/* readonly */ CProxy_states_occ states_occ_proxy;
/* readonly */ CProxy_states_unocc states_unocc_proxy;


// Entry point of Charm++ application
Main::Main(CkArgMsg* msg) {

  GWBSE *gwbse = GWBSE::get();
  char input_name[10000];
  Config gw_configure;
    
    
  sprintf (input_name,"minjung_test_file");
  gw_configure.readConfig (input_name, gwbse);
    
  gwbse->state_class_out();

  // We are done with msg so delete it.
  delete msg;

  // Display some info about this execution for the user
  CkPrintf("Running GWBSE_V0 using %d processors.\n\n",CkNumPes());

  // Set the mainProxy readonly to point to a
  //   proxy for the Main chare object (this
  //   chare object).
  mainProxy = thisProxy;

  // Create the array of state chare objects.
  int nocc = gwbse->gwbseopts.nocc;  
  states_occ_proxy = CProxy_states_occ::ckNew(nocc);
  int nunocc = gwbse->gwbseopts.nunocc;  
  states_unocc_proxy = CProxy_states_unocc::ckNew(nunocc);


  // Invoke the "beamoutMyState()" entry method on the first
  //   element of the state array of chare objects.
#ifdef _ROUND_ROBIN_
  states_occ_proxy[0].beamoutMyState(1,1);
  states_unocc_proxy[0].beamoutMyState(1,1);
#else
  for(int i=0; i<nocc; i++){states_occ_proxy[i].beamoutMyState(i,1);}
  for(int i=0; i<nunocc; i++){states_unocc_proxy[i].beamoutMyState(i,1);}//parameter marshalling
#endif

  
}


// Constructor needed for chare object migration (ignore for now)
// NOTE: This constructor does not need to appear in the ".ci" file
Main::Main(CkMigrateMessage* msg) { }


// When called, the "done()" entry method will cause the program
//   to exit.
void Main::done() {
  CkExit();
}


#include "main.def.h"
