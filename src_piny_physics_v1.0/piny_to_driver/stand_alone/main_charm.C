//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                         PI_MD:                                           
//             The future of simulation technology                          
//             ------------------------------------                         
//                      Module: main                                        
//                                                                          
// This program performs MD on a classical potential energy surface (PES),  
// minimization on a classical PES,                                         
// MD on a mixed classical-density functional PES,                          
// minimization on a mixed classical-density functional PES,                
// PIMD on a classical PES,                                                 
// centroid minimization on a classical PES,                                
// PIMD on a mixed classical-density functional PES,                        
// centroid minimization on a mixed classical-density functional PES.       
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#include "standard_include.h"

#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_mdintegrate.h"
#include "../class_defs/allclass_mdatoms.h"
#include "../class_defs/allclass_mdinter.h"
#include "../class_defs/allclass_mdintra.h"

#include "../class_defs/ATOMS/class_mdclatoms_pos.h"
#include "../class_defs/INTEGRATE/class_mdtherm_pos.h"

#include "../proto_defs/proto_main_entry.h"
#include "../proto_defs/proto_main_local.h"
#include "../proto_defs/proto_parse_entry.h"
#include "../proto_defs/proto_coords_local.h"

#include "../class_defs/main_charm.h"

#include "../class_defs/AtomPosInit.h"

// class to test whether pupped variables made it to other nodes
#define TEST_MY_PUPS_NOT
#ifdef TEST_MY_PUPS
#include "../../main/test/test_parse.h"
#endif

// Defined as globals because we want 1 copy per node
// not one copy per object which would be silly given
// the classes, by construction, do not contain data.
// May have to copy out some variables to objects to 
// avoid memory conflicts but that can be handled on a
// case by case basis in the constructor of each object.
  MDINTEGRATE  readonly_mdintegrate;
  MDATOMS      readonly_mdatoms;
  MDINTER      readonly_mdinter;
  MDINTRA      readonly_mdintra;
  GENERAL_DATA readonly_general_data;
  CP           readonly_cp; 

  int nElements;


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// charm mainchare : class main : charm++ has its own main function. 
//                                You no longer control it.
//                                By default the mainchare is invoked
//                                on one processor only.
//==========================================================================

OATOM_MAINCHARE::OATOM_MAINCHARE (CkArgMsg *m)

//==========================================================================
    {// begin routine  
//==========================================================================
//   Local Variables 

  int pi_beads,ip,natm_tot;
  int iextended_on,num_nhc,len_nhc;

//-------------------------------------------------------------------------
// Nifty way to unglobally access the global variables.
// Easy to implement unglobalize versions later if necessary.

  MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();
  MDATOMS      *mdatoms      = MDATOMS::get();
  MDINTER      *mdinter      = MDINTER::get();
  MDINTRA      *mdintra      = MDINTRA::get();
  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();

//-------------------------------------------------------------------------
// local classes to read in the data for leanMD

  MDCLATOMS_POS *mdclatoms_pos;
  MDTHERM_POS    therm_class;
  MDTHERM_POS   *therm_bead;

//=======================================================================
//  I)             Check for input file 

  if(m->argc < 2) {
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("No input file specified\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }//endif

//=======================================================================
// II)            Invoke User Interface on 1 processor

  parse(mdintegrate,mdatoms,mdinter,mdintra,general_data,cp,m->argv[1]);

//=======================================================================
// III) output state of classes on 1 processor


//#define DEBUG_PARSE
//#ifdef DEBUG_PARSE
  mdintegrate->state_class_out(); 
  mdatoms->state_class_out(); 
  mdinter->state_class_out();
  mdintra->state_class_out(); 
  general_data->state_class_out();
  cp->state_class_out(); 
//#endif

//=======================================================================
// IV) Allocate classes, read the coordinates for leanMD
/*

  pi_beads      = (mdatoms->mdclatoms_info.pi_beads);
  natm_tot      = (mdatoms->mdclatoms_info.natm_tot);
  iextended_on  = (mdintegrate->mdtherm_info.iextended_on);

  mdclatoms_pos = (MDCLATOMS_POS *)
                   cmalloc(pi_beads*sizeof(MDCLATOMS_POS),"main.C")-1;

  for(ip=1;ip<=pi_beads;ip++){
   mdclatoms_pos[ip].natm_tot = natm_tot;
   mdclatoms_pos[ip].x = (double *)cmalloc(natm_tot*sizeof(double),"main.C")-1;
   mdclatoms_pos[ip].y = (double *)cmalloc(natm_tot*sizeof(double),"main.C")-1;
   mdclatoms_pos[ip].z = (double *)cmalloc(natm_tot*sizeof(double),"main.C")-1;
   mdclatoms_pos[ip].vx =(double *)cmalloc(natm_tot*sizeof(double),"main.C")-1;
   mdclatoms_pos[ip].vy =(double *)cmalloc(natm_tot*sizeof(double),"main.C")-1;
   mdclatoms_pos[ip].vz =(double *)cmalloc(natm_tot*sizeof(double),"main.C")-1;
  }// endfor

  therm_bead    = (MDTHERM_POS *)
                   cmalloc(pi_beads*sizeof(MDTHERM_POS),"main.C")-1;

  if(iextended_on==1){
    num_nhc = (mdintegrate->mdtherm_info.num_nhc);
    len_nhc = (mdintegrate->mdtherm_info.len_nhc);
    therm_class.num_nhc = num_nhc;
    therm_class.len_nhc = len_nhc;
    therm_class.x_nhc   = cmall_mat(1,len_nhc,1,num_nhc,"main.C");
    therm_class.v_nhc   = cmall_mat(1,len_nhc,1,num_nhc,"main.C");
  }// endif

  read_coord(mdintegrate,mdatoms,mdinter,mdintra,general_data,cp,
             mdclatoms_pos,&therm_class,therm_bead);
*/
//  AtomPosInit atomPosInit;
//=======================================================================
// V) Test pups  : invokes function that fires up a chare array
//              and looks at variables on different procs to check
//              the pupping. Define TEST_MY_PUPS above so function 
//              prototype is included along with function.

#ifdef TEST_MY_PUPS
   test_parse_invoke();
#endif

//=======================================================================
// VI) Control the simulations : Parallelism invoked within calls below.
//                               readonly classes automatically communicated
//                               to all NODES by pup magic.

#define CP_OFF
#ifdef CP_ON
  control_cp(mdintegrate,mdatoms,mdinter,mdintra,general_data,cp);
#endif

//=======================================================================
// VI)  Exit, stage left  : Actually you hang because well
//                          charm++ doesn't have a soft-exit.
//                          If you call CkExit, the pups don't pup.

  fflush(stdout);
  fflush(stderr);

  nElements=15;

  CkPrintf("Running Hello on %d processors for %d elements\n",
             CkNumPes(),nElements);
                                                                                
  CProxy_Hello arr = CProxy_Hello::ckNew(nElements);
                                                                                
  arr[0].SayHi(17);

//----------------------------------------------------------------------
   }//end main
//==========================================================================

class Hello : public CBase_Hello
{
public:
  Hello()
  {
    CkPrintf("Hello %d created\n",thisIndex);
  }
                                                                                
  Hello(CkMigrateMessage *m) {}
                                                                                
  void SayHi(int hiNo)
  {
    CkPrintf("Hi[%d] from element %d\n",hiNo,thisIndex);

/*
    readonly_mdintegrate.state_class_out();
    readonly_mdatoms.state_class_out();
    readonly_mdinter.state_class_out();
    readonly_mdintra.state_class_out();
    readonly_general_data.state_class_out();
    readonly_cp.state_class_out();
*/
    if (thisIndex < nElements-1)
      //Pass the hello on:
      thisProxy[thisIndex+1].SayHi(hiNo+1);
    else
      //We've been around once-- we're done.
      // mainProxy.done();
      CkExit ();
  }
};

#include "../charm_defs/main_charm.def.h"
