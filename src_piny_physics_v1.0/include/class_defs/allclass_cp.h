//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                         PI_MD:                                           
//             The future of simulation technology                          
//             ------------------------------------                         
//                   Structures: typedefs_cp.h                              
//                                                                          
// The include file with the typedefs of all the cp structures              
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#include "../class_defs/CP/class_cpopts.h"
#include "../class_defs/CP/class_cpatom_maps.h"
#include "../class_defs/CP/class_cpcoeffs_info.h"
#include "../class_defs/CP/class_cptherm_info.h"
#include "../class_defs/CP/class_cpconstrnt.h"
#include "../class_defs/CP/class_cpewald.h"
#include "../class_defs/CP/class_cppseudo.h"
#include "../class_defs/CP/class_cpdual_pme.h"
#include "../class_defs/CP/class_cpatom_pme.h"
#include "../class_defs/CP/class_cpvel_samp.h"
#include "../class_defs/CP/class_cpylm_cons.h"
#include "../class_defs/CP/class_gen_wave.h"

//==========================================================================

#ifndef _CP_
#define _CP_

#ifdef CHARM_ON
// declare an instance of the class, called readonly_name, as extern or global
class CP; extern CP readonly_cp;
#endif

class CP {

 public:
  int num_pup;
  CPOPTS        cpopts;
  CPCOEFFS_INFO cpcoeffs_info;
  CPTHERM_INFO  cptherm_info;
  CPATOM_MAPS   cpatom_maps;
  CPCONSTRNT    cpconstrnt;
  CPEWALD       cpewald;
  CPPSEUDO      cppseudo; 
  CPDUAL_PME    cpdual_pme;
  CPATOM_PME    cpatom_pme;
  CPVEL_SAMP    cpvel_samp;
  CPYLM_CONS    cpylm_cons;
// GEN_WAVE      cpgen_wave; // not yet

//-------------------------------------------------------------------------
  CP(){num_pup=0;};
 ~CP(){};

#ifdef CHARM_ON
  static CP *get(){
    return &readonly_cp;  // return the pointer of the global instance
  }
#endif

//-------------------------------------------------------------------------
  void state_class_out(){
   PRINTF("\n");
   PRINT_LINE_STAR;
   PRINTF("CP state_class_out \n");
   cpopts.state_class_out();
   cpcoeffs_info.state_class_out();
   cptherm_info.state_class_out();
   cpatom_maps.state_class_out();
   cpconstrnt.state_class_out();
   cpewald.state_class_out();
   cppseudo.state_class_out();
   cpdual_pme.state_class_out();
   cpatom_pme.state_class_out();
   cpvel_samp.state_class_out();
   cpylm_cons.state_class_out();
   PRINT_LINE_STAR;
   PRINTF("\n");
  }// end routine

//-------------------------------------------------------------------------
#ifdef PUP_ON
  void pup(PUP::er &p){
    p | num_pup;
    cpopts.pup(p);
    cpcoeffs_info.pup(p);
    cptherm_info.pup(p);
    cpatom_maps.pup(p);
    cpconstrnt.pup(p);
    cpewald.pup(p);
    cppseudo.pup(p);
    cpdual_pme.pup(p);
    cpatom_pme.pup(p);
    cpvel_samp.pup(p);
    cpylm_cons.pup(p);
    num_pup++;
    if(num_pup==2){
     PUP_PRINTF("\n");
     PUP_PRINT_LINE_STAR;
     PUP_PRINTF("CP PUP COMPLETED\n");
     PUP_PRINT_LINE_STAR;
     PUP_PRINTF("\n");
    }//endif
  }
#endif

//-------------------------------------------------------------------------
  }; // end class CP
//==========================================================================

#ifdef PUP_ON
PUPmarshall(CP);
#endif

#endif

//==========================================================================



