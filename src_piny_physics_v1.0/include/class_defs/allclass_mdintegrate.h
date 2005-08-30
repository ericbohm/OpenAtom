//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
//                   Classes allclass_mdintegrate.h
//
// The include file all of the classical atom integration
//
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
   
#include "../class_defs/INTEGRATE/class_mdtherm_info.h"
#include "../class_defs/INTEGRATE/class_mdvel_samp.h"
#include "../class_defs/INTEGRATE/class_mdbaro.h"
#include "../class_defs/INTEGRATE/class_mdpar_rahman.h"

//==========================================================================

#ifndef _MDINTEGRATE_
#define _MDINTEGRATE_

#ifdef CHARM_ON
// declare an instance of the class, called readonly_name, as extern or global
class MDINTEGRATE; extern MDINTEGRATE readonly_mdintegrate;
#endif

class MDINTEGRATE {

 public:
  int num_pup;
  MDTHERM_INFO mdtherm_info;
  MDTHERM_INFO mdtherm_info_bead;
  MDVEL_SAMP   mdvel_samp;
  MDBARO       mdbaro;
  MDPAR_RAHMAN mdpar_rahman;

//------------------------------------------------------------------------

  MDINTEGRATE(){num_pup = 0;};
 ~MDINTEGRATE(){};

#ifdef CHARM_ON
  static MDINTEGRATE *get(){
    return &readonly_mdintegrate;  // return the pointer of the global instance
  }
#endif

//------------------------------------------------------------------------

  void state_class_out(){
    PRINTF("\n");
    PRINT_LINE_STAR;
    mdtherm_info.state_class_out(); 
    mdtherm_info_bead.state_class_out(); 
    mdvel_samp.state_class_out(); 
    mdbaro.state_class_out(); 
    mdpar_rahman.state_class_out();
    PRINT_LINE_STAR;
    PRINTF("\n");
  }// end routine

//------------------------------------------------------------------------
#ifdef PUP_ON
  void pup(PUP::er &p){
    p | num_pup;
    mdtherm_info.pup(p);
    mdtherm_info_bead.pup(p);
    mdvel_samp.pup(p);
    mdbaro.pup(p);
    mdpar_rahman.pup(p);
    num_pup++;
    if(num_pup==2){
     PUP_PRINTF("\n");
     PUP_PRINT_LINE_STAR;
     PUP_PRINTF("INTEGRATE PUP COMPLETED\n");
     PUP_PRINT_LINE_STAR;
     PUP_PRINTF("\n");
   }//endif
  }// end pack/unpack
#endif
//------------------------------------------------------------------------

};/* end class definition */

#ifdef PUP_ON
PUPmarshall(MDINTEGRATE);
#endif

#endif
//==========================================================================
