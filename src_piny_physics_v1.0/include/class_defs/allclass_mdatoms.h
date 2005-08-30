//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
//                   Classes allclass_mdatoms.h
//
// The include file all of the classical atom relevant classes
//
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
                                                                                
#include "../class_defs/ATOMS/class_mdatom_maps.h"
#include "../class_defs/ATOMS/class_mdclatoms_info.h"
#include "../class_defs/ATOMS/class_mdclatoms_pimd.h"

//==========================================================================

#ifndef _MDATOMS_
#define _MDATOMS_

#ifdef CHARM_ON
// declare an instance of the class, called readonly_name, as extern or global
class MDATOMS; extern MDATOMS readonly_mdatoms;  
#endif

class MDATOMS {

 public:
  int num_pup;
  MDATOM_MAPS mdatom_maps;
  MDCLATOMS_INFO mdclatoms_info;
  MDCLATOMS_PIMD mdclatoms_pimd;

//--------------------------------------------------------------------------
  MDATOMS(){num_pup=0;};
 ~MDATOMS(){};

#ifdef CHARM_ON
  static MDATOMS *get(){
    return &readonly_mdatoms;  // return the pointer of the global instance
  }
#endif

//--------------------------------------------------------------------------
  void state_class_out(){
   PRINTF("\n");
   PRINT_LINE_STAR;
   PRINTF("MDATOMS state_class_out\n");
   mdatom_maps.state_class_out();
   mdclatoms_info.state_class_out();
   mdclatoms_pimd.state_class_out();
   PRINT_LINE_STAR;
   PRINTF("\n");
  }// end routine

//--------------------------------------------------------------------------
#ifdef PUP_ON
  void pup(PUP::er &p){
    p | num_pup;
    mdatom_maps.pup(p);
    mdclatoms_info.pup(p);
    mdclatoms_pimd.pup(p);
    num_pup++;
    if(num_pup==2){
     PUP_PRINTF("\n");
     PUP_PRINT_LINE_STAR;
     PUP_PRINTF("ATOMS PUP COMPLETED\n");
     PUP_PRINT_LINE_STAR;
     PUP_PRINTF("\n");
   }//endif
  }/* end pack/unpack */
#endif

 };/* end class definition */
//==========================================================================

#ifdef PUP_ON
PUPmarshall(MDATOMS);
#endif
                                                                             
#endif
//==========================================================================
