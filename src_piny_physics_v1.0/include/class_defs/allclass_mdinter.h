//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
//                   Classes allclass_mdinter.h
//
// The include file all of the classical atom intermolecular interactions
//
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#include "../class_defs/INTER/class_mdbrnch_root.h"
#include "../class_defs/INTER/class_mdenergy_ctrl.h"
#include "../class_defs/INTER/class_mdinteract.h"
#include "../class_defs/INTER/class_mdlnklist.h"
#include "../class_defs/INTER/class_mdnbr_list.h"
#include "../class_defs/INTER/class_mdpart_mesh.h"
#include "../class_defs/INTER/class_mdsurface.h"
#include "../class_defs/INTER/class_mdverlist.h"

//==========================================================================

#ifndef _MDINTER_
#define _MDINTER_

#ifdef CHARM_ON
// declare an instance of the class, called readonly_name, as extern or global
class MDINTER; extern MDINTER readonly_mdinter;
#endif


class MDINTER {

  public:
    int num_pup;

    MDBRNCH_ROOT   mdbrnch_root;
    MDENERGY_CTRL  mdenergy_ctrl;
    MDINTERACT     mdinteract;
    MDLNKLIST      mdlnklist;
    MDVERLIST      mdverlist;
    MDNBR_LIST     mdnbr_list;
    MDPART_MESH    mdpart_mesh;
    MDSURFACE      mdsurface;

    MDINTER(){num_pup=0;};
    ~MDINTER(){};

#ifdef CHARM_ON
    static MDINTER *get(){
      return &readonly_mdinter;  // return the pointer of the global instance
    }
#endif

    void state_class_out(){
      PRINTF("\n");
      PRINT_LINE_STAR;
      PRINTF("MDINTER state_class_out\n");
      mdbrnch_root.state_class_out();  
      mdenergy_ctrl.state_class_out(); 
      mdinteract.state_class_out();
      mdlnklist.state_class_out();
      mdverlist.state_class_out();
      mdnbr_list.state_class_out();
      mdpart_mesh.state_class_out();
      mdsurface.state_class_out();
      PRINT_LINE_STAR;
      PRINTF("\n");
    }//end routine


#ifdef PUP_ON
    void pup(PUP::er &p){
      p | num_pup;
      mdbrnch_root.pup(p);
      mdenergy_ctrl.pup(p);
      mdinteract.pup(p);
      mdlnklist.pup(p);
      mdverlist.pup(p);
      mdnbr_list.pup(p);
      mdpart_mesh.pup(p);
      mdsurface.pup(p);
      num_pup++;
      if(num_pup==2){
        PUP_PRINTF("\n");
        PUP_PRINT_LINE_STAR;
        PUP_PRINTF("INTER PUP COMPLETED\n");
        PUP_PRINT_LINE_STAR;
        PUP_PRINTF("\n");
      }//endif
    }// end pack/unpack
#endif


};/* end class definition */

#ifdef PUP_ON
PUPmarshall(MDINTER);
#endif

#endif
