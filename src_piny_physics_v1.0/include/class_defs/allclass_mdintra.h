//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
//                   Classes allclass_mdintra.h
//
// The include file all of the classical bonded (intramolecular) 
//     interaction classes
//
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#include "../class_defs/INTRA/class_mdbend_bnd.h"
#include "../class_defs/INTRA/class_mdbend_free.h"
#include "../class_defs/INTRA/class_mdbend.h"
#include "../class_defs/INTRA/class_mdbond_free.h"
#include "../class_defs/INTRA/class_mdbond.h"
#include "../class_defs/INTRA/class_mdconstrnt.h"
#include "../class_defs/INTRA/class_mdecor.h"
#include "../class_defs/INTRA/class_mdexcl.h"
#include "../class_defs/INTRA/class_mdgrp_bond_con.h"
#include "../class_defs/INTRA/class_mdgrp_bond_watts.h"
#include "../class_defs/INTRA/class_mdonfo.h"
#include "../class_defs/INTRA/class_mdrbar_sig_free.h"
#include "../class_defs/INTRA/class_mdtors_free.h"
#include "../class_defs/INTRA/class_mdtors.h"
#include "../class_defs/INTRA/class_mdghost_atoms.h"

//==========================================================================

#ifndef _MDINTRA_
#define _MDINTRA_

#ifdef CHARM_ON
// declare an instance of the class, called readonly_name, as extern or global
class MDINTRA; extern MDINTRA readonly_mdintra;
#endif

class MDINTRA {

  public:
    int num_pup;
    MDBOND mdbond;
    MDBEND mdbend;
    MDTORS mdtors;
    MDONFO mdonfo;
    MDECOR mdecor;
    MDEXCL mdexcl;
    MDBOND_FREE mdbond_free;
    MDBEND_FREE mdbend_free;
    MDTORS_FREE mdtors_free;
    MDBEND_BND mdbend_bnd;
    MDGRP_BOND_CON mdgrp_bond_con;
    MDGRP_BOND_WATTS mdgrp_bond_watts;
    MDCONSTRNT mdconstrnt;
    MDRBAR_SIG_FREE mdrbar_sig_free;
    MDGHOST_ATOMS mdghost_atoms;

    MDINTRA(){num_pup=0;};
    ~MDINTRA(){};

#ifdef CHARM_ON
    static MDINTRA *get(){
      return &readonly_mdintra;  // return the pointer of the global instance
    }
#endif

    void state_class_out(){
      PRINTF("\n");
      PRINT_LINE_STAR;
      PRINTF("MDINTRA state_class_out\n");
      mdbond.state_class_out(); 
      mdbend.state_class_out(); 
      mdtors.state_class_out(); 
      mdonfo.state_class_out(); 
      mdecor.state_class_out(); 
      mdexcl.state_class_out(); 
      mdbond_free.state_class_out();
      mdbend_free.state_class_out();
      mdtors_free.state_class_out();
      mdbend_bnd.state_class_out(); 
      mdgrp_bond_con.state_class_out(); 
      mdgrp_bond_watts.state_class_out();
      mdconstrnt.state_class_out(); 
      mdrbar_sig_free.state_class_out();
      mdghost_atoms.state_class_out();
      PRINT_LINE_STAR;
      PRINTF("\n");
    }//end routine

#ifdef PUP_ON
    void pup(PUP::er &p){
      p | num_pup;
      mdbond.pup(p);
      mdbend.pup(p);
      mdtors.pup(p);
      mdonfo.pup(p);
      mdecor.pup(p);
      mdexcl.pup(p);
      mdbond_free.pup(p);
      mdbend_free.pup(p);
      mdtors_free.pup(p);
      mdbend_bnd.pup(p);
      mdgrp_bond_con.pup(p);
      mdgrp_bond_watts.pup(p);
      mdconstrnt.pup(p);
      mdrbar_sig_free.pup(p);
      mdghost_atoms.pup(p);
      num_pup++;
      if(num_pup==2){
        PUP_PRINTF("\n");
        PUP_PRINT_LINE_STAR;
        PUP_PRINTF("INTRA PUP COMPLETED\n");
        PUP_PRINT_LINE_STAR;
        PUP_PRINTF("\n");
      }//endif
    }// end pack/unpack
#endif


};/* end class definition */

#ifdef PUP_ON
PUPmarshall(MDINTRA);
#endif


#endif
