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

#include "../include/class_defs/Interface_ctrl.h"

// Defined as globals because we want 1 copy per node
// not one copy per object which would be silly given
// the classes, by construction, do not contain data.
// May have to copy out some variables to objects to 
// avoid memory conflicts but that can be handled on a
// case by case basis in the contructor of each object.
MDINTEGRATE  readonly_mdintegrate;
MDATOMS      readonly_mdatoms;
MDINTER      readonly_mdinter;
MDINTRA      readonly_mdintra;
GENERAL_DATA readonly_general_data;
CP           readonly_cp; 


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//
//==========================================================================

Interface_ctrl::Interface_ctrl (char* fname, CkCallback c)

  //==========================================================================
{// begin routine  
  //==========================================================================
  //   Local Variables 

  int iii;

  //-------------------------------------------------------------------------
  // Nifty way to unglobally access the global variables.
  // Easy to implement unglobalize versions later if necessary.

  MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();
  MDATOMS      *mdatoms      = MDATOMS::get();
  MDINTER      *mdinter      = MDINTER::get();
  MDINTRA      *mdintra      = MDINTRA::get();
  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();

  //=======================================================================
  // I)            Invoke User Interface on 1 processor

  PRINTF("\n"); 
  PRINT_LINE_STAR
    PRINTF("Invoking PINY Interface\n");
  PRINT_LINE_DASH
    PRINTF("\n"); 

  parse(mdintegrate,mdatoms,mdinter,mdintra,general_data,cp,fname);

  //=======================================================================
  // II) output state of classes on 1 processor

#ifdef DEBUG_PARSE
  mdintegrate->state_class_out(); 
  mdatoms->state_class_out(); 
  mdinter->state_class_out();
  mdintra->state_class_out(); 
  general_data->state_class_out();
  cp->state_class_out(); 
#endif

  //=======================================================================
  // IV)  Exit, stage left 

  // inform caller that initialization is done
  c.send ();
  fflush(stdout);
  fflush(stderr);

  PRINTF("\n"); 
  PRINT_LINE_DASH
    PRINTF("Completed Piny Parameter setup\n");
  PRINT_LINE_STAR 
    PRINTF("\n");

  //----------------------------------------------------------------------
}//end main
//==========================================================================

