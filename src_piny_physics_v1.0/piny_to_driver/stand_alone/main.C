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


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

int main (int argc, char *argv[])

  //==========================================================================
{// begin routine  
  //==========================================================================
  //   Local Variables 

  int pi_beads,ip,natm_tot;
  int iextended_on,num_nhc,len_nhc;

  MDINTEGRATE    mdintegrate;
  MDATOMS        mdatoms;
  MDINTER        mdinter;
  MDINTRA        mdintra;
  GENERAL_DATA   general_data;
  CP             cp; 

  // classes to read in position data for leanMD
  MDCLATOMS_POS *mdclatoms_pos;
  MDTHERM_POS    therm_class;
  MDTHERM_POS   *therm_bead;


  //=======================================================================
  //  I)             Check for input file                                  

  if(argc < 2) {
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("No input file specified\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }//endif

  //=======================================================================
  // II)            Invoke User Interface                                  

  parse(&mdintegrate,&mdatoms,&mdinter,&mdintra,&general_data,&cp,argv[1]);

  //=======================================================================
  // III) output state of classes 

#define DEBUG_PARSE
#ifdef DEBUG_PARSE
  mdintegrate.state_class_out(); 
  mdatoms.state_class_out(); 
  mdinter.state_class_out();
  mdintra.state_class_out(); 
  general_data.state_class_out();
  cp.state_class_out(); 
#endif

  //=======================================================================
  // IV) Allocate classes, read the coordinates 

  pi_beads      = (mdatoms.mdclatoms_info.pi_beads);
  natm_tot      = (mdatoms.mdclatoms_info.natm_tot);
  iextended_on  = (mdintegrate.mdtherm_info.iextended_on);

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
    num_nhc = (mdintegrate.mdtherm_info.num_nhc);
    len_nhc = (mdintegrate.mdtherm_info.len_nhc);
    therm_class.num_nhc = num_nhc;
    therm_class.len_nhc = len_nhc;
    therm_class.x_nhc   = cmall_mat(1,len_nhc,1,num_nhc,"main.C");
    therm_class.v_nhc   = cmall_mat(1,len_nhc,1,num_nhc,"main.C");
  }// endif

  read_coord(&mdintegrate,&mdatoms,&mdinter,&mdintra,&general_data,&cp,
      mdclatoms_pos,&therm_class,therm_bead);

  //=======================================================================
  // V)  Exit, stage left 

  fflush(stdout);
  fflush(stderr);
  exit(0); 
  return 0;

  //----------------------------------------------------------------------
}//end routine
//==========================================================================




