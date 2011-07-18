//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================


//============================================================================
#include "standard_include.h"
#include "../include/class_defs/allclass_mdintegrate.h"
#include "../include/class_defs/allclass_mdatoms.h"
#include "../include/class_defs/allclass_mdinter.h"
#include "../include/class_defs/allclass_mdintra.h"
#include "../include/class_defs/allclass_gen.h"
#include "../include/class_defs/allclass_cp.h"
#include "../include/class_defs/typedefs_par.h"
#include "../include/class_defs/PINY_INIT/PhysicsAtomPosInit.h"
#include "../include/class_defs/ATOM_OPERATIONS/class_vx_smpl.h"
#include "../proto_defs/proto_coords_entry.h"
#include "../proto_defs/proto_coords_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../../../include/Atoms.h"
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
PhysicsAtomPosInit::PhysicsAtomPosInit (int ibead_in , int itemper_in){
//============================================================================
// Get the readonly structures ready to go

  MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();
  MDATOMS      *mdatoms      = MDATOMS::get();
  MDINTER      *mdinter      = MDINTER::get();
  MDINTRA      *mdintra      = MDINTRA::get();
  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();

 MDTHERM_INFO *mdtherm_info = &(mdintegrate->mdtherm_info);
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_mdatoms.h"

//============================================================================
// Copy out some useful variables and error check

  pi_beads_true = (mdclatoms_info->pi_beads); //input controled by directory
  pi_beads      = 1;
  ntemper       = gensimopts->ntemper;
  natm_tot      = (mdclatoms_info->natm_tot);
  natm_nl       = (cppseudo->nonlocal.natm);
  iextended_on  = (mdtherm_info->iextended_on);
  kT            = (genstatepoint->t_ext)/BOLTZ;
  cp_min_opt    = (gensimopts->cp_wave_min+gensimopts->cp_min);
  istart_typ    = gensimopts->istart;
  num_nhc       = (mdtherm_info->num_nhc);
  len_nhc       = (mdtherm_info->len_nhc);
  cp_wave_opt   = (gensimopts->cp_wave);
  isokin_opt    = mdtherm_info->isokin_opt;
  ibead         = ibead_in;
  itemper       = itemper_in;

  if( (ibead>=pi_beads_true) || (ibead < 0) || (itemper >=ntemper) || (itemper < 0) ){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("Bad bead or temperer index to DriverAtomInit : %d %d : %d %d\n",
            ibead,pi_beads_true,itemper,ntemper);
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
  }//endif

  if(iextended_on==1){
    if(num_nhc != 3*natm_tot){
       PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       PRINTF("Presently only massive thermstatting is supported.\n");
       PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       EXIT(1);
    }//endif
  }//endif
  if(len_nhc>5){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("The maximum allowed Atm NHC len is 5.\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

//============================================================================
// Malloc the piny data structure : Lazy man way. Use piny input routines.

  mdclatoms_pos = (MDCLATOMS_POS *)
                cmalloc(pi_beads*sizeof(MDCLATOMS_POS),"PhysicsAtomPosInit.C")-1;

  for(int ip=1;ip<=pi_beads;ip++){
   mdclatoms_pos[ip].natm_tot = natm_tot;
   mdclatoms_pos[ip].x = (double *)
                       cmalloc(natm_tot*sizeof(double),"PhysicsAtomPosInit.C")-1;
   mdclatoms_pos[ip].y = (double *)
                       cmalloc(natm_tot*sizeof(double),"PhysicsAtomPosInit.C")-1;
   mdclatoms_pos[ip].z = (double *)
                       cmalloc(natm_tot*sizeof(double),"PhysicsAtomPosInit.C")-1;
   mdclatoms_pos[ip].vx =(double *)
                       cmalloc(natm_tot*sizeof(double),"PhysicsAtomPosInit.C")-1;
   mdclatoms_pos[ip].vy =(double *)
                       cmalloc(natm_tot*sizeof(double),"PhysicsAtomPosInit.C")-1;
   mdclatoms_pos[ip].vz =(double *)
                       cmalloc(natm_tot*sizeof(double),"PhysicsAtomPosInit.C")-1;
  }// endfor

  therm_bead    = (MDTHERM_POS *)
                   cmalloc(pi_beads*sizeof(MDTHERM_POS),"PhysicsAtomPosInit.C")-1;

  if(iextended_on==1 && num_nhc>0){
    therm_class.num_nhc = num_nhc;
    therm_class.len_nhc = len_nhc;
    therm_class.x_nhc   = cmall_mat(1,len_nhc,1,num_nhc,"PhysicsAtomPosInit.C");
    therm_class.v_nhc   = cmall_mat(1,len_nhc,1,num_nhc,"PhysicsAtomPosInit.C");
    mass_nhc            = cmall_mat(1,len_nhc,1,num_nhc,"PhysicsAtomPosInit.C");
    double **mass_nhc_in= mdtherm_info->mass_nhc;
    for(int j=1;j<=len_nhc;j++){
    for(int i=1;i<=num_nhc;i++){
      mass_nhc[j][i] = mass_nhc_in[j][i];
    }}
  }// endif

//============================================================================
//  Fill the Piny Data Structures : Resampl atom velocities if necessary

  read_coord(mdintegrate,mdatoms,mdinter,mdintra,general_data,cp,
             mdclatoms_pos,&therm_class,therm_bead,ibead,itemper);

  if(istart_typ<3){
    for(int ip=1;ip<=pi_beads;ip++){
      double *vx   = mdclatoms_pos[ip].vx;
      double *vy   = mdclatoms_pos[ip].vy;
      double *vz   = mdclatoms_pos[ip].vz;
      double *mass = mdclatoms_info->mass;
      VX_SMPL::ctrlSamplAtomVel(natm_tot,&vx[1],&vy[1],&vz[1],&mass[1]);
    }//endif
  }//endif

//----------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PhysicsAtomPosInit::DriverAtomInit(int natm_in,Atom *atoms,AtomNHC *atomsNHC,
                                        int ibead_in, int itemper_in){
//============================================================================
// Local pointers and Error checking

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP *cp                     = CP::get();
  MDATOMS      *mdatoms      = MDATOMS::get();
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"
#include "../include/class_defs/allclass_strip_mdatoms.h"
  int istart = gensimopts->istart;

  if(natm_tot != natm_in || ibead != ibead_in || itemper != itemper_in){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("Bad input to DriverAtomInit %d %d : %d %d : %d %d\n",
  	     natm_tot,natm_in,ibead,ibead_in,itemper,itemper_in);
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
  }//endif

//============================================================================
// Atom Initialization 

  double *x  = mdclatoms_pos[1].x;
  double *y  = mdclatoms_pos[1].y;
  double *z  = mdclatoms_pos[1].z;
  double *vx = mdclatoms_pos[1].vx;
  double *vy = mdclatoms_pos[1].vy;
  double *vz = mdclatoms_pos[1].vz;
  double *q  = mdclatoms_info->q;
  double *m  = mdclatoms_info->mass;

  for(int i=0,j=1;i<natm_tot;i++,j++){
    atoms[i].x  =  x[j];
    atoms[i].y  =  y[j];
    atoms[i].z  =  z[j];
    atoms[i].vx = vx[j];
    atoms[i].vy = vy[j];
    atoms[i].vz = vz[j];
    atoms[i].q  = q[j];
    atoms[i].m  = m[j];
    atoms[i].fx = 0.0;
    atoms[i].fy = 0.0;
    atoms[i].fz = 0.0;
  }//endfor

//============================================================================
// AtomNHC Initialization 

  for(int i=0;i<natm_tot;i++){
    atomsNHC[i].len_nhc = len_nhc;
    atomsNHC[i].kT      = kT;
    atomsNHC[i].posKT   = 0.0;
    for(int j=0;j<len_nhc;j++){
      atomsNHC[i].m[j] =1.0;
      atomsNHC[i].vx[j]=0.0;
      atomsNHC[i].vy[j]=0.0;
      atomsNHC[i].vz[j]=0.0;
      atomsNHC[i].fx[j]=0.0;
      atomsNHC[i].fy[j]=0.0;
      atomsNHC[i].fz[j]=0.0;
    }//endif
  }//endfor

  if(cp_min_opt==0 && iextended_on==1  && num_nhc>0){
    for(int i=0;i<natm_tot;i++){
      int iii = 3*i+1;
      for(int j=0,j1=1;j<len_nhc;j++,j1++){
        atomsNHC[i].m[j]  = mass_nhc[j1][iii];
      }//endfor
    }//endfor
    if(istart>=4){
      double **v_nhc = therm_class.v_nhc;
      for(int i=0;i<natm_tot;i++){
        int iii = 3*i+1;
        for(int j=0,j1=1;j<len_nhc;j++,j1++){
          atomsNHC[i].m[j]  = mass_nhc[j1][iii];
          atomsNHC[i].vx[j] = v_nhc[j1][iii];
          atomsNHC[i].vy[j] = v_nhc[j1][(iii+1)];
          atomsNHC[i].vz[j] = v_nhc[j1][(iii+2)];
        }//endfor
      }//endfor
    }else{
      VX_SMPL::ctrlSamplAtomNhcVel(natm_tot,atomsNHC);
    }//endif
  }//endif

 

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
PhysicsAtomPosInit::~PhysicsAtomPosInit (){
//============================================================================

  if(iextended_on==1  && num_nhc>0){
    cfree_mat(therm_class.x_nhc,1,len_nhc,1,num_nhc);
    cfree_mat(therm_class.v_nhc,1,len_nhc,1,num_nhc);
    cfree_mat(mass_nhc,1,len_nhc,1,num_nhc);
  }// endif

  cfree (&(therm_bead[1]),"~AtomPosInit()");

  for(int ip=1;ip<=pi_beads;ip++){
   cfree (&(mdclatoms_pos[ip].x[1]),"~AtomPosInit()");
   cfree (&(mdclatoms_pos[ip].y[1]),"~AtomPosInit()");
   cfree (&(mdclatoms_pos[ip].z[1]),"~AtomPosInit()");
   cfree (&(mdclatoms_pos[ip].vx[1]),"~AtomPosInit()");
   cfree (&(mdclatoms_pos[ip].vy[1]),"~AtomPosInit()");
   cfree (&(mdclatoms_pos[ip].vz[1]),"~AtomPosInit()");
  }// endfor

  cfree (&(mdclatoms_pos[1]),"~AtomPosInit()");

//---------------------------------------------------------------------------
  }//end routine
//============================================================================
