//=================================================================== 
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=================================================================== 
#include "../../include/class_defs/ATOM_OPERATIONS/class_vx_smpl.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_mdintegrate.h"
#include "../class_defs/allclass_mdintra.h"
#include "../class_defs/allclass_mdatoms.h"
//=================================================================== 


//=================================================================== 
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=================================================================== 
// Some list of velocites
//=================================================================== 
void VX_SMPL::ctrlSamplAtomVel(int natm,double *vx,double *vy,
    double *vz,double *mass)
  //=================================================================== 
{//begin routine 
  //=================================================================== 

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  MDATOMS      *mdatoms      = MDATOMS::get();
  MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();
  MDINTRA      *mdintra      = MDINTRA::get();

#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_mdintegrate.h"
#include "../class_defs/allclass_strip_mdintra.h"

  int iconstrnt      = mdconstrnt->iconstrnt;
  int nghost_tot     = mdghost_atoms->nghost_tot;
  int *ighost_map    = mdghost_atoms->ighost_map;
  int nfreeze        = mdconstrnt->nfreeze;
  int *freeze_map    = mdconstrnt->freeze_map;
  int npt_f          = genensopts->npt_f;
  int npt_i          = genensopts->npt_i;
  int nve            = genensopts->nve;
  int nvt            = genensopts->nvt;
  double *t_ext      = mdclatoms_info->text_atm;

  //=================================================================== 
  // 0) Write to screen                                                 

  PRINT_LINE_STAR;
  PRINTF("Sampling atomic velocities\n");
  PRINT_LINE_DASH;printf("\n");

  //=================================================================== 
  // I) Sample atom velocities : velocities passed to start at vx[0]

  sampl3DVelMultiT(natm,vx,vy,vz,mass,&t_ext[1],
      &(mdvel_samp->iseed),&(mdvel_samp->iseed2),
      &(mdvel_samp->qseed));

  //==================================================================== 
  // II) Project onto surface of constraint                              

  if(iconstrnt==1){
    PRINTF("Correction to velocities due to projection\n");
    PRINTF("on surface on constraint not implemented in vx_smpl\n");
  }//endif

  //==================================================================== 
  // III) Zero velocities of ghost atoms if any                          

  if(nghost_tot > 0) {
    for(int ighost=1;ighost <= nghost_tot;ighost++){
      int igloc = ighost_map[ighost]-1;
      vx[igloc] = 0.0;
      vy[igloc] = 0.0;
      vz[igloc] = 0.0;
    }//endfor 
  }//endif 

  //==================================================================== 
  // III) Zero velocities of freeze atoms if any                          

  if(nfreeze > 0) {
    for(int i=1;i <= nfreeze;i++){
      int igloc = freeze_map[i]-1;
      vx[igloc] = 0.0;
      vy[igloc] = 0.0;
      vz[igloc] = 0.0;
    }//endfor 
  }//endif 

  //==================================================================== 
  // III) Write to screen                                                

  PRINT_LINE_DASH;
  PRINTF("Atomic velocity sampling complete\n");
  PRINT_LINE_STAR;printf("\n");

  //--------------------------------------------------------------------
}//end routine 
//==================================================================== 


//=================================================================== 
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=================================================================== 
// Atom Velocities
//=================================================================== 
void VX_SMPL::ctrlSamplAtomVel(int natm, Atom *atoms){
  //=================================================================== 

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  MDATOMS      *mdatoms      = MDATOMS::get();
  MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();
  MDINTRA      *mdintra      = MDINTRA::get();

#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_mdintegrate.h"
#include "../class_defs/allclass_strip_mdintra.h"

  int iconstrnt   = mdconstrnt->iconstrnt;
  int nghost_tot  = mdghost_atoms->nghost_tot;
  int *ighost_map = mdghost_atoms->ighost_map;
  int nfreeze     = mdconstrnt->nfreeze;
  int *freeze_map = mdconstrnt->freeze_map;
  int npt_f       = genensopts->npt_f;
  int npt_i       = genensopts->npt_i;
  int nve         = genensopts->nve;
  int nvt         = genensopts->nvt;
  double *t_ext   = mdclatoms_info->text_atm;

  //=================================================================== 
  // 0) Write to screen                                                 

  PRINT_LINE_STAR;
  PRINTF("Sampling atomic velocities\n");
  PRINT_LINE_DASH;printf("\n");

  //=================================================================== 
  // I) Sample atom velocities : atoms start at 0 : atoms[0].vx 

  sampl3DVelMultiT(natm,atoms,&t_ext[1],
      &(mdvel_samp->iseed),&(mdvel_samp->iseed2),
      &(mdvel_samp->qseed));

  //==================================================================== 
  // II) Project onto surface of constraint                              

  if(iconstrnt==1){
    PRINTF("Correction to velocities due to projection\n");
    PRINTF("on surface on constraint not implemented in vx_smpl\n");
  }//endif

  //==================================================================== 
  // III) Zero velocities of ghost atoms if any                          

  if(nghost_tot > 0) {
    for(int ighost=1;ighost <= nghost_tot;ighost++){
      int igloc = ighost_map[ighost]-1;
      atoms[igloc].vx = 0.0;
      atoms[igloc].vy = 0.0;
      atoms[igloc].vz = 0.0;
    }//endfor 
  }//endif 

  //==================================================================== 
  // III) Zero velocities of freeze atoms if any                          

  if(nfreeze > 0) {
    for(int i=1;i <= nfreeze;i++){
      int igloc = freeze_map[i]-1;
      atoms[igloc].vx = 0.0;
      atoms[igloc].vy = 0.0;
      atoms[igloc].vz = 0.0;
    }//endfor 
  }//endif 

  //==================================================================== 
  // III) Write to screen                                                

  PRINT_LINE_DASH;
  PRINTF("Atomic velocity sampling complete\n");
  PRINT_LINE_STAR;printf("\n");

  //--------------------------------------------------------------------
}//end routine 
//==================================================================== 


//=================================================================== 
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================
// NHC Velocities
//=================================================================== 
void VX_SMPL::ctrlSamplAtomNhcVel(int natm,AtomNHC *atomsNHC){

  MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();
#include "../class_defs/allclass_strip_mdintegrate.h"

  for(int i=0;i<natm;i++){
    double kT    = atomsNHC[i].kT;
    int len_nhc  = atomsNHC[i].len_nhc;
    double *vx   = atomsNHC[i].vx;
    double *vy   = atomsNHC[i].vy;
    double *vz   = atomsNHC[i].vz;
    double *mass = atomsNHC[i].m;
    sampl3DVelOneT(len_nhc,vx,vy,vz,mass,kT,&(mdvel_samp->iseed),
        &(mdvel_samp->iseed2),&(mdvel_samp->qseed));
  }//endfor

}//end routine
//=================================================================== 


//==================================================================== 
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==================================================================== 
// Velocities  
//==================================================================== 
void VX_SMPL::sampl3DVelMultiT(int natm, double* vx, double* vy,
    double* vz, double* mass, double* t_ext,
    long* iseed, long* iseed2, double* qseed)
  //=================================================================== 
{//begin routine 
  //===================================================================
  // Sample unit gaussian

  double *temp = new double[natm];
  double *ptemp = temp-1;

  gaussran(natm,iseed,iseed2,qseed,ptemp);
  for(int i=0;i<natm;i++){vx[i] = temp[i];}
  gaussran(natm,iseed,iseed2,qseed,ptemp);
  for(int i=0;i<natm;i++){vy[i] = temp[i];}
  gaussran(natm,iseed,iseed2,qseed,ptemp);
  for(int i=0;i<natm;i++){vz[i] = temp[i];}

  delete [] temp;

  //===================================================================
  // Apply the width

  for(int i=0;i<natm;i++){
    double start_temp = t_ext[i];
    double width = sqrt(start_temp/(mass[i]*BOLTZ));
    vx[i] *= width;
    vy[i] *= width;
    vz[i] *= width;
  }//endfor 

  //------------------------------------------------------------------
} //end routine 
//=================================================================== 


//==================================================================== 
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==================================================================== 
// Velocities of Atoms
//==================================================================== 
void VX_SMPL::sampl3DVelMultiT(int natm, Atom *atoms, double* t_ext,
    long* iseed, long* iseed2, double* qseed)
  //=================================================================== 
{//begin routine 
  //===================================================================
  // Sample unit gaussian

  double *temp = new double[natm];
  double *ptemp = temp-1;

  gaussran(natm,iseed,iseed2,qseed,ptemp);
  for(int i=0;i<natm;i++){atoms[i].vx = temp[i];}
  gaussran(natm,iseed,iseed2,qseed,ptemp);
  for(int i=0;i<natm;i++){atoms[i].vy = temp[i];}
  gaussran(natm,iseed,iseed2,qseed,ptemp);
  for(int i=0;i<natm;i++){atoms[i].vz = temp[i];}

  delete [] temp;

  //===================================================================
  // Apply the width 

  for(int i=0;i<natm;i++){
    double start_temp = t_ext[i];
    double width = sqrt(start_temp/(atoms[i].m*BOLTZ));
    atoms[i].vx *= width;
    atoms[i].vy *= width;
    atoms[i].vz *= width;
  }//endfor 

  //------------------------------------------------------------------
} //end routine 
//=================================================================== 


//==================================================================== 
// Velocities  
//==================================================================== 
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==================================================================== 
void VX_SMPL::sampl3DVelOneT(int natm, double* vx, double* vy,
    double* vz,double* mass,double  kT,
    long* iseed,long* iseed2, double* qseed)
  //=================================================================== 
{//begin routine 
  //===================================================================
  // Sample unit gaussian

  double *temp = new double[natm];
  double *ptemp = temp-1;

  gaussran(natm,iseed,iseed2,qseed,ptemp);
  for(int i=0;i<natm;i++){vx[i] = temp[i]; }
  gaussran(natm,iseed,iseed2,qseed,ptemp);
  for(int i=0;i<natm;i++){vy[i] = temp[i]; }
  gaussran(natm,iseed,iseed2,qseed,ptemp);
  for(int i=0;i<natm;i++){vz[i] = temp[i]; }

  delete [] temp;

  //===================================================================
  // Apply the width

  for(int i=0;i<natm;i++){
    double width = sqrt(kT/mass[i]); //kT has boltz
    vx[i] *= width;
    vy[i] *= width;
    vz[i] *= width;
  }//endfor 

  //------------------------------------------------------------------
} //end routine 
//=================================================================== 



