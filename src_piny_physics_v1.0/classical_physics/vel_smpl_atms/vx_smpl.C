#include "../../include/class_defs/vx_smpl.h"

//=================================================================== 
//=================================================================== 

vx_smpl::vx_smpl(double*        velx,
                 double*        vely,
                 double*        velz,    
                 double*        mass,    
                 double*        text_atm, 
                 int            natm,      
                 GENENSOPTS*    genensopts,  
                 MDVEL_SAMP*    vel_samp,    
                 MDGHOST_ATOMS* ghost_atoms,
                 MDCONSTRNT*    mdconstrnt)
//=================================================================== 
{//begin routine 
//=================================================================== 

 //-------------------  Local pointers ----------------------

  int iconstrnt      = mdconstrnt->iconstrnt;

  int nghost_tot     = ghost_atoms->nghost_tot;
  int *ighost_map    = ghost_atoms->ighost_map;

  int nfreeze        = mdconstrnt->nfreeze;
  int *freeze_map    = mdconstrnt->freeze_map;

  int npt_f          = genensopts->npt_f;
  int npt_i          = genensopts->npt_i;
  int nve            = genensopts->nve;
  int nvt            = genensopts->nvt;

//=================================================================== 
// 0) Write to screen                                                 

    PRINT_LINE_STAR;
    PRINTF("Sampling atomic velocities\n");
    PRINT_LINE_DASH;printf("\n");

//=================================================================== 
// I) Sample atom velocities                                           

   vx_smpl::sampl_vx(velx,vely,velz,mass,text_atm,natm,
         	           &(vel_samp->iseed),
                     &(vel_samp->iseed2),
                     &(vel_samp->qseed));

//==================================================================== 
//==================================================================== 

// II) Project onto surface of constraint                              

  if(iconstrnt==1){
     int iproj_vel = 1;

     if(iproj_vel == 1) {
     if(npt_f==1 || npt_i == 1){
       PRINTF("Correction to velocities due to projection onto surface of constraint has \n");
       PRINTF(" not implemented for npt_i or npt_f \n");
       exit(1);
     }

     if(nve  ==1 || nvt==1){
       PRINTF("Correction to velocities due to projection onto surface of constraint has \n");
       PRINTF(" not implemented for nve or nvt \n");
       //       proj_vel( );
     }

    }// endif iproj_vel  
  }//endif iconstrnt 
 
//==================================================================== 
// III) Zero velocities of ghost atoms if any                          

  if(nghost_tot > 0) {
    int ighost;
   for(ighost=1;ighost <= nghost_tot;ighost++){
    int igloc = ighost_map[ighost];
     velx[igloc] = 0.0;
     vely[igloc] = 0.0;
     velz[igloc] = 0.0;
   }//endfor 
  }//endif 

//==================================================================== 
// III) Zero velocities of freeze atoms if any                          

  if(nfreeze > 0) {
    int i;
   for(i=1;i <= nfreeze;i++){
     int igloc = freeze_map[i];
     velx[igloc] = 0.0;
     vely[igloc] = 0.0;
     velz[igloc] = 0.0;
   }//endfor 
  }//endif 

//==================================================================== 
// III) Write to screen                                                

    PRINT_LINE_DASH;
    PRINTF("Atomic velocity sampling complete\n");
    PRINT_LINE_STAR;printf("\n");

//==================================================================== 
   }//end routine 
//==================================================================== 


//==================================================================== 
// Particle Velocities  
//==================================================================== 

void vx_smpl::sampl_vx(double* velx,              // x-comp of atom velocities
                       double* vely,              // y-comp of atom velocities
                       double* velz,              // z-comp of atom velocities
                       double* mass,
                       double* text_atm,
                       int     natm_tot,
                       int*    iseed,
                       int*    iseed2,
                       double* qseed)
// =================================================================== 
{//begin routine 
// ===================================================================

   int iatm;

#define COMPARING_TO_OLD_PINY
#ifdef  COMPARING_TO_OLD_PINY

   double *temp = (new double[natm_tot]) - 1;

   gaussran(natm_tot,iseed,iseed2,qseed,temp);
   int i;
   for(i=1; i <= natm_tot; i++){ velx[i] = temp[i]; }
   gaussran(natm_tot,iseed,iseed2,qseed,temp);
   for(i=1; i <= natm_tot; i++){ vely[i] = temp[i]; }
   gaussran(natm_tot,iseed,iseed2,qseed,temp);
   for(i=1; i <= natm_tot; i++){ velz[i] = temp[i]; }

   delete [] (temp+1);

#else
   for(iatm=1;iatm<= natm_tot;iatm++){
     gaussran(1,iseed,iseed2,qseed,&(velx[iatm]));
     gaussran(1,iseed,iseed2,qseed,&(vely[iatm]));
     gaussran(1,iseed,iseed2,qseed,&(velz[iatm]));
   }//endfor 
#endif

   for(iatm=1; iatm <= natm_tot; iatm++){

    double start_temp = text_atm[iatm];
    double width = sqrt(start_temp/(mass[iatm]*BOLTZ));
    velx[iatm] *= width;
    vely[iatm] *= width;
    velz[iatm] *= width;

   }//endfor 

//=================================================================== 
 } //end routine 
//=================================================================== 


