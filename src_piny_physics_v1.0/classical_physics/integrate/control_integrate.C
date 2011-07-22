//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

#include "standard_include.h"
#include "../../../include/debug_flags.h"
#include "../../../include/Atoms.h"
#include "../class_defs/ATOM_OPERATIONS/class_atomintegrate.h"
#include "../class_defs/ATOM_OPERATIONS/class_atomoutput.h"
#include "../class_defs/allclass_mdintegrate.h"
#include "../class_defs/allclass_gen.h"

//============================================================================



//============================================================================
//  Atom Integration controller :
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::ctrl_atom_integrate(int itime,int natm,int len_nhc,
                        int cp_min_opt, int cp_wave_opt, int iextended_on,
                        Atom *atoms,AtomNHC *atomsNHC,int myid,
                        double *eKinetic,double *eKineticNhc,double *potNhc,
                        int *iwrite_atm,int output_on,
			int natmNow,int natmStr,int natmEnd,int mybead)
//============================================================================
   {//begin routine 
//============================================================================
// Local Variables and Pointers : Verbose output

   GENERAL_DATA *general_data = GENERAL_DATA::get();
#include "../class_defs/allclass_strip_gen.h"
    int pi_beads        = gensimopts->pi_beads;

#ifdef GJM_DBG_ATMS
   PRINTF("GJM_DBG : Inside integrate %d\n",myid);  
#endif
//============================================================================
// (I) Evolve to the last 1/2 step of previous step : only changes velocities
 
   if(cp_min_opt==0 && cp_wave_opt==0){
     integrate_2nd_half_step(itime,natm,len_nhc,iextended_on,atoms,atomsNHC,
                             eKinetic,eKineticNhc,potNhc,natmNow,natmStr,natmEnd);
   }//endif
 
//============================================================================
// (II) Save the end step velocities and positions

   if(cp_min_opt==0 && cp_wave_opt==0){
     if(pi_beads==1){
      for(int i=natmStr;i<natmEnd;i++){
        atoms[i].xold = atoms[i].x;  atoms[i].vxold = atoms[i].vx;
        atoms[i].yold = atoms[i].y;  atoms[i].vyold = atoms[i].vy;
        atoms[i].zold = atoms[i].z;  atoms[i].vzold = atoms[i].vz;
      }//endif
     }else{
      for(int i=natmStr;i<natmEnd;i++){
        atoms[i].vxold = atoms[i].vx;
        atoms[i].vyold = atoms[i].vy;
        atoms[i].vzold = atoms[i].vz;
      }//endif
     }//endfor
   }//endfor

//============================================================================
// (III) Evolve to first 1/2 step of the present step

   if(cp_min_opt==0 && cp_wave_opt==0){
     integrate_1st_half_step(natm,len_nhc,iextended_on,atoms,atomsNHC,
                             natmNow,natmStr,natmEnd); 
   }//endif

//============================================================================
// (III) Debug output

#ifdef _CP_DEBUG_ATM_FORC_
  if(myid==0){
      FILE *fp = fopen("atom_forc.out","w");
      for(int i=0;i<natm;i++){
        fprintf(fp,"%d %g %g %g\n",i+1,atoms[i].fx,atoms[i].fy,atoms[i].fz);
      }//endfor
      fclose(fp);
  }//endif
#endif

//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//  Atom Integration : 2nd half time step
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::integrate_2nd_half_step(int itime,int natm,int len_nhc,
                            int iextended_on,Atom *atoms,AtomNHC *atomsNHC,
                            double *eKinetic,double *eKineticNhc,double *potNhc,
                            int natmNow,int natmStr,int natmEnd)
//============================================================================
   {//begin routine 
//============================================================================

   MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();
   MDTHERM_INFO *mdtherm_info = &(mdintegrate->mdtherm_info);
   int isokin_opt = mdtherm_info->isokin_opt;
 
   (*eKinetic)    = 0.0;
   (*eKineticNhc) = 0.0;
   (*potNhc)      = 0.0;
   switch(iextended_on){
     case 0 : integrate_nve_2nd_half(itime,natm,atoms,
                                     eKinetic,natmNow,natmStr,natmEnd);
              break;
     case 1 : if(isokin_opt==0){
                integrate_nvt_2nd_half(itime,natm,len_nhc,atoms,atomsNHC,
                                       eKinetic,eKineticNhc,potNhc,natmNow,natmStr,natmEnd); 
              }else{
                integrate_isonvt_2nd_half(itime,natm,len_nhc,atoms,atomsNHC,
                                       eKinetic,eKineticNhc,potNhc,natmNow,natmStr,natmEnd); 
	      }//endif
              break;
   }//endif

//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//  Atom Integration : 2nd half time step
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::integrate_1st_half_step(int natm,int len_nhc,int iextended_on,
                                            Atom *atoms,AtomNHC *atomsNHC,
                                            int natmNow,int natmStr,int natmEnd)
//============================================================================
   {//begin routine 
//============================================================================

   MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();
   MDTHERM_INFO *mdtherm_info = &(mdintegrate->mdtherm_info);
   int isokin_opt = mdtherm_info->isokin_opt;

   switch(iextended_on){
     case 0 : integrate_nve_1st_half(natm,atoms,natmNow,natmStr,natmEnd); break;
     case 1 : if(isokin_opt==0){
               integrate_nvt_1st_half(natm,len_nhc,atoms,atomsNHC,
                                      natmNow,natmStr,natmEnd); 
              }else{
 	       integrate_isonvt_1st_half(natm,len_nhc,atoms,atomsNHC,
                                         natmNow,natmStr,natmEnd);
              }//endif
              break;
   }//endif

//---------------------------------------------------------------------------
   }//end routine
//============================================================================

