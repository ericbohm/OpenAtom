//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

#include "standard_include.h"
#include "../../../include/debug_flags.h"
#include "../../../include/Atoms.h"
#include "../class_defs/ATOM_OPERATIONS/class_atomintegrate.h"
#include "../class_defs/ATOM_OPERATIONS/class_atomoutput.h"
#include "../class_defs/allclass_mdintegrate.h"

//============================================================================



//============================================================================
//  Atom Integration controller : also invokes output
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::ctrl_atom_integrate(int itime,int natm,int len_nhc,
                        int cp_min_opt, int cp_wave_opt, int iextended_on,
                        Atom *atoms,AtomNHC *atomsNHC,int myid,
                        double *eKinetic,double *eKineticNhc,double *potNhc,
                        int *iwrite_atm,int output_on)
//============================================================================
   {//begin routine 
//============================================================================
// Local Variables and Pointers : Verbose output

#ifdef GJM_DBG_ATMS
   PRINTF("GJM_DBG : Inside integrate %d\n",myid);  
#endif
   int pi_beads = 1;

//============================================================================
// (I) Evolve to the last 1/2 step of previous step
 
   if(cp_min_opt==0 && cp_wave_opt==0){
     integrate_2nd_half_step(itime,natm,len_nhc,iextended_on,atoms,atomsNHC,
                             eKinetic,eKineticNhc,potNhc); 
   }//endif
 
//============================================================================
// (II) Invoke atom output to dump and config files

   ATOMOUTPUT::ctrl_piny_output(itime,natm,len_nhc,pi_beads,myid,atoms,atomsNHC,
                                iwrite_atm,output_on);

//============================================================================
// (III) Evolve to first 1/2 step of the present step

   if(cp_min_opt==0 && cp_wave_opt==0){
     integrate_1st_half_step(natm,len_nhc,iextended_on,atoms,atomsNHC); 
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
                            double *eKinetic,double *eKineticNhc,double *potNhc)
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
     case 0 : integrate_nve_2nd_half(itime,natm,atoms,eKinetic);
              break;
     case 1 : if(isokin_opt==0){
                integrate_nvt_2nd_half(itime,natm,len_nhc,atoms,atomsNHC,
                                        eKinetic,eKineticNhc,potNhc); 
              }else{
                integrate_isonvt_2nd_half(itime,natm,len_nhc,atoms,atomsNHC,
                                        eKinetic,eKineticNhc,potNhc); 
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
                                            Atom *atoms,AtomNHC *atomsNHC)
//============================================================================
   {//begin routine 
//============================================================================

   MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();
   MDTHERM_INFO *mdtherm_info = &(mdintegrate->mdtherm_info);
   int isokin_opt = mdtherm_info->isokin_opt;

   switch(iextended_on){
     case 0 : integrate_nve_1st_half(natm,atoms); break;
     case 1 : if(isokin_opt==0){
                integrate_nvt_1st_half(natm,len_nhc,atoms,atomsNHC); 
              }else{
                integrate_isonvt_1st_half(natm,len_nhc,atoms,atomsNHC); 
              }//endif
              break;
   }//endif

//---------------------------------------------------------------------------
   }//end routine
//============================================================================
