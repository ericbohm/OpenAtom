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
    int cp_min_opt, int cp_wave_opt, int cp_bomd_opt, int tol_reached, int iextended_on,
    Atom *atoms,AtomNHC *atomsNHC,int myid,
    double *eKinetic,double *eKineticNhc,double *potNhc,
    int *iwrite_atm,int output_on,
    int natmNow,int natmStr,int natmEnd,int mybead, bool switchMoveNow, 
    double new_t_ext, double old_t_ext)
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
  // (0) Check whether to move
  int move_atoms = 0;
  if (cp_min_opt == 0 && cp_wave_opt == 0) { move_atoms = 1; }
  if (cp_bomd_opt == 1 && tol_reached == 1) { move_atoms = 1; }
  if (switchMoveNow && move_atoms==0) { 
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("move_atoms must not be zero during a temper switch\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
  }
  //============================================================================
  // (I) Evolve to the last 1/2 step of previous step : only changes velocities

#ifdef _NAN_CHECK_
  for(int i=0; i<natm ; i++)
    {
     CkAssert(finite(atoms[i].x));
     CkAssert(finite(atoms[i].y));
     CkAssert(finite(atoms[i].z));
    }
#endif
  if(move_atoms){
    integrate_2nd_half_step(itime,natm,len_nhc,iextended_on,atoms,atomsNHC,
        eKinetic,eKineticNhc,potNhc,natmNow,natmStr,natmEnd);
  }//endif

#ifdef _NAN_CHECK_
  for(int i=0; i<natm ; i++)
    {
     CkAssert(finite(atoms[i].x));
     CkAssert(finite(atoms[i].y));
     CkAssert(finite(atoms[i].z));
    }
#endif

  //============================================================================
  // (II) Save the end step velocities and positions : pos already saved for PIMD
  
  if(move_atoms){
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

#ifdef _NAN_CHECK_
  for(int i=0; i<natm ; i++)
    {
     CkAssert(finite(atoms[i].x));
     CkAssert(finite(atoms[i].y));
     CkAssert(finite(atoms[i].z));
    }
#endif
  //============================================================================
  // (IIb) implement parallel tempering
  //  1. rescale particle velocities 
  //  2. rescale extended system velocities 
  //  3. modify all extended system parameters that depend on kT
  //     to reflect the new temperature 
  double fact=sqrt(new_t_ext/old_t_ext);
  if(switchMoveNow)
    {
      for(int i=natmStr;i<natmEnd;i++){
        atoms[i].vx *= fact;
        atoms[i].vy *= fact; 
        atoms[i].vz *= fact;
      }
      if(iextended_on)
	{
	  
	}
    }

  //============================================================================
  // (III) Evolve to first 1/2 step of the present step
#ifdef _NAN_CHECK_
  for(int i=0; i<natm ; i++)
    {
     CkAssert(finite(atoms[i].x));
     CkAssert(finite(atoms[i].y));
     CkAssert(finite(atoms[i].z));
    }
#endif


  if(move_atoms){
    integrate_1st_half_step(natm,len_nhc,iextended_on,atoms,atomsNHC,
			    natmNow,natmStr,natmEnd, switchMoveNow, fact); 
  }//endif
#ifdef _NAN_CHECK_
  for(int i=0; i<natm ; i++)
    {
     CkAssert(finite(atoms[i].x));
     CkAssert(finite(atoms[i].y));
     CkAssert(finite(atoms[i].z));
    }
#endif

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
    int natmNow,int natmStr,int natmEnd, bool switchMoveNow, double fact)
  //============================================================================
{//begin routine 
  //============================================================================

  MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();
  MDTHERM_INFO *mdtherm_info = &(mdintegrate->mdtherm_info);
  int isokin_opt = mdtherm_info->isokin_opt;
#ifdef _NAN_CHECK_
  CkPrintf("iextended_on is %d isokin_opt %d\n",iextended_on, isokin_opt);
  for(int i=0; i<natm ; i++)
    {
     CkAssert(finite(atoms[i].x));
     CkAssert(finite(atoms[i].y));
     CkAssert(finite(atoms[i].z));
     CkAssert(finite(atoms[i].vx));
     CkAssert(finite(atoms[i].vy));
     CkAssert(finite(atoms[i].vz));
    }
#endif


  switch(iextended_on){
    case 0 : integrate_nve_1st_half(natm,atoms,natmNow,natmStr,natmEnd); break;
    case 1 : if(isokin_opt==0){
               integrate_nvt_1st_half(natm,len_nhc,atoms,atomsNHC,
				      natmNow,natmStr,natmEnd, switchMoveNow, fact); 
             }else{
               integrate_isonvt_1st_half(natm,len_nhc,atoms,atomsNHC,
					 natmNow,natmStr,natmEnd, switchMoveNow, fact);
             }//endif
             break;
  }//endif
#ifdef _NAN_CHECK_
  CkPrintf("iextended_on is %d isokin_opt %d\n",iextended_on, isokin_opt);
  for(int i=0; i<natm ; i++)
    {
     CkAssert(finite(atoms[i].x));
     CkAssert(finite(atoms[i].y));
     CkAssert(finite(atoms[i].z));
     CkAssert(finite(atoms[i].vx));
     CkAssert(finite(atoms[i].vy));
     CkAssert(finite(atoms[i].vz));
    }
#endif

  //---------------------------------------------------------------------------
}//end routine
//============================================================================

