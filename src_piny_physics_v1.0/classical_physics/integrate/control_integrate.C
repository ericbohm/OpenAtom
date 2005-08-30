#include "standard_include.h"

#include "../../../include/debug_flags.h"
#include "../../../include/Atoms.h"
#include "../class_defs/ATOM_OPERATIONS/class_atomintegrate.h"

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// 
//============================================================================
  
void ATOMINTEGRATE::ATOM_integrate(int natm, Atom *atoms, int myid)

//============================================================================
   {//begin routine 
//============================================================================

#ifdef GJM_DBG_ATMS
   PRINTF("GJM_DBG : Inside integrate %d\n",myid);  
#endif


#ifdef _CP_DEBUG_ATM_FORC_

  if(myid==0){

      FILE *fp = fopen("atom_forc.out","w");
      for(int i=0;i<natm;i++){
        fprintf(fp,"%d %g %g %g\n",i+1,atoms[i].fx,atoms[i].fy,atoms[i].fz);
      }//endfor
      fclose(fp);
//      CkExit();

  }//endif

#endif

//---------------------------------------------------------------------------
   }//end routine
//============================================================================
