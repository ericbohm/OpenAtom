#include "standard_include.h"
#include "ckcomplex.h"

#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"

#include "../class_defs/CP_OPERATIONS/class_cpintegrate.h"


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// CG minimization
//============================================================================

void CPINTEGRATE::CP_integrate_dyn(const int ncoef, const int istate,int iteration,
              complex *forces,complex *vpsi,complex *psi,
              const int *k_x, const int *k_y, const int *k_z,
              double *cmass, double *fictEke_ret)

//============================================================================
   { /* Begin Function */
//----------------------------------------------------------------------------
// Local Variables

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP *cp = CP::get();
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

   double dt   = gentimeinfo->dt;
   double dt2  = dt*0.5;

//============================================================================
// Update the velcoties to full step (if not the intial step).
// Compute the Fictious Kinetic Energy at time t (e.g. t=0 on initial step)

   if(iteration>1){
      for(int i = 0; i < ncoef; i++){
        vpsi[i] += forces[i]*(dt2/cmass[i]);     
      }//endfor
      // apply nhc
   }//endif

   double fictEke = 0.0;
   for(int i = 0; i < ncoef; i++){
     fictEke += vpsi[i].getMagSqr()*(0.5*cmass[i]);
   }//endfor
   (*fictEke_ret) = fictEke;

//============================================================================
// Update the velocities to half step (t+dt/2): 
// Update the positions to full step

   // apply nhc
   for(int i = 0; i < ncoef; i++){
     vpsi[i] += forces[i]*(dt2/cmass[i]);     
   }//endfor

   for(int i = 0; i < ncoef; i++){
     psi[i] += vpsi[i]*dt;
   }//endfor

//============================================================================
  } /* End function */
//============================================================================
