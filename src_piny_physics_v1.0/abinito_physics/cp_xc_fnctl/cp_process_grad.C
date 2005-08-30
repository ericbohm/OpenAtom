#include "standard_include.h"
#include "ckcomplex.h"

#include "../class_defs/CP_OPERATIONS/class_cpxcfnctls.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"


//============================================================================
/* The following function takes the g-space density and multiplies by
*  the vector ig, in order to create a g-space representation of the
*  gradient of the density
*  
*         rho(r) = sum_g rho(g) exp(ig.r)
*         grad rho(r) = sum_g ig rho(g) exp(ig.r)
*   
*  The g-space density is currently held in gradientCorrection and the
*  result of the multiplication is placed in rhoIGX, rhoIGY, rhoIGZ
*/
//============================================================================

void
CPXCFNCTS::CP_div_rho_gspace_calc(
             complex* gradientCorrection,
             const int *k_x, const int *k_y, const int *k_z, 
             int size, complex* rhoIGX, complex* rhoIGY, complex* rhoIGZ) 

//============================================================================
{/* Begin function */
//============================================================================
// Local variables

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

    double gx,gy,gz;

    double *hmati    = gencell->hmati;
    double ecut      = cpcoeffs_info->ecut_psi; // KS-state cutoff in Ryd
    double tpi       = 2.0*M_PI; 

//============================================================================
// Loop over the g-vectors I have, construct g-vectors, compute ig*rho(g)

#ifdef GJM_DEBUG_SIZE
    PRINTF(" %d : coefs in div_rho_gspace_calc\n",size);
#endif
    for(int i = 0; i < size; i++){

       gx = tpi*(k_x[i]*hmati[1] + k_y[i]*hmati[2] + k_z[i]*hmati[3]);
       gy = tpi*(k_x[i]*hmati[4] + k_y[i]*hmati[5] + k_z[i]*hmati[6]);
       gz = tpi*(k_x[i]*hmati[7] + k_y[i]*hmati[8] + k_z[i]*hmati[9]);
      
 
       rhoIGX[i] = (gradientCorrection[i].multiplyByi())*gx;
       rhoIGY[i] = (gradientCorrection[i].multiplyByi())*gy;
       rhoIGZ[i] = (gradientCorrection[i].multiplyByi())*gz;
 
      //if(fabs(gradientCorrection[i].re) > 1e-15)
      //	ckout<<gradientCorrection[i].re<<" "<<gradientCorrection[i].im<<" "<<k_x[i]<<" "<<k_y[i]<<" "<<k_z[i]<<endl;
    }

//============================================================================
}/* end function */
//============================================================================




//============================================================================
/* The following function takes the g-space representation of
*  df/d|grad rho(r)| multiplied by grad rho(r)/|grad rho(r)| and
*  computes a dot product with ig.  The result is stored in gradientCorrection
*  The input is stored in rhoIGX, rhoIGY, rhoIGZ.
*/
//============================================================================

void
CPXCFNCTS::CP_white_byrd_gspace_calc(
	                     complex *rhoIGX, complex *rhoIGY, complex *rhoIGZ, 
                             int *k_x, int *k_y, int *k_z, const int size,
                             const int nf1,const int nf2,const int nf3,
                             complex *gradientCorrection) 
  
  //============================================================================
{/* Begin function */
  //============================================================================
  // Local variables

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

  double *hmati    = gencell->hmati;
  double ecut      = cpcoeffs_info->ecut_psi; // KS-state cutoff in Ryd
  double tpi       = 2.0*M_PI; 
  
  double gx,gy,gz,g2;
  complex tmp;
  double sc; 

  sc =  1.0/((double)(nf1*nf2*nf3));

  //============================================================================
  // Loop over the g-vectors I have, construct g-vectors, compute dot product with ig
  
#ifdef GJM_DEBUG_SIZE
  PRINTF(" %d : CP_white_byrd_gspace_calc\n",size);
#endif

  for(int i = 0; i < size; i++){
      
    gx = tpi*(k_x[i]*hmati[1] + k_y[i]*hmati[2] + k_z[i]*hmati[3]);
    gy = tpi*(k_x[i]*hmati[4] + k_y[i]*hmati[5] + k_z[i]*hmati[6]);
    gz = tpi*(k_x[i]*hmati[7] + k_y[i]*hmati[8] + k_z[i]*hmati[9]);
    g2 = gx*gx+gy*gy+gz*gz;
    
    if(g2 <= 4.0*ecut){  /* in Rydberg */

      // -,+ gives PINY convention
      tmp = rhoIGX[i]*gx + rhoIGY[i]*gy + rhoIGZ[i]*gz;
      gradientCorrection[i] = tmp.multiplyByi(); 

#ifdef CHECK_INPUT_WHITE
      if(k_y[i]==0){
	FILE *fp=fopen("white_fft_input.out","a");
	
        fprintf(fp,"input FFT3 : gx %g gy %g gz %g %d %d %d : %g %g \n",
               gx,gy,gz,k_x[i], k_y[i],k_z[i],
               gradientCorrection[i].re,gradientCorrection[i].im);
	fclose(fp);

      }/*endif*/
#endif

    }else{

      gradientCorrection[i] = complex(0.0,0.0);

    }/*endif*/

      
  }/*end for*/


  
  //============================================================================
}/* end function */
//============================================================================





