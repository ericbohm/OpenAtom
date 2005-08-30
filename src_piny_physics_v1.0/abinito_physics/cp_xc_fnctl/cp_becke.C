#include "standard_include.h"
#include "ckcomplex.h"

#include "../class_defs/CP_OPERATIONS/class_cpxcfnctls.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void 
CPXCFNCTS::becke_gcx_lda(double rho,double g_rho2,
                         double *fn,double *df_dn,double *df_dgn,
                         double beta)

/*==========================================================================*/
/*         Begin Routine                                                    */
{/*Begin Function*/
/*=======================================================================*/
/*local variables */
   double gmag;
   double rho43,rho13,x,root,asinch;
   double bterm,bx;
   double bpterm,bpx,bx_by_x;

/*   static double beta; */
   static double sixb;
   static double twelveb;
   static double fth;

/*   beta = 0.0042; */
   sixb = 6.0*beta;
   twelveb = 12.0*beta;
   fth = 4.0/3.0;


/*=======================================================================*/
/* calculate b(x)   the Becke function */

   rho    *= 0.5;
   g_rho2 *= 0.25;

   rho43 = pow(rho,(4.0/3.0));
   rho13 = pow(rho,(1.0/3.0));
   gmag = sqrt(g_rho2);
   x = gmag/rho43;
   root = sqrt(1.0 + x*x);
   asinch = log(x + root);
   bterm = 1.0 + 6.0*beta*x*asinch;
   bx = x*x/bterm;
   bx_by_x = x/bterm;

/*=======================================================================*/
/* calculate b'(x)  */

   bpterm = (asinch + x/root);
   bpx = 2.0*bx_by_x - sixb*(bx_by_x)*(bx_by_x)*bpterm;

/*=======================================================================*/
/* calculate function and its derivatives wrt */
/* density and the gradient of the density    */

   *fn = -2.0*beta*rho43*bx;
   *df_dn = -beta*(fth*rho13*(bx - x*bpx));
   *df_dgn = -beta*bpx;

/*==========================================================================*/
  }/*end function*/
/*==========================================================================*/
