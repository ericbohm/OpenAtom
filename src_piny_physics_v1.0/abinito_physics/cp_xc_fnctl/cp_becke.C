#include "standard_include.h"
#include "ckcomplex.h"

#include "../class_defs/CP_OPERATIONS/class_cpxcfnctls.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void CPXCFNCTS::becke_gcx_lda(double rho,double g_rho2,
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
   
   rho43   = pow(rho,(4.0/3.0));
   rho13   = pow(rho,(1.0/3.0));
   gmag    = sqrt(g_rho2);
   x       = gmag/rho43;
   root    = sqrt(1.0 + x*x);
   asinch  = log(x + root);
   bterm   = 1.0 + 6.0*beta*x*asinch;
   bx      = x*x/bterm;
   bx_by_x = x/bterm;

/*=======================================================================*/
/* calculate b'(x)  */

   bpterm = (asinch + x/root);
   bpx    = 2.0*bx_by_x - sixb*(bx_by_x)*(bx_by_x)*bpterm;

/*=======================================================================*/
/* calculate function and its derivatives wrt */
/* density and the gradient of the density    */

   fn[0]     = -2.0*beta*rho43*bx;
   df_dn[0]  = -beta*(fth*rho13*(bx - x*bpx));
   df_dgn[0] = -beta*bpx;

/*==========================================================================*/
  }/*end function*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void CPXCFNCTS::lyp_gcc(double rho,double g_rho2,double *fn,
                        double *df_dn,double *df_dgn)

/*==========================================================================*/
/*         Begin Routine                                                    */
  {/*Begin Routine*/
/*=======================================================================*/
/*  local variables   */
   double rho2,drho;
   double denom,denom2,cexp,omega,domega,delta,ddelta;
   double arho,darho;
   double dprho;
   double rho13,rho113,rho43;
   double term1,brack1,brack2,brack3,term2;
   double dterm1,dterm2;

/* static variables */
   static double a = 0.049180;
   static double b = 0.1320;
   static double c = 0.25330;
   static double d = 0.3490;
   static double pi = M_PI;
   static double pow13;
   static double pow113;
   static double pow43;
   static double pow53;
   static double pow83;
   static double cf;

/*=======================================================================*/
/* I) Get square root density coefficient                                */

   pow13 = -1.0/3.0;
   pow113 = 11.0*pow13;
   pow43 = 4.0*pow13;
   pow53 = -5.0*pow13;
   pow83 = -8.0*pow13;
   cf = 0.30*pow((3.0*pi*pi),(2.0/3.0));
   drho = sqrt(g_rho2);

/*=======================================================================*/
/*  II) evaluate some grouped quantities                                 */

   rho13 = pow(rho,pow13);
   rho113 = pow(rho,pow113);
   rho43 = rho13/rho;
   denom = 1.0 + d*rho13;
   denom2 = denom*denom;
   cexp = exp(-c*rho13);

   omega = (cexp*rho113)/denom;
   brack1 = (d*rho43 + denom*(c*rho43 - 11.0/rho))/3.0;
   domega = omega*brack1/denom;

   delta = c*rho13 + d*rho13/denom;
   brack1 = (d*rho13 - denom)*d*rho43/denom2;
   ddelta = (-c*rho43 + brack1)/3.0;

/*=======================================================================*/
/* III) get all the brackets                                             */

   rho2 = rho*rho;
   term1 = -a*rho/denom;
   brack1 = ((5.0 - 7.0/6.0*delta)/3.0)*g_rho2;
   brack2 = 4.0*cf*pow(rho,pow83) + brack1;
   brack3 = (0.250*brack2 - 11.0*g_rho2/24.0)*rho2;
   term2 = -a*b*omega*brack3;

   arho = 0.250*brack2 - 11.0*g_rho2/24.0;
   brack1 = -7.0/18.0*ddelta*g_rho2;
   darho = 0.250*(32.0*cf*(pow(rho,pow53))/3.0 + brack1);

   dprho = 2.0*arho*rho + darho*rho2;
   dterm1 = -a*(1.0 - pow43*d*rho13)/denom2;
   dterm2 = -a*b*(domega*brack3 + omega*dprho);

/*=======================================================================*/
/* IV) evaluate the function                                             */

   *fn = term1 + term2;
/*=======================================================================*/
/*  V) evaluate the derivatives of the funcion                           */

    *df_dn = dterm1 + dterm2;
    brack1 = 0.50*((5.0 - 7.0/6.0*delta)/3.0) - 11.0/12.0;
    *df_dgn = -a*b*omega*brack1*rho2*drho;

/*==============================================================================*/
   }/*end routine*/
/*==============================================================================*/
