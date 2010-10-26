#include "standard_include.h"
#include "ckcomplex.h"

#include "../class_defs/CP_OPERATIONS/class_cpxcfnctls.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void CPXCFNCTS::pbe_gcx_lda(double rho,double g_rho2,
                            double *fn,double *df_dn,double *df_dgn)

/*==========================================================================*/
/*         Begin Routine                                                    */
{/*Begin Function*/
/*=======================================================================*/
/*local variables */

    double s,fx,dfx_dn,dfx_dgn;
    double denom,denom2;
    double gmag;
    double kf;
    double rhoth;
    double rhofth;

    /* Constants */

    static double pi = M_PI;
    static double pow1 = 1.0/3.0;
    static double cfact;
    static double aa = 3.0*M_PI*M_PI;
    static double mu = 0.21951;
    static double kappa = 0.804;
    static double fth = 4.0/3.0;
    static double sth = 7.0/3.0;

    /*==============================================================================*/
    /* Calculate the F_X function */

    cfact = -(0.75/M_PI)*pow((3.0*M_PI*M_PI),pow1);
    gmag = sqrt(g_rho2);
    rhoth = pow(rho,pow1);
    kf = aa*rhoth;
    s = gmag/(2.0*kf*rho);
    denom = 1.0 + mu*s*s/kappa;
    fx = kappa - kappa/denom;

    /*==============================================================================*/
    /* Calculate the functional */

    *fn = cfact*rho*rhoth*fx;

    /*==============================================================================*/
    /* Calculate the derivative of the F_X function wrt density */

    denom2 = denom*denom;
    dfx_dn = -8.0*mu*s*s/(3.0*denom2*rho);


    /*==============================================================================*/
    /* Calculate the derivative of the functional wrt density */

    rhofth = rho*rhoth;
    *df_dn = fth*cfact*rhoth*fx + cfact*rhofth*dfx_dn;

    /*==============================================================================*/
    /* Calculate the derivative of the F_X function wrt |grad n| */

    dfx_dgn = mu*s/(denom2*aa*rhofth);

    /*==============================================================================*/
    /* Calculate the derivative of the functional wrt |grad n| */

    *df_dgn = cfact*rhofth*dfx_dgn;

/*==============================================================================*/
    }/*end routine*/
/*==============================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void CPXCFNCTS::pbe_gcc_lda(double rho,double g_rho2,double *fn,double *df_dn,double *df_dgn)

/*==========================================================================*/
/*         Begin Routine
 */
{/*Begin Routine*/
  /*=======================================================================*/
  /*local variables */

  double ec,dec_dn,rr;
  double rs,drs;
  double srs,rstt,rs2,rs4;
  double gmag;
  double logterm;
  double preterm;
  double kf,ks,t,t2,t4;
  double rhoth;
  double expf,Afun,A2;
  double num,denom;
  double Hfun,Hlogterm;
  double rat;
  double Bfun;

  double dt_dn;
  double dA_dn;
  double dP_dn;
  double dB_dn;
  double dH_dn;

  double dt_dgn;
  double dP_dgn;
  double dB_dgn;
  double dH_dgn;

  /* static variables for lda functional         */
  static double a = 0.03109070;
  static double aa = 3.0*M_PI*M_PI;
  static double alpha1 = 0.213700;
  static double beta1  = 7.59570;
  static double beta2  = 3.58760;
  static double beta3  = 1.63820;
  static double beta4  = 0.492940;

  /* static variables for gc part */

  static double gamma = 0.031091;
  static double beta = 0.066725;
  static double pow1 = 1.0/3.0;
  static double bet_by_gam;

  /*==============================================================================*/
  /* First do e_c and its derivative  */

  rs = pow((3.0/(4.0*M_PI*rho)),pow1);
  rs2 = rs*rs;
  rs4 = rs2*rs2;
  rstt = pow(rs,1.5);
  srs = sqrt(rs);
  rr = 2.0*a*(beta1*srs + beta2*rs + beta3*rstt + beta4*rs2);
  logterm = log(1.0 + 1.0/rr);
  preterm = 1.0 + alpha1*rs;
  ec = -2.0*a*preterm*logterm;


  drs = 2.0*a*(0.5*beta1/srs + beta2 + 1.5*beta3*srs + 2.0*beta4*rs);
  dec_dn = (8.0*M_PI*a*rs4/9.0)*(alpha1*logterm - preterm*drs/(rr*(1.0 + rr)));


  /*==============================================================================*/
  /* Calculate the H function  */

  bet_by_gam = beta/gamma;
  gmag = sqrt(g_rho2);
  rhoth = pow(rho,pow1);
  kf = aa*rhoth;
  ks = sqrt(4.0*kf/M_PI);
  t = gmag/(2.0*ks*rho);

  /* Afun */

  expf = exp(-ec/gamma) - 1.0;
  Afun = bet_by_gam/expf;

  /* rest of H */

  t2 = t*t;
  t4 = t2*t2;
  A2 = Afun*Afun;
  num = 1.0 * Afun*t2;
  denom = 1.0 + Afun*t2 + A2*t4;
  rat = num/denom;
  Bfun = t2*bet_by_gam*rat;
  Hlogterm = log(1.0 + Bfun);
  Hfun = gamma*Hlogterm;

  /*==============================================================================*/
  /* Calculate the functional  */

  *fn = rho*Hfun;

  /*==============================================================================*/
  /* Calculate the derivative of Hfun wrt density  */

  dt_dn = -7.0*t/(6.0*rho);
  dA_dn = bet_by_gam*(expf + 1.0)*dec_dn/(expf*expf*gamma);
  dP_dn = (t2*dA_dn + 2.0*Afun*t*dt_dn)/denom
    - (num/(denom*denom))*(t2*(1.0 + 2.0*Afun*t2)*dA_dn
			   + 2.0*Afun*t*(1.0 + 2.0*Afun*t2)*dt_dn);

  dB_dn = bet_by_gam*(2.0*t*dt_dn*rat + t2*dP_dn);
  dH_dn = gamma*dB_dn/(1.0 + Bfun);

  /*==============================================================================*/
  /* Calculate the derivative of the functional wrt density  */

  *df_dn = Hfun + rho*dH_dn;

  /*==============================================================================*/
  /*  Calculate the derivative of the Hfun wrt |grad n|  */

  dt_dgn = 1.0/(2.0*ks*rho);
  dP_dgn = 2.0*Afun*t*dt_dgn/denom
    - num*2.0*Afun*t*(1.0 + 2.0*Afun*t2)*dt_dgn/(denom*denom);
  dB_dgn = bet_by_gam*(2.0*t*dt_dgn*rat + t2*dP_dgn);
  dH_dgn = gamma*dB_dgn/(1.0 + Bfun);

  /*==============================================================================*/
  /* Calculate the derivative of the functional wrt density  */

  *df_dgn = rho*dH_dgn;

/*==============================================================================*/
   }/*end routine*/
/*==============================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/




