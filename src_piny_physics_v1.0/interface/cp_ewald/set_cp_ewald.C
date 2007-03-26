/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_ewald                                    */
/*                                                                          */
/* This reads in and sets up the electron-atom interaction pseudopotential  */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_mdintegrate.h"
#include "../class_defs/allclass_mdatoms.h"
#include "../class_defs/allclass_mdinter.h"
#include "../proto_defs/proto_cp_ewald_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"

#define CP_EMAGIC 0.00050


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_cpmass(int ncoef,int *kastore,int *kbstore,int *kcstore,
                double *cmass,double *hmati,
                double *cmass_tau_def,double cmass_cut_def,
                int *icmass_unif)

/*==========================================================================*/
/*            Begin subprogram:                                          */
    {/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations:                               */

  double aka,akb,akc,xk,yk,zk;              /* Num: k-vector components */
  double tpi,enow,cmass0;                /* Num: Some useful consts  */
  int i,iii;                              /* Num: For loop counter    */

/*=======================================================================*/
/*  A) Set up CP masses                                                  */

  *icmass_unif = 1;
  tpi = 2.0*M_PI;
  *cmass_tau_def /= TIME_CONV;
  cmass0 = 4.0*(*cmass_tau_def)*(*cmass_tau_def)*CP_EMAGIC;
  for(i=1;i<=ncoef-1;i++) {
    aka = (double)(kastore[i]);
    akb = (double)(kbstore[i]);
    akc = (double)(kcstore[i]);
    xk = (aka*hmati[1]+akb*hmati[2]+akc*hmati[3])*tpi;
    yk = (aka*hmati[4]+akb*hmati[5]+akc*hmati[6])*tpi;
    zk = (aka*hmati[7]+akb*hmati[8]+akc*hmati[9])*tpi;
    enow = (xk*xk+yk*yk+zk*zk)*0.5;
    if(enow > (0.5*cmass_cut_def)) {
      cmass[i] =  cmass0*enow/(0.5*cmass_cut_def);
      *icmass_unif = 0;
    } else {
      cmass[i] =  cmass0;
    } /* endif */
  } /* endfor */

  cmass[ncoef] = 0.5*cmass0;


/*-------------------------------------------------------------------------*/
} /* end routine */
/*=======================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void check_kvec(int nktot,int kastore[],int kbstore[],int kcstore[],
                 int nktot_sm,
                 int kastore_sm[],int kbstore_sm[],int kcstore_sm[])

/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations:                               */
      int i,icount;                           /* Num: Counters              */
/*==========================================================================*/
      icount = 0;

      for(i=1;i<=nktot;i++) {
       if((kastore[i] == kastore_sm[icount+1]) &&
          (kbstore[i] == kbstore_sm[icount+1]) &&
          (kcstore[i] == kcstore_sm[icount+1])) {
        icount++;
       } /* endif */
      } /* endfor */
      if(icount != nktot_sm) {
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("Mismatch in large and small kvectors\n");
        PRINTF("%d vs %d\n",icount,nktot_sm);
        PRINTF("Contact technical support\n");
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
      } /* endif */

/*-------------------------------------------------------------------------*/
     } /* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void calc_cutoff(int kmax_ewd, double *ecut,double *ecut_cp,int cp_on,
                 int *kmax_cp, int *kmaxv, double *hmatik, double deth, 
                 int *nfft, int fft_opt)

/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations:                               */

  int iii,na,nb,nc;
  double rtwoth,tpi,rvol23;
  double d1,d2,d3;
  double try1,try2,try3;
  double temp1,temp2,temp3;
  int itemp1,itemp2,itemp3;
  double frac1,frac2,frac3;

/*==========================================================================*/
/* II) Calculate Ewald Cutoff */

   rtwoth = -(2./3.);
   rvol23 = pow(deth,rtwoth);
   tpi = M_PI * 2.0;
   *ecut = M_PI * .5 * M_PI * (double) (kmax_ewd * kmax_ewd)*rvol23;

/*==========================================================================*/
/* III) Compare Ewald to CP cutoff  */

   if (cp_on == 1) {
     if (*ecut_cp * .5 < *ecut) {
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
	PRINTF("Ewald cutoff greater than cp cutoff %g vs %g \n",
                                     (*ecut)*2.,*ecut_cp);
        PRINTF("Therefore, you must lower the maximum k-vector\n");
        PRINTF("required by your Ewald sum, perhaps change alp_ewd\n");
        PRINTF("and try again\n");
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
     }/*endif*/
     d1 = *ecut, d2 = *ecut_cp * .5;
     *ecut = MAX(d1,d2);
   }/*endif*/
   *ecut_cp = *ecut;

/*==========================================================================*/
/* III) Adjust shape of reciprocal space                                    */

   d1 = hmatik[1];  d2 = hmatik[4];   d3 = hmatik[7];
   try1 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
   d1 = hmatik[2];  d2 = hmatik[5];   d3 = hmatik[8];
   try2 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
   d1 = hmatik[3];  d2 = hmatik[6];   d3 = hmatik[9];
   try3 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
   temp1 = sqrt(*ecut * .5) / (M_PI * try1);
   temp2 = sqrt(*ecut * .5) / (M_PI * try2);
   temp3 = sqrt(*ecut * .5) / (M_PI * try3);
   if(cp_on==1){
      itemp1     = (int)temp1;
      itemp2     = (int)temp2;
      itemp3     = (int)temp3;
      frac1      = temp1-(double)itemp1;
      frac2      = temp2-(double)itemp2;
      frac3      = temp3-(double)itemp3;
      if(frac1>0.99)itemp1++;
      if(frac2>0.99)itemp2++;
      if(frac3>0.99)itemp3++;
      kmax_cp[1] = itemp1;
      kmax_cp[2] = itemp2;
      kmax_cp[3] = itemp3;
      radixme(kmax_cp[1],kmax_cp[2],kmax_cp[3],&na,&nb,&nc,fft_opt);
      kmaxv[1] = 2*kmax_cp[1];
      kmaxv[2] = 2*kmax_cp[2];
      kmaxv[3] = 2*kmax_cp[3];
      nfft[1]  = na;
      nfft[2]  = nb;
      nfft[3]  = nc;
   }else{
     itemp1     = (int)(2.0*temp1);
     itemp2     = (int)(2.0*temp2);
     itemp3     = (int)(2.0*temp3);
     frac1      = temp1-(double)itemp1;
     frac2      = temp2-(double)itemp2;
     frac3      = temp3-(double)itemp3;
     if(frac1>0.99)itemp1++;
     if(frac2>0.99)itemp2++;
     if(frac3>0.99)itemp3++;
     kmaxv[1] = itemp1;
     kmaxv[2] = itemp2;
     kmaxv[3] = itemp3;
    }/*endif*/

/*-------------------------------------------------------------------------*/
     } /* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Count the number of k vectors on large grid */
/*==========================================================================*/

void countkvec3d(int *nktot,double ecut,int *kmaxv,double *hmatik, 
                 double *gmin_spl,double *gmin_true,double *gmax_spl)
/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

  int iii,icount;
  int i, kbmin, kcmin, kbmax, kcmax, kamax, ka, kb, kc;
  double xk, yk, zk;
  double tryme;
  double tpi;
  double aka, akb, akc,g;

/*==========================================================================*/
/* Count the kvectors */

  tpi = 2.0*M_PI;
  icount = 0;
  kamax = kmaxv[1];

  gmin_spl[0] = 1.0e10;
  gmax_spl[0] = 0.0;

/*********************************/

   for (i = 1; i <= kamax; ++i) {
    aka = (double) i;
    xk = aka * hmatik[1] * tpi;
    yk = aka * hmatik[4] * tpi;
    zk = aka * hmatik[7] * tpi;
    tryme = (xk * xk + yk * yk + zk * zk) * .5;
    if (tryme > ecut * 4.) {break;}
  }

  kamax = i - 1;

/***********************************/

  for (ka = 0; ka <= kamax; ++ka) {
    aka = (double) ka;
    kbmin = -kmaxv[2];
    if (ka == 0) {
      kbmin = 0;
    }
    for (i = kbmin; i <= 0; ++i) {
      akb = (double) i;
      xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
      yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
      zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
      tryme = (xk * xk + yk * yk + zk * zk) * .5;
      if (tryme <= ecut * 4.) {break;}
    }

/*********************************/

    kbmin = i;
    for (i = 1; i <= kmaxv[2]; ++i) {
      akb = (double) i;
      xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
      yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
      zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
      tryme = (xk * xk + yk * yk + zk * zk) * .5;
      if (tryme > ecut * 4.) {break;}
    }
    
    kbmax = i - 1;
    for (kb = kbmin; kb <= kbmax; ++kb) {
      
      akb = (double) kb;
      kcmin = -kmaxv[3];
      if (ka == 0 && kb == 0) {
	kcmin = 1;
      }
      for (i = kcmin; i <= 0; ++i) {
	akc = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	tryme = (xk * xk + yk * yk + zk * zk) * .5;
	if (tryme <= ecut * 4.) {break;}
      }
/*********************************/

      kcmin = i;
      for (i = 1; i <= kmaxv[3]; ++i) {
	akc = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	tryme = (xk * xk + yk * yk + zk * zk) * .5;
	if (tryme > ecut * 4.) {break;}
      }//endfor : kc

      kcmax = i - 1;
      akc = (double) kcmin;
      for (kc = kcmin; kc <= kcmax; ++kc) {
	akc = (double) kc;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
        g  = sqrt(xk * xk + yk * yk + zk * zk);
        gmin_spl[0] = MIN(gmin_spl[0],g);
        gmax_spl[0] = MAX(gmax_spl[0],g);
	++icount;
      }//endfor : kc

    }//endfor : kb
  }//end for :ka
  *nktot = icount;

  gmin_true[0] = gmin_spl[0];
  gmin_spl[0] *= 0.75;
  gmax_spl[0] *= (4.0/3.0);

/*--------------------------------------------------------------------------*/
  } /* countkvec3d */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* This subroutine determines the allowed spherically truncated             */
/* half space k (.i.e. g) vectors given a cutoff. it is used by             */
/* both the ewald and cp modules.                                           */
/* Sets up the k vectors for a given cutoff and shape                       */
/*==========================================================================*/

void setkvec3d(int nktot,double ecut,int *kmaxv,double *hmatik,
		int *kastore, int *kbstore, int *kcstore, 
                int *ibreak1, int *ibreak2, int cp_on, 
                double *gmin_spl, double *gmin_true,double *gmax_spl)

/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

    int iii;
    int i1, i2, i3;

    int i, kbmin, kcmin, kbmax, kcmax, kamax, ka, kb, kc;
    double xk, yk, zk;
    int icount;
    double aka, akb, akc;
    double tpi;
    double tryme;

/*==========================================================================*/
/* Count the K-vectors */

    tpi = 2.0*M_PI;
    for(i=1;i<=nktot;i++){
      ibreak1[i] = 0;
      ibreak2[i] = 0;
    }/*endfor*/
    icount = 0;
    (*gmin_spl) = 1.0e10;
    (*gmax_spl) = 0.0;

/*=============================*/

    kamax = kmaxv[1];
    i1 = kamax;
    for (ka = 0; ka <= i1; ++ka) {
      aka = (double) ka;
      kbmin = -kmaxv[2];
      if (ka == 0) {
	kbmin = 0;
      }
      for (i = kbmin; i <= 0; ++i) {
	akb = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
	tryme = (xk * xk + yk * yk + zk * zk) * .5;
	if (tryme <= ecut * 4.) {
	  break;
	}
      }
      kbmin = i;
      i2 = kmaxv[2];
      for (i = 1; i <= i2; ++i) {
	akb = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
	tryme = (xk * xk + yk * yk + zk * zk) * .5;
	if (tryme > ecut * 4.) {
	  break;
	}
      }

/*=============================*/

      kbmax = i - 1;
      i2 = kbmax;
      for (kb = kbmin; kb <= i2; ++kb) {
	akb = (double) kb;
	kcmin = -kmaxv[3];
	if (ka == 0 && kb == 0) {
	  kcmin = 1;
	}
	for (i = kcmin; i <= 0; ++i) {
	  akc = (double) i;
	  xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	  yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	  zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	  tryme = (xk * xk + yk * yk + zk * zk) * .5;
	  if (tryme <= ecut * 4.) {
	    break;
	  }
	}

/*=============================*/

	kcmin = i;
	i3 = kmaxv[3];
	for (i = 1; i <= i3; ++i) {
	  akc = (double) i;
	  xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	  yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	  zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	  tryme = (xk * xk + yk * yk + zk * zk) * .5;
	  if (tryme > ecut * 4.) {
	    break;
	  }
	}
	kcmax = i - 1;
	i3 = kcmax;
	for (kc = kcmin; kc <= i3; ++kc) {
	  ++icount;
	  aka = (double) ka;
	  akb = (double) kb;
	  akc = (double) kc;
	  xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	  yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	  zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	  kastore[icount] = ka;
	  kbstore[icount] = kb;
	  kcstore[icount] = kc;
	  tryme = sqrt(xk * xk + yk * yk + zk * zk);
	  (*gmin_spl) = MIN(tryme,(*gmin_spl));
	  (*gmax_spl) = MAX(tryme,(*gmax_spl));
	  if (kc == kcmin) {
	    ibreak1[icount] = 1;
	  }
	  if (kc < kcmax) {
	    ibreak2[icount] = 1;
	  }
	}
      }
    }
    if(nktot!=icount){
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("Mismatch number of kvectors\n");
        PRINTF("%d vs %d\n",icount,nktot);
        PRINTF("Contact technical support\n");
        PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        FFLUSH(stdout);
        EXIT(1);
    }
    if (cp_on == 1) {
      ++icount;
      kastore[icount] = 0;
      kbstore[icount] = 0;
      kcstore[icount] = 0;
      (*gmin_true) = (*gmin_spl);
      (*gmin_spl) *= 0.75;
      (*gmax_spl) *= (4.0/3.0);
    }
/*--------------------------------------------------------------------------*/
} /* setkvec3d */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* This subroutine determines the allowed spherically truncated             */
/* half space k (.i.e. g) vectors given a cutoff. It is used by             */
/* cp modules. Sets up the k vectors for a given cutoff and shape           */
/*==========================================================================*/

void setkvec3d_sm(int nktot,double ecut,int *kmax_cp,double *hmatik,
                  int *kastore, int *kbstore, int *kcstore, 
		  int *ibreak1, int *ibreak2, double *gmin, double *gmax)

/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

  int i1, i2, i3;

  int i, kbmin, kcmin, kbmax, kcmax, kamax, ka, kb, kc;
  double xk, yk, zk,g;
  int icount;
  double aka, akb, akc;
  double tpi, tryme;
  tpi = M_PI * 2.;

/*==========================================================================*/
/* SETUP THE KVECTORS */
  
  for(i=1;i<=nktot;i++){
      ibreak1[i] = 0;
      ibreak2[i] = 0;
  }

  *gmin = 10000.0;
  *gmax = .0;
  icount = 0;

  i1 = kmax_cp[1];
  for (i = 1; i <= i1; ++i) {
    aka = (double) i;
    xk = aka * hmatik[1] * tpi;
    yk = aka * hmatik[4] * tpi;
    zk = aka * hmatik[7] * tpi;
    tryme = (xk * xk + yk * yk + zk * zk) * .5;
    if (tryme > ecut) {
      break;
    }
  }

  kamax = i - 1;
  i1 = kamax;
  for (ka = 0; ka <= i1; ++ka) {
    aka = (double) ka;
    kbmin = -kmax_cp[2];
    if (ka == 0) {
      kbmin = 0;
    }
    for (i = kbmin; i <= 0; ++i) {
      akb = (double) i;
      xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
      yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
      zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
      tryme = (xk * xk + yk * yk + zk * zk) * .5;
      if (tryme <= ecut) {
	break;
      }
    }
    kbmin = i;
    i2 = kmax_cp[2];
    for (i = 1; i <= i2; ++i) {
      akb = (double) i;
      xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
      yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
      zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
      tryme = (xk * xk + yk * yk + zk * zk) * .5;
      if (tryme > ecut) {
	break;
      }
    }

    kbmax = i - 1;
    i2 = kbmax;
    for (kb = kbmin; kb <= i2; ++kb) {
      akb = (double) kb;
      kcmin = -kmax_cp[3];
      if (ka == 0 && kb == 0) {
	kcmin = 1;
      }
      for (i = kcmin; i <= 0; ++i) {
	akc = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	tryme = (xk * xk + yk * yk + zk * zk) * .5;
	if (tryme <= ecut) {
	  break;
	}
      }
      
      kcmin = i;
      i3 = kmax_cp[3];
      for (i = 1; i <= i3; ++i) {
	akc = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	tryme = (xk * xk + yk * yk + zk * zk) * .5;
	if (tryme > ecut) {
	  break;
	}
      }
      
      kcmax = i - 1;
      i3 = kcmax;
      for (kc = kcmin; kc <= i3; ++kc) {
	++icount;
	aka = (double) ka;
	akb = (double) kb;
	akc = (double) kc;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
        g  = sqrt(xk*xk+yk*yk+zk*zk);
        *gmin = MIN(*gmin,g);
        *gmax = MAX(*gmax,g);
	kastore[icount] = ka;
	kbstore[icount] = kb;
	kcstore[icount] = kc;
	if (kc == kcmin) {
	  ibreak1[icount] = 1;
	}
	if (kc < kcmax) {
	  ibreak2[icount] = 1;
	}
      }
    }
  }
  if(nktot!=icount){
       PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       PRINTF("Mismatch number of small kvectors\n");
       PRINTF("%d vs %d\n",icount,nktot);
       PRINTF("Contact technical support\n");
       PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       FFLUSH(stdout);
       EXIT(1);
  }

  ++icount;
  kastore[icount] = 0;
  kbstore[icount] = 0;
  kcstore[icount] = 0;

/*--------------------------------------------------------------------------*/
    }/* setkvec3d_sm */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Determines which k vectors of the full set belong to respa               */
/*==========================================================================*/

void setkvec3d_res(int kmax_res, double *hmatik, 
		  int *kastore, int *kbstore, int *kcstore, int *ibreak3, 
		  int nktot, int nktot_res)

/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

  int i1;
  
  double rvol23, akmax2_res;
  int ka, kb, kc;
  double xk, yk, zk;
  int icount,jcount;
  double aka, akb, akc;
  double tryme,voli,temp;

/*==========================================================================*/
  /* START THE ROUTINE */

  voli =  getdeth(hmatik);
  temp =  (2./3.);
  rvol23 = pow(voli,temp);
  akmax2_res = (double) ((kmax_res) * (kmax_res));
  /* SETUP THE KVECTORS */
  i1 = nktot;
  jcount = 0;
  for (icount = 1; icount <= i1; ++icount) {
    ka = kastore[icount];
    kb = kbstore[icount];
    kc = kcstore[icount];
    aka = (double) ka;
    akb = (double) kb;
    akc = (double) kc;
    xk = aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3];
    yk = aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6];
    zk = aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9];
    tryme = (xk * xk + yk * yk + zk * zk)/rvol23;
    ibreak3[icount] = 0;
    if (tryme <= akmax2_res) {
      ibreak3[icount] = 1;
      jcount +=1;
    }
  }
  if(jcount!=nktot_res){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Error setting break3 in set_recip.c\n");
    PRINTF("%d k-vectors found versus %d expected\n",jcount,nktot_res);
    PRINTF("Your respa kmax is %d\n",kmax_res);
    PRINTF("*Contact technical support*\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }
/*--------------------------------------------------------------------------*/
} /* setkvec3d_res */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* This routine makes sure that the k1,k2,k3 for the fft                    */
/*  satisfy the radix condition                                             */
/*==========================================================================*/

void radixme(int kmax1, int kmax2, int kmax3, int *n1, int *n2, int *n3,
             int fft_opt)

/*==========================================================================*/
/* Calculate the quantity to be radicized */
/* it written with the plus one to get rid */
/* of the stupid annoying edge vectors. */
/* the factor of 4 appears because the kmax's are the */
/* maximum k vector along a direction. the normal fft */
/* grid is therefore (2*(kmax1+1))(2*(kmax2+1))(2*(kmax3+1)). */
/* the additional factor of two comes from the fact the density */
/* needs to be defined on twice as fine an fft grid */
/* (4*(kmax1+1))(4*(kmax2+1))(4*(kmax3+1)). */
/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

  int i1,i2,i3;
  int i, k1, k2, k3;
  int kk1, kk2, kk3;

  int nrad;
  int krad[181]; 
  int nrad_in = 180;

/*==========================================================================*/

  set_fftsizes(nrad_in,&nrad,krad,fft_opt); /* radix conditions */

  kk1 = 2*(2*kmax1 + 1);
  kk2 = 2*(2*kmax2 + 1);
  kk3 = 2*(2*kmax3 + 1);

  i1=0;
  for(i=nrad;i>=1;i--){
    if(krad[i]>=kk1){k1 = krad[i]; i1=i;}
  }//endfor
  if(i1==0){PRINTF("Bad Radix\n");EXIT(1);}

  i2=0;
  for(i=nrad;i>=1;i--){
    if(krad[i]>=kk2){k2 = krad[i]; i2=i;}
  }//endfor
  if(i2==0){PRINTF("Bad Radix\n");EXIT(1);}

  i3=0;
  for(i=nrad;i>=1;i--){
    if(krad[i]>=kk3){k3 = krad[i]; i3=i;}
  }//endfor
  if(i3==0){PRINTF("Bad Radix\n");EXIT(1);}

  PRINTF("   Radix data : k1 %d kmax1 %d kk1 %d\n",k1,kmax1,kk1);
  PRINTF("   Radix data : k2 %d kmax2 %d kk2 %d\n",k2,kmax2,kk2);
  PRINTF("   Radix data : k3 %d kmax3 %d kk3 %d\n",k3,kmax3,kk3);

  n1[0]    = k1;
  n2[0]    = k2;
  n3[0]    = k3;

/*--------------------------------------------------------------------------*/
  } /* radixme */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Count number of k vectors on small grid */
/*==========================================================================*/

void countkvec3d_sm(int *nktot, double ecut, int *kmax_cp, double *hmatik ,
                    double *gmin, double *gmax)

/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

  int i1, i2, i3;
  int i, kbmin, kcmin, kbmax, kcmax, kamax, ka, kb, kc;
  double xk, yk, zk;
  int icount;
  double aka, akb, akc,g;
  double tpi, tryme;

/*==========================================================================*/

  tpi = M_PI * 2.;
  gmin[0] = 1.0e10;
  gmax[0] = 0.0;

  icount = 0;
  i1 = kmax_cp[1];
  for (i = 1; i <= i1; ++i) {
    aka = (double) i;
    xk = aka * hmatik[1] * tpi;
    yk = aka * hmatik[4] * tpi;
    zk = aka * hmatik[7] * tpi;
    tryme = (xk * xk + yk * yk + zk * zk) * .5;
    if (tryme > ecut) {break;}
  }
  kamax = i - 1;
  i1 = kamax;
  for (ka = 0; ka <= i1; ++ka) {
    aka = (double) ka;
    kbmin = -kmax_cp[2];
    if (ka == 0) {
      kbmin = 0;
    }
    for (i = kbmin; i <= 0; ++i) {
      akb = (double) i;
      xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
      yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
      zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
      tryme = (xk * xk + yk * yk + zk * zk) * .5;
      if (tryme <= ecut) {break;}
    }

    kbmin = i;
    i2 = kmax_cp[2];
    for (i = 1; i <= i2; ++i) {
      akb = (double) i;
      xk = (aka * hmatik[1] + akb * hmatik[2]) * tpi;
      yk = (aka * hmatik[4] + akb * hmatik[5]) * tpi;
      zk = (aka * hmatik[7] + akb * hmatik[8]) * tpi;
      tryme = (xk * xk + yk * yk + zk * zk) * .5;
      if (tryme > ecut) {break;}
    }

    kbmax = i - 1;
    i2 = kbmax;
    for (kb = kbmin; kb <= i2; ++kb) {
      akb = (double) kb;
      kcmin = -kmax_cp[3];
      if (ka == 0 && kb == 0) {
	kcmin = 1;
      }
      for (i = kcmin; i <= 0; ++i) {
	akc = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	tryme = (xk * xk + yk * yk + zk * zk) * .5;
	if (tryme <= ecut) {break;}
      }

      kcmin = i;
      i3 = kmax_cp[3];
      for (i = 1; i <= i3; ++i) {
	akc = (double) i;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
	tryme = (xk * xk + yk * yk + zk * zk) * .5;
	if (tryme > ecut) {break;}
      }

      kcmax = i - 1;
      i3 = kcmax;
      for (kc = kcmin; kc <= i3; ++kc) {
	akc = (double) kc;
	xk = (aka * hmatik[1] + akb * hmatik[2] + akc * hmatik[3]) * tpi;
	yk = (aka * hmatik[4] + akb * hmatik[5] + akc * hmatik[6]) * tpi;
	zk = (aka * hmatik[7] + akb * hmatik[8] + akc * hmatik[9]) * tpi;
        g  = sqrt(xk * xk + yk * yk + zk * zk);
        gmin[0] = MIN(gmin[0],g);
        gmax[0] = MAX(gmax[0],g);
	++icount;
      }
    }
  }
  *nktot = icount;

/*--------------------------------------------------------------------------*/
  } /* countkvec3d_sm */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Set the PME maps   */
/*==========================================================================*/

void set_pme_wght(int nktot,int *kastore,int *kbstore,int *kcstore,
                  int nkf1,int nkf2,int nkf3,
                  int ncoef_proc,int ncoef_use,int icoef_off,
                  int pme_b_opt,double *bfact_r,double *bfact_i,
                  double *bweight_tot, int n_interp,
                  double *aj,double *rn,double *rn1)

/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

  size_t dim_k;
  int i,ka,kb,kc,kap,kbp,kcp;
  int ngrid_a,ngrid_b,ngrid_c;
  double tmp_a_r,tmp_a_i;
  double tmp_b_r,tmp_b_i;
  double tmp_ab_r,tmp_ab_i;
  double tmp_c_r,tmp_c_i;
  double *bden_a_r,*bden_a_i,*bweight_a;
  double *bden_b_r,*bden_b_i,*bweight_b;
  double *bden_c_r,*bden_c_i,*bweight_c;
  double *uk,*mn_k;
  int *map_a,*map_b,*map_c;
  
/*==========================================================================*/
/* I) Spherical Map */

  dim_k = (size_t) nktot;
  map_a         = (int *) cmalloc(dim_k*sizeof(int),"set_cp_ewald")-1;
  map_b         = (int *) cmalloc(dim_k*sizeof(int),"set_cp_ewald")-1;
  map_c         = (int *) cmalloc(dim_k*sizeof(int),"set_cp_ewald")-1;
  for(i=1;i<=nktot;i++){
    ka = kastore[i];
    kb = kbstore[i];
    kc = kcstore[i];
    if (kc < 0) {
         kcp = kc + nkf3 + 1;
    } else {
         kcp = kc + 1;
    }/*endif*/
    if (kb < 0) {
         kbp = kb + nkf2 + 1;
    } else {
         kbp = kb + 1;
    }/*endif*/
    if (ka < 0) {
         kap = ka + nkf1 + 1;
    } else {
         kap = ka + 1;
    }/*endif*/
    map_a[i]   = kap;
    map_b[i]   = kbp;
    map_c[i]   = kcp;
  }/*endfor*/

/*==========================================================================*/
/* V) Calculate bweight on the spherically cutoff grid                      */

/*--------------------------------------------------------------------------*/
/*     A) Malloc memory and define constants                                */
   ngrid_a = nkf1;
   ngrid_b = nkf2;
   ngrid_c = nkf3;

   dim_k = (size_t) (ngrid_a);
   bden_a_r  = (double *)cmalloc(dim_k*sizeof(double),"set_cp_ewald")-1;
   bden_a_i  = (double *)cmalloc(dim_k*sizeof(double),"set_cp_ewald")-1;
   bweight_a = (double *)cmalloc(dim_k*sizeof(double),"set_cp_ewald")-1;

   dim_k = (size_t) (ngrid_b);
   bden_b_r  = (double *)cmalloc(dim_k*sizeof(double),"set_cp_ewald")-1;
   bden_b_i  = (double *)cmalloc(dim_k*sizeof(double),"set_cp_ewald")-1;
   bweight_b = (double *)cmalloc(dim_k*sizeof(double),"set_cp_ewald")-1;

   dim_k = (size_t) (ngrid_c);
   bden_c_r  = (double *)cmalloc(dim_k*sizeof(double),"set_cp_ewald")-1;
   bden_c_i  = (double *)cmalloc(dim_k*sizeof(double),"set_cp_ewald")-1;
   bweight_c = (double *)cmalloc(dim_k*sizeof(double),"set_cp_ewald")-1;

   dim_k = (size_t) (n_interp);
   uk   = (double *)cmalloc(dim_k*sizeof(double),"set_cp_ewald")-1;
   mn_k = (double *)cmalloc(dim_k*sizeof(double),"set_cp_ewald")-1;

/*--------------------------------------------------------------------------*/
/*     B) Construct the weighting Function                                  */

   get_bspline_wght1d(n_interp,ngrid_a,aj,rn,rn1,mn_k,uk,
                      bden_a_r,bden_a_i,bweight_a);
   get_bspline_wght1d(n_interp,ngrid_b,aj,rn,rn1,mn_k,uk,
                      bden_b_r,bden_b_i,bweight_b);
   get_bspline_wght1d(n_interp,ngrid_c,aj,rn,rn1,mn_k,uk,
                      bden_c_r,bden_c_i,bweight_c);
   if(pme_b_opt > 0){
    for(i=1;i <= nktot; ++i) {
     bweight_tot[i] = bweight_a[map_a[i]]
                     *bweight_b[map_b[i]]
                     *bweight_c[map_c[i]];
    }/*endfor*/
   }/*endif*/

   if(pme_b_opt == 0 || pme_b_opt ==2){
    for(i=1;i <=ncoef_use; ++i){
      tmp_a_r = bden_a_r[map_a[i+icoef_off]];
      tmp_a_i = bden_a_i[map_a[i+icoef_off]];
      tmp_b_r = bden_b_r[map_b[i+icoef_off]];
      tmp_b_i = bden_b_i[map_b[i+icoef_off]];
      tmp_c_r = bden_c_r[map_c[i+icoef_off]];
      tmp_c_i = bden_c_i[map_c[i+icoef_off]];
      tmp_ab_r = tmp_a_r*tmp_b_r - tmp_a_i*tmp_b_i;
      tmp_ab_i = tmp_a_i*tmp_b_r + tmp_a_r*tmp_b_i;
      bfact_r[i] = tmp_ab_r*tmp_c_r - tmp_ab_i*tmp_c_i;
      bfact_i[i] = tmp_ab_i*tmp_c_r + tmp_ab_r*tmp_c_i;
    }/*endfor*/
    if(ncoef_proc > ncoef_use){
      bfact_r[(ncoef_proc)] = bden_a_r[1]*bden_b_r[1]*bden_c_r[1];
      bfact_i[(ncoef_proc)] = 0.0;
    }
  }/*endif*/

/*========================================================================*/
/* Free memory */

  cfree(&map_a[1],"set_pme_wght");
  cfree(&map_b[1],"set_pme_wght");
  cfree(&map_c[1],"set_pme_wght");
  cfree(&bden_a_r[1],"set_pme_wght"); 
  cfree(&bden_a_i[1],"set_pme_wght"); 
  cfree(&bweight_a[1],"set_pme_wght");
  cfree(&bden_b_r[1],"set_pme_wght"); 
  cfree(&bden_b_i[1],"set_pme_wght"); 
  cfree(&bweight_b[1],"set_pme_wght");
  cfree(&bden_c_r[1],"set_pme_wght"); 
  cfree(&bden_c_i[1],"set_pme_wght"); 
  cfree(&bweight_c[1],"set_pme_wght");
  cfree(&uk[1],"set_pme_wght"); 
  cfree(&mn_k[1],"set_pme_wght"); 
/*--------------------------------------------------------------------------*/
} /* set_pme_map */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Calculate the B spline weighting function   */
/*==========================================================================*/

void get_bspline_wght1d(int n_interp,int ngrid,double *aj,double *rn,
                        double *rn1,double *mn_k,double *uk,
                        double *bden_r,double *bden_i,double *bweight)

/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

  int k,k1,n,m,j;
  double arg;
  double tpi_n,mn_k_tmp;
  double bnum_real,bnum_imag;
  double bden_real,bden_imag,denom;
  double tmp_real,tmp_imag;
  double grid;

/*==========================================================================*/
/* I) Get B spline coefficients                                         */

   grid  = (double) ngrid;
   for(j=1;j<=n_interp;j++){
     aj[j] = (double) (j-1);
     rn[j] = (double) (j);
     if(j > 1){rn1[j] = 1.0/((double)(j-1));}
   }/*endfor*/
   rn1[1] = 0.0;
   mn_k[1] = 1.0; 
   uk[1]   = 1.0;
   for(k=2;k<=n_interp;k++){
     uk[k] = (double) (k);
     mn_k[k] = 0.0; 
   }/*endfor*/
   for(n=3;n<=n_interp;n++){
     for(k=n;k>=2;k--){
       k1 = k-1;
       mn_k_tmp  = (uk[k]*mn_k[k]+(rn[n]-uk[k])*mn_k[k1])*rn1[n];
       mn_k[k] = mn_k_tmp;
     }/*endfor*/
     mn_k[1] = uk[1]*mn_k[1]*rn1[n];
   }/*endfor*/

/*==========================================================================*/
/* II) Transform coeffs                                                     */

/*-----------------------------------------------*/
/*    i)Perform transform using 1D Slow FT */

   tpi_n = 2.0*M_PI/grid;
   for(m=1;m<=ngrid;m++){
     bden_r[m] = 0.0;
     bden_i[m] = 0.0;
     for(k=1;k<=n_interp-1;k++){
       arg = tpi_n*((double)((k-1)*(m-1)));
       bden_r[m] += cos(arg)*mn_k[k];
       bden_i[m] += sin(arg)*mn_k[k];
     }/*endfor*/
   }/*endfor*/

/*==========================================================================*/
/* III) Make separable bweights                                             */

   tpi_n = 2.0*M_PI*((double) (n_interp-1))/grid;
   for(m=1;m<=ngrid;m++){

     arg = tpi_n*((double)(m-1));
     bnum_real  = cos(arg);
     bnum_imag  = sin(arg);
     bden_real  = bden_r[m];
     bden_imag  = bden_i[m];

     denom      = bden_real*bden_real + bden_imag*bden_imag;
     tmp_real   = (bnum_real*bden_real + bnum_imag*bden_imag)/denom;
     tmp_imag   = (bnum_imag*bden_real - bnum_real*bden_imag)/denom;
     bweight[m] =  tmp_real*tmp_real + tmp_imag*tmp_imag;

     bden_r[m] = tmp_real;
     bden_i[m] = tmp_imag;

   }/*endfor*/

/*--------------------------------------------------------------------------*/
     }/*end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void init_nonlocal_ees(int *kmax,double ecut,PSNONLOCAL *psnonlocal,int fft_opt)

/*==========================================================================*/
   {/*begin routine */
/*==========================================================================*/
/* Local variables */

   int i;  

   int n_interp = psnonlocal->n_interp;
   double scale = psnonlocal->fft_size_scale;

   int nk1      = kmax[1]+1;
   int nk2      = kmax[2]+1;
   int nk3      = kmax[3]+1;

   int nkf1t,nkf2t,nkf3t;
   int nkf1,nkf2,nkf3,nrad;

   int krad[181];
   int nrad_in=180;

/*==========================================================================*/
/* (I) Compute the grid size : */

   set_fftsizes(nrad_in,&nrad,krad,fft_opt); /* radix conditions */

   nkf1 = (int)( 2.0*scale*((double)nk1) );
   nkf2 = (int)( 2.0*scale*((double)nk2) );
   nkf3 = (int)( 2.0*scale*((double)nk3) );

   if(nkf1 < 2*nk1){nkf1++;}
   if(nkf2 < 2*nk2){nkf2++;}
   if(nkf3 < 2*nk3){nkf3++;}

   nkf1t = 0;  nkf2t = 0;  nkf3t = 0;
   for(i=1;i<=nrad;i++){if(krad[i]>=nkf1){nkf1t=krad[i];break;}}
   for(i=1;i<=nrad;i++){if(krad[i]>=nkf2){nkf2t=krad[i];break;}}
   for(i=1;i<=nrad;i++){if(krad[i]>=nkf3){nkf3t=krad[i];break;}}

   if( (nkf1t==0) || (nkf2t==0) || (nkf3t==0) ){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("3D-FFT size required too large %d %d %d\n",nkf1,nkf2,nkf3);
     PRINTF("for the hard coded radix conditions. Time to update!\n");
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");    
     FFLUSH(stdout);
     EXIT(1);
   }/*endif*/

   nkf1 = nkf1t;
   nkf2 = nkf2t;
   nkf3 = nkf3t;

   if( (nkf1<n_interp) || (nkf2<n_interp) || (nkf3<n_interp) ){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("The PME n_interp parameter > number of grid points \n");
     PRINTF("This is not allowed\n");      
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     FFLUSH(stdout);
     EXIT(1);      
   }/*endif*/

 
/*==========================================================================*/
/* (III) Pack the structure                                                 */

   psnonlocal->nk1       = nk1;
   psnonlocal->nk2       = nk2;
   psnonlocal->nk3       = nk3;
   psnonlocal->ngrid_a   = nkf1;
   psnonlocal->ngrid_b   = nkf2;
   psnonlocal->ngrid_c   = nkf3;
   psnonlocal->nfft      = nkf1*nkf2;
   psnonlocal->ecut      = ecut;
   psnonlocal->n_interp2 = n_interp*n_interp;

   PRINTF("Your Non-local EES g-space grid is: (-%d,%d) x (-%d,%d) x (-%d,%d)\n",
                  kmax[1],kmax[1],kmax[2],kmax[2],kmax[3],kmax[3]);
   PRINTF("Your Non-local EES r-space grid is : %d x %d x %d\n",
                                       nkf1,nkf2,nkf3);

/*--------------------------------------------------------------------------*/
   }/*end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void init_eext_ees(int *kmax,CPPSEUDO *cppseudo,int fft_opt)

/*==========================================================================*/
   {/*begin routine */
/*==========================================================================*/
/* Local variables */

   int i;  

   int n_interp = cppseudo->n_interp_ps;
   double scale = cppseudo->fft_size_scale_ps;

   int nk1      = 2*kmax[1]+1;
   int nk2      = 2*kmax[2]+1;
   int nk3      = 2*kmax[3]+1;

   int nkf1t,nkf2t,nkf3t;
   int nkf1,nkf2,nkf3,nrad;

   int krad[181];
   int nrad_in=180;

/*==========================================================================*/
/* (I) Compute the grid size : */

   set_fftsizes(nrad_in,&nrad,krad,fft_opt); /* radix conditions */

   nkf1 = (int)( 2.0*scale*((double)nk1) );
   nkf2 = (int)( 2.0*scale*((double)nk2) );
   nkf3 = (int)( 2.0*scale*((double)nk3) );

   if(nkf1 < 2*nk1){nkf1++;}
   if(nkf2 < 2*nk2){nkf2++;}
   if(nkf3 < 2*nk3){nkf3++;}

   nkf1t = 0;  nkf2t = 0;  nkf3t = 0;
   for(i=1;i<=nrad;i++){if(krad[i]>=nkf1){nkf1t=krad[i];break;}}
   for(i=1;i<=nrad;i++){if(krad[i]>=nkf2){nkf2t=krad[i];break;}}
   for(i=1;i<=nrad;i++){if(krad[i]>=nkf3){nkf3t=krad[i];break;}}

   if( (nkf1t==0) || (nkf2t==0) || (nkf3t==0) ){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("3D-FFT size required too large %d %d %d\n",nkf1,nkf2,nkf3);
     PRINTF("for the hard coded radix conditions. Time to update!\n");
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");    
     FFLUSH(stdout);
     EXIT(1);
   }/*endif*/

   nkf1 = nkf1t;
   nkf2 = nkf2t;
   nkf3 = nkf3t;

   if( (nkf1<n_interp) || (nkf2<n_interp) || (nkf3<n_interp) ){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("The PME n_interp parameter > number of grid points \n");
     PRINTF("This is not allowed\n");      
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     FFLUSH(stdout);
     EXIT(1);      
   }/*endif*/
 
/*==========================================================================*/
/* (III) Pack the structure                                                 */

   cppseudo->ngrid_eext_a = nkf1;
   cppseudo->ngrid_eext_b = nkf2;
   cppseudo->ngrid_eext_c = nkf3;
   cppseudo->nka_eext     = 2*kmax[1];
   cppseudo->nkb_eext     = 2*kmax[2];
   cppseudo->nkc_eext     = 2*kmax[3];

   PRINTF("Your Eext EES g-space grid is: (-%d,%d) x (-%d,%d) x (-%d,%d)\n",
             2*kmax[1],2*kmax[1],2*kmax[2],2*kmax[2],2*kmax[3],2*kmax[3]);
   PRINTF("Your Eext EES r-space grid is : %d x %d x %d\n",nkf1,nkf2,nkf3);

/*--------------------------------------------------------------------------*/
   }/*end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_fftsizes(int nrad_in,int *nrad_ret, int *krad, int fft_opt)

/*==========================================================================*/
  {/*begin routine */
/*-------------------------------------------------------------------------*/

  int nrad=179;
 (*nrad_ret)  = nrad;

  if(nrad_in<nrad){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Internal Error in hardcoded radix size array.\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }/*endif*/

  krad[1]   = 4;
  krad[2]   = 6;
  krad[3]   = 8;
  krad[4]   = 10;
  krad[5]   = 12;
  krad[6]   = 14;
  krad[7]   = 16;
  krad[8]   = 18;
  krad[9]   = 20;
  krad[10]  = 22;
  krad[11]  = 24;
  krad[12]  = 28;
  krad[13]  = 30;
  krad[14]  = 32;
  krad[15]  = 36;
  krad[16]  = 40;
  krad[17]  = 42;
  krad[18]  = 44;
  krad[19]  = 48;
  krad[20]  = 56;
  krad[21]  = 60;
  krad[22]  = 64;
  krad[23]  = 66;
  krad[24]  = 70;
  krad[25]  = 72;
  krad[26]  = 80;
  krad[27]  = 84;
  krad[28]  = 88;
  krad[29]  = 90;
  krad[30]  = 96;
  krad[31]  = 112;
  //  if(fft_opt==0){krad[31]  = 100;}
  krad[32]  = 112;
  krad[33]  = 120;
  krad[34]  = 128;
  krad[35]  = 128;
  krad[36]  = 132;
  krad[37]  = 140;
  krad[38]  = 144;
  krad[39]  = 154;
  krad[40]  = 160;
  krad[41]  = 168;
  krad[42]  = 176;
  krad[43]  = 180;
  krad[44]  = 192;
  krad[45]  = 198;
  krad[46]  = 210;
  krad[47]  = 220;
  krad[48]  = 224;
  krad[49]  = 240;
  krad[50]  = 252;
  krad[51]  = 256;
  krad[52]  = 264;
  krad[53]  = 280;
  krad[54]  = 288;
  krad[55]  = 308;
  krad[56]  = 320;
  krad[57]  = 330;
  krad[58]  = 336;
  krad[59]  = 352;
  krad[60]  = 360;
  krad[61]  = 384;
  krad[62]  = 396;
  krad[63]  = 420;
  krad[64]  = 440;
  krad[65]  = 448;
  krad[66]  = 462;
  krad[67]  = 480;
  krad[68]  = 504;
  krad[69]  = 512;
  krad[70]  = 528;
  krad[71]  = 560;
  krad[72]  = 576;
  krad[73]  = 616;
  krad[74]  = 630;
  krad[75]  = 640;
  krad[76]  = 660;
  krad[77]  = 672;
  krad[78]  = 704;
  krad[79]  = 720;
  krad[80]  = 768;
  krad[81]  = 770;
  krad[82]  = 792;
  krad[83]  = 840;
  krad[84]  = 880;
  krad[85]  = 896;
  krad[86]  = 924;
  krad[87]  = 960;
  krad[88]  = 990;
  krad[89]  = 1008;
  krad[90]  = 1024;
  krad[91]  = 1056;
  krad[92]  = 1120;
  krad[93]  = 1152;
  krad[94]  = 1232;
  krad[95]  = 1260;
  krad[96]  = 1280;
  krad[97]  = 1320;
  krad[98]  = 1344;
  krad[99]  = 1386;
  krad[100] = 1408;
  krad[101] = 1440;
  krad[102] = 1536;
  krad[103] = 1540;
  krad[104] = 1584;
  krad[105] = 1680;
  krad[106] = 1760;
  krad[107] = 1792;
  krad[108] = 1848;
  krad[109] = 1920;
  krad[110] = 1980;
  krad[111] = 2016;
  krad[112] = 2048;
  krad[113] = 2112;
  krad[114] = 2240;
  krad[115] = 2304;
  krad[116] = 2310;
  krad[117] = 2464;
  krad[118] = 2520;
  krad[119] = 2560;
  krad[120] = 2640;
  krad[121] = 2688;
  krad[122] = 2772;
  krad[123] = 2816;
  krad[124] = 2880;
  krad[125] = 3072;
  krad[126] = 3080;
  krad[127] = 3168;
  krad[128] = 3360;
  krad[129] = 3520;
  krad[130] = 3584;
  krad[131] = 3696;
  krad[132] = 3840;
  krad[133] = 3960;
  krad[134] = 4032;
  krad[135] = 4096;
  krad[136] = 4224;
  krad[137] = 4480;
  krad[138] = 4608;
  krad[139] = 4620;
  krad[140] = 4928;
  krad[141] = 5040;
  krad[142] = 5120;
  krad[143] = 5280;
  krad[144] = 5376;
  krad[145] = 5544;
  krad[146] = 5632;
  krad[147] = 5760;
  krad[148] = 6144;
  krad[149] = 6160;
  krad[150] = 6336;
  krad[151] = 6720;
  krad[152] = 6930;
  krad[153] = 7040;
  krad[154] = 7168;
  krad[155] = 7392;
  krad[156] = 7680;
  krad[157] = 7920;
  krad[158] = 8064;
  krad[159] = 8192;
  krad[160] = 8448;
  krad[161] = 8960;
  krad[162] = 9216;
  krad[163] = 9240;
  krad[164] = 9856;
  krad[165] = 10080;
  krad[166] = 10240;
  krad[167] = 10560;
  krad[168] = 10752;
  krad[169] = 11088;
  krad[170] = 11264;
  krad[171] = 11520;
  krad[172] = 12288;
  krad[173] = 12320;
  krad[174] = 12672;
  krad[175] = 13440;
  krad[176] = 13860;
  krad[177] = 14080;
  krad[178] = 14336;
  krad[179] = 14784;

/*--------------------------------------------------------------------------*/
   }/*end routine */
/*==========================================================================*/


