#include "standard_include.h"
#include "ckcomplex.h"

#include "../../../include/debug_flags.h"
#include "../../../include/Atoms.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_mdatoms.h"

#include "../class_defs/CP_OPERATIONS/class_cplocal.h"

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void CPLOCAL::CP_hart_eext_calc(
           const int ncoef, complex *chunk,const int natm, Atom *atoms,
           complex *result, double *ehart_ret,double *eext_ret,double *ewd_ret,
           const int *k_x, const int *k_y, const int *k_z, int mydoublePack)

//============================================================================
// Function:  Hartree and External potentials
//
// NOTE FOR RAMKUMAR:  INVERSE BOX MATRIX (hmati) AND VOLUME (vol) 
//                     MUST BE PASSED IN AND ehart_ret AND eext_ret
//                     MUST BE SENT OUT.  FOR CUBIC SYSTEMS, HMATI 
//                     JUST 1/L ON ITS DIAGONAL, BUT ONE SHOULD ALLOW
//                     FOR A GENERAL 3x3 MATRIX (hmat) AND ITS INVERSE (hmati).
//                     I ALSO ASSUME result IS ZEROED SOMEWHERE SO THAT I
//                     CAN ACCUMULATE IT.
//                     FINALLY, THE ATOMIC COORDINATES (x,y,z) AND THEIR CHARGES (q)
//                     NEED TO BE PASSED IN
//============================================================================
{ /* Begin Function */
//----------------------------------------------------------------------------
// Local Variables and local pointers

  MDATOMS      *mdatoms      = MDATOMS::get();
  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();

#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_gen.h"

/*--------------------------------------------*/
/*         Local Pointer declarations         */

  /*------------------*/
  /* Atom information */
  int natm_piny     = mdclatoms_info->natm_tot; 
  double *q         = mdclatoms_info->q;
  int natm_typ      = mdatom_maps->natm_typ;
  int *iatm_typ     = mdatom_maps->iatm_atm_typ;

  /*--------------------------------*/
  /* Cell and pressure information  */
  int iperd             = gencell->iperd;
  double *hmat          = gencell->hmat;
  double *hmati         = gencell->hmati;
  double vol            = gencell->vol;
  double rvol           = 1.0/gencell->vol;

  /*----------------------------------------------*/
  /* Pseudo-potential and Reduced Periodicity info*/

  int nsplin_g         = cppseudo->nsplin_g;
  int n_rad_max        = cppseudo->n_rad_max;
  double dg_spl        = cppseudo->dg_spl;
  double gmin_spl      = cppseudo->gmin_spl;
  double *vps0         = cppseudo->vps0;
  double *vps1         = cppseudo->vps1;
  double *vps2         = cppseudo->vps2;
  double *vps3         = cppseudo->vps3;
  double *gzvps        = cppseudo->gzvps;
  int n_ang_max        = cppseudo->n_ang_max;
  int *loc_opt         = cppseudo->loc_opt;

  /*---------------------------------*/
  /* Ewald and ewald scr information */

  double alp_ewd    = genewald->alp_ewd;

  /*---------------------------------*/
  /* Local Pointers from leanCP      */

  // -----------------------------------------------------------
  // Local variables 
   double ehart=0.0,eext=0.0,EwdEnergy=0.0;
   double gx,gy,gz,g2,g;
   double HartreeFact,EwdFact;
   complex vext,sewd;
   double falp2;
   double wght,wght_now;
   complex s;
   double pi  = M_PI;
   double tpi = 2.0*M_PI;
   double fpi = 4.0*M_PI;


   double ecut4 = 8.0*cpcoeffs_info->ecut; // convert to Ryd extra factor of 2.0

   int izero;

//----------------------------------------------------------------------------
// Allocate temporary memory for break points

   int *index_atm  = new int[(natm+1)];
   double *vtemp   = new double[natm];
   complex *ei_inc = new complex[natm];
   complex *h      = new complex[natm];

//----------------------------------------------------------------------------
// Expressions for Hartree and external energies:
//
//  ehart = (1/vol) sum_{gx,gy,gz} (4*pi/g**2)|rho_g|**2  (excluding (0,0,0))
//
//  eext  = (1/vol) sum_{gx,gy,gz} (rho_g)^* S_g V_g
//
//         where S_g = sum_{I=1}^{natm} q_I exp[i(gx*x[I] + gy*y[I] + gz*z[I])]
//
//               V_g = exp(-g^2/4*alpha**2)/(g**2)

//----------------------------------------------------------------------------
// Some useful constants

   falp2 = 4.0*alp_ewd*alp_ewd;
   wght = 1.0;  if(mydoublePack==1){wght = 2.0;}

//----------------------------------------------------------------------------
// Initialize the results = vks(g) from Hartree and Exchange Correlation

#ifdef GJM_DEBUG_SIZE
   PRINTF(" %d : coefs in CP_hart_eext_calc\n",ncoef);
#endif
   memset(result,0,sizeof(complex)*ncoef);

//----------------------------------------------------------------------------
// Set up variables for break point calculations (helpful vectors!)

#ifdef GJM_DEBUG_SIZE
   if(ncoef>10){
     for(int i=0;i<10;i++){
       PRINTF("%d %d %d\n",k_x[i],k_y[i],k_z[i]);
     }
   }
#endif

   for(int iatm = 0; iatm < natm; iatm++){
      double arg_tmp = tpi*(hmati[1]*atoms[iatm].x + hmati[4]*atoms[iatm].y 
                          + hmati[7]*atoms[iatm].z);
      ei_inc[iatm] = complex(cos(arg_tmp),sin(arg_tmp));
   } /* endfor */


   for(int itype=1;itype<=natm_typ;itype++){
      index_atm[itype] =  (itype-1)*nsplin_g*(n_ang_max+1)*n_rad_max
                       +  loc_opt[itype]*nsplin_g*n_rad_max;
      //      CkPrintf("index_atm %d %d %d %d %d %d  %12.8lg \n",
      // itype,index_atm[itype],nsplin_g,n_ang_max,n_rad_max,loc_opt[itype],gzvps[itype]);

    }/*endfor*/
   //   CkPrintf("gmin_spl %lg dg_spl %lg \n",gmin_spl,dg_spl);


//----------------------------------------------------------------------------
// Begin loop over FFT grid.

   izero             = -10;
   double e0         = 0.0;
   double count      = 0.0;
   double count_slow = 0.0;
   int igo           = 0;
   int kx_old        = k_x[0];
   int ky_old        = k_y[0];
   int kz_old        = k_z[0];

#ifdef _CP_DEBUG_VKS_HART_EEXT_
   FILE *fp = fopen("vks_hart_eext_stuff.out","a+");
#endif

   for(int i = 0; i < ncoef; i++){/* Note that the (0,0,0) term is excluded! */

//----------------------------------------------------------------------------
// I.  Construct the reciprocal space vectors and their square magnitudes

     gx = tpi*(k_x[i]*hmati[1] + k_y[i]*hmati[2] + k_z[i]*hmati[3]);
     gy = tpi*(k_x[i]*hmati[4] + k_y[i]*hmati[5] + k_z[i]*hmati[6]);
     gz = tpi*(k_x[i]*hmati[7] + k_y[i]*hmati[8] + k_z[i]*hmati[9]);

     g2 = gx*gx + gy*gy + gz*gz;

     if(g2==0){izero=i;}
     if(g2 <= ecut4 && g2 != 0){

       wght_now = (k_x[i]==0 ? 1.0 : wght);
//----------------------------------------------------------------------------
// II. Use these to construct the Hartree energy and its potential

       HartreeFact = fpi/(g2*vol);
       ehart      += HartreeFact*chunk[i].getMagSqr()*wght_now;

#ifndef GLENN_DEBUG_ON
       result[i] += chunk[i]*HartreeFact;
#endif
       count+=1.0;

//----------------------------------------------------------------------------
// III. Get the structure factor:  If condition chosen to ensure we get in first time

       if(kx_old != k_x[i]-1 || ky_old != k_y[i] || kz_old != k_z[i] || igo==0) {
         for(int iatm = 0; iatm < natm; iatm++){
           double arg = atoms[iatm].x*gx + atoms[iatm].y*gy + atoms[iatm].z*gz;
           h[iatm] = complex(cos(arg),sin(arg));
         } /* endfor */
         count_slow+=1.0;
       }/* endif */

       g = sqrt(g2);
       CP_get_vpsnow(index_atm,nsplin_g,gmin_spl,dg_spl,g,
                 vps0,vps1,vps2,vps3,vtemp,iatm_typ,natm_typ,natm,vol); 

       sewd = complex(0.0,0.0);
       vext  = complex(0.0,0.0);
       for(int iatm = 0; iatm < natm; iatm++){
         vext  += h[iatm]*vtemp[iatm];
         sewd  += h[iatm]*atoms[iatm].q;
       }//endfor

#ifdef CHECK_STRUCT_FACT
       CkPrintf("%d %d %d %.12g %.12g\n",k_x[i],k_y[i],k_z[i],s.re,s.im);
       if(k_x[i] == 1 && k_y[i] == 0 && k_z[i] == 0){
 	 CkPrintf("%.12g %.12g %.12g\n",gx,gy,gz);
         for(int iatm = 0; iatm < natm; iatm++){
  	  CkPrintf("%d %.12g %.12g %.12g\n",
             iatm,atoms[iatm].x,atoms[iatm].y,atoms[iatm].z);
          CkPrintf("%d %.12g\n",iatm,atoms[iatm].q);
         }
       }
#endif
        
//----------------------------------------------------------------------------
// IV. Use structure factor to evaluate external potential

       EwdFact    = fpi*exp(-g2/falp2)/(g2*vol);

       EwdEnergy += EwdFact*sewd.getMagSqr()*wght_now;

       double sumr    = sewd.re*EwdFact*wght_now;
       double sumi    = sewd.im*EwdFact*wght_now;
       double chunk_r = chunk[i].re*wght_now;
       double chunk_i = chunk[i].im*wght_now;

       for(int iatm=0; iatm < natm; iatm++){
         double chunk_temp_r = (chunk_r*vtemp[iatm]
                             +  sumr*atoms[iatm].q);
         double chunk_temp_i = (-chunk_i*vtemp[iatm]
                             +   sumi*atoms[iatm].q);
         double srx = (gx*chunk_temp_r);
         double sry = (gy*chunk_temp_r);
         double srz = (gz*chunk_temp_r);
         double six = (gx*chunk_temp_i);
         double siy = (gy*chunk_temp_i);
         double siz = (gz*chunk_temp_i);

         atoms[iatm].fx += (srx*h[iatm].im - six*h[iatm].re);
         atoms[iatm].fy += (sry*h[iatm].im - siy*h[iatm].re);
         atoms[iatm].fz += (srz*h[iatm].im - siz*h[iatm].re); 
       }/*endfor*/

#ifndef GLENN_DEBUG_ON
       result[i] += vext.conj(); // minus (PINY CONV) but plus (THIS CODE CONV)
#endif

       eext += (chunk[i]*vext).re*wght_now;

       for(int iatm=0; iatm < natm; iatm++){h[iatm] = h[iatm]*ei_inc[iatm];}

#ifdef _CP_DEBUG_VKS_HART_EEXT_
       if(k_x[i] == 0 && k_y[i] == 1 && k_z[i] == 4){
         fprintf(fp,"0 1 4 : %g %g : %g %g : %g %g\n",
          chunk[i].re,chunk[i].im,vext.re,vext.im,result[i].re,result[i].im);
       }//endif
       if(k_x[i] == 0 && k_y[i] == 4 && k_z[i] == 1){
         fprintf(fp,"0 4 1 : %g %g : %g %g : %g %g\n",
          chunk[i].re,chunk[i].im,vext.re,vext.im,result[i].re,result[i].im);
       }//endif
       if(k_x[i] == 1 && k_y[i] == 0 && k_z[i] == 4){
         fprintf(fp,"1 0 4 : %g %g : %g %g : %g %g\n",
          chunk[i].re,chunk[i].im,vext.re,vext.im,result[i].re,result[i].im);
       }//endif
       if(k_x[i] == 1 && k_y[i] == 4 && k_z[i] == 0){
         fprintf(fp,"1 4 0 : %g %g : %g %g : %g %g\n",
          chunk[i].re,chunk[i].im,vext.re,vext.im,result[i].re,result[i].im);
       }//endif
       if(k_x[i] == 4 && k_y[i] == 1 && k_z[i] == 0){
         fprintf(fp,"4 1 0 : %g %g : %g %g : %g %g\n",
          chunk[i].re,chunk[i].im,vext.re,vext.im,result[i].re,result[i].im);
       }//endif
       if(k_x[i] == 4 && k_y[i] == 0 && k_z[i] == 1){
         fprintf(fp,"4 0 1 : %g %g : %g %g : %g %g\n",
          chunk[i].re,chunk[i].im,vext.re,vext.im,result[i].re,result[i].im);
       }//endif
#endif
       kx_old = k_x[i];
       ky_old = k_y[i];
       kz_old = k_z[i];
       igo    = 1;
     }else{
       igo    = 0;
     }// endif : ecut

   }/* endfor */

//----------------------------------------------------------------------------
// Deal with g=0 : double pack weight is 1 just like non-double pack

   if(izero>=0){
     vext.re = 0.0;
     sewd.re = 0.0;
     for(int iatm=0;iatm< natm; iatm++){
        vtemp[iatm] = gzvps[iatm_typ[(iatm+1)]]/vol;
        sewd.re += atoms[iatm].q;   
        vext.re += vtemp[iatm];
     }/*endfor*/
     int i = izero;
#ifndef GLENN_DEBUG_ON
      result[i].re += vext.re;
#endif
      e0 = vext.re*chunk[i].re;
      eext += e0;

   }//endif

#ifdef GJM_DEBUG_SIZE
   PRINTF("Hart eext : tot %g slow %g\n",count,count_slow);
#endif

//----------------------------------------------------------------------------
// VI. Return values
//FACTOR 2 FIX LATER WITH LESS COMPUTATION

   *ehart_ret = ehart/2.0;
   *eext_ret  = eext;
   *ewd_ret   = EwdEnergy/2.0;


//----------------------------------------------------------------------------
// VII. Free temporary memory

   delete [] index_atm;
   delete [] vtemp;
   delete [] ei_inc;
   delete [] h;

#ifdef _CP_DEBUG_VKS_HART_EEXT_
   fclose(fp);
#endif

   
//============================================================================
} /* End function */
//============================================================================

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void CPLOCAL::CP_get_vpsnow(int *index_atm,int nsplin_g,
               double gmin_spl,double  dg_spl,double g,
               double *vps0,double *vps1,double *vps2,double *vps3,
               double *vtemp,int *iatm_typ,int natm_typ,int npart,double vol)

/*========================================================================*/
 {/*begin routine*/
/*========================================================================*/
/*        Local variable declarations            */

  int ipart,itype,iii;
  int index_now;
  double vtemp_atyp[200];
  double h,h0;
  double partem1,partem2,partem3;
  double partem4;

/*==========================================================================*/
/* Loop over atom types to calculate pseudopotential                        */


  for(itype=1;itype<=natm_typ;itype++){
    iii = (int)((g-gmin_spl)/dg_spl + 1);
    iii = MIN(iii,nsplin_g);
    iii = MAX(iii,1);
    h0  = (double)(iii-1)*dg_spl+gmin_spl;
    h   = g-h0;
    index_now = index_atm[itype] + iii;
    partem1   = vps0[index_now];
    partem2   = vps1[index_now];
    partem3   = vps2[index_now];
    partem4   = vps3[index_now];

    vtemp_atyp[itype] = (((partem4*h+partem3)*h+partem2)*h + partem1)/vol;
  }/*endfor*/

  for(ipart=1;ipart<=npart;ipart++){
    vtemp[(ipart-1)] = vtemp_atyp[iatm_typ[ipart]];  // changed vtemp to start at 0
  }/*endfor*/

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/
