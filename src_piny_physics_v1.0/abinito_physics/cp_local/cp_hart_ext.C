#include "standard_include.h"
#include "ckcomplex.h"

#include "../../../include/debug_flags.h"
#include "../../../include/Atoms.h"
#include "../../../include/eesDataClass.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_mdatoms.h"

#include "../class_defs/CP_OPERATIONS/class_cplocal.h"
#include "../proto_defs/proto_cp_ewald_local.h"

//#define _CP_DEBUG_VKS_HART_EEXT_

//============================================================================
//  N^2 method  : Invoked from gchare
//============================================================================
//  N^2 method  : Invoked from gchare
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPLOCAL::CP_hart_eext_calc(int ncoef, complex *rho,int natm, Atom *atoms,
                                complex *vks, double *ehart_ret,double *eext_ret,
                                double *ewd_ret,int *k_x, int *k_y, int *k_z, 
                                int index)
//============================================================================
// Function:  Hartree and External potentials
//
// NOTE FOR RAMKUMAR:  INVERSE BOX MATRIX (hmati) AND VOLUME (vol) 
//                     MUST BE PASSED IN AND ehart_ret AND eext_ret
//                     MUST BE SENT OUT.  FOR CUBIC SYSTEMS, HMATI 
//                     JUST 1/L ON ITS DIAGONAL, BUT ONE SHOULD ALLOW
//                     FOR A GENERAL 3x3 MATRIX (hmat) AND ITS INVERSE (hmati).
//                     I ALSO ASSUME vks IS ZEROED SOMEWHERE SO THAT I
//                     CAN ACCUMULATE IT.
//                     FINALLY, THE ATOMIC COORDINATES (x,y,z) AND THEIR CHARGES (q)
//                     NEED TO BE PASSED IN
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
//
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
  int iperd         = gencell->iperd;
  double *hmat      = gencell->hmat;
  double *hmati     = gencell->hmati;
  double vol        = gencell->vol;
  double rvol       = 1.0/gencell->vol;

  /*----------------------------------------------*/
  /* Pseudo-potential and Ewald stuff             */

  int nsplin_g      = cppseudo->nsplin_g;
  int n_rad_max     = cppseudo->n_rad_max;
  double dg_spl     = cppseudo->dg_spl;
  double gmin_spl   = cppseudo->gmin_spl;
  double *vps0      = cppseudo->vps0;
  double *vps1      = cppseudo->vps1;
  double *vps2      = cppseudo->vps2;
  double *vps3      = cppseudo->vps3;
  double *gzvps     = cppseudo->gzvps;
  int n_ang_max     = cppseudo->n_ang_max;
  int *loc_opt      = cppseudo->loc_opt;
  double ecut4      = 8.0*cpcoeffs_info->ecut; // convert to Ryd extra factor of 2.0
  double alp_ewd    = genewald->alp_ewd;
  double pi         = M_PI;
  double tpi        = 2.0*M_PI;
  double fpi        = 4.0*M_PI;

// -----------------------------------------------------------
// Local variables 

   double gx,gy,gz,g2,g;
   double HartreeFact,EwdFact;
   complex s;
   complex vext,sewd;

//============================================================================
// Debug output

#ifdef _CP_DEBUG_VKS_HART_EEXT_
   char myFileName[100];
   sprintf(myFileName, "Vext_Gspace_%d.out", index);
   FILE *fp = fopen(myFileName,"w");
   sprintf(myFileName, "Ewd_Gspace_%d.out", index);
   FILE *fp2 = fopen(myFileName,"w");
   double *fxt = new double [natm];
   double *fyt = new double [natm];
   double *fzt = new double [natm];
   double *fxt_ewd = new double [natm];
   double *fyt_ewd = new double [natm];
   double *fzt_ewd = new double [natm];
   for(int i =0;i<natm;i++){
     fxt[i] = 0.0;
     fyt[i] = 0.0;
     fzt[i] = 0.0;
     fxt_ewd[i] = 0.0;
     fyt_ewd[i] = 0.0;
     fzt_ewd[i] = 0.0;
   }//endforp
#endif

#ifdef GJM_DEBUG_SIZE
   PRINTF(" %d : coefs in CP_hart_eext_calc\n",ncoef);
#endif

//============================================================================
// Initialize vks(g) from Hartree and Exchange Correlation
//                   and  useful constants

   double falp2     = 4.0*alp_ewd*alp_ewd;
   double wght_now  = 1.0;
   double wght      = 2.0;
   double ehart     = 0.0;
   double eext      = 0.0;
   double EwdEnergy = 0.0;
   memset(vks,0,sizeof(complex)*ncoef);

//============================================================================
// Set up variables for break point calculations (helpful vectors!)

   int *index_atm  = new int[(natm+1)];
   double *vtemp   = new double[natm];
   complex *ei_inc = new complex[natm];
   complex *h      = new complex[natm];
   int izero         = -10;
   int igo           = 0;
   double count      = 0.0;
   double count_slow = 0.0;
   int kx_old        = k_x[0];
   int ky_old        = k_y[0];
   int kz_old        = k_z[0];

   // kx moves fastest through memory : see CP_Rho_GSpacePlane::computeK
   for(int iatm = 0; iatm < natm; iatm++){
      double arg_tmp = tpi*(hmati[3]*atoms[iatm].x + hmati[6]*atoms[iatm].y 
                          + hmati[9]*atoms[iatm].z);
      ei_inc[iatm] = complex(cos(arg_tmp),sin(arg_tmp));
   }/* endfor */

   for(int itype=1;itype<=natm_typ;itype++){
      index_atm[itype] =  (itype-1)*nsplin_g*(n_ang_max+1)*n_rad_max
                       +  loc_opt[itype]*nsplin_g*n_rad_max;
   }/*endfor*/

//============================================================================
// Begin loop over FFT grid.

   for(int i = 0; i < ncoef; i++){/* Note that the (0,0,0) term is excluded! */

#ifdef CMK_VERSION_BLUEGENE
      int nfreq_cmi_update=4; 
      if((i+1)%nfreq_cmi_update==0){CmiNetworkProgress();}
#endif

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
       ehart      += HartreeFact*rho[i].getMagSqr()*wght_now;

       vks[i] += rho[i]*HartreeFact;
       count+=1.0;

  //----------------------------------------------------------------------------
  // III. Get the structure factor:  If condition chosen to ensure we get in first time

       if(kx_old != k_x[i] || ky_old != k_y[i] || kz_old != k_z[i]-1 || igo==0) {
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
       double vsave = 0;
       for(int iatm = 0; iatm < natm; iatm++){
#ifdef _CP_DEBUG_VKS_HART_EEXT_
         if(iatm_typ[(iatm+1)]==2){
           vsave  = vtemp[iatm];
           vext  += h[iatm]*vtemp[iatm];
         }//endif
#else
         vext  += h[iatm]*vtemp[iatm];
#endif
         sewd  += h[iatm]*atoms[iatm].q;
       }//endfor

  //----------------------------------------------------------------------------
  // IV. Use structure factor to evaluate external potential

       EwdFact    = fpi*exp(-g2/falp2)/(g2*vol);
       EwdEnergy += EwdFact*sewd.getMagSqr()*wght_now;

       double sumr   = sewd.re*EwdFact*wght_now;
       double sumi   = sewd.im*EwdFact*wght_now;
       double rho_r  = rho[i].re*wght_now;
       double rho_i  = rho[i].im*wght_now;
       for(int iatm=0; iatm < natm; iatm++){
         double rho_temp_r = ( rho_r*vtemp[iatm] + sumr*atoms[iatm].q);
         double rho_temp_i = (-rho_i*vtemp[iatm] + sumi*atoms[iatm].q);
         double srx        = (gx*rho_temp_r);
         double sry        = (gy*rho_temp_r);
         double srz        = (gz*rho_temp_r);
         double six        = (gx*rho_temp_i);
         double siy        = (gy*rho_temp_i);
         double siz        = (gz*rho_temp_i);
         atoms[iatm].fx   += (srx*h[iatm].im - six*h[iatm].re);
         atoms[iatm].fy   += (sry*h[iatm].im - siy*h[iatm].re);
         atoms[iatm].fz   += (srz*h[iatm].im - siz*h[iatm].re); 
#ifdef _CP_DEBUG_VKS_HART_EEXT_
        if(iatm_typ[(iatm+1)]==2){
 	   rho_temp_r = ( rho_r*vtemp[iatm] );
   	   rho_temp_i = (-rho_i*vtemp[iatm] );
           srx = (gx*rho_temp_r);
           sry = (gy*rho_temp_r);
           srz = (gz*rho_temp_r);
           six = (gx*rho_temp_i);
           siy = (gy*rho_temp_i);
           siz = (gz*rho_temp_i);
           fxt[iatm] +=  (srx*h[iatm].im - six*h[iatm].re);
           fyt[iatm] +=  (sry*h[iatm].im - siy*h[iatm].re);
           fzt[iatm] +=  (srz*h[iatm].im - siz*h[iatm].re);
        }//endif
        rho_temp_r = sumr*atoms[iatm].q;
        rho_temp_i = sumi*atoms[iatm].q;
        srx = (gx*rho_temp_r);
        sry = (gy*rho_temp_r);
        srz = (gz*rho_temp_r);
        six = (gx*rho_temp_i);
        siy = (gy*rho_temp_i);
        siz = (gz*rho_temp_i);
        fxt_ewd[iatm] +=  (srx*h[iatm].im - six*h[iatm].re);
        fyt_ewd[iatm] +=  (sry*h[iatm].im - siy*h[iatm].re);
        fzt_ewd[iatm] +=  (srz*h[iatm].im - siz*h[iatm].re);
#endif
       }/*endfor*/

       vks[i] += vext.conj(); 
       eext += (rho[i]*vext).re*wght_now;
       for(int iatm=0; iatm < natm; iatm++){h[iatm] = h[iatm]*ei_inc[iatm];}

#ifdef _CP_DEBUG_VKS_HART_EEXT_
       fprintf(fp,"%d %d %d : %g %g : %g %g : %g %g : %g %g : %g\n",k_x[i],k_y[i],k_z[i],
	       rho[i].re,rho[i].im,vext.re,vext.im,vks[i].re,vks[i].im,eext,ehart,vsave);
       fprintf(fp2,"%d %d %d : %g %g : %g %g\n",k_x[i],k_y[i],k_z[i],
  	          sewd.re,sewd.im,sewd.getMagSqr(),EwdFact);
#endif
       kx_old = k_x[i];
       ky_old = k_y[i];
       kz_old = k_z[i];
       igo    = 1;
     }else{
       igo    = 0;
     }// endif : ecut

   }/* endfor : icoef*/

//============================================================================
// Deal with g=0 : double pack weight is 1 just like non-double pack

   if(izero>=0){
     vext.re = 0.0;   vext.im = 0.0;
     sewd.re = 0.0;   sewd.im = 0.0;
     for(int iatm=0;iatm< natm; iatm++){
        sewd.re += atoms[iatm].q;   
        vtemp[iatm] = gzvps[iatm_typ[(iatm+1)]]/vol;
#ifdef _CP_DEBUG_VKS_HART_EEXT_
        if(iatm_typ[(iatm+1)]==2){
          vext.re += vtemp[iatm];
	}//endif
#else
        vext.re += vtemp[iatm];
#endif
     }/*endfor*/
     int i = izero;
     vks[i].re += vext.re;
     vks[i].im  = 0.0;
     eext += (vext.re*rho[i].re);
#ifdef _CP_DEBUG_VKS_HART_EEXT_
     fprintf(fp,"0 0 0 : %g %g : %g %g : %g %g\n",
		rho[i].re,rho[i].im,vext.re,vext.im,vks[i].re,vks[i].im);
#endif

   }//endif

//============================================================================
// VI. Return values and memory clean up

   *ehart_ret = ehart/2.0;
   *eext_ret  = eext;
   *ewd_ret   = EwdEnergy/2.0;
   delete [] index_atm;
   delete [] vtemp;
   delete [] ei_inc;
   delete [] h;

//============================================================================
// VIII. Debug output Control

#ifdef GJM_DEBUG_SIZE
   PRINTF("Hart eext : tot %g slow %g\n",count,count_slow);
#endif

#ifdef _CP_DEBUG_VKS_HART_EEXT_
   fclose(fp);
   fclose(fp2);

   sprintf(myFileName, "fatm_%d.out",index);
   fp = fopen(myFileName,"w");
   for(int i=0;i<natm;i++){
     fprintf(fp,"%d %g %g %g\n",i,fxt[i],fyt[i],fzt[i]);
   }//endfor
   fclose(fp);   
   delete [] fxt;
   delete [] fyt;
   delete [] fzt;
   sprintf(myFileName, "fatm_%d.out.ewd",index);
   fp = fopen(myFileName,"w");
   for(int i=0;i<natm;i++){
     fprintf(fp,"%d %g %g %g\n",i,fxt_ewd[i],fyt_ewd[i],fzt_ewd[i]);
   }//endfor
   fclose(fp);   
   delete [] fxt_ewd;
   delete [] fyt_ewd;
   delete [] fzt_ewd;
#endif
   
//============================================================================
   }/* End function */
//============================================================================


//============================================================================
// Invoked local to class
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
  double vtemp_atyp[201];
  double h,h0;
  double partem1,partem2,partem3;
  double partem4;
  if(natm_typ>200){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("vemp_atyp hard coded too small in cp_hart_ext.C\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout); EXIT(1);
  }//endif

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
  }//endfor

  for(ipart=1;ipart<=npart;ipart++){
    vtemp[(ipart-1)] = vtemp_atyp[iatm_typ[ipart]];  // changed vtemp to start at 0
  }//endfor

//--------------------------------------------------------------------------
 }//end routine
//==========================================================================


//==========================================================================
// Fetch local pseudo ees parameters
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void CPLOCAL::getEesPrms(int *ngrid_a, int *ngrid_b, int *ngrid_c,
                            int *n_interp, int *natm)
//==========================================================================
  { // begin routine 
//==========================================================================
  CP           *cp           = CP::get();
  MDATOMS      *mdatoms      = MDATOMS::get();
#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_mdatoms.h"

  ngrid_a[0]  = cppseudo->ngrid_eext_a;
  ngrid_b[0]  = cppseudo->ngrid_eext_b;
  ngrid_c[0]  = cppseudo->ngrid_eext_c;
  n_interp[0] = cppseudo->n_interp_ps;
  natm[0]     = mdclatoms_info->natm_tot;

  }//end routine
//==========================================================================


//==========================================================================
// Compute the g-space weight for the k-vectors given.
//--------------------------------------------------------------------------
// The weight only depends on k-vectors given.
// The EESgroup should call the routine ONCE for each collection it is assigned.
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void CPLOCAL::eesSetEesWghtGgrp(int ncoef, int *ka_in, int *kb_in, int *kc_in,
                                double *b_re, double *b_im, 
                                int nkf1,int nkf2,int nkf3,int n_interp)
//==========================================================================
  { // begin routine 
  CP           *cp           = CP::get();
#include "../class_defs/allclass_strip_cp.h"
  int nkf1_t = cppseudo->ngrid_eext_a;
  int nkf2_t = cppseudo->ngrid_eext_b;
  int nkf3_t = cppseudo->ngrid_eext_c;

  if(nkf1_t!=nkf1 || nkf2_t!=nkf2 || nkf3_t!=nkf3){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Incorrect Eext FFT size %d %d %d vs %d %d %d",
            nkf1,nkf2,nkf3,nkf1_t,nkf2_t,nkf3_t);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout); EXIT(1);
  }//endif

//==========================================================================
// I) Calculate bweight for this set of points

//--------------------------------------------------------------------------
//   A) Malloc memory and define constants                                

   int dim_k;
   int ngrid_a = nkf1;
   int ngrid_b = nkf2;
   int ngrid_c = nkf3;

   dim_k = ngrid_a;
   double *bden_a_r  = (double *)cmalloc(dim_k*sizeof(double),"set_ees_wgt")-1;
   double *bden_a_i  = (double *)cmalloc(dim_k*sizeof(double),"set_ees_wgt")-1;
   double *bweight_a = (double *)cmalloc(dim_k*sizeof(double),"set_ees_wgt")-1;

   dim_k = ngrid_b;
   double *bden_b_r  = (double *)cmalloc(dim_k*sizeof(double),"set_ees_wgt")-1;
   double *bden_b_i  = (double *)cmalloc(dim_k*sizeof(double),"set_ees_wgt")-1;
   double *bweight_b = (double *)cmalloc(dim_k*sizeof(double),"set_ees_wgt")-1;

   dim_k = ngrid_c;
   double *bden_c_r  = (double *)cmalloc(dim_k*sizeof(double),"set_ees_wgt")-1;
   double *bden_c_i  = (double *)cmalloc(dim_k*sizeof(double),"set_ees_wgt")-1;
   double *bweight_c = (double *)cmalloc(dim_k*sizeof(double),"set_ees_wgt")-1;

   dim_k = n_interp;
   double *aj   = (double *)cmalloc(dim_k*sizeof(double),"set_ees_wgt")-1;
   double *rn   = (double *)cmalloc(dim_k*sizeof(double),"set_ees_wgt")-1;
   double *rn1  = (double *)cmalloc(dim_k*sizeof(double),"set_ees_wgt")-1;
   double *mn_k = (double *)cmalloc(dim_k*sizeof(double),"set_ees_wgt")-1;
   double *uk   = (double *)cmalloc(dim_k*sizeof(double),"set_ees_wgt")-1;

//--------------------------------------------------------------------------
//   B) Construct the weighting Function                                  

   get_bspline_wght1d(n_interp,ngrid_a,aj,rn,rn1,mn_k,uk,
                      bden_a_r,bden_a_i,bweight_a);
   get_bspline_wght1d(n_interp,ngrid_b,aj,rn,rn1,mn_k,uk,
                      bden_b_r,bden_b_i,bweight_b);
   get_bspline_wght1d(n_interp,ngrid_c,aj,rn,rn1,mn_k,uk,
                      bden_c_r,bden_c_i,bweight_c);

   for(int i=1;i <=ncoef; ++i){
     int ka = ka_in[(i-1)];
     int kb = kb_in[(i-1)];
     int kc = kc_in[(i-1)];
     int kap,kbp,kcp;
     if (kc < 0) {kcp = kc + nkf3 + 1;} else {kcp = kc + 1;}
     if (kb < 0) {kbp = kb + nkf2 + 1;} else {kbp = kb + 1;}
     if (ka < 0) {kap = ka + nkf1 + 1;} else {kap = ka + 1;}
     double tmp_a_r,tmp_a_i;
     tmp_a_r     = bden_a_r[kap];
     tmp_a_i     = bden_a_i[kap];
     double tmp_b_r,tmp_b_i;
     tmp_b_r     = bden_b_r[kbp];
     tmp_b_i     = bden_b_i[kbp];
     double tmp_c_r,tmp_c_i;
     tmp_c_r     = bden_c_r[kcp];
     tmp_c_i     = bden_c_i[kcp];
     double tmp_ab_r,tmp_ab_i;
     tmp_ab_r    = tmp_a_r*tmp_b_r - tmp_a_i*tmp_b_i;
     tmp_ab_i    = tmp_a_i*tmp_b_r + tmp_a_r*tmp_b_i;
     b_re[(i-1)] = tmp_ab_r*tmp_c_r - tmp_ab_i*tmp_c_i;
     b_im[(i-1)] = tmp_ab_i*tmp_c_r + tmp_ab_r*tmp_c_i;
     if(ka==0&&kb==0&&kc==0){ //safety
       b_re[(i-1)] = bden_a_r[1]*bden_b_r[1]*bden_c_r[1];
       b_im[(i-1)] = 0.0;
    }//endfor
  }//endif

//==========================================================================

  cfree(&bden_a_r[1],"set_ees_wgt"); 
  cfree(&bden_a_i[1],"set_ees_wgt"); 
  cfree(&bweight_a[1],"set_ees_wgt");
  cfree(&bden_b_r[1],"set_ees_wgt"); 
  cfree(&bden_b_i[1],"set_ees_wgt"); 
  cfree(&bweight_b[1],"set_ees_wgt");
  cfree(&bden_c_r[1],"set_ees_wgt"); 
  cfree(&bden_c_i[1],"set_ees_wgt"); 
  cfree(&bweight_c[1],"set_ees_wgt");
  cfree(&aj[1],"set_ees_wgt"); 
  cfree(&rn[1],"set_ees_wgt"); 
  cfree(&rn1[1],"set_ees_wgt"); 
  cfree(&uk[1],"set_ees_wgt"); 
  cfree(&mn_k[1],"set_ees_wgt"); 

//--------------------------------------------------------------------------
  }//end routine
//==========================================================================


//==========================================================================
// At every time step, EesGroup with allowed planes Eext planes calls this routine
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void CPLOCAL::eesAtmBsplineRgrp(Atom *atoms, int *allowed_planes, 
                                RHORHARTDATA *RhoRHartData)
//==========================================================================
  {// begin routine 
//==========================================================================

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();
  MDATOMS      *mdatoms      = MDATOMS::get();

#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_mdatoms.h"
  PSNONLOCAL *non_local = &(cppseudo->nonlocal);

  int i,j,ja,jb,jc,n;
  int j1,j2,jj,ia,ib,ic,kkk;
  double grid_a,grid_b,grid_c;
  double atemp,btemp,ctemp;
  double mn_a_tmp,mn_b_tmp,mn_c_tmp;
  double x,y,z;

  double *hmati  = gencell->hmati;
  int natm       = mdclatoms_info->natm_tot;

  int n_interp   = cppseudo->n_interp_ps;
  int ngrid_a    = cppseudo->ngrid_eext_a;
  int ngrid_b    = cppseudo->ngrid_eext_b;
  int ngrid_c    = cppseudo->ngrid_eext_c;

  double *aj     = non_local->aj;
  double *rn     = non_local->rn;
  double *rn1    = non_local->rn1;
  int *iatemp    = non_local->iatemp;
  int *ibtemp    = non_local->ibtemp;
  int *ictemp    = non_local->ictemp;
  double *frac_a = non_local->frac_a;
  double *frac_b = non_local->frac_b;
  double *frac_c = non_local->frac_c;
  int **igrid_a  = non_local->igrid_a;
  int **igrid_b  = non_local->igrid_b;
  int **igrid_c  = non_local->igrid_c;
  double **mn_a  = non_local->mn_a;
  double **mn_b  = non_local->mn_b;
  double **mn_c  = non_local->mn_c;
  double **ua  = non_local->ua;
  double **ub  = non_local->ub;
  double **uc  = non_local->uc;
  double **dmn_a = non_local->dmn_a;
  double **dmn_b = non_local->dmn_b;
  double **dmn_c = non_local->dmn_c;
  double cpu1,cpu2;

//==========================================================================
// 0) Useful Constants 

   for(j=1;j<=n_interp;j++){
     aj[j] = (double) (j-1);
     rn[j] = (double) (j);
     if(j > 1){rn1[j] = 1.0/((double)(j-1));}
   }//endfor
   rn1[1] = 0.0;

   grid_a    = (double) ngrid_a;
   grid_b    = (double) ngrid_b;
   grid_c    = (double) ngrid_c;

//==========================================================================
// I) scaled coordinates                                                   

#ifdef TIME_BSPLINE
   cputime(&cpu1);
#endif

   for(i=0;i<natm;i++){
     x = atoms[i].x;
     y = atoms[i].y;
     z = atoms[i].z;
     atemp = x*hmati[1] + y*hmati[4] + z*hmati[7];
     btemp = x*hmati[2] + y*hmati[5] + z*hmati[8];
     ctemp = x*hmati[3] + y*hmati[6] + z*hmati[9];
     atemp -= NINT((atemp-0.5));
     btemp -= NINT((btemp-0.5));
     ctemp -= NINT((ctemp-0.5));
     atemp = (atemp==1.0 ? 0.0 : atemp);
     btemp = (btemp==1.0 ? 0.0 : btemp);
     ctemp = (ctemp==1.0 ? 0.0 : ctemp);
     atemp = atemp*grid_a;
     btemp = btemp*grid_b;
     ctemp = ctemp*grid_c;
     iatemp[i] = (int) (atemp);
     ibtemp[i] = (int) (btemp);
     ictemp[i] = (int) (ctemp);
     frac_a[i] = atemp - (double) (iatemp[i]);
     frac_b[i] = btemp - (double) (ibtemp[i]);
     frac_c[i] = ctemp - (double) (ictemp[i]);
   }//endfor
#ifdef CMK_VERSION_BLUEGENE
    CmiNetworkProgress();
#endif

//==========================================================================
// II) Using current fraction, find the grid points on which M_n is non-zero 

   for(j=1;j<=n_interp;j++){
     for(i=0;i<natm;i++){
       ua[j][i] = frac_a[i] + aj[j];
       ub[j][i] = frac_b[i] + aj[j];
       uc[j][i] = frac_c[i] + aj[j];
       j2       = j-2;
       ia       = iatemp[i] - j2;
       ib       = ibtemp[i] - j2;
       ic       = ictemp[i] - j2;
       ia       = (ia>0 ? ia:ngrid_a+ia);
       ib       = (ib>0 ? ib:ngrid_b+ib);
       ic       = (ic>0 ? ic:ngrid_c+ic);
       ia       = (ia<=ngrid_a ? ia:ia-ngrid_a);
       ib       = (ib<=ngrid_b ? ib:ib-ngrid_b);
       ic       = (ic<=ngrid_c ? ic:ic-ngrid_c);
       igrid_a[j][i] = ia-1;
       igrid_b[j][i] = (ib - 1)*(ngrid_a+2);
       igrid_c[j][i] = (ic - 1);  // use to assign to planes
     }//endfor
#ifdef CMK_VERSION_BLUEGENE
    CmiNetworkProgress();
#endif
   }//endfor
  

//==========================================================================
// III) Initialize M2 and get the Mn's using the recursion relation         
//    Note: M_n is defined on 0 to n. Since frac between 0 and 1, the       
//          calculation is performed in an order that takes advantage of it 

   for(i=0;i<natm;i++){
     mn_a[1][i] = 1.0 - fabs(ua[1][i]-1.0);
     mn_b[1][i] = 1.0 - fabs(ub[1][i]-1.0);
     mn_c[1][i] = 1.0 - fabs(uc[1][i]-1.0);
     mn_a[2][i] = 1.0 - fabs(ua[2][i]-1.0);
     mn_b[2][i] = 1.0 - fabs(ub[2][i]-1.0);
     mn_c[2][i] = 1.0 - fabs(uc[2][i]-1.0);
   }//endfor
#ifdef CMK_VERSION_BLUEGENE
   CmiNetworkProgress();
#endif
   for(j=3;j<=n_interp;j++){
     for(i=0;i<natm;i++){
       mn_a[j][i]   = 0.0;
       mn_b[j][i]   = 0.0;
       mn_c[j][i]   = 0.0;
     }//endfor
   }//endfor
#ifdef CMK_VERSION_BLUEGENE
   CmiNetworkProgress();
#endif
   for(n=3;n<=n_interp;n++){

     for(j=n;j>=2;j--){
       j1 = j-1;
       for(i=0;i<natm;i++){
         mn_a_tmp  = (ua[j][i]*mn_a[j][i]+(rn[n]-ua[j][i])*mn_a[j1][i])*rn1[n];
         mn_b_tmp  = (ub[j][i]*mn_b[j][i]+(rn[n]-ub[j][i])*mn_b[j1][i])*rn1[n];
         mn_c_tmp  = (uc[j][i]*mn_c[j][i]+(rn[n]-uc[j][i])*mn_c[j1][i])*rn1[n];
         mn_a[j][i] = mn_a_tmp;
         mn_b[j][i] = mn_b_tmp;
         mn_c[j][i] = mn_c_tmp;
       }//end for: i
#ifdef CMK_VERSION_BLUEGENE
    CmiNetworkProgress();
#endif
     }//end for: j
     for(i=0;i<natm;i++){
       mn_a[1][i] = ua[1][i]*mn_a[1][i]*rn1[n];
       mn_b[1][i] = ub[1][i]*mn_b[1][i]*rn1[n];
       mn_c[1][i] = uc[1][i]*mn_c[1][i]*rn1[n];
     }//endfor 
#ifdef CMK_VERSION_BLUEGENE
    CmiNetworkProgress();
#endif

     if(n==(n_interp-1)){
       for(i=0;i<natm;i++){
         dmn_a[1][i] = mn_a[1][i];
         dmn_b[1][i] = mn_b[1][i];
         dmn_c[1][i] = mn_c[1][i];
       }//endfor
       for(j=2;j<=n_interp;j++){    
         j1 = j-1;
         for(i=0;i<natm;i++){
           dmn_a[j][i] = mn_a[j][i] - mn_a[j1][i];
           dmn_b[j][i] = mn_b[j][i] - mn_b[j1][i];
           dmn_c[j][i] = mn_c[j][i] - mn_c[j1][i];
         }//endfor: i
#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif
       }//endfor : j
     }//endif : get derivative

   }//end for: n

//==========================================================================
// IV) put together the separable bits and store 

   // Zero array mapping plane index, j, to interpolation index, jc.
   int    *plane_index;
   int    **igrid; 
   double **mn,**dmn_x,**dmn_y,**dmn_z;

   // Zero array mapping plane index, j, to interpolation index, jc.
   for(j=0;j<ngrid_c;j++){
     if(allowed_planes[j]==1){
       plane_index = RhoRHartData[j].plane_index;
       for(i=0;i<natm;i++){plane_index[i] = 0;}
     }//endif
   }//endfor

   for(i=0;i<natm;i++){
   for(jc=1;jc<=n_interp;jc++){
     int ip = igrid_c[jc][i];
     if(allowed_planes[ip]==1){ // if the group wants this plane
       plane_index = RhoRHartData[ip].plane_index;
       igrid       = RhoRHartData[ip].igrid;
       mn          = RhoRHartData[ip].mn;
       dmn_x       = RhoRHartData[ip].dmn_x;
       dmn_y       = RhoRHartData[ip].dmn_y;
       dmn_z       = RhoRHartData[ip].dmn_z;
       jj = 1;
       for(jb=1;jb<=n_interp;jb++){
         for(ja=1,j=jj;ja<=n_interp;ja++,j++){
           igrid[i][j] = igrid_a[ja][i]+igrid_b[jb][i]; //plane index only
           atemp           = dmn_a[ja][i]* mn_b[jb][i]* mn_c[jc][i]*grid_a;
           btemp           =  mn_a[ja][i]*dmn_b[jb][i]* mn_c[jc][i]*grid_b;
           ctemp           =  mn_a[ja][i]* mn_b[jb][i]*dmn_c[jc][i]*grid_c;
           mn[i][j]    =  mn_a[ja][i]* mn_b[jb][i]* mn_c[jc][i];
           dmn_x[i][j] = atemp*hmati[1]+btemp*hmati[2]+ctemp*hmati[3];
           dmn_y[i][j] = atemp*hmati[4]+btemp*hmati[5]+ctemp*hmati[6];
           dmn_z[i][j] = atemp*hmati[7]+btemp*hmati[8]+ctemp*hmati[9];
         }//endfor : ja
         jj += n_interp;
       }//endfor : jb
       plane_index[i] = jc; // each jc is a different plane
#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif
     }//endif : allowed
   }}//endfor : iatm and jc


#ifdef DEBUG_GJM_BSPLINE
   for(i=0;i<natm;i++){
    printf("frac[%d] : %g %g %g\n",i,frac_a[i],frac_b[i],frac_c[i]);
    for(j=1;j<=n_interp;j++){
      printf("mn[%d][%d] : %g %g %g\n",j,i,mn_a[j][i],mn_b[j][i],mn_c[j][i]);
    }//endfor
    printf("\n"); scanf("%d",&kkk);
   }//endfor
#endif

#ifdef TIME_BSPLINE
   cputime(&cpu2);
   printf("Bspline timing %g\n",cpu2-cpu1);
#endif

//--------------------------------------------------------------------------
  }//end routine
//==========================================================================



//==========================================================================
// Use the B-spline coefs to generate sfAtmTypR (real space SF for spec. type)
//     Invoke by RchareEesSF which controls type loop
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void CPLOCAL::eesPackGridRchare(int natm, int ityp, double *sfAtmTypR, int iplane,
                                int **igrid, double **mn, int *plane_index)
//==========================================================================
  {// begin routine 
//==========================================================================

  MDATOMS      *mdatoms      = MDATOMS::get();
  CP           *cp           = CP::get();
#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_mdatoms.h"

  int ngrid_a   = cppseudo->ngrid_eext_a;
  int ngrid_b   = cppseudo->ngrid_eext_b;
  int n_interp  = cppseudo->n_interp_ps;
  int n_interp2 = n_interp*n_interp;
  int *iatm_typ = mdatom_maps->iatm_atm_typ;

//==========================================================================

  for(int i=0;i<(ngrid_a+2)*ngrid_b;i++){sfAtmTypR[i]=0.0;}

  for(int iatm=0;iatm<natm;iatm++){
    if(iatm_typ[(iatm+1)]==ityp){
      int jc = plane_index[iatm];  // interpolation pt of plane
      if(jc>0){
        for(int j=1;j<=n_interp2;j++){
          int ind         = igrid[iatm][j]; // index of pt in the plane
          sfAtmTypR[ind] += mn[iatm][j]; // contribute to zmatrix
        }//endfor
#ifdef CMK_VERSION_BLUEGENE
        CmiNetworkProgress();
#endif
      }//endif
    }//endif
  }//endfor

//==========================================================================
// sfatmtypr is ready to be ffted and used to generate eext

//--------------------------------------------------------------------------
  }//end routine
//==========================================================================


//==========================================================================
// Invoked from g-chare : SfatmTypG is result of fft of sfatmtypr
//                        sfAtomTotG is charged weighted SF that is accumulated
//                        during the `ityp' control loop along with vks.
//                        Generate eext and hartree(ityp==1 only)
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void CPLOCAL::eesHartEextGchare(int ncoef, int ityp, complex *rho, complex *vks, 
                                complex *sfAtmTypG, complex *sfAtmTotG,
                                double *b_re, double *b_im,
                                double *ehart_ret,double *eext_ret,
                                int *k_x,int *k_y,int *k_z,int index)
//==========================================================================
  {// begin routine 
//==========================================================================

  MDATOMS      *mdatoms      = MDATOMS::get();
  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();

#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_mdatoms.h"
#include "../class_defs/allclass_strip_gen.h"

  //---------------------------------------------
  // Atom look up information : we are already sorted by atm type
  int natm_typ_now     = 1;
  int natm_now         = 1;
  int iatm_typ_now[2];
  int index_atm[2];
  double vtemp[2];

  //---------------------------------------------
  // Cell and pressure information  
  int iperd             = gencell->iperd;
  double *hmat          = gencell->hmat;
  double *hmati         = gencell->hmati;
  double vol            = gencell->vol;
  double rvol           = 1.0/gencell->vol;

  //----------------------------------------------
  // Local pseudo-potential information
  int nsplin_g         = cppseudo->nsplin_g;
  int n_rad_max        = cppseudo->n_rad_max;
  int n_ang_max        = cppseudo->n_ang_max;
  double dg_spl        = cppseudo->dg_spl;
  double gmin_spl      = cppseudo->gmin_spl;
  int *loc_opt         = cppseudo->loc_opt;
  double *vps0         = cppseudo->vps0;
  double *vps1         = cppseudo->vps1;
  double *vps2         = cppseudo->vps2;
  double *vps3         = cppseudo->vps3;
  double *gzvps        = cppseudo->gzvps;
  double *q_typ        = cppseudo->q_typ;
  double ecut4         = 8.0*cpcoeffs_info->ecut; // in Ryd; extra factor of 2.0
  double pi            = M_PI;
  double tpi           = 2.0*M_PI;
  double fpi           = 4.0*M_PI;
  double wght          = 2.0;

  // -------------------------------------------
  // Local variables 

  double gx,gy,gz,g2,g;
  double HartreeFact;
  double wght_now;
  double temp_r,temp_i;
  complex vext;

//==========================================================================
// DEBUG output

#ifdef _CP_DEBUG_VKS_HART_EEXT_
   printf("Qtyp[%d] %g\n",ityp,q_typ[ityp]);
   char myFileName[100];
   sprintf(myFileName, "Vext_Gspace_%d.out.ees.%d", index,ityp);
   FILE *fp = fopen(myFileName,"w");
#endif

//==========================================================================
// Initialize atom information

  index_atm[1]    =  (ityp-1)*nsplin_g*(n_ang_max+1)*n_rad_max
                  +  loc_opt[ityp]*nsplin_g*n_rad_max;
  iatm_typ_now[1] = 1;

//==========================================================================
// Loop over g-space and do it up!

   double ehart = 0.0;
   double eext  = 0.0;
   int izero    = -10;  // deal with g=0

   for(int i = 0; i < ncoef; i++){

   //----------------------------------------------------------------------------
#ifdef CMK_VERSION_BLUEGENE
     int nfreq_cmi_update = 100; 
     if( ((i+1)%nfreq_cmi_update)==0 ){CmiNetworkProgress();}
#endif
   //----------------------------------------------------------------------------
   // I.  Construct the reciprocal space vectors and their square magnitudes
     gx = tpi*(k_x[i]*hmati[1] + k_y[i]*hmati[2] + k_z[i]*hmati[3]);
     gy = tpi*(k_x[i]*hmati[4] + k_y[i]*hmati[5] + k_z[i]*hmati[6]);
     gz = tpi*(k_x[i]*hmati[7] + k_y[i]*hmati[8] + k_z[i]*hmati[9]);
     g2 = gx*gx + gy*gy + gz*gz;
     g  = sqrt(g2);
     if(g2==0){izero=i;}
     if(g2 <= ecut4 && g2 != 0){
       wght_now = (k_x[i]==0 ? 1.0 : wght);
   //----------------------------------------------------------------------------
   // II. Use these to construct the Hartree energy and its potential
       if(ityp==1){
          HartreeFact = fpi/(g2*vol);
          ehart      += HartreeFact*rho[i].getMagSqr()*wght_now;
          vks[i]     += rho[i]*HartreeFact;
       }/*endif*/
   //----------------------------------------------------------------------------
   // III. Get the structure factor and pseudopotential interaction
      //------------------------------------------------------
      // a)Atm SF : apply ees g-space weight : add to the total SF
       temp_r           = (sfAtmTypG[i].re*b_re[i]-sfAtmTypG[i].im*b_im[i]);
       temp_i           = (sfAtmTypG[i].re*b_im[i]+sfAtmTypG[i].im*b_re[i]);
       sfAtmTotG[i].re += (sfAtmTypG[i].re*q_typ[ityp]); // bweight added later
       sfAtmTotG[i].im += (sfAtmTypG[i].im*q_typ[ityp]); // bweight added later
      //------------------------------------------------------
      // b)Atm-electron interaction
       CP_get_vpsnow(index_atm,nsplin_g,gmin_spl,dg_spl,g,vps0,vps1,vps2,vps3,
                     vtemp,iatm_typ_now,natm_typ_now,natm_now,vol); 
       vext.re = temp_r*vtemp[0];
       vext.im = temp_i*vtemp[0];
       eext   += (rho[i]*vext).re*wght_now;
       //------------------------------------------------------
       // c)electron force
       vks[i] += vext.conj();
#ifdef _CP_DEBUG_VKS_HART_EEXT_
       fprintf(fp,"%d %d %d : %g %g : %g %g : %g %g : %g %g : %g\n",
               k_x[i],k_y[i],k_z[i],
	       rho[i].re,rho[i].im,vext.re,vext.im,vks[i].re,vks[i].im,eext,ehart,vtemp[0]);
#endif
       //------------------------------------------------------
       // d)atom force : set up backtransform
       temp_r          =  (rho[i].re*b_re[i]-rho[i].im*b_im[i]);
       temp_i          = -(rho[i].re*b_im[i]+rho[i].im*b_re[i]);
       sfAtmTypG[i].re = vtemp[0]*temp_r;
       sfAtmTypG[i].im = vtemp[0]*temp_i;
     }else{
       if(g2!=0.0){
         sfAtmTypG[i].re = 0.0;
         sfAtmTypG[i].im = 0.0;
         sfAtmTotG[i].re = 0.0;
         sfAtmTotG[i].im = 0.0;
         vks[i].re       = 0.0;
         vks[i].im       = 0.0;
       }//endif
     }// endif : ecut

   }// endfor : coefs

//=======================================================================
// Deal with g=0 : double pack weight is 1 just like non-double pack

   if(izero>=0){
     int i = izero;
    //------------------------------------------------------
    // Atm SF : apply ees gspace weight : add to the total SF
     temp_r           = sfAtmTypG[i].re*b_re[i];
     sfAtmTotG[i].re += (sfAtmTypG[i].re*q_typ[ityp]);
     sfAtmTotG[i].im  = 0.0;
    //------------------------------------------------------
    // Atm-electron interaction
     vtemp[1]     = gzvps[ityp]/vol;
     vext.re      = vtemp[1]*temp_r;
     vext.im      = 0.0;
     eext        += (rho[i].re*vext.re);
    //------------------------------------------------------
    // electron force
     vks[i].re   += vext.re;
     vks[i].im    = 0.0;
    //------------------------------------------------------
    // atom force : set up backtransform
     temp_r          = rho[i].re*b_re[i];
     sfAtmTypG[i].re = vtemp[1]*temp_r;
     sfAtmTypG[i].im = 0.0;
#ifdef _CP_DEBUG_VKS_HART_EEXT_
       fprintf(fp,"0 0 0 : %g %g : %g %g : %g %g : %g %g \n",
	       rho[i].re,rho[i].im,vext.re,vext.im,vks[i].re,vks[i].im,eext,ehart);
#endif
   }//endif

//==========================================================================
// Set return values

   if(ityp==1){ehart_ret[0] = ehart/2.0;}
   eext_ret[0] += eext;

//==========================================================================
// Debug output

#ifdef _CP_DEBUG_VKS_HART_EEXT_
   fclose(fp);
#endif

//==========================================================================
// sfAtmTypG is now inversed FFT to get atom forces
// if itype loop is done, sfAtmTotG is complete and Ewald energy is computed.

//--------------------------------------------------------------------------
  }//end routine
//==========================================================================


//==========================================================================
// Invoked from g-chare : sfAtomTotG is charged weighted SF e.g. summed over 
//                        all atm types. It is used to generate ewald energy.
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void CPLOCAL::eesEwaldGchare(int ncoef, complex *sfAtmTotG,
                             double *b_re, double *b_im,double *ewd_ret,
                             int *k_x, int *k_y, int *k_z,int index)
//==========================================================================
  {// begin routine 
//==========================================================================
  GENERAL_DATA *general_data = GENERAL_DATA::get();

#include "../class_defs/allclass_strip_gen.h"

 // Cell information
  double *hmati   = gencell->hmati;
  double vol      = gencell->vol;

 // Ewald and g-space constants
  double alp_ewd  = genewald->alp_ewd;
  double tpi      = 2.0*M_PI;
  double falp2    = 4.0*alp_ewd*alp_ewd;
  double fpi      = 4.0*M_PI;
  double wght_now = 2.0;
  double wght     = 2.0;

 // Local variables
  double aka,akb,akc;
  double xk,yk,zk,g2,preg,smag,bmag,prep;
  double temp_r,temp_i;

//==========================================================================
// Debug output

#ifdef _CP_DEBUG_VKS_HART_EEXT_
   char myFileName[100];
   sprintf(myFileName, "Ewd_Gspace_%d.out.ees", index);
   FILE *fp = fopen(myFileName,"w");
#endif

//==========================================================================
// Do the Ewald sum

  double vnow = 0.0;

  for(int i=0;i<ncoef;i++) {
   //---------------------------------------
   // Make Sameer Happy with some progress
#ifdef CMK_VERSION_BLUEGENE
    int nfreq_cmi_update = 100; 
    if( ((i+1)%nfreq_cmi_update)==0 ){CmiNetworkProgress();}
#endif
   //---------------------------------------
   // Compute kvectors form recip lattice puppies
    aka = tpi*( (double) k_x[i] );
    akb = tpi*( (double) k_y[i] );
    akc = tpi*( (double) k_z[i] );
    xk  = (aka*hmati[1] + akb*hmati[2] + akc*hmati[3]);
    yk  = (aka*hmati[4] + akb*hmati[5] + akc*hmati[6]);
    zk  = (aka*hmati[7] + akb*hmati[8] + akc*hmati[9]);
    g2  = xk*xk + yk*yk + zk*zk;
    if(g2!=0.0){
     //---------------------------------------
     // Compute Ewald energy
      wght_now = (k_x[i]==0 ? 1.0 : wght);
      preg     = fpi*exp(-g2/falp2)/(g2*vol);
      smag     = sfAtmTotG[i].getMagSqr();
      bmag     = (b_re[i]*b_re[i]+b_im[i]*b_im[i]);
      vnow    += (bmag*smag*preg*wght_now);
      temp_r   = (sfAtmTotG[i].re*b_re[i]-sfAtmTotG[i].im*b_im[i]);
      temp_i   = (sfAtmTotG[i].re*b_im[i]+sfAtmTotG[i].im*b_re[i]);
#ifdef _CP_DEBUG_VKS_HART_EEXT_
      fprintf(fp,"%d %d %d : %g %g : %g %g : %g %g\n",k_x[i],k_y[i],k_z[i],
  	          temp_r,temp_i,b_re[i],b_im[i],sfAtmTotG[i].getMagSqr()*bmag,preg);
#endif
     //---------------------------------------
     // Ewald force : Setup inverse FFT
      prep     = preg*bmag;      // there is 1/2 below canceling 2x
      sfAtmTotG[i].re *= prep;
      sfAtmTotG[i].im *= prep;
    }else{
      sfAtmTotG[i].re = 0.0;
      sfAtmTotG[i].im = 0.0;
    }//endif
  }//endfor

//==========================================================================
// Energy return values

  ewd_ret[0] = vnow/2.0;

//==========================================================================
// Debug me baby

#ifdef _CP_DEBUG_VKS_HART_EEXT_
  fclose(fp);
#endif

//==========================================================================
// sfAtmTotG is ready to be FFTed and used to generate atom forces

//==========================================================================
  }//end routine
//==========================================================================


//==========================================================================
// Invoked from R-chare : SfatmtypR is back FFT of SFatmTypG (flag==0)
//                        or SFatmTotG (flag==1)
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void CPLOCAL::eesAtmForceRchare(int natm, Atom *atoms,int ityp, 
                int **igrid, double **dmn_x, double **dmn_y, double **dmn_z, 
                int *plane_index, double *sfAtmTypR, int iplane,int flag)
//==========================================================================
  {// begin routine 
//==========================================================================

  MDATOMS      *mdatoms      = MDATOMS::get();
  CP           *cp           = CP::get();
#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_mdatoms.h"

  int n_interp  = cppseudo->n_interp_ps;
  int n_interp2 = n_interp*n_interp;
  int *iatm_typ = mdatom_maps->iatm_atm_typ;
  double *q     = mdclatoms_info->q;

//==========================================================================
// Setup some debug output

#ifdef _CP_DEBUG_VKS_HART_EEXT_
   double *fxt = new double [natm];
   double *fyt = new double [natm];
   double *fzt = new double [natm];
   for(int i =0;i<natm;i++){
     fxt[i] = 0.0;
     fyt[i] = 0.0;
     fzt[i] = 0.0;
   }//endfor
#endif

//==========================================================================
// Compute the Eext forces from the back transformed ityp SF   : flag=0
// Compute the Ewald forces from the back transformed total SF : flag=1

  for(int iatm=0;iatm<natm;iatm++){
    if(iatm_typ[(iatm+1)]==ityp || flag==1){
      int jc = plane_index[iatm];  // interpolation pt of plane
      if(jc>0){
        double qnow = (flag==1 ? q[(iatm+1)] : 1.0);
        for(int j=1;j<=n_interp2;j++){
          int ind         = igrid[iatm][j]; // index of pt in the plane
          double p        = sfAtmTypR[ind]*qnow;
          atoms[iatm].fx -= (dmn_x[iatm][j]*p);
          atoms[iatm].fy -= (dmn_y[iatm][j]*p);
          atoms[iatm].fz -= (dmn_z[iatm][j]*p);
#ifdef _CP_DEBUG_VKS_HART_EEXT_
          fxt[iatm] -=(dmn_x[iatm][j]*p);
          fyt[iatm] -=(dmn_y[iatm][j]*p);
          fzt[iatm] -=(dmn_z[iatm][j]*p);
#endif
        }//endfor
#ifdef CMK_VERSION_BLUEGENE
        CmiNetworkProgress();
#endif
      }//endif
    }//endif
  }//endfor

//==========================================================================
// Debug output

#ifdef _CP_DEBUG_VKS_HART_EEXT_
   char myFileName[100];
   if(flag==0){
     sprintf(myFileName, "fatm_%d.out.ees.%d",iplane,ityp);
   }else{
     sprintf(myFileName, "fatm_%d.out.ees.ewd",iplane);
   }//endif
   FILE *fp = fopen(myFileName,"w");
     for(int i=0;i<natm;i++){
       fprintf(fp,"%d %g %g %g\n",i,fxt[i],fyt[i],fzt[i]);
     }//endfor
   fclose(fp);   
   delete []fxt;
   delete []fyt;
   delete []fzt;
#endif

//==========================================================================
// Done with the present atom type or if flag==1 we're totally done!

//==========================================================================
  }//end routine
//==========================================================================
