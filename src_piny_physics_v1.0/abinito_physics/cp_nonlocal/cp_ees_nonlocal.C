#include "standard_include.h"
#include "ckcomplex.h"

#include "../../../include/Atoms.h"    
#include "../../../include/eesDataClass.h"    
#include "../../../src_mathlib/mathlib.h"    
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_gen.h"

#include "../class_defs/CP_OPERATIONS/class_cpnonlocal.h"
#include "../proto_defs/proto_cp_ewald_local.h"

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void CPNONLOCAL::getEesPrms(int *ngrid_a, int *ngrid_b, int *ngrid_c,
                            int *n_interp, int *natm)
//==========================================================================
  { // begin routine 
//==========================================================================
  CP           *cp           = CP::get();
#include "../class_defs/allclass_strip_cp.h"
  PSNONLOCAL *nonlocal = &(cppseudo->nonlocal);

  ngrid_a[0]  = nonlocal->ngrid_a;
  ngrid_b[0]  = nonlocal->ngrid_b;
  ngrid_c[0]  = nonlocal->ngrid_c;
  n_interp[0] = nonlocal->n_interp;
  natm[0]     = nonlocal->natm;

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
void CPNONLOCAL::eesSetEesWghtGgrp(int ncoef, int *ka_in, int *kb_in, int *kc_in,
                                   double *b_re, double *b_im, 
                                   int nkf1,int nkf2,int nkf3,int n_interp)
//==========================================================================
  { // begin routine 
//==========================================================================
  CP           *cp           = CP::get();
#include "../class_defs/allclass_strip_cp.h"
  int nkf1_t = cppseudo->nonlocal.ngrid_a;
  int nkf2_t = cppseudo->nonlocal.ngrid_b;
  int nkf3_t = cppseudo->nonlocal.ngrid_c;

  if(nkf1_t!=nkf1 || nkf2_t!=nkf2 || nkf3_t!=nkf3){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Incorrect Nonlocal FFT size %d %d %d vs %d %d %d",
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
    }//endif
  }//endfor

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
// Compute the g-space spline coefs for the k-vectors given.
//--------------------------------------------------------------------------
// The spline coefs only depends on k-vectors given.
// The EESgroup should call the routine ONCE for each collection it is assigned.
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void CPNONLOCAL::eesSplProjectorGgrp(int ncoef, int *ka, int *kb, int *kc,
                                     double *h_gspl,int *ind_gspl)
//==========================================================================
  { // begin routine 
//==========================================================================

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();

#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_gen.h"
  int i,iii;
  double aka,akb,akc;
  double xk,yk,zk;
  double g2,g,h0,h;
  double tpi = 2.0*M_PI;

  double *hmati   = gencell->hmati;
  int nsplin_g    = cppseudo->nsplin_g;
  double gmin_spl = cppseudo->gmin_spl;
  double dg_spl   = cppseudo->dg_spl;

//==========================================================================
// Spline look up stuff : never changes. Done starting with 0 convention

  for(i=0;i<ncoef;i++){
    aka = ((double)(ka[i]));
    akb = ((double)(kb[i]));
    akc = ((double)(kc[i]));
    xk  = (aka*hmati[1]+akb*hmati[2]+akc*hmati[3])*tpi;
    yk  = (aka*hmati[4]+akb*hmati[5]+akc*hmati[6])*tpi;
    zk  = (aka*hmati[7]+akb*hmati[8]+akc*hmati[9])*tpi;
    g2  = xk*xk+yk*yk+zk*zk;
    g   = sqrt(g2);
    iii = (g-gmin_spl)/dg_spl + 1;  //look up array convention in piny, starts at 1
    iii = MIN(iii,nsplin_g);
    iii = MAX(iii,1);
    h0  = ((double)(iii-1))*dg_spl+gmin_spl;
    h   = g-h0;
    h_gspl[i]   = h;
    ind_gspl[i] = iii;
    if(g==0.0){
      h_gspl[i]   = 0.0;
      ind_gspl[i] = 1;
    }//endif
  }//endfor

//--------------------------------------------------------------------------
  }//end routine
//==========================================================================


//==========================================================================
// Compute in real space, the B-Spline coefficients of the atoms
//--------------------------------------------------------------------------
// Should be invoked by EESgroup members with real space planes at each time step:
// EESGroup should provide a list of ALL allowed real space planes.
//    allowed_planes[ip] = 1, if I want the Bspline coefs on plane ip
//    allowed_planes[ip] = 0, if I DONT want the Bspline coefs on plane ip
// where ip=0 ... nfftc or ngrid_c. Taking allowed_planes[ip]=1 for all ip
// decreases efficiency a tiny bit. The atoms are a group, too, so this 
// group based approach is totally cool.
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void CPNONLOCAL::eesAtmBsplineRgrp(Atom *atoms, int *allowed_planes, RPPDATA *RPPData)
//==========================================================================
  {// begin routine 
//==========================================================================

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();

#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_gen.h"
  PSNONLOCAL *nonlocal = &(cppseudo->nonlocal);

  int i,j,ja,jb,jc,n;
  int j1,j2,jj,ia,ib,ic,kkk;
  double grid_a,grid_b,grid_c;
  double atemp,btemp,ctemp;
  double mn_a_tmp,mn_b_tmp,mn_c_tmp;
  double x,y,z;

  double *hmati  = gencell->hmati;
  int natm       = nonlocal->natm;
  int n_interp   = nonlocal->n_interp;
  int ngrid_a    = nonlocal->ngrid_a;
  int ngrid_b    = nonlocal->ngrid_b;
  int ngrid_c    = nonlocal->ngrid_c;
  int *map_nl    = nonlocal->map_nl;

  double *aj     = nonlocal->aj;
  double *rn     = nonlocal->rn;
  double *rn1    = nonlocal->rn1;
  int *iatemp    = nonlocal->iatemp;
  int *ibtemp    = nonlocal->ibtemp;
  int *ictemp    = nonlocal->ictemp;
  double *frac_a = nonlocal->frac_a;
  double *frac_b = nonlocal->frac_b;
  double *frac_c = nonlocal->frac_c;
  int **igrid_a  = nonlocal->igrid_a;
  int **igrid_b  = nonlocal->igrid_b;
  int **igrid_c  = nonlocal->igrid_c;
  double **mn_a  = nonlocal->mn_a;
  double **mn_b  = nonlocal->mn_b;
  double **mn_c  = nonlocal->mn_c;
  double **ua  = nonlocal->ua;
  double **ub  = nonlocal->ub;
  double **uc  = nonlocal->uc;
  double **dmn_a = nonlocal->dmn_a;
  double **dmn_b = nonlocal->dmn_b;
  double **dmn_c = nonlocal->dmn_c;

  double cpu1,cpu2;

//==========================================================================
// 0) Useful Constants 

   for(j=1;j<=n_interp;j++){
     aj[j] = (double) (j-1);
     rn[j] = (double) (j);
     if(j > 1){rn1[j] = 1.0/((double)(j-1));}
   }//endfor
   rn1[1] = 0.0;

   grid_a   = (double) ngrid_a;
   grid_b   = (double) ngrid_b;
   grid_c   = (double) ngrid_c;

//==========================================================================
// I) scaled coordinates                                                   

#ifdef TIME_BSPLINE
   cputime(&cpu1);
#endif

   for(i=0;i<natm;i++){
     int k = map_nl[(i+1)]-1;
     x = atoms[k].x;
     y = atoms[k].y;
     z = atoms[k].z;
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
   }//endfor
#ifdef CMK_VERSION_BLUEGENE
        CmiNetworkProgress();
#endif
  
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

   int    *plane_index;
   int    **igrid; 
   double **mn,**dmn_x,**dmn_y,**dmn_z;

   // Zero array mapping plane index, j, to interpolation index, jc.
   for(j=0;j<ngrid_c;j++){
     if(allowed_planes[j]==1){
       plane_index = RPPData[j].plane_index;
       for(i=0;i<natm;i++){plane_index[i] = 0;}
     }//endif
   }//endfor

#ifdef CMK_VERSION_BLUEGENE
        CmiNetworkProgress();
#endif

   for(i=0;i<natm;i++){
   for(jc=1;jc<=n_interp;jc++){
     int ip = igrid_c[jc][i];
     if(allowed_planes[ip]==1){ // if the group wants this plane
       plane_index = RPPData[ip].plane_index;
       igrid       = RPPData[ip].igrid;
       mn          = RPPData[ip].mn;
       dmn_x       = RPPData[ip].dmn_x;
       dmn_y       = RPPData[ip].dmn_y;
       dmn_z       = RPPData[ip].dmn_z;
       jj = 1;
       for(jb=1;jb<=n_interp;jb++){
         for(ja=1,j=jj;ja<=n_interp;ja++,j++){
           igrid[i][j] = igrid_a[ja][i]+igrid_b[jb][i]; //plane index only
           atemp       = dmn_a[ja][i]* mn_b[jb][i]* mn_c[jc][i]*grid_a;
           btemp       =  mn_a[ja][i]*dmn_b[jb][i]* mn_c[jc][i]*grid_b;
           ctemp       =  mn_a[ja][i]* mn_b[jb][i]*dmn_c[jc][i]*grid_c;
           mn[i][j]    =  mn_a[ja][i]* mn_b[jb][i]* mn_c[jc][i];
           dmn_x[i][j] = atemp*hmati[1]+btemp*hmati[2]+ctemp*hmati[3];
           dmn_y[i][j] = atemp*hmati[4]+btemp*hmati[5]+ctemp*hmati[6];
           dmn_z[i][j] = atemp*hmati[7]+btemp*hmati[8]+ctemp*hmati[9];
         }//endfor : ja
         jj += n_interp;
       }//endfor : jb
       plane_index[i] = jc; // Each jc is a different plane.
                            // Ex : plane_index[3] = 4
                            //      For the 3rd atom, 
                            //      the 4th c-interpolation pt
                            //      is on plane number ip=25
#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif
     }//endif : allowed
   }}//endfor : iatm and jc

#ifdef DEBUG_GJM_BSPLINE
   for(i=0;i<natm;i++){
   PRINTF("frac[%d] : %g %g %g\n",i,frac_a[i],frac_b[i],frac_c[i]);
   for(j=1;j<=n_interp;j++){
     PRINTF("mn[%d][%d] : %g %g %g\n",j,i,mn_a[j][i],mn_b[j][i],mn_c[j][i]);
   }
   PRINTF("\n"); scanf("%d",&kkk);
   }
#endif

#ifdef TIME_BSPLINE
   cputime(&cpu2);
   PRINTF("Bspline timing %g\n",cpu2-cpu1);
#endif

//--------------------------------------------------------------------------
  }//end routine
//==========================================================================



//==========================================================================
// Every g-chare must invoke this puppy for at every time step for every
// iteration of the non-local loop (iter_nl).
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void CPNONLOCAL::eesProjGchare(int ncoef, complex *psi,int *ka,int *kb, int *kc,
                     int ihave_g0, int ind_g0, int iter_nl,
                     double *d_re, double *d_im, double *dyp_re, double *dyp_im,
                     complex *projPsiG, int *ind_gspl, double *h_gspl)
//==========================================================================
   {//begin routine
//==========================================================================
  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();

#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_gen.h"
  PSNONLOCAL *nonlocal = &(cppseudo->nonlocal);

  double *hmati   = gencell->hmati;
  int nsplin_g    = cppseudo->nsplin_g;
  int n_ang_max   =  cppseudo->n_ang_max;

  int *natm_typ_lang  = cppseudo->nonlocal.natm_typ_lang;
  int *natm_lang      = cppseudo->nonlocal.natm_lang;
  int *iatm_str_lang  = cppseudo->nonlocal.iatm_str_lang;
  int **iatm_typ_lang = cppseudo->nonlocal.iatm_typ_lang;
  int *lang_v         = cppseudo->nonlocal.lang_v;
  int *mang_v         = cppseudo->nonlocal.mang_v;
  int *ityp_v         = cppseudo->nonlocal.ityp_v;

  double *vps0   = cppseudo->vps0;
  double *vps1   = cppseudo->vps1;
  double *vps2   = cppseudo->vps2;
  double *vps3   = cppseudo->vps3;
  double *gzvps0 =  cppseudo->gzvps0;

  int lang     = lang_v[iter_nl];  // piny indices
  int mang     = mang_v[iter_nl];
  int ityp     = ityp_v[iter_nl];
  int lang1    = lang+1;

  int iatm_typ = iatm_typ_lang[ityp][lang1]; // true atom type      
  int natm     = natm_lang[iatm_typ];        // # atms of this type 
  int iatm_str = iatm_str_lang[iatm_typ];    // where atms begin    

  int nfreq=100;

//==========================================================================
// The dyp_re,dyp_im,psi, must be memory on each chare!! It is reused below!! 
// Different states can have different iter_nl leading to different dyp's hence
// they cannot be cached without synchonizing.
// The d_re,d_im are from the cache and only depend on g-plane index.

  eesYlmOnD(lang,mang,ncoef,ka,kb,kc,dyp_re,dyp_im,d_re,d_im,hmati);

  int ind_off = (ityp-1)*nsplin_g*(n_ang_max+1) + lang*nsplin_g;
  for(int i=0;i<ncoef;i++){
    int ind_now = ind_off + ind_gspl[i];
    double h    = h_gspl[i];
    double v0   = vps0[ind_now];
    double v1   = vps1[ind_now];
    double v2   = vps2[ind_now];
    double v3   = vps3[ind_now];
    double vnow = ((v3*h+v2)*h+v1)*h+v0;
    dyp_re[i]  *= vnow;
    dyp_im[i]  *= vnow;
#ifdef CMK_VERSION_BLUEGENE
    if(i%nfreq==0){CmiNetworkProgress();}
#endif
  }//endfor

  int i = ind_g0;
  if(lang==0 && ihave_g0==1){
    double vnow = gzvps0[ityp];
    dyp_re[i] *= vnow;
    dyp_im[i]  = 0.0;
  }//endif
  if(lang!=0 && ihave_g0==1){
    dyp_re[i] = 0.0;
    dyp_im[i] = 0.0;
  }//endif

//==========================================================================
// Using your psi and your dyp create the puppy to be 3DFFTed to real space

  // congj(psi * struct_fact)
  for(int ig=0;ig<ncoef;ig++){ 
     projPsiG[ig].re =  dyp_re[ig]*psi[ig].re-dyp_im[ig]*psi[ig].im;
     projPsiG[ig].im =-(dyp_im[ig]*psi[ig].re+dyp_re[ig]*psi[ig].im);
#ifdef CMK_VERSION_BLUEGENE
    if(ig%nfreq==0){CmiNetworkProgress();}
#endif
  }//endfor
#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif

//-------------------------------------------------------------------------
   }//end routine
//==========================================================================


//==========================================================================
// Only invoked by class function eesProjGchare
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void CPNONLOCAL::eesYlmOnD(int lang,int mang,int ncoef,int *ka,int *kb,int *kc,
                     double *dy_re,double *dy_im, double *d_re, double *d_im,
                     double *hmati)
//==========================================================================
  {// begin routine 
//==========================================================================
  CP           *cp           = CP::get();
#include "../class_defs/allclass_strip_cp.h"

  int i;
  int m_ang_abs;
  double sgn_m_ang;
  double y00;
  double y10;
  double y11;
  double y20;
  double y21;
  double y22;
  double y30;
  double y31;
  double y32;
  double y33;
  double aka,akb,akc;
  double xk,yk,zk,g2,g,gs;
  double ctheta,stheta,cphi,sphi;

  double tpi        = 2.0*M_PI;
  double rt_fpi     = cpylm_cons->rt_fpi;
  double rt_thrfpi  = cpylm_cons->rt_thrfpi;
  double rt_threpi  = cpylm_cons->rt_threpi;

  int nfreq = 100;

//==========================================================================
// Spherical harmonic times the ees g-space weight, d                       
//   Note : since E = |S|^2 for imaginary Ylm = Yl*exp(im*phi)              
//          can use cos(m*phi) projector and sin(m*phi) projector           
//          as two real projectors                                          
//          | \int cos + i \int sin|^2 = |\int cos|^2 + |\int sin |^2       
//          when the states are real as they are here.                      

  switch(lang){  
    //-----------------------------------------------------------
    case 0:
      for(i=0;i<ncoef;i++){
        y00      = rt_fpi;
        dy_re[i] = d_re[i]*y00;
        dy_im[i] = d_im[i]*y00;
#ifdef CMK_VERSION_BLUEGENE
        if(i%nfreq==0){CmiNetworkProgress();}
#endif
      }//endfor
    break;
    //-----------------------------------------------------------
    case 1:
      switch(mang){  
        case 0:
          for(i=0;i<ncoef;i++){
            aka = (double)(ka[i]);
            akb = (double)(kb[i]);
            akc = (double)(kc[i]);
            xk  = (aka*hmati[1]+akb*hmati[2]+akc*hmati[3])*tpi;
            yk  = (aka*hmati[4]+akb*hmati[5]+akc*hmati[6])*tpi;
            zk  = (aka*hmati[7]+akb*hmati[8]+akc*hmati[9])*tpi;
            g2  = xk*xk+yk*yk+zk*zk;
            g   = sqrt(g2);
            if(g!=0.0){
              ctheta   = zk/g;
              y10      =  rt_thrfpi*ctheta;
              dy_re[i] = d_re[i]*y10;
              dy_im[i] = d_im[i]*y10;
	    }else{
              dy_re[i] = 0.0;
              dy_im[i] = 0.0;
	    }
#ifdef CMK_VERSION_BLUEGENE
           if(i%nfreq==0){CmiNetworkProgress();}
#endif
	  }//endfor
        break;
        case 1:
          for(i=0;i<ncoef;i++){
            aka = (double)(ka[i]);
            akb = (double)(kb[i]);
            akc = (double)(kc[i]);
            xk  = (aka*hmati[1]+akb*hmati[2]+akc*hmati[3])*tpi;
            yk  = (aka*hmati[4]+akb*hmati[5]+akc*hmati[6])*tpi;
            zk  = (aka*hmati[7]+akb*hmati[8]+akc*hmati[9])*tpi;
            g2  = xk*xk+yk*yk+zk*zk;
            g   = sqrt(g2);
            if(g!=0.0){
              gs  = sqrt(xk*xk + yk*yk);
              stheta = gs/g;
              cphi   = (gs==0.0 ? 1.0 : xk/gs);
              y11    =  rt_threpi*stheta*cphi;
              dy_re[i] = d_re[i]*y11;
              dy_im[i] = d_re[i]*y11;
	    }else{
              dy_re[i] = 0.0;
              dy_im[i] = 0.0;
	    }//endif
#ifdef CMK_VERSION_BLUEGENE
            if(i%nfreq==0){CmiNetworkProgress();}
#endif
          }//endfor
        break;
        case -1:
          for(i=0;i<ncoef;i++){
            aka = (double)(ka[i]);
            akb = (double)(kb[i]);
            akc = (double)(kc[i]);
            xk  = (aka*hmati[1]+akb*hmati[2]+akc*hmati[3])*tpi;
            yk  = (aka*hmati[4]+akb*hmati[5]+akc*hmati[6])*tpi;
            zk  = (aka*hmati[7]+akb*hmati[8]+akc*hmati[9])*tpi;
            g2  = xk*xk+yk*yk+zk*zk;
            g   = sqrt(g2);
            if(g!=0.0){
              gs  = sqrt(xk*xk + yk*yk);
              stheta = gs/g;
              sphi   = (gs==0.0 ? 0.0 : yk/gs);
              y11    =  rt_threpi*stheta*sphi;
              dy_re[i] = d_re[i]*y11;
              dy_im[i] = d_re[i]*y11;
	    }else{
              dy_re[i] = 0.0;
              dy_im[i] = 0.0;
	    }//endif
#ifdef CMK_VERSION_BLUEGENE
            if(i%nfreq==0){CmiNetworkProgress();}
#endif
          }//endfor
        break;
      }//end : switch m : l = 0
    break;
    //-----------------------------------------------------------
    case 2:
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("No l=2 ees projectors yet\n"); 
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout); EXIT(1);
    break;
    //-----------------------------------------------------------
    case 3:
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("No l=3 ees projectors yet\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout); EXIT(1);
    break;
  }//end swithc : l

#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();
#endif

//==========================================================================
  }//end routine
//==========================================================================


//==========================================================================
// Invoked by all real space nl chares : Compute your part of the current zmat.
//  zmat must be large enough so that all NL interations can be stored at the
//  same time. Otherwise, we have to synchronize! zmat[natm_nl*nlcnt]
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void CPNONLOCAL::eesZmatRchare(double *projPsiR, int iter_nl, double *zmat, 
                               int **igrid, double **mn,int *plane_index,
                               int plane){
//==========================================================================
  CP           *cp           = CP::get();
  GENERAL_DATA *general_data = GENERAL_DATA::get();
#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_gen.h"

  double vol         = gencell->vol;

  int n_ang_max  = cppseudo->n_ang_max;
  int n_interp   = cppseudo->n_interp_ps;
  int n_interp2  = n_interp*n_interp;

  int *natm_typ_lang = cppseudo->nonlocal.natm_typ_lang;
  int *natm_lang     = cppseudo->nonlocal.natm_lang;
  int *iatm_str_lang = cppseudo->nonlocal.iatm_str_lang;
  int **iatm_typ_lang = cppseudo->nonlocal.iatm_typ_lang;
  int *lang_v        = cppseudo->nonlocal.lang_v;
  int *mang_v        = cppseudo->nonlocal.mang_v;
  int *ityp_v        = cppseudo->nonlocal.ityp_v;
  double *vpsnorm    = cppseudo->vpsnorm; 

  int lang           = lang_v[iter_nl];
  int mang           = mang_v[iter_nl];
  int ityp           = ityp_v[iter_nl];
  int lang1          = lang+1;

  int iatm_typ       = iatm_typ_lang[ityp][lang1]; // true atom type      
  int natm           = natm_lang[iatm_typ];        // # atms of this type 
  int iatm_str       = iatm_str_lang[iatm_typ];    // where atms begin    

  int ind_now        = (ityp-1)*(n_ang_max+1) + lang + 1;
  double vnorm_now   = vpsnorm[ind_now];

//==========================================================================
// The projPsiR should be FFT3D(projPsiG)
// Mr. zmat should be passed in with the correct offset for interation counter.

   for(int jatm=0;jatm<natm;jatm++){ // atms of this type
     int iatm   = iatm_str+jatm-1;   // non-local atom index
     int jc     = plane_index[iatm]; // interpolation 
     zmat[jatm] = 0.0; 
     if(jc>0){
       for(int j=1;j<=n_interp2;j++){
         int ind     = igrid[iatm][j];             // index of pt in the plane
         double p    = projPsiR[ind]*mn[iatm][j]; // projection operator
         zmat[jatm] += p;                          // contribute to zmatrix
       }//endfor : B-spline interp to get Zmat contributions
#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif
     }//endif
   }//endfor : atoms of this type

//==========================================================================
// zmat must now be reduced over all planes of the current state.
// We malloc enough zmat memory to contain iterations. The msg should contain
// the iteration counter. Really the FFTs guarantee you won't get into trouble 
// but thats OK. Saftey first.
//==========================================================================
  }//end routine
//==========================================================================


//==========================================================================
// Every Real Space Chare array invokes me for each non-local loop iteration
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void CPNONLOCAL::eesEnergyAtmForcRchare(int iter_nl, double *cp_enl_tot, double *zmat, 
                int **igrid,double **dmn_x,double **dmn_y,double **dmn_z,
		double *projPsiR, int *plane_index, int plane, Atom *atoms)
//==========================================================================
   {//Begin Routine 
//==========================================================================
  CP           *cp           = CP::get();
  GENERAL_DATA *general_data = GENERAL_DATA::get();
#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_gen.h"
  PSNONLOCAL *nonlocal = &(cppseudo->nonlocal);

  double vol         = gencell->vol;

  int n_ang_max  =  cppseudo->n_ang_max;
  int n_interp   = cppseudo->n_interp_ps;
  int n_interp2  = n_interp*n_interp;

  int **iatm_typ_lang = cppseudo->nonlocal.iatm_typ_lang;
  int *natm_typ_lang = cppseudo->nonlocal.natm_typ_lang;
  int *natm_lang     = cppseudo->nonlocal.natm_lang;
  int *iatm_str_lang = cppseudo->nonlocal.iatm_str_lang;
  int *lang_v        = cppseudo->nonlocal.lang_v;
  int *mang_v        = cppseudo->nonlocal.mang_v;
  int *ityp_v        = cppseudo->nonlocal.ityp_v;
  double *vpsnorm    = cppseudo->vpsnorm; 
  int *map_nl    = nonlocal->map_nl;

  int lang           = lang_v[iter_nl];
  int mang           = mang_v[iter_nl];
  int ityp           = ityp_v[iter_nl];
  int lang1          = lang+1;

  int iatm_typ       = iatm_typ_lang[ityp][lang1]; // true atom type      
  int natm           = natm_lang[iatm_typ];        // # atms of this type 
  int iatm_str       = iatm_str_lang[iatm_typ];    // where atms begin    

  int ind_now        = (ityp-1)*(n_ang_max+1) + lang + 1;
  double vnorm_now   = vpsnorm[ind_now];

//==========================================================================
// Compute cp_enl whether you are the reduction plane or not.
// Its cheap, so why not? Then Piny code does not care what plane it is given.

   double cp_enl = 0.0; 
   for(int jatm=0;jatm<natm;jatm++){ // atms of this type
     cp_enl += ((zmat[jatm]*zmat[jatm])*(vnorm_now/vol));
   }//endfor
   cp_enl_tot[0] += cp_enl;

// When all non-local interations are done, the reduction plane, ONLY,
// contributes to a reduction over all states to get the total 
// non-local energy, ENL. 
//==========================================================================
// Atom forces

   for(int jatm=0;jatm<natm;jatm++){// loop over all atms of this type
     int iatm    = iatm_str+jatm-1;    // index of atm in non-local atm list
     int katm    = map_nl[(iatm+1)]-1; // index of atm in full atom list
     int jc      = plane_index[iatm];  // interpolation ind to plane ind mapping
     zmat[jatm] *= (2.0*vnorm_now/vol); 
     if(jc>0){ // jc starts at 1 in piny-like fashion
       for(int j=1;j<=n_interp2;j++){
         int ind   = igrid[iatm][j];            // index of pt in the plane
         double pz = projPsiR[ind]*zmat[jatm];  // projector*psi
         atoms[katm].fx -= pz*dmn_x[iatm][j];   // forces
         atoms[katm].fy -= pz*dmn_y[iatm][j];
         atoms[katm].fz -= pz*dmn_z[iatm][j];
       }//endfor : B spline interpolation for fatm
#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif
     }//endif : atom is interpolated on this plane
   }//endfor : iatm

//--------------------------------------------------------------------------
  }//end routine
//==========================================================================


//==========================================================================
// Every Real Space Chare array invokes me for each non-local loop iteration
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void CPNONLOCAL::eesPsiForcRchare(int iter_nl, double *zmat, 
 	             int **igrid,double **mn,double *projPsiR, 
                     int *plane_index,int plane)
//==========================================================================
   {// begin routine
//==========================================================================

  CP           *cp           = CP::get();
  GENERAL_DATA *general_data = GENERAL_DATA::get();
#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_gen.h"

  double vol         = gencell->vol;

  int nsplin_g   = cppseudo->nsplin_g;
  int n_ang_max  = cppseudo->n_ang_max;
  int n_interp   = cppseudo->n_interp_ps;
  int n_interp2  = n_interp*n_interp;

  int *natm_typ_lang = cppseudo->nonlocal.natm_typ_lang;
  int *natm_lang     = cppseudo->nonlocal.natm_lang;
  int *iatm_str_lang = cppseudo->nonlocal.iatm_str_lang;
  int **iatm_typ_lang = cppseudo->nonlocal.iatm_typ_lang;
  int *lang_v        = cppseudo->nonlocal.lang_v;
  int *mang_v        = cppseudo->nonlocal.mang_v;
  int *ityp_v        = cppseudo->nonlocal.ityp_v;
  int ngrid_a        = cppseudo->nonlocal.ngrid_a;
  int ngrid_b        = cppseudo->nonlocal.ngrid_b;
  double *vpsnorm    = cppseudo->vpsnorm; 

  int lang           = lang_v[iter_nl];
  int mang           = mang_v[iter_nl];
  int ityp           = ityp_v[iter_nl];
  int lang1          = lang+1;

  int iatm_typ       = iatm_typ_lang[ityp][lang1]; // true atom type      
  int natm           = natm_lang[iatm_typ];        // # atms of this type 
  int iatm_str       = iatm_str_lang[iatm_typ];    // where atms begin    

  int ind_now        = (ityp-1)*(n_ang_max+1) + lang + 1;
  double vnorm_now   = vpsnorm[ind_now];

//==========================================================================
// The zmat had better be reduced over all planes of this state.

   for(int i=0;i<(ngrid_a+2)*ngrid_b;i++){projPsiR[i]=0.0;}

   for(int jatm=0;jatm<natm;jatm++){// atms of this type
     int iatm = iatm_str+jatm-1;   // non-local atom index
     int jc   = plane_index[iatm]; // interpolation index to plane index mapping
     if(jc>0){ // jc starts at 1 in piny-like fashion
       for(int j=1;j<=n_interp2;j++){
         int ind        = igrid[iatm][j];         // index of pt in the plane
         double q       = zmat[jatm]*mn[iatm][j]; // contrib of this atm to psi force
         projPsiR[ind] += q;                      // add contrib into total
       }//endfor : B spline interpolation for fcoef
#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif
     }//endif
   }//endfor : iatm

// projPSiR now goeth forth to FFT3D land
//-------------------------------------------------------------------------
  }//end routine
//==========================================================================



//==========================================================================
// Every G-Space Chare array invokes me for each non-local loop iteration
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void CPNONLOCAL::eesPsiForcGspace(int ncoef, int ihave_g0, int ind_g0,
    			          complex *projPsiG,complex *fPsiG, 
                                  double *dyp_re,double *dyp_im)
//==========================================================================
   {//Begin Routine 
//==========================================================================
// projPsiG is FFT3Dinv(projPsiR)

   int nfreq = 100;

   for(int ig=0;ig<ncoef;ig++){ 
     double qtilde_re = projPsiG[ig].re;
     double qtilde_im = projPsiG[ig].im;
     fPsiG[ig].re    -= 2.0*(dyp_re[ig]*qtilde_re - dyp_im[ig]*qtilde_im);
     fPsiG[ig].im    += 2.0*(dyp_im[ig]*qtilde_re + dyp_re[ig]*qtilde_im);
#ifdef CMK_VERSION_BLUEGENE
     if(ig%nfreq==0){CmiNetworkProgress();}
#endif
   }//endfor

   if(ihave_g0==1){
     int ig           = ind_g0;
     double qtilde_re = projPsiG[ig].re;
     fPsiG[ig].re    -= 2.0*dyp_re[ig]*qtilde_re;
     fPsiG[ig].im     = 0.0;
   }//endif

#ifdef CMK_VERSION_BLUEGENE
   CmiNetworkProgress();
#endif
   // Done with this iteration

//-------------------------------------------------------------------------
  }//end routine
//==========================================================================
