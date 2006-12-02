#include "standard_include.h"
#include "ckcomplex.h"

#include "../../../include/Atoms.h"    
#include "../../../include/eesDataClass.h"    
#include "../../../src_mathlib/mathlib.h"    
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_gen.h"

#include "../class_defs/CP_OPERATIONS/class_cpnonlocal.h"
#include "../proto_defs/proto_cp_ewald_local.h"

//#define _CP_DEBUG_EES_NONLOCAL_

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
  double h0,h;
  double tpi = 2.0*M_PI;
  double g2,g;

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
void CPNONLOCAL::eesAtmBsplineRgrp(FastAtoms *atoms, int *allowed_planes, RPPDATA *RPPData)
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

  int n_interp2  = n_interp*n_interp;

  double *aj     = nonlocal->aj;
  double *rn     = nonlocal->rn;
  double *rn1    = nonlocal->rn1;
  int *index_a   = nonlocal->index_a;
  int *index_b   = nonlocal->index_b;
  int *igrid_at  = nonlocal->igrid_at;
  int *igrid_bt  = nonlocal->igrid_bt;
  int *igo       = nonlocal->iatemp;
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
  double **ua    = nonlocal->ua;
  double **ub    = nonlocal->ub;
  double **uc    = nonlocal->uc;
  double **dmn_a = nonlocal->dmn_a;
  double **dmn_b = nonlocal->dmn_b;
  double **dmn_c = nonlocal->dmn_c;

  double *xatm   = atoms->x;
  double *yatm   = atoms->y;
  double *zatm   = atoms->z;

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
     x     = xatm[k];
     y     = yatm[k];
     z     = zatm[k];
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
#ifdef CMK_VERSION_BLUEGENE
     if(j%100==0)
       CmiNetworkProgress();
#endif
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
       igrid       = RPPData[j].igrid;
       for(i=0;i<natm;i++){
          plane_index[i] = 0;
          for(jc=1;jc<=n_interp*n_interp;jc++){igrid[i][jc]=-2;}
       }//endfor
     }//endif
   }//endfor

#ifdef CMK_VERSION_BLUEGENE
        CmiNetworkProgress();
#endif

   for(i=0;i<natm;i++){
     igo[i] = 0;
     for(jc=1;jc<=n_interp;jc++){
       int ip = igrid_c[jc][i];
       if(allowed_planes[ip]==1){igo[i]=1; break;}
     }//endfor
   }//endfor

   for(i=0;i<natm;i++){
    if(igo[i]==1){   
     for(j=1;j<=n_interp;j++){
       igrid_at[j] = igrid_a[j][i]; index_a[j] = j;  
       igrid_bt[j] = igrid_b[j][i]; index_b[j] = j;
     }//endfor
#define _CP_CACHE_BETTER_OFF
#ifdef _CP_CACHE_BETTER_
     sort_commence_piny(n_interp,igrid_at,index_a);
     sort_commence_piny(n_interp,igrid_bt,index_b);
#endif
     for(jc=1;jc<=n_interp;jc++){
       int ip = igrid_c[jc][i];
       if(allowed_planes[ip]==1){
        plane_index = RPPData[ip].plane_index;
        igrid       = RPPData[ip].igrid;
        mn          = RPPData[ip].mn;
        dmn_x       = RPPData[ip].dmn_x;
        dmn_y       = RPPData[ip].dmn_y;
        dmn_z       = RPPData[ip].dmn_z;
        plane_index[i] = jc; // Each jc is a different plane. Ex : plane_index[3] = 4
                            //      For the 3rd atom, the 4th c-interpolation pt
                            //      is on plane number ip=25
        jj = 1;
        for(jb=1;jb<=n_interp;jb++){
         for(ja=1,j=jj;ja<=n_interp;ja++,j++){
           igrid[i][j] = igrid_at[ja]+igrid_bt[jb]; //plane index only
           int jaa     = index_a[ja];
           int jbb     = index_b[jb];
           atemp       = dmn_a[jaa][i]* mn_b[jbb][i]* mn_c[jc][i]*grid_a;
           btemp       =  mn_a[jaa][i]*dmn_b[jbb][i]* mn_c[jc][i]*grid_b;
           ctemp       =  mn_a[jaa][i]* mn_b[jbb][i]*dmn_c[jc][i]*grid_c;
           mn[i][j]    =  mn_a[jaa][i]* mn_b[jbb][i]* mn_c[jc][i];
           dmn_x[i][j] = atemp*hmati[1]+btemp*hmati[2]+ctemp*hmati[3];
           dmn_y[i][j] = atemp*hmati[4]+btemp*hmati[5]+ctemp*hmati[6];
           dmn_z[i][j] = atemp*hmati[7]+btemp*hmati[8]+ctemp*hmati[9];
         }//endfor : ja
         jj += n_interp;
        }//endfor : jb
#define _CP_BRK_BETTER_OFF // must flip cp_ees_nonlocal.h too
#ifdef _CP_BRK_BETTER_
        int *nBreakJ  = RPPData[ip].nBreakJ;
        int **sBreakJ = RPPData[ip].sBreakJ;
        nBreakJ[i]    = 1;   // length natm
        sBreakJ[i][1] = 1;   // length [natm][n_interp2+1]
        for(int jj=2;jj<=n_interp2;jj++){
          if(igrid[i][jj]!=igrid[i][(jj-1)]+1){
	    nBreakJ[i]            += 1;
            sBreakJ[i][nBreakJ[i]] = jj;
	  }//endif
	}//endfor
        sBreakJ[i][(nBreakJ[i]+1)] = n_interp2+1;
#endif
#define _CP_TEST_BREAK_OFF
#ifdef _CP_TEST_BREAK_
        int num = 0;
        for(int ib=1;ib<=nBreakJ[i];ib++){
          num += sBreakJ[i][(ib+1)]-sBreakJ[i][ib];
          if(sBreakJ[i][(ib+1)]!=num+1){PRINTF("Help.1 %d\n",num); EXIT(1);}
        }//endif
        if(num!=n_interp2){PRINTF("Help %d\n",num); EXIT(1);}
#endif
#ifdef CMK_VERSION_BLUEGENE
        CmiNetworkProgress();
#endif
       }//endif : this jc is cool
      }//endfor jc
    }//endif : allowed
   }//endfor : iatm

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
                     complex *projPsiG, int *ind_gspl, double *h_gspl,
                     int istate,int ichare)
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
  double *gzvps0 = cppseudo->gzvps0;

  // piny indices : start at 1
  int lang     = lang_v[iter_nl]; // l-channel 
  int mang     = mang_v[iter_nl]; // m-channel
  int ityp     = ityp_v[iter_nl]; // the 3rd type is this lm channel

  int lang1    = lang+1;
  int iatm_typ = iatm_typ_lang[ityp][lang1]; // true atom type      
  int natm     = natm_lang[iatm_typ];        // # atms of this type 
  int iatm_str = iatm_str_lang[iatm_typ];    // where atms begin    

  int nfreq=100;

//==========================================================================
// Set up the debugging stuff

#ifdef _CP_DEBUG_EES_NONLOCAL_
  char myFileName[1000];
  sprintf(myFileName, "proj_Gpsi_%d.out.%d.%d",istate,ichare,iter_nl);
  FILE *fp;
  if(istate==4){fp = fopen(myFileName,"w");}
#endif  

//==========================================================================
// The dyp_re,dyp_im,psi, must be memory on each chare!! It is reused below!! 
// Different states can have different iter_nl leading to different dyp's hence
// they cannot be cached without synchonizing.
// The d_re,d_im are from the cache and only depend on g-plane index.

  eesYlmOnD(lang,mang,ncoef,ka,kb,kc,dyp_re,dyp_im,d_re,d_im,hmati);

  int ii = (ihave_g0==1 ? ind_g0 : ncoef);
  if(lang==0 && ihave_g0==1){
    double vnow = gzvps0[iatm_typ];
    dyp_re[ii] *= vnow;
    dyp_im[ii]  = 0.0;
  }//endif
  if(lang!=0 && ihave_g0==1){
    dyp_re[ii] = 0.0;
    dyp_im[ii] = 0.0;
  }//endif

  int ind_off = (iatm_typ-1)*nsplin_g*(n_ang_max+1) + lang*nsplin_g;
  for(int i=0;i<ii;i++){
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
#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();
#endif
  for(int i=ii+1;i<ncoef;i++){
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

//==========================================================================
// Using your psi and your dyp create the puppy to be 3DFFTed to real space

  // congj(psi * struct_fact)
  for(int ig=0;ig<ncoef;ig++){ 
     projPsiG[ig].re =  dyp_re[ig]*psi[ig].re-dyp_im[ig]*psi[ig].im;
     projPsiG[ig].im =-(dyp_im[ig]*psi[ig].re+dyp_re[ig]*psi[ig].im);
#ifdef _CP_DEBUG_EES_NONLOCAL_
     if(istate==4){
      fprintf(fp,"%d %d %d : %g %g\n",ka[ig],kb[ig],kc[ig],projPsiG[ig].re,projPsiG[ig].im);
     }/*endif*/
#endif
#ifdef CMK_VERSION_BLUEGENE
    if(ig%nfreq==0){CmiNetworkProgress();}
#endif
  }//endfor

#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif

//==========================================================================
// close any debug files

#ifdef _CP_DEBUG_EES_NONLOCAL_
  if(istate==4){fclose(fp);}
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
                               int state,int plane){
//==========================================================================
  CP           *cp           = CP::get();
  GENERAL_DATA *general_data = GENERAL_DATA::get();
#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_gen.h"

  double vol         = gencell->vol;

  int n_ang_max  = cppseudo->n_ang_max;
  int n_interp   = cppseudo->n_interp_ps;
  int n_interp2  = n_interp*n_interp;

  int *natm_typ_lang  = cppseudo->nonlocal.natm_typ_lang;
  int *natm_lang      = cppseudo->nonlocal.natm_lang;
  int *iatm_str_lang  = cppseudo->nonlocal.iatm_str_lang;
  int **iatm_typ_lang = cppseudo->nonlocal.iatm_typ_lang;
  int *lang_v         = cppseudo->nonlocal.lang_v;
  int *mang_v         = cppseudo->nonlocal.mang_v;
  int *ityp_v         = cppseudo->nonlocal.ityp_v;
  double *vpsnorm     = cppseudo->vpsnorm; 

  int ngrida          = cppseudo->nonlocal.ngrid_a;
  int ngridb          = cppseudo->nonlocal.ngrid_b;

  int lang            = lang_v[iter_nl];
  int mang            = mang_v[iter_nl];
  int ityp            = ityp_v[iter_nl];
  int lang1           = lang+1;

  int iatm_typ        = iatm_typ_lang[ityp][lang1]; // true atom type      
  int natm            = natm_lang[iatm_typ];        // # atms of this type 
  int iatm_str        = iatm_str_lang[iatm_typ];    // where atms begin    

//==========================================================================
// A little debugging for you

#ifdef _CP_DEBUG_EES_NONLOCAL
  if(state==4){
    char myFileName[1000];
    sprintf(myFileName, "proj_Rpsi_%d.out.%d.%d",state,plane,iter_nl);
    FILE *fp = fopen(myFileName,"w");
    int ic = 0;
    for(int j=0;j<ngridb;j++){
     for(int i=0;i<ngrida;i++){
       fprintf(fp,"%d %g\n",i,projPsiR[ic]);
       ic++;
     }//endfor
     ic+=2;
    }//endfor
    fclose(fp);
  }//endif : state=0
#endif  

//==========================================================================
// The projPsiR should be FFT3D(projPsiG)
// Mr. zmat should be passed in with the correct offset for interation counter.

   int nroll = 5; // you can't check this without modifying the code below
   int nrem  = (n_interp2 % nroll);
   int jstrt = (n_interp2-nrem+1);
   int jend  = (n_interp2-nrem);

   bzero(zmat,sizeof(double)*natm);
   for(int jatm=0;jatm<natm;jatm++){ // atms of this type
     int iatm   = iatm_str+jatm-1;   // non-local atom index
     int jc     = plane_index[iatm]; // interpolation 
     if(jc>0){
       for(int j0=1,j1=2,j2=3,j3=4,j4=5;j0<=jend;
           j0+=nroll,j1+=nroll,j2+=nroll,j3+=nroll,j4+=nroll){
         double p0   = projPsiR[igrid[iatm][j0]]*mn[iatm][j0];  // projection operator
         double p1   = projPsiR[igrid[iatm][j1]]*mn[iatm][j1]; 
         double p2   = projPsiR[igrid[iatm][j2]]*mn[iatm][j2]; 
         double p3   = projPsiR[igrid[iatm][j3]]*mn[iatm][j3]; 
         double p4   = projPsiR[igrid[iatm][j4]]*mn[iatm][j4]; 
         zmat[jatm] += (p0+p1+p2+p3+p4);         // add to zmatrix
       }//endfor : B-spline interp to get Zmat contributions
       for(int j=jstrt;j<=n_interp2;j++){
         double p   = projPsiR[igrid[iatm][j]]*mn[iatm][j]; // projection operator
         zmat[jatm] += p;                                   // add to zmatrix
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
// A little debugging for you

#ifdef _CP_DEBUG_EES_NONLOCAL_
  if(state==4){
    char myFileName2[1000];
    sprintf(myFileName2, "proj_zmat_%d.out.%d.%d",state,plane,iter_nl);
    FILE *fp = fopen(myFileName2,"w");
     for(int jatm=0;jatm<natm;jatm++){ // atms of this type
       fprintf(fp,"%d %g\n",jatm,zmat[jatm]);
     }//endfor
    fclose(fp);
  }//endif : state=0
#endif  

//==========================================================================
  }//end routine
//==========================================================================


//==========================================================================
// Every Real Space Chare array invokes me for each non-local loop iteration
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void CPNONLOCAL::eesEnergyAtmForcRchare(int iter_nl, double *cp_enl_tot, double *zmat, 
                int **igrid,double **mn,double **dmn_x,double **dmn_y,double **dmn_z,
		double *projPsiR, double *projPsiRScr, int *plane_index, 
		int *nBreak, int **sBreak,
                int plane, int state, FastAtoms *atoms)
//==========================================================================
   {//Begin Routine 
//==========================================================================
  CP           *cp           = CP::get();
  GENERAL_DATA *general_data = GENERAL_DATA::get();
#include "../class_defs/allclass_strip_cp.h"
#include "../class_defs/allclass_strip_gen.h"
  PSNONLOCAL *nonlocal = &(cppseudo->nonlocal);

  double vol     = gencell->vol;

  double *vpsnorm = cppseudo->vpsnorm; 
  int n_ang_max   = cppseudo->n_ang_max;
  int n_interp    = cppseudo->n_interp_ps;
  int n_interp2   = n_interp*n_interp;

  int ngrid_a         = nonlocal->ngrid_a;
  int ngrid_b         = nonlocal->ngrid_b;
  int ngrid_c         = nonlocal->ngrid_c;
  int **iatm_typ_lang = nonlocal->iatm_typ_lang;
  int *natm_typ_lang  = nonlocal->natm_typ_lang;
  int *natm_lang      = nonlocal->natm_lang;
  int *iatm_str_lang  = nonlocal->iatm_str_lang;
  int *lang_v         = nonlocal->lang_v;
  int *mang_v         = nonlocal->mang_v;
  int *ityp_v         = nonlocal->ityp_v;
  int *map_nl         = nonlocal->map_nl;

  double *fx         = atoms->fx;
  double *fy         = atoms->fy;
  double *fz         = atoms->fz;

  int lang           = lang_v[iter_nl];
  int mang           = mang_v[iter_nl];
  int ityp           = ityp_v[iter_nl];
  int lang1          = lang+1;

  int iatm_typ       = iatm_typ_lang[ityp][lang1]; // true atom type      
  int natm           = natm_lang[iatm_typ];        // # atms of this type 
  int iatm_str       = iatm_str_lang[iatm_typ];    // where atms begin    

  int ind_now        = (iatm_typ-1)*(n_ang_max+1) + lang + 1;
  double vnorm_now   = vpsnorm[ind_now];
  double vnormVol    = vnorm_now/vol;

  int nroll = 5; // you can't check this without modifying the code below
  int nrem,jstrt,jend;

//==========================================================================
// Set up some debuggin stuff

#ifdef _CP_DEBUG_EES_NONLOCAL_
   double *fxt = new double [natm];
   double *fyt = new double [natm];
   double *fzt = new double [natm];
   for(int i =0;i<natm;i++){fxt[i]=0.;fyt[i]=0.;fzt[i]=0.;}
#endif

//==========================================================================
// Compute cp_enl whether you are the reduction plane or not.
// Its cheap, so why not? Then Piny code does not care what plane it is given.
// When all non-local interations are done, the reduction plane, ONLY,
// contributes to a reduction over all states to get the total 
// non-local energy, ENL. 

   nrem  = (natm % nroll);
   jstrt = (natm-nrem);
   jend  = (natm-nrem);

   double cp_enl = 0.0; 
   for(int j=0;j<jend;j+=nroll){
     cp_enl += (zmat[j]    *zmat[j]     + zmat[(j+1)]*zmat[(j+1)]
               +zmat[(j+2)]*zmat[(j+2)] + zmat[(j+3)]*zmat[(j+3)]
               +zmat[(j+4)]*zmat[(j+4)]);
   }//endfor
   for(int j=jstrt;j<natm;j++){cp_enl += (zmat[j]*zmat[j]);}
   cp_enl_tot[0] += (cp_enl*vnormVol);
#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif
//==========================================================================
// Atom forces

   nrem  = (n_interp2 % nroll);
   jstrt = (n_interp2-nrem+1);
   jend  = (n_interp2-nrem);

   //   bzero(projPsiRScr,(ngrid_a+2)*ngrid_b*sizeof(double));
   for(int z=0; z<(ngrid_a+2)*ngrid_b;z++)
     projPsiRScr[z]=0.0;

//#define _DEBUG_STUFF_   
#ifdef  _DEBUG_STUFF_   
    char myFileName[1000];
    sprintf(myFileName, "Rproj_%d.out.%d.%d",state,plane,iter_nl);
    FILE *fp = fopen(myFileName,"w");
#define _UNROLL_OFF_
#endif

   for(int jatm=0;jatm<natm;jatm++){// loop over all atms of this type
     int iatm = iatm_str+jatm-1;    // index of atm in non-local atm list
     int jc   = plane_index[iatm];  // interpolation ind to plane ind mapping
     if(jc>0){ // jc starts at 1 in piny-like fashion
       int katm    = map_nl[(iatm+1)]-1;  // index of atm in full atom list
       zmat[jatm] *= (2.0*vnormVol); 
       double fxx=0.0,fyy=0.0,fzz=0.0;
#ifndef _CP_BRK_BETTER_
#ifndef _UNROLL_OFF_
       for(int j=1,j1=2,j2=3,j3=4,j4=5;j<=jend;
           j+=nroll,j1+=nroll,j2+=nroll,j3+=nroll,j4+=nroll){
         double pz0 = projPsiR[igrid[iatm][j]];   // ProjPsi
         double pz1 = projPsiR[igrid[iatm][j1]];     
         double pz2 = projPsiR[igrid[iatm][j2]];     
         double pz3 = projPsiR[igrid[iatm][j3]];     
         double pz4 = projPsiR[igrid[iatm][j4]];     
         double q0  = zmat[jatm]*mn[iatm][j];     // contrib of this atm to psi force
         double q1  = zmat[jatm]*mn[iatm][j1];
         double q2  = zmat[jatm]*mn[iatm][j2];
         double q3  = zmat[jatm]*mn[iatm][j3];
         double q4  = zmat[jatm]*mn[iatm][j4];
         fxx       += (pz0*dmn_x[iatm][j]  + pz1*dmn_x[iatm][j1] // forces_x
		       +pz2*dmn_x[iatm][j2] + pz3*dmn_x[iatm][j3]
		       +pz4*dmn_x[iatm][j4]);
         fyy       += (pz0*dmn_y[iatm][j]  + pz1*dmn_y[iatm][j1] // forces_y
		       +pz2*dmn_y[iatm][j2] + pz3*dmn_y[iatm][j3]
		       +pz4*dmn_y[iatm][j4]);
         fzz       += (pz0*dmn_z[iatm][j]  + pz1*dmn_z[iatm][j1] // forces_z
		       +pz2*dmn_z[iatm][j2] + pz3*dmn_z[iatm][j3]
		       +pz4*dmn_z[iatm][j4]);

         projPsiRScr[igrid[iatm][j]]  += q0;       // add contrib into total
         projPsiRScr[igrid[iatm][j1]] += q1;  
         projPsiRScr[igrid[iatm][j2]] += q2;  
         projPsiRScr[igrid[iatm][j3]] += q3;  
         projPsiRScr[igrid[iatm][j4]] += q4;  
       }//endfor
       for(int j=jstrt;j<=n_interp2;j++){
#else //unroll is on
       for(int j=1;j<=n_interp2;j++){
#endif //unroll is off
         double pz = projPsiR[igrid[iatm][j]]; // psi
         double q  = zmat[jatm]*mn[iatm][j];  
         fxx      += (pz*dmn_x[iatm][j]);
         fyy      += (pz*dmn_y[iatm][j]);
         fzz      += (pz*dmn_z[iatm][j]);
         projPsiRScr[igrid[iatm][j]] += q;        // add contrib into total
#ifdef  _DEBUG_STUFF_   
         fprintf(fp,"%d %d %d %d %d %g %g %g %g %g %d\n",jatm,iatm,katm,jc,igrid[iatm][j],
                                     pz,q,zmat[jatm],mn[iatm][j],projPsiRScr[igrid[iatm][j]], (ngrid_a+2)*ngrid_b);

#endif
       }//endfor
#else
       for(int ib=1;ib<=nBreak[iatm];ib++){
         int koff = igrid[iatm][sBreak[iatm][ib]];
         int joff = sBreak[iatm][ib];
         int nn   = sBreak[iatm][(ib+1)]-joff;  
         nrem     = (nn % nroll);
         jend     = nn-nrem+joff-1;
         jstrt    = jend+1;
         int kstrt= jstrt-joff+koff;
         for(int j=joff,k=koff;j<=jend;j+=5,k+=5){
           double p0 = projPsiR[k];           // psi
           double p1 = projPsiR[k+1];           // psi
           double p2 = projPsiR[k+2];           // psi
           double p3 = projPsiR[k+3];           // psi
           double p4 = projPsiR[k+4];           // psi
           fxx      += (p0*dmn_x[iatm][j]   + p1*dmn_x[iatm][j+1] + p2*dmn_x[iatm][j+2] 
			+p3*dmn_x[iatm][j+3] + p4*dmn_x[iatm][j+4]);
           fyy      += (p0*dmn_y[iatm][j]   + p1*dmn_y[iatm][j+1] + p2*dmn_y[iatm][j+2] 
			+p3*dmn_y[iatm][j+3] + p4*dmn_y[iatm][j+4]);
           fzz      += (p0*dmn_z[iatm][j]   + p1*dmn_z[iatm][j+1] + p2*dmn_z[iatm][j+2] 
			+p3*dmn_z[iatm][j+3] + p4*dmn_z[iatm][j+4]);
           double q0  = zmat[jatm]*mn[iatm][j];  
           double q1  = zmat[jatm]*mn[iatm][j+1];  
           double q2  = zmat[jatm]*mn[iatm][j+2];  
           double q3  = zmat[jatm]*mn[iatm][j+3];  
           double q4  = zmat[jatm]*mn[iatm][j+4];  
           projPsiRScr[k]   += q0;               // add contrib into total
           projPsiRScr[k+1] += q1;
           projPsiRScr[k+2] += q2;
           projPsiRScr[k+3] += q3;
           projPsiRScr[k+4] += q4;
	 }//endfor
         for(int j=jstrt,k=kstrt;j<sBreak[iatm][(ib+1)];j++,k++){
           double p = projPsiR[k];           // psi
           fxx      += (p*dmn_x[iatm][j]);
           fyy      += (p*dmn_y[iatm][j]);
           fzz      += (p*dmn_z[iatm][j]);
           double q  = zmat[jatm]*mn[iatm][j];  
           projPsiRScr[k] += q;               // add contrib into total
	 }//endfor
       }//endfor
#endif
       fx[katm]  -= (fxx*zmat[jatm]); // finish up
       fy[katm]  -= (fyy*zmat[jatm]);
       fz[katm]  -= (fzz*zmat[jatm]);
#ifdef _CP_DEBUG_EES_NONLOCAL_
       for(int j=1;j<=n_interp2;j++){
         double pz  = projPsiR[igrid[iatm][j]]*zmat[jatm]; // projector*psi
         fxt[jatm] -= pz*dmn_x[iatm][j];                   // forces 
         fyt[jatm] -= pz*dmn_y[iatm][j];
         fzt[jatm] -= pz*dmn_z[iatm][j];
       }//endif
#endif
#ifdef CMK_VERSION_BLUEGENE
       //       CmiNetworkProgress();
#endif
     }//endif : atom is interpolated on this plane
   }//endfor : iatm

//==========================================================================
// output the debugging stuff

#ifdef _CP_DEBUG_EES_NONLOCAL_
    char myFileName[1000];
    sprintf(myFileName, "fatm_Rproj_%d.out.%d.%d",state,plane,iter_nl);
    FILE *fp = fopen(myFileName,"w");
     for(int jatm=0;jatm<natm;jatm++){// loop over all atms of this type
       int iatm  = iatm_str+jatm-1;    // index of atm in non-local atm list
       int katm  = map_nl[(iatm+1)]-1; // index of atm in full atom list
       fprintf(fp,"%d %.10g %.10g %.10g\n",katm,fxt[jatm],fyt[jatm],fzt[jatm]);
     }//endfor
    fclose(fp);
    delete [] fxt;
    delete [] fyt;
    delete [] fzt;
#endif  

#ifdef  _DEBUG_STUFF_   
    fclose(fp);
#endif


//--------------------------------------------------------------------------
  }//end routine
//==========================================================================


//==========================================================================
// Every G-Space Chare array invokes me for each non-local loop iteration
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void CPNONLOCAL::eesPsiForcGspace(int ncoef, int ihave_g0, int ind_g0,int nkx0,
    			          complex *projPsiG,complex *fPsiG, 
                                  double *dyp_re,double *dyp_im,
                                  int *ka, int *kb,int *kc,
                                  int istate,int ichare,int iter_nl)
//==========================================================================
   {//Begin Routine 
//==========================================================================
// Setup some debugging stuff

#ifdef _CP_DEBUG_EES_NONLOCAL_
  char myFileName[1000];
  sprintf(myFileName, "proj_fpsiG_%d.out.%d.%d",istate,ichare,iter_nl);
  FILE *fp;
  if(istate==0){fp = fopen(myFileName,"w");}
#endif

//==========================================================================
// G=0 : Compute Psi forces : projPsiG is FFT3Dinv(projPsiR)

  if(ihave_g0==1){
    int ig           = ind_g0;
    projPsiG[ig].im  = 0.0;
    dyp_im[ig]       = 0.0;
  }//endif

//==========================================================================
// Gx=0 : Compute Psi forces : projPsiG is FFT3Dinv(projPsiR)
//        The sum over iterations leads to a sum reduction here.

  int nfreq   = 100;
  for(int ig=0;ig<nkx0;ig++){ 
    fPsiG[ig].re    -= (dyp_re[ig]*projPsiG[ig].re - dyp_im[ig]*projPsiG[ig].im);
    fPsiG[ig].im    += (dyp_im[ig]*projPsiG[ig].re + dyp_re[ig]*projPsiG[ig].im);
#ifdef _CP_DEBUG_EES_NONLOCAL_
    if(istate==0){fprintf(fp,"%d %d %d : %g %g\n",ka[ig],kb[ig],kc[ig],
                 fPsiG[ig].re,fPsiG[ig].im);}
#endif
#ifdef CMK_VERSION_BLUEGENE
    if(ig%nfreq==0){CmiNetworkProgress();}
#endif
  }//endfor

#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();
#endif

//==========================================================================
// Gx!=0 : Compute Psi forces : projPsiG is FFT3Dinv(projPsiR)

  double wght = 2.0;
  for(int ig=nkx0;ig<ncoef;ig++){ 
    fPsiG[ig].re    -= wght*(dyp_re[ig]*projPsiG[ig].re - dyp_im[ig]*projPsiG[ig].im);
    fPsiG[ig].im    += wght*(dyp_im[ig]*projPsiG[ig].re + dyp_re[ig]*projPsiG[ig].im);
#ifdef _CP_DEBUG_EES_NONLOCAL_
    if(istate==0){fprintf(fp,"%d %d %d : %g %g\n",ka[ig],kb[ig],kc[ig],
                 fPsiG[ig].re,fPsiG[ig].im);}
#endif
#ifdef CMK_VERSION_BLUEGENE
    if(ig%nfreq==0){CmiNetworkProgress();}
#endif
  }//endfor

#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();
#endif

//==========================================================================
// Close the files if debugging

#ifdef _CP_DEBUG_EES_NONLOCAL_
  if(istate==0){fclose(fp);}
#endif

//==========================================================================
// Done with this iteration

//-------------------------------------------------------------------------
  }//end routine
//==========================================================================


//==============================================================================
//  A generic routine to set kvectors and indices from runs
//==============================================================================
void CPNONLOCAL::genericSetKvector(int numPoints, int *k_x, int *k_y, int *k_z,
                        double *g, double *g2,int numRuns, RunDescriptor *runs, 
                        GCHAREPKG *gCharePkg, int checkFill, 
                        int ngrid_a, int ngrid_b, int ngrid_c){
//======================================================================

  GENERAL_DATA *general_data = GENERAL_DATA::get();
#include "../class_defs/allclass_strip_gen.h"
  double *hmati   = gencell->hmati;

//======================================================================
// Construct the k-vectors
  
  int dataCovered = 0;
  for (int r = 0; r < numRuns; r++) { // 2*number of lines z
    int x, y, z;
    x = runs[r].x;
    if (x > ngrid_a/2){x -= ngrid_a;}
    y = runs[r].y;
    if (y > ngrid_b/2){y -= ngrid_b;}
    z = runs[r].z;
    if (z > ngrid_c/2){z -= ngrid_c;}
    for(int i = 0; i < runs[r].length; i++) { //pts in lines of z
      k_x[dataCovered] = x;
      k_y[dataCovered] = y;
      k_z[dataCovered] = (z+i);
      dataCovered++;
    }//endfor
  }//endfor

  CkAssert(dataCovered == numPoints);

//======================================================================
// give me some g's

  double tpi = 2.0*M_PI;
  for(int i=0;i<numPoints;i++){
    double aka = ((double)(k_x[i]));
    double akb = ((double)(k_y[i]));
    double akc = ((double)(k_z[i]));
    double xk  = (aka*hmati[1]+akb*hmati[2]+akc*hmati[3])*tpi;
    double yk  = (aka*hmati[4]+akb*hmati[5]+akc*hmati[6])*tpi;
    double zk  = (aka*hmati[7]+akb*hmati[8]+akc*hmati[9])*tpi;
    g2[i]      = xk*xk+yk*yk+zk*zk;
    g[i]       = sqrt(g2[i]);
  }//endfor

//======================================================================
// Find pts with k_x==0 then check the layout : kx=0 first

  int ihave_g000 =  0;
  int ind_g000   = -1;
  int ihave_kx0  = 0;
  int nkx0       = 0;
  int nkx0_uni   = 0;
  int nkx0_red   = 0;
  int nkx0_zero  = 0;
  int kx0_strt   = 0;
  for(int i=0;i<numPoints;i++){
    if(k_x[i]==0 && k_y[i]>0){nkx0_uni++;}
    if(k_x[i]==0 && k_y[i]<0){nkx0_red++;}
    if(k_x[i]==0 && k_y[i]==0 && k_z[i]>=0){nkx0_uni++;}
    if(k_x[i]==0 && k_y[i]==0 && k_z[i]<0){nkx0_red++;}
    if(k_x[i]==0 && k_y[i]==0 && k_z[i]==0){nkx0_zero++;ihave_g000=1;ind_g000=i;}
    if(k_x[i]==0){
      if(ihave_kx0==0){kx0_strt=i;}
      ihave_kx0=1;
      nkx0++;
    }//endif
  }//endif
  int kx0_end = kx0_strt + nkx0;

  if(checkFill==1){
   if(kx0_strt!=0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("kx=0 should be stored first | kx_srt !=0\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
   }//endif
  
   if(nkx0!=nkx0_uni+nkx0_red){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Incorrect count of redundant guys\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
   }//endif

   for(int i=0;i<nkx0;i++){  
    if(k_x[i]!=0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("kx should be stored consecutively and first\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
   }//endfor

   for(int i=0;i<nkx0_red;i++){  
    if(k_y[i]>0 || (k_y[i]==0 && k_z[i]>=0)){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("ky <0 should be stored first\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
     }//endif
   }//endfor

   for(int i=nkx0_red;i<nkx0_uni;i++){  
    if(k_y[i]<0 || (k_y[i]==0 && k_z[i]<0)){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("ky <0 should be stored first\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
   }//endfor

  }//endif : check the ordering : states only

//==============================================================================
// Set the return values

  gCharePkg->ihave_g000 = ihave_g000;
  gCharePkg->ind_g000   = ind_g000;
  gCharePkg->ihave_kx0  = ihave_kx0;
  gCharePkg->nkx0       = nkx0;
  gCharePkg->nkx0_uni   = nkx0_uni;
  gCharePkg->nkx0_red   = nkx0_red;
  gCharePkg->nkx0_zero  = nkx0_zero;
  gCharePkg->kx0_strt   = kx0_strt;
  gCharePkg->kx0_end    = kx0_end;

//==============================================================================
  }//end routine
//==============================================================================
