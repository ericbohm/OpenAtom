#include "PMEGSpaceData.h"

complex* PMEGSpaceData::getQGrid () {
  return &(qgrid_gspace[1]);
}

void PMEGSpaceData::initialize (int sy_, int nPlanes_, double alp_ewd_, 
                                double vol_, double* hmati_, int n_interp_) {

  sy      = sy_;
  nPlanes = nPlanes_;
  
  double* hmati = hmati_-1;
  
  alp_ewd    = alp_ewd_;
  vol        = vol_;
  n_interp   = n_interp_;

  ecut_recip = 2.0*M_PI*M_PI*alp_ewd*alp_ewd;
  
  count_recip(hmati);
  ka = ((int *) malloc(nktot*sizeof(int)))-1;
  kb = ((int *) malloc(nktot*sizeof(int)))-1;
  kc = ((int *) malloc(nktot*sizeof(int)))-1;
  fill_recip(hmati);

  ngrid_a  = (int)(3.44*ka_max); if( (ngrid_a % 2) == 1){ngrid_a++;}
  ngrid_b  = (int)(3.44*kb_max); if( (ngrid_b % 2) == 1){ngrid_b++;}
  ngrid_c  = (int)(3.44*kc_max); if( (ngrid_c % 2) == 1){ngrid_c++;}

  qgridSize = ngrid_a*ngrid_b*nPlanes;
  qgrid_gspace = ((complex*)malloc(qgridSize*sizeof(complex)))-1;
  
  memset (&(qgrid_gspace[1]), 0, sizeof (complex)*qgridSize);
  
  aj     = (double *)malloc(n_interp*sizeof(double))-1; // size [interpOrder]
  rn     = (double *)malloc(n_interp*sizeof(double))-1; // size [interpOrder]
  rn1    = (double *)malloc(n_interp*sizeof(double))-1; // size [interpOrder]
  
  for(int j=1;j<=n_interp;j++){
    aj[j] = (double) (j-1);
    rn[j] = (double) (j);
    if(j > 1){rn1[j] = 1.0/((double)(j-1));}
  }//endfor
  rn1[1] = 0.0;  
  
  int pme_b_opt = 1;
  double* bfact_r = NULL;  // dummy argument 
  double* bfact_i = NULL;  // dummy argument 
  
  bfact_r = ((double *)malloc(nktot*sizeof(double)))-1;
  bfact_i = ((double *)malloc(nktot*sizeof(double)))-1;

  bweight_tot = ((double *)malloc(nktot*sizeof(double)))-1;
  
  set_pme_wght (pme_b_opt, bfact_r, bfact_i);
    
  delete [] &(bfact_r[1]);
  delete [] &(bfact_i[1]);  
}

void PMEGSpaceData::pup (PUP::er& p) {
  p|sy;
  p|nPlanes;
  p|qgridSize;
    
  p|alp_ewd;
  p|vol;
  p|ecut_recip;

  p|n_interp;

  p|ngrid_a;
  p|ngrid_b;
  p|ngrid_c;

  p|ka_max;
  p|kb_max;
  p|kc_max;
  
  p|nktot;

  if (0 < nktot) {
    pup1d_int(p, &ka, nktot); // size [nktot]
    pup1d_int(p, &kb, nktot); // size [nktot]
    pup1d_int(p, &kc, nktot); // size [nktot]
    
    pup1d_dbl(p, &bweight_tot, nktot); // size [nktot]
  }

  if (0 < qgridSize) {
    if (p.isUnpacking()) {
      qgrid_gspace = ((complex*)malloc(qgridSize*sizeof(complex)))-1;
    }
    p(&(qgrid_gspace[1]),qgridSize);
  }
  
  if (0 < n_interp) {
    pup1d_dbl(p, &aj, n_interp); // size [n_interp]
    pup1d_dbl(p, &rn, n_interp); // size [n_interp]
    pup1d_dbl(p, &rn1, n_interp); // size [n_interp]
  }
}

void PMEGSpaceData::clear() {
  if (0 < nktot) {
    free(&(ka[1])); // size [nktot]
    free(&(kb[1])); // size [nktot]
    free(&(kc[1])); // size [nktot]
    
    free(&(bweight_tot[1])); // size [nktot]

    ka = NULL;
    kb = NULL;
    kc = NULL;
    bweight_tot = NULL;
    
    nktot = 0;
  }

  if (0 < qgridSize) {
    free (&(qgrid_gspace[1]));
    
    qgrid_gspace = NULL;
    qgridSize = 0;
  }
  
  if (0 < n_interp) {
    free (&(aj[1]));
    free (&(rn[1]));
    free (&(rn1[1]));
    
    aj  = NULL;
    rn  = NULL;
    rn1 = NULL;
    
    n_interp = 0;
  }
  
  sy = 0;
  nPlanes = 0;
    
  alp_ewd = 0;
  vol = 0;
  ecut_recip = 0;

  ngrid_a = 0;
  ngrid_b = 0;
  ngrid_c = 0;

  ka_max = 0;
  kb_max = 0;
  kc_max = 0;
}

void PMEGSpaceData::count_recip(double *hmati)

//=========================================================================
   {//begin routine 
//=========================================================================

   int iii;
   int iga,igb,igc;
   int nbmin,nbmax;
   int namin,namax;

   double ecut_scal;
   double sumsq_a,sumsq_b,sumsq_c,sdot_bc;
   double gb,gc;
   double dgx,dgy,dgz,dg2;
   double aa,bb,cc,des;

//=========================================================================
// calculate maximum k-vector

   ecut_scal = ecut_recip*0.5/(M_PI*M_PI);

   sumsq_a = hmati[1]*hmati[1]+hmati[4]*hmati[4]+hmati[7]*hmati[7];
   ka_max  = (int) (sqrt(ecut_scal/sumsq_a));

   sumsq_b = hmati[2]*hmati[2]+hmati[5]*hmati[5]+hmati[8]*hmati[8];
   kb_max  = (int) (sqrt(ecut_scal/sumsq_b));

   sumsq_c = hmati[3]*hmati[3]+hmati[6]*hmati[6]+hmati[9]*hmati[9];
   kc_max  = (int) (sqrt(ecut_scal/sumsq_c));

   sdot_bc = 2.0*(hmati[2]*hmati[3]+hmati[5]*hmati[6]+hmati[8]*hmati[9]);

//=========================================================================

   iii = 0;
   for(igc=0;igc<=kc_max;igc++){
     gc  = (double)igc;
     aa  = sumsq_b;
     bb  = gc*sdot_bc;
     cc  = gc*gc*sumsq_c-ecut_scal;
     des = bb*bb-4.0*aa*cc;
     if(des>=0.0){
       nbmax = (int)( (-bb+sqrt(des))/(2.0*aa) );
       nbmin = (int)( (-bb-sqrt(des))/(2.0*aa) );
     }else{
       nbmax = 0;
       nbmin = 0;
     }//endif
     if(igc==0){nbmin=0;}
     for(igb=nbmin;igb<=nbmax;igb++){
       gb  = (double)igb;
       dgx = gb*hmati[2] + gc*hmati[3];
       dgy = gb*hmati[5] + gc*hmati[6];
       dgz = gb*hmati[8] + gc*hmati[9];
       dg2 = dgx*dgx+dgy*dgy+dgz*dgz;
       aa  = sumsq_a;
       bb  = 2.0*(hmati[1]*dgx+hmati[4]*dgy+hmati[7]*dgz);
       cc  = dg2-ecut_scal;
       des = bb*bb-4.0*aa*cc;
       if(des>=0.0){
         namax = (int) ( (-bb+sqrt(des))/(2.0*aa) );
         namin = (int) ( (-bb-sqrt(des))/(2.0*aa) );
       }else{
         namax = 0;
         namin = 0;
       }//endif
       if(igc==0 && igb==0){namin=0;}
       for(iga=namin;iga<=namax;iga++){iii++;}
     }// endfor : igb
   }// endfor : igc

   nktot = iii;

//=========================================================================

//----------------------------------------------------------------------- 
   }//end routine
//=======================================================================

//==========================================================================
// Count the number of lattice vectors
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

void PMEGSpaceData::fill_recip(double *hmati)

//=========================================================================
   {//begin routine 
//=========================================================================

   int iii;
   int iga,igb,igc;
   int nbmin,nbmax;
   int namin,namax;

   double ecut_scal;
   double sumsq_a,sumsq_b,sumsq_c,sdot_bc;
   double gb,gc;
   double dgx,dgy,dgz,dg2;
   double aa,bb,cc,des;

//=========================================================================
// calculate maximum k-vector

    ecut_scal = ecut_recip*0.5/(M_PI*M_PI);

    sumsq_a = hmati[1]*hmati[1]+hmati[4]*hmati[4]+hmati[7]*hmati[7];
    sumsq_b = hmati[2]*hmati[2]+hmati[5]*hmati[5]+hmati[8]*hmati[8];
    sumsq_c = hmati[3]*hmati[3]+hmati[6]*hmati[6]+hmati[9]*hmati[9];
    sdot_bc = 2.0*(hmati[2]*hmati[3]+hmati[5]*hmati[6]+hmati[8]*hmati[9]);

//=========================================================================

    iii = 0;
    for(igc=0;igc<=kc_max;igc++){
      gc  = (double)igc;
      aa  = sumsq_b;
      bb  = gc*sdot_bc;
      cc  = gc*gc*sumsq_c-ecut_scal;
      des = bb*bb-4.0*aa*cc;
      if(des>=0.0){
       nbmax = (int) ( (-bb+sqrt(des))/(2.0*aa) );
       nbmin = (int) ( (-bb-sqrt(des))/(2.0*aa) );
      }else{
       nbmax = 0;
       nbmin = 0;
      }//endif
      if(igc==0){nbmin=0;}
      for(igb=nbmin;igb<=nbmax;igb++){
        gb  = (double)igb;
        dgx = gb*hmati[2] + gc*hmati[3];
        dgy = gb*hmati[5] + gc*hmati[6];
        dgz = gb*hmati[8] + gc*hmati[9];
        dg2 = dgx*dgx+dgy*dgy+dgz*dgz;
        aa  = sumsq_a;
        bb  = 2.0*(hmati[1]*dgx+hmati[4]*dgy+hmati[7]*dgz);
        cc  = dg2-ecut_scal;
        des = bb*bb-4.0*aa*cc;
        if(des>=0.0){
         namax = (int) ( (-bb+sqrt(des))/(2.0*aa) );
         namin = (int) ( (-bb-sqrt(des))/(2.0*aa) );
        }else{
         namax = 0;
         namin = 0;
        }//endif
        if(igc==0 && igb==0){namin=0;}
        for(iga=namin;iga<=namax;iga++){
          iii++;
          ka[iii] = iga;  kb[iii] = igb; kc[iii] = igc;
        }//endfor : iga
      }// endfor : igb
    }// endfor : igc

//----------------------------------------------------------------------- 
   }//end routine
//=======================================================================

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Set the PME weight 3D  */
/*==========================================================================*/

void PMEGSpaceData::set_pme_wght(int pme_b_opt,double *bfact_r,double *bfact_i)

/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

  int* kastore = ka;
  int* kbstore = kb;
  int* kcstore = kc;
  
  int nkf1 = ngrid_a;
  int nkf2 = ngrid_b;
  int nkf3 = ngrid_c;
  
  size_t dim_k;
  int i,k_a,k_b,k_c,kap,kbp,kcp;
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
  map_a         = (int *) malloc(dim_k*sizeof(int))-1;
  map_b         = (int *) malloc(dim_k*sizeof(int))-1;
  map_c         = (int *) malloc(dim_k*sizeof(int))-1;
  for(i=1;i<=nktot;i++){
    k_a = kastore[i];
    k_b = kbstore[i];
    k_c = kcstore[i];
    if (k_c < 0) {
         kcp = k_c + nkf3 + 1;
    } else {
         kcp = k_c + 1;
    }/*endif*/
    if (k_b < 0) {
         kbp = k_b + nkf2 + 1;
    } else {
         kbp = k_b + 1;
    }/*endif*/
    if (k_a < 0) {
         kap = k_a + nkf1 + 1;
    } else {
         kap = k_a + 1;
    }/*endif*/
    map_a[i]   = kap;
    map_b[i]   = kbp;
    map_c[i]   = kcp;
  }/*endfor*/

/*==========================================================================*/
/* V) Calculate bweight on the spherically cutoff grid                      */

/*--------------------------------------------------------------------------*/
/*     A) Malloc memory and define constants                                */


   dim_k = (size_t) (ngrid_a);
   bden_a_r  = (double *)malloc(dim_k*sizeof(double))-1;
   bden_a_i  = (double *)malloc(dim_k*sizeof(double))-1;
   bweight_a = (double *)malloc(dim_k*sizeof(double))-1;

   dim_k = (size_t) (ngrid_b);
   bden_b_r  = (double *)malloc(dim_k*sizeof(double))-1;
   bden_b_i  = (double *)malloc(dim_k*sizeof(double))-1;
   bweight_b = (double *)malloc(dim_k*sizeof(double))-1;

   dim_k = (size_t) (ngrid_c);
   bden_c_r  = (double *)malloc(dim_k*sizeof(double))-1;
   bden_c_i  = (double *)malloc(dim_k*sizeof(double))-1;
   bweight_c = (double *)malloc(dim_k*sizeof(double))-1;

   dim_k = (size_t) (n_interp);
   uk   = (double *)malloc(dim_k*sizeof(double))-1;
   mn_k = (double *)malloc(dim_k*sizeof(double))-1;

/*--------------------------------------------------------------------------*/
/*     B) Construct the weighting Function                                  */

   get_bspline_wght1d(ngrid_a,mn_k,uk,
                      bden_a_r,bden_a_i,bweight_a);
   get_bspline_wght1d(ngrid_b,mn_k,uk,
                      bden_b_r,bden_b_i,bweight_b);
   get_bspline_wght1d(ngrid_c,mn_k,uk,
                      bden_c_r,bden_c_i,bweight_c);

   if(pme_b_opt > 0){
    for(i=1;i <= nktot; ++i) {
     bweight_tot[i] = bweight_a[map_a[i]]
                     *bweight_b[map_b[i]]
                     *bweight_c[map_c[i]];
    }/*endfor*/
   }/*endif*/

   if(pme_b_opt == 0 || pme_b_opt ==2){
    for(i=1;i <=nktot; ++i){
      tmp_a_r = bden_a_r[map_a[i]];
      tmp_a_i = bden_a_i[map_a[i]];
      tmp_b_r = bden_b_r[map_b[i]];
      tmp_b_i = bden_b_i[map_b[i]];
      tmp_c_r = bden_c_r[map_c[i]];
      tmp_c_i = bden_c_i[map_c[i]];
      tmp_ab_r = tmp_a_r*tmp_b_r - tmp_a_i*tmp_b_i;
      tmp_ab_i = tmp_a_i*tmp_b_r + tmp_a_r*tmp_b_i;
      bfact_r[i] = tmp_ab_r*tmp_c_r - tmp_ab_i*tmp_c_i;
      bfact_i[i] = tmp_ab_i*tmp_c_r + tmp_ab_r*tmp_c_i;
    }/*endfor*/
  }/*endif*/
  
#ifdef PINY_PME_DEBUG
  for (int i=1; i<=nktot; i++) {
    CmiPrintf ("bweight_tot[%d]=%e\n",i,bweight_tot[i]);
  }
#endif
/*========================================================================*/
/* Free memory */

  free(&map_a[1]);
  free(&map_b[1]);
  free(&map_c[1]);
  free(&bden_a_r[1]); 
  free(&bden_a_i[1]); 
  free(&bweight_a[1]);
  free(&bden_b_r[1]); 
  free(&bden_b_i[1]); 
  free(&bweight_b[1]);
  free(&bden_c_r[1]); 
  free(&bden_c_i[1]); 
  free(&bweight_c[1]);
  free(&uk[1]); 
  free(&mn_k[1]); 
/*--------------------------------------------------------------------------*/
} /* set_pme_wght */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Calculate the B spline weighting function   */
/*==========================================================================*/

void PMEGSpaceData::get_bspline_wght1d(int ngrid,
                                     double *mn_k,
                                     double *uk,
                                     double *bden_r,
                                     double *bden_i,
                                     double *bweight)

/*==========================================================================*/
/*       Begin routine */
{/*begin routine */
/*==========================================================================*/
/* Local variables */

  int k,k1,n,m;
  double arg;
  double tpi_n,mn_k_tmp;
  double bnum_real,bnum_imag;
  double bden_real,bden_imag,denom;
  double tmp_real,tmp_imag;
  double grid;

/*==========================================================================*/
/* I) Get B spline coefficients                                         */

   grid  = (double) ngrid;
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

void PMEGSpaceData::compute_recip_energy (double* hmati, StepOutput& out) {
//==========================================================================
// VI) Compute the potential energy on the spherically cutoff grid.         

   const double falp2 = 4.0*alp_ewd*alp_ewd;
   const double tpi   = 2.0*M_PI;
   const double pivol = vol/(4.0*M_PI);

   double pvten_tmp [10];
   
   for(int i=1;i<=9;i++){pvten_tmp[i]=0.0;}
   double pot_pme = 0.0;

// in parallel my part of reciprocal space
// in scalar its all of it.
   for(int i=1;i <= nktot; i++) {
     const double aka  = tpi*( (double) ka[i] );
     const double akb  = tpi*( (double) kb[i] );
     const double akc  = tpi*( (double) kc[i] );
     const double xk   = (aka*hmati[1] + akb*hmati[2] + akc*hmati[3]);
     const double yk   = (aka*hmati[4] + akb*hmati[5] + akc*hmati[6]);
     const double zk   = (aka*hmati[7] + akb*hmati[8] + akc*hmati[9]);
     const double g2   = xk*xk + yk*yk + zk*zk;
     if(0.5*g2 < ecut_recip && g2 != 0){
       const double preg = exp((-g2/falp2))/(g2*pivol);
       const double smag = (qgrid_gspace[i].getMagSqr())*bweight_tot[i];
       const double prep = -2.0*preg*smag*((g2/falp2)+1.0)/g2;
       pot_pme      += (smag*preg);
       pvten_tmp[1] += prep*xk*xk;
       pvten_tmp[5] += prep*yk*yk;
       pvten_tmp[9] += prep*zk*zk;
       pvten_tmp[2] += prep*xk*yk;
       pvten_tmp[3] += prep*xk*zk;
       pvten_tmp[6] += prep*yk*zk;
       qgrid_gspace[i].re     *= (preg*bweight_tot[i]);
       qgrid_gspace[i].im     *= (preg*bweight_tot[i]);
     }else{
       qgrid_gspace[i].re     = 0.0;
       qgrid_gspace[i].im     = 0.0;
     }//endif
   }//endfor

   pvten_tmp[4]  = pvten_tmp[2];
   pvten_tmp[7]  = pvten_tmp[3];
   pvten_tmp[8]  = pvten_tmp[6];
   pvten_tmp[1] += pot_pme;
   pvten_tmp[5] += pot_pme;
   pvten_tmp[9] += pot_pme;

   for(int i=1;i<=9;i++){
     out.pt.pvten[i-1] += pvten_tmp[i];
   }//endfor
   
   out.vpme += pot_pme;
}

PMEGSpaceData::~PMEGSpaceData () {
  clear();
}

void PMEGSpaceData::calcGridSize (double* _hmati, double _alpha_ewald,
                                  int& _ngrid_a, int& _ngrid_b, int& _ngrid_c) {
  double _ecut_recip = 2.0*M_PI*M_PI*_alpha_ewald*_alpha_ewald;
  double* hmati = _hmati-1;
  int _ka_max, _kb_max, _kc_max;
  
  double ecut_scal;
  double sumsq_a,sumsq_b,sumsq_c;

//=========================================================================
// calculate maximum k-vector

   ecut_scal = _ecut_recip*0.5/(M_PI*M_PI);

   sumsq_a = hmati[1]*hmati[1]+hmati[4]*hmati[4]+hmati[7]*hmati[7];
   _ka_max  = (int) (sqrt(ecut_scal/sumsq_a));

   sumsq_b = hmati[2]*hmati[2]+hmati[5]*hmati[5]+hmati[8]*hmati[8];
   _kb_max  = (int) (sqrt(ecut_scal/sumsq_b));

   sumsq_c = hmati[3]*hmati[3]+hmati[6]*hmati[6]+hmati[9]*hmati[9];
   _kc_max  = (int) (sqrt(ecut_scal/sumsq_c));

//=========================================================================
  
   _ngrid_a  = (int)(3.44*_ka_max); if((_ngrid_a % 2) == 1){_ngrid_a++;}
   _ngrid_b  = (int)(3.44*_kb_max); if((_ngrid_b % 2) == 1){_ngrid_b++;}
   _ngrid_c  = (int)(3.44*_kc_max); if((_ngrid_c % 2) == 1){_ngrid_c++;}
//========================================================================
}
