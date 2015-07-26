/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_ewald                                    */
/*                                                                          */
/* This reads in and sets up the k-space                                    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../class_defs/typedefs_par.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_mdinter.h"
#include "../class_defs/allclass_mdintra.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_cp_ewald_entry.h"
#include "../proto_defs/proto_cp_ewald_corrs.h"
#include "../proto_defs/proto_friend_lib_entry.h"

#if 0 && GLENN_PERIODIC_CORRECTION
//=========================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=========================================================================================
void setput_nd_eext_corrs(int nktot, vector< gridPoints > myPoints, double *perdCorr){
  //=========================================================================================
  // Strip out data and then check to see if you belong here

  GENERAL_DATA *general_data = GENERAL_DATA::get();
#include "../class_defs/allclass_strip_gen.h"

  int ngaMax    = genewald->nka_max;
  int ngbMax    = genewald->nkb_max;
  int ngcMax    = genewald->nkc_max;
  int iperd     = gencell->iperd;
  double *hmat  = gencell->hmat;
  double *hmati = gencell->hmati;

  // No corrections in 3D
  if(iperd<0 || iperd>2){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Incorrect Periodicity in setup_nd_eext_corrs %d\n",iperd);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }//endif

  int kamax = 0;  int kcmax = 0;  int kbmax = 0;
  for(int i=0;i<nktot;i++){
    kamax = MAX(kamax,abs(myPoints[i].d1));
    kbmax = MAX(kbmax,abs(myPoints[i].d2));
    kcmax = MAX(kcmax,abs(myPoints[i].d3));
  }//endfor

  if(kamax>ngaMax || kbmax>ngbMax || kcmax>ngcMax){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Incorrect kspace in setup_nd_eext_corrs %d:%d %d:%d %d:%d \n",
        kamax,ngaMax,kbmax,ngbMax,kcmax,ngcMax);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }//endkif

  //=========================================================================================
  // Create the periodicity correction to the kernel at specfied k-vectors

  //-----------------------------------------------------------------------------------------
  // Cluster and Wire corrections computed via Gaussian quadrature numerical integration
  if(iperd!=2){
    int nquad,M;
    double wmax, wmin;
    double *anode,*weight,**data;
    create_perd_corr_data(ngaMax,ngbMax,ngcMax,&M,&nquad,&anode,&weight,&data,&wmax,&wmin);

    if(iperd==0){
      create_clus_corr(nktot,myPoints,hmat,hmati,nquad,anode,weight,data,wmax,wmin,perdCorr);
    }//endif

    if(iperd==1){
      create_wire_corr(nktot,myPoints,hmat,hmati,nquad,anode,weight,data,wmax,wmin,perdCorr);
    }//endif

    free(&anode[1]);
    free(&weight[1]);
    cfree_mat(data,1,nquad,0,M+1);    
  }//endif

  //-------------------------------------------------
  // Surface correction has a simple analytical form
  if(iperd==2){create_surf_corr(nktot,myPoints,hmat,hmati,perdCorr);}

  //=========================================================================================
}//end routine
//=========================================================================================


//=========================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=========================================================================================
void setput_nd_ewd_corrs(GENEWALD *genewald, GENCELL *gencell){
  //=========================================================================================
  // Get the parameters

  double ecut4 = 2.0*genewald->ecut; // convert to Ryd explains the 2.0
  int kamax    = genewald->nka_max;
  int kbmax    = genewald->nkb_max;
  int kcmax    = genewald->nkc_max;
  int iperd    = genewald->iperd;

  int nka_fix  = genewald->nfix_para;
  int nkb_fix  = genewald->nfix_para;
  int nkc_fix  = genewald->nfix_para;

  int kcadd    = (genewald->nfix_perp)*kcmax;
  int kbadd    = (genewald->nfix_perp)*kbmax;

  double *hmat  = gencell->hmat;
  double *hmati = gencell->hmati;

  int nktot;
  int *kastr,*kbstr,*kcstr;

  //=========================================================================================
  // Checks

  if(iperd<1 || iperd>2){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Incorrect Periodicity in setup_nd_ewd_corrs %d\n",iperd);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }//endif

  if(nka_fix>kamax || nkb_fix >kbmax){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("nka_fix %d nkb_fix %d greater than %d %d %d\n",nka_fix,nkb_fix,kamax,kbmax);
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    FFLUSH(stdout);
    EXIT(1);
  }//endif

  if(iperd==1){
    if(nkc_fix>kcmax){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("nkc_fix %d greater than %d\n",nkc_fix,kcmax);
      PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      FFLUSH(stdout);
      EXIT(1);
    }//endif
  }//endif

  //=========================================================================================
  // Set up k-space

  set_large_kvectors(ecut4,hmati,&kastr,&kbstr,&kcstr,&nktot,kamax,kbmax,kcmax);

  //=========================================================================================
  // Find k-space regions where you need to add a little more juice

  //--------------------------------------------------------
  // When ka<= nka_max and |kb|<=nkb_max find kcmax and kcmin
  // so more kc can be added

  int **kcmax_str = cmall_int_mat(0,nka_fix,0,2*nkb_fix,"set_perd_corrs.C");
  int **kcmin_str = cmall_int_mat(0,nka_fix,0,2*nkb_fix,"set_perd_corrs.C");
  for(int i=0;i<nktot;i++){
    int ka    = kastr[i];
    int kb    = abs(kbstr[i]);
    if(ka<=nka_fix && kb<=nkb_fix){
      int kboff = kbstr[i]+nkb_fix;
      kcmax_str[ka][kboff] = MAX(kcmax_str[ka][kboff],kcstr[i]);
      kcmin_str[ka][kboff] = MIN(kcmin_str[ka][kboff],kcstr[i]);
    }//endif
  }//endfor

  //--------------------------------------------------------
  // When ka<= nka_max and |kc|<=nkc_max find kbmax and kbmin
  // so more kb can be added

  int **kbmax_str = NULL;
  int **kbmin_str = NULL;
  if(iperd==1){
    kbmax_str = cmall_int_mat(0,nka_fix,0,2*nkc_fix,"set_perd_corrs.C");
    kbmin_str = cmall_int_mat(0,nka_fix,0,2*nkc_fix,"set_perd_corrs.C");
    for(int i=0;i<nktot;i++){
      int ka    = kastr[i];
      int kc    = abs(kcstr[i]);
      if(ka<=nka_fix && kc<=nkc_fix){
        int kcoff = kcstr[i]+nkc_fix;
        kbmax_str[ka][kcoff] = MAX(kbmax_str[ka][kcoff],kbstr[i]);
        kbmin_str[ka][kcoff] = MIN(kbmin_str[ka][kcoff],kbstr[i]);
      }//endif
    }//endfor
  }//endif

  //=======================================================================
  // Add kc when ka and kb are small

  // Count the corrections along c
  int ic=0;
  for(int ka=0;ka<=nka_fix;ka++){
    int kbmin = (ka==0 ? 0 : -nkb_fix);
    for(int kb=kbmin;kb<=nkb_fix;kb++){
      int kboff = kb+nkb_fix;
      int jj  = (ka==0 && kb==0 ? 1 : 2);
      for(int ii=1;ii<=jj;ii++){
        int kcmin,kcmax;
        if(ii==1){kcmin=kcmax_str[ka][kboff];kcmax=kcmin+kcadd;}
        if(ii==2){kcmax=kcmin_str[ka][kboff];kcmin=kcmax-kcadd;}
        for(int kc=kcmin;kc<=kcmax;kc++){ic++;}
      }//endfor
    }//endfor
  }//endfor
  int ncorr_c = ic;

  // Store the corrections along c
  int *ka_corr_c = (int *)cmalloc(ncorr_c*sizeof(int),"set_perd_corrs.C");
  int *kb_corr_c = (int *)cmalloc(ncorr_c*sizeof(int),"set_perd_corrs.C");
  int *kc_corr_c = (int *)cmalloc(ncorr_c*sizeof(int),"set_perd_corrs.C");
  double *kernel_corr_c = (double *)cmalloc(ncorr_c*sizeof(double),"set_perd_corrs.C");
  ic = 0;
  for(int ka=0;ka<=nka_fix;ka++){
    int kbmin    = (ka==0 ? 0 : -nkb_fix);
    for(int kb=kbmin;kb<=nkb_fix;kb++){
      int kboff = kb+nkb_fix;
      int jj  = (ka==0 && kb==0 ? 1 : 2);
      for(int ii=1;ii<=jj;ii++){
        int kcmin,kcmax;
        if(ii==1){kcmin=kcmax_str[ka][kboff];kcmax=kcmin+kcadd;}
        if(ii==2){kcmax=kcmin_str[ka][kboff];kcmin=kcmax-kcadd;}
        for(int kc=kcmin;kc<=kcmax;kc++){
          ka_corr_c[ic] = ka;
          kb_corr_c[ic] = kb;
          kc_corr_c[ic] = kc;
          ic++;
        }//endfor
      }//endfor
    }//endfor
  }//endfor

  //============================================================================
  // kb corrections : note kb is never 0, so kc range is unrestricted

  int ncorr_b           = 0;
  int *ka_corr_b        = NULL;
  int *kb_corr_b        = NULL;
  int *kc_corr_b        = NULL;
  double *kernel_corr_b = NULL;

  if(iperd==1){
    // Count the corrections along b
    int ic=0;
    for(int ka=0;ka<=nka_fix;ka++){
      int jj  = (ka==0 ? 1 : 2);
      for(int kc=-nkc_fix;kc<=nkc_fix;kc++){
        int kcoff = kc+nkc_fix;
        for(int ii=1;ii<=jj;ii++){
          int kbmin,kbmax;
          if(ii==1){kbmin=kbmax_str[ka][kcoff];kbmax=kbmin+kbadd;}
          if(ii==2){kbmax=kbmin_str[ka][kcoff];kbmin=kbmax-kbadd;}
          for(int kb=kbmin;kb<=kbmax;kb++){ic++;}
        }//endfor
      }//endfor
    }//endfor
    ncorr_b = ic;

    // Store the corrections along b
    ka_corr_b = (int *)cmalloc(ncorr_b*sizeof(int),"set_perd_corrs.C");
    kb_corr_b = (int *)cmalloc(ncorr_b*sizeof(int),"set_perd_corrs.C");
    kc_corr_b = (int *)cmalloc(ncorr_b*sizeof(int),"set_perd_corrs.C");
    kernel_corr_b = (double *)cmalloc(ncorr_b*sizeof(double),"set_perd_corrs.C");
    ic=0;
    for(int ka=0;ka<=nka_fix;ka++){
      int jj  = (ka==0 ? 1 : 2);
      for(int kc=-nkc_fix;kc<=nkc_fix;kc++){
        int kcoff = kc+nkc_fix;
        for(int ii=1;ii<=jj;ii++){
          int kbmin,kbmax;
          if(ii==1){kbmin=kbmax_str[ka][kcoff];kbmax=kbmin+kbadd;}
          if(ii==2){kbmax=kbmin_str[ka][kcoff];kbmin=kbmax-kbadd;}
          for(int kb=kbmin;kb<=kbmax;kb++){
            ka_corr_b[ic] = ka;
            kb_corr_b[ic] = kb;
            kc_corr_b[ic] = kc;
            ic++;
          }//endfor :kb
        }//endfor : plus/minus
      }//endfor : kc
    }//endfor : ka

  }//endif

  //============================================================================
  // Compute kernel fixes

  // In 1D 
  if(iperd==1){
    int ngaMax = nka_fix;
    int ngbMax = kbmax+kbadd;
    int ngcMax = kcmax+kcadd;
    int nquad,M;
    double wmax, wmin;
    double *anode,*weight,**data;
    create_perd_corr_data(ngaMax,ngbMax,ngcMax,&M,&nquad,&anode,&weight,&data,
        &wmax,&wmin);
    create_wire_corr(ncorr_c,ka_corr_c,kb_corr_c,kc_corr_c,hmat,hmati,nquad,anode,weight,data,
        wmax,wmin,kernel_corr_c);
    create_wire_corr(ncorr_b,ka_corr_b,kb_corr_b,kc_corr_b,hmat,hmati,nquad,anode,weight,data,
        wmax,wmin,kernel_corr_b);
    free(&anode[1]);
    free(&weight[1]);
    cfree_mat(data,1,nquad,0,M+1);    
  }//endif

  // In 2D 
  if(iperd==2){
    create_surf_corr(ncorr_c,ka_corr_c,kb_corr_c,kc_corr_c,hmat,hmati,kernel_corr_c);
  }//endif

  //============================================================================
  // Clean the memory; Store the goodies

  free(kastr);
  free(kbstr);
  free(kcstr);

  cfree_int_mat(kcmax_str,0,nka_fix,0,2*nkb_fix);
  cfree_int_mat(kcmin_str,0,nka_fix,0,2*nkb_fix);
  if(iperd==1){
    cfree_int_mat(kbmax_str,0,nka_fix,0,2*nkc_fix);
    cfree_int_mat(kbmin_str,0,nka_fix,0,2*nkc_fix);
  }//endif

  genewald->ncorr_c       = ncorr_c;
  genewald->ncorr_b       = ncorr_b;
  genewald->ka_corr_c     = ka_corr_c;
  genewald->kb_corr_c     = kb_corr_c;
  genewald->kc_corr_c     = kc_corr_c;
  genewald->ka_corr_b     = ka_corr_b;
  genewald->kb_corr_b     = kb_corr_b;
  genewald->kc_corr_b     = kc_corr_b;
  genewald->kernel_corr_c = kernel_corr_c;
  genewald->kernel_corr_b = kernel_corr_b;

  //============================================================================
}//end routine
//============================================================================


//============================================================================
// Generate the kernel for Poisson's equations for clusters
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void create_clus_corr(int nktot, vector< gridPoints> & myPoints, double *hmat, double *hmati, 
    int nquad, double *anode, double *weight, double **data,
    double wmax, double wmin,double *kernel_corr){
  //============================================================================
  // Check the box

  if(hmat[5]!=hmat[9] || hmat[1]!=hmat[9] || hmat[2]!=0 || hmat[3]!=0 || 
      hmat[4]!=0 || hmat[6]!=0 || hmat[7]!=0 || hmat[8]!=0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("The box for clusters must be cubic with Lx=Ly=Lz\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

  //============================================================================
  // Precomputed g=0 term

  double L  = hmat[1];
  double L2 = L*L;

  double gzero = 0.7576021548; // well defined : computed off line numerically
  gzero       *= (L2*M_PI);   // \int_0^infty [erf(0.5*t^{-1/2})]^3

  //============================================================================
  // Create the kernel

  double pre = L2/sqrt(M_PI);
  for(int i=0;i<nktot;i++){
    double gx = 2.0*M_PI*((double) myPoints[i].d1);
    double gy = 2.0*M_PI*((double) myPoints[i].d2);
    double gz = 2.0*M_PI*((double) myPoints[i].d3);
    double g2 = gx*gx+gy*gy+gz*gz;
    double g  = sqrt(g2);
    int ind_x = abs(myPoints[i].d1);
    int ind_y = abs(myPoints[i].d2);
    int ind_z = abs(myPoints[i].d3);
    if(g==0.0){kernel_corr[i]=gzero;}
    if(g>0.0){
      kernel_corr[i]=0.0;
      for(int k=1;k<=nquad;k++){
        double x   = anode[k]*sqrt(anode[k]); // missing L^3 factors out with data[]^3
        kernel_corr[i] += ( (weight[k]/x)*data[k][ind_x]*data[k][ind_y]*data[k][ind_z] );
      }//endfor
      double gx2 = gx*gx; double gy2 = gy*gy; double gz2 = gz*gz;
      double wmax32 = sqrt(wmax)*wmax;
      double wmax52 = wmax32*wmax;
      double wmax72 = wmax52*wmax;
      //--------------------------------
      // 1 guy is non-zero 
      if(myPoints[i].d1!=0&&myPoints[i].d2==0&&myPoints[i].d3==0){
        kernel_corr[i] -= (2.0*cos(gx*0.5)/(gx2*wmax32))
          *(2.0/3.0+(12.0/(5.0*gx2)-1.0/6.0)/wmax);
      }//endif
      if(myPoints[i].d1==0&&myPoints[i].d2!=0&&myPoints[i].d3==0){
        kernel_corr[i] -= (2.0*cos(gy*0.5)/(gy2*wmax32))
          *(2.0/3.0+(12.0/(5.0*gy2)-1.0/6.0)/wmax);
      }//endif
      if(myPoints[i].d1==0&&myPoints[i].d2==0&&myPoints[i].d3!=0){
        kernel_corr[i] -= (2.0*cos(gz*0.5)/(gz2*wmax32))
          *(2.0/3.0+(12.0/(5.0*gz2)-1.0/6.0)/wmax);
      }//endif
      //--------------------------------
      // 1 guy is zero
      if(myPoints[i].d1==0&&myPoints[i].d2!=0&&myPoints[i].d3!=0){
        kernel_corr[i] += (8.0/5.0)*cos(gy*0.5)*cos(gz*0.5)/(gy2*gz2*wmax52);
      }//endif
      if(myPoints[i].d1!=0&&myPoints[i].d2==0&&myPoints[i].d3!=0){
        kernel_corr[i] += (8.0/5.0)*cos(gx*0.5)*cos(gz*0.5)/(gx2*gz2*wmax52);
      }//endif
      if(myPoints[i].d1!=0&&myPoints[i].d2!=0&&myPoints[i].d3==0){
        kernel_corr[i] += (8.0/5.0)*cos(gx*0.5)*cos(gy*0.5)/(gx2*gy2*wmax52);
      }//endif
      //--------------------------------
      // No one is zero
      if(myPoints[i].d1!=0&&myPoints[i].d2!=0&&myPoints[i].d3!=0){
        kernel_corr[i] -= (16.0/7.0)*cos(gx*0.5)*cos(gy*0.5)*cos(gz*0.5)/(gx2*gy2*gz2*wmax72);
      }//endif
      //--------------------------------
      // Complete kernel
      kernel_corr[i] *=pre;
      kernel_corr[i] -= (4.0*M_PI*L2/g2)*exp(-g2*0.25*wmin);
    }//endif
  }//endfor

  //----------------------------------------------------------------------------
}//end routine
//============================================================================


//============================================================================
// Generate the kernel for Poisson's equations for wires
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void create_wire_corr(int nktot, vector< gridPoints> &myPoints, double *hmat, double *hmati, 
    int nquad, double *anode, double *weight,double **data,
    double wmax, double wmin,double *kernel_corr){
  //============================================================================
  // Check the box

  if(hmat[5]!=hmat[9] || hmat[2]!=0 || hmat[3]!=0 || 
      hmat[4]!=0 || hmat[6]!=0 || hmat[7]!=0 || hmat[8]!=0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("The box for wires must be orthorhombic with Ly=Lz\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

  if(2.0*M_PI*hmat[5]/hmat[1]<1 || 2.0*M_PI*hmat[9]/hmat[1]<1){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("The box for wires has too small an aspect ratio.\n");
    PRINTF("Please make the faux b/c-axises longer\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

  //============================================================================
  // Precomputed g=0 term

  double L  = hmat[5];
  double L2 = L*L;

  double gzero = 0.491831811; // has a logarithmic sigularity subtracted out
  gzero       *= (L2*M_PI);   //  \int_0^B dt [erf(0.5/sqrt(t))]^2 - log(B)/pi lim B->infty
  // computed off-line numerically

  //============================================================================
  // Create the kernel

  double rat = L/hmat[1];
  double pre = L2;
  for(int i=0;i<nktot;i++){
    double gx = 2.0*M_PI*rat*((double) myPoints[i].d1);
    double gy = 2.0*M_PI*((double) myPoints[i].d2);
    double gz = 2.0*M_PI*((double) myPoints[i].d3);
    double g2 = gx*gx+gy*gy+gz*gz;
    double g  = sqrt(g2);
    int ind_y = abs(myPoints[i].d2);
    int ind_z = abs(myPoints[i].d3);
    if(g==0.0){kernel_corr[i]=gzero;}
    if(fabs(gx)<200.0 && g>0.0){
      kernel_corr[i]=0.0;
      for(int k=1;k<=nquad;k++){
        double x = anode[k]; // missing L^2 factors out with data[]^2
        kernel_corr[i] += ((weight[k]/x)*exp(-gx*gx*anode[k]*0.25)*
            data[k][ind_y]*data[k][ind_z]);
      }//endfor
      //---------------------------------------------------
      // long range correction in the absense of gx damping (gx==0).
      if(myPoints[i].d1==0){
        double gz2 = gz*gz; double gy2 = gy*gy;
        //----------------------------------------
        // 1 guy = 0 1 guy !=0
        if(myPoints[i].d2==0 && myPoints[i].d3!=0){
          kernel_corr[i] += cos(0.5*gz)/(gz2*wmax)*( (1.0/3.0-6.0/gz2)/wmax-2.0 );
        }//endif
        if(myPoints[i].d2!=0 && myPoints[i].d3==0){
          kernel_corr[i] += cos(0.5*gy)/(gy2*wmax)*( (1.0/3.0-6.0/gy2)/wmax-2.0 );
        }//endif
        //----------------------------------------
        // both non-zero
        if(myPoints[i].d2!=0 && myPoints[i].d3!=0){
          kernel_corr[i] += 2.0*cos(0.5*gy)*cos(0.5*gz)/(gy2*gz2*wmax*wmax);
        }//endif
      }//endif
      //--------------------------------
      // Complete kernel
      kernel_corr[i] *=pre;
      kernel_corr[i] -= (4.0*M_PI*L2/g2)*exp(-g2*0.25*wmin);
    }//endif
    if(fabs(gx)>=200){kernel_corr[i]=0.0;}
  }//endfor

  //----------------------------------------------------------------------------
}//end routine
//============================================================================


//============================================================================
// Generate the kernel for Poisson's equations for surfaces
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void create_surf_corr(int nktot, vector< gridPoints> & myPoints, double *hmat, double *hmati, 
    double *kernel_corr){
  //============================================================================
  // Check the box

  if(hmat[3]!=0.0 || hmat[6]!=0 || hmat[7]!=0.0 || hmat[8]!=0.0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("The box for surfaces must have the c-axis perp to the surface\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

  //============================================================================
  // Compute the kernel correction analytically

  double gzero = -0.5*M_PI*hmat[9]*hmat[9];  // analytically with 
  // singularity taken out
  for(int i=0;i<nktot;i++){
    double aka = 2.0*M_PI*( (double) myPoints[i].d1 );
    double akb = 2.0*M_PI*( (double) myPoints[i].d2 );
    double akc = 2.0*M_PI*( (double) myPoints[i].d3 );
    double xk  = aka*hmati[1] + akb*hmati[2];
    double yk  = aka*hmati[4] + akb*hmati[5];
    double zk  = akc*hmati[9];
    double g2  = xk*xk + yk*yk + zk*zk;
    double gs  = sqrt(xk*xk+yk*yk);
    if(myPoints[i].d1==0 && myPoints[i].d2==0 && myPoints[i].d3==0){
      kernel_corr[i] = gzero;
    }else{
      kernel_corr[i] = -(4.0*M_PI/g2)*cos(zk*hmat[9]*0.5)*exp(-gs*hmat[9]*0.5);
    }//endif
  }//endfor

  //----------------------------------------------------------------------------
}//end routine
//============================================================================


//============================================================================
// Generate the kernel for Poisson's equations for surfaces
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void create_surf_corr_dummy(int nktot, vector< gridPoints > & myPoints, double *hmat, double *hmati, 
    int nquad, double *anode, double *weight,double **data,
    double wmax, double wmin,double *kernel_corr){
  //============================================================================
  // Check the box

  if(hmat[3]!=0.0 || hmat[6]!=0 || hmat[7]!=0.0 || hmat[8]!=0.0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("The box for surfaces must have the c-axis perp to the surface\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

  if(2.0*M_PI*hmati[1]*hmat[9]< 1 || 2.0*M_PI*hmati[1]*hmat[9]< 1){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("The box for surfaces has too small an aspect ratio.\n");
    PRINTF("Please make the faux c-axis longer\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

  //============================================================================
  // Precomputed g=0 term

  double L  = hmat[9];
  double L2 = L*L;

  double gzero = -0.5*M_PI*L2;  // analytically with singularity taken out

  //============================================================================
  // Create the kernel

  double pre = L2*sqrt(M_PI);
  for(int i=0;i<nktot;i++){
    double aka = 2.0*M_PI*( (double) myPoints[i].d1 );
    double akb = 2.0*M_PI*( (double) myPoints[i].d2 );
    double akc = 2.0*M_PI*( (double) myPoints[i].d3 );
    double xk  = (aka*hmati[1] + akb*hmati[2])*hmat[9];
    double yk  = (aka*hmati[4] + akb*hmati[5])*hmat[9];
    double zk  = akc;
    double g2  = (xk*xk + yk*yk + zk*zk);
    double gs  = sqrt(xk*xk+yk*yk);
    double g   = sqrt(g2);
    int ind_z  = abs(myPoints[i].d3);
    if(g==0.0){kernel_corr[i]=gzero;}
    if(gs<200.0 && g>0.0){
      kernel_corr[i]=0.0;
      for(int k=1;k<=nquad;k++){
        double x = sqrt(anode[k]); // missing L factors out with data[]
        kernel_corr[i] += ((weight[k]/x)*exp(-gs*gs*anode[k]*0.25)*
            data[k][ind_z]);
      }//endfor
      if(myPoints[i].d1==0 && myPoints[i].d2==0){
        kernel_corr[i] -= 4.0*cos(zk*0.5)/(zk*zk*sqrt(wmax));
        kernel_corr[i] += (1.0/3.0)*cos(zk*0.5)*(1.0/(zk*zk)-24.0/(zk*zk*zk*zk))
          /(sqrt(wmax)*wmax);
      }//endif
      kernel_corr[i] *=pre;
      kernel_corr[i] -= (4.0*M_PI*L2/g2)*exp(-g2*0.25*wmin);
    }//endif
    if(gs>=200.0){kernel_corr[i]=0.0;}
  }//endfor

  //----------------------------------------------------------------------------
}//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void create_perd_corr_data(int ngxMax,int ngyMax, int ngzMax,int *M_ret,
    int *nquad_ret, double **anode_ret, double **weight_ret,
    double ***data_ret, double *wmax_ret, double *wmin_ret){
  //============================================================================
  // Collect the weights and nodes

  int nquad       =  2048;
  double *anode   = (double *)cmalloc(nquad*sizeof(double),"set_perd_corrs.C")-1;
  double *weight  = (double *)cmalloc(nquad*sizeof(double),"set_perd_corrs.C")-1;

#include "../proto_defs/gauss_2048.h"

  // 0 < |gs| < 200
  double wmax  = 80;
  double wmin  = 0.01;
  double scale = (wmax-wmin)/2.0;
  double shift = (wmax+wmin)/2.0;

  for(int k=1;k<=nquad;k++){
    anode[k]   = scale*anode[k]+shift;
    weight[k]  = scale*weight[k];
  }//endfor

  //============================================================================
  // Compute nquad ffts for low g and high g stuff

  int M = MAX(ngxMax,ngyMax);
  M = MAX(M,ngzMax);

  int N = 8192;               //Make a nice quadrature
  N     = MAX(N,8*ngxMax);
  N     = MAX(N,8*ngyMax);
  N     = MAX(N,8*ngzMax);

  double *workfft  = (double *)cmalloc(8*N*sizeof(double),"set_perd_corrs.C")-1;
  double *datafft  = (double *)cmalloc(2*N*sizeof(double),"set_perd_corrs.C")-1;
  double **data    = cmall_mat(1,nquad,0,M+1,"set_perd_corrs.C");
  DCFFTI_GENERIC(&N,&workfft[1]);

  double dl = 1.0/(double) N;
  double dN = 1.0/((double)N);  // L factors out (see below)
  for(int k=1;k<=nquad;k++){
    for(int i=0,j=1;i<N;i++,j+=2){
      double x       = dN*((double)i)-0.5;
      datafft[j]     = exp(-x*x/anode[k]);
      datafft[(j+1)] = 0.0;
    }//endfor
    DCFFTF_GENERIC(&N,&datafft[1],&workfft[1]);
    for(int i=0,j=1;i<=M;i++,j+=2){
      double g   = 2.0*M_PI*((double) i);
      double ss  = dl*cos(0.5*g); //missing L in dl factors out in computation
      data[k][i] = datafft[j]*ss;
    }//endfor
  }//endfor

  //============================================================================
  // Clean up and return 

  free(&datafft[1]);
  free(&workfft[1]);

  *M_ret      = M;
  *nquad_ret  = nquad;
  *anode_ret  = anode;
  *weight_ret = weight;
  *data_ret   = data;
  *wmax_ret   = wmax;
  *wmin_ret   = wmin;

  //============================================================================
}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void set_large_kvectors(double ecut4, double *hmati, int **kx_ret, int **ky_ret, 
    int **kz_ret, int *nPacked_ret,
    int ka_max, int kb_max, int kc_max)

  //============================================================================
{//begin routine
  //============================================================================
  // count the k-vectors : preserve nice symmetry for non-cubic boxes even
  //                       though you need more kvectors   

  int iii;
  double tpi  = 2.0*M_PI;

  int nPacked = 0;
  for(int ka=0;ka<=ka_max;ka++){
    int kb_strt = -kb_max;
    if(ka==0)kb_strt = 0;
    for(int kb=kb_strt;kb<=kb_max;kb++){
      int kc_strt = -kc_max;
      if(ka==0 && kb==0)kc_strt = 0;
      for(int kc=kc_strt;kc<=kc_max;kc++){
        double gx = tpi*(ka*hmati[1] + kb*hmati[2] + kc*hmati[3]);
        double gy = tpi*(ka*hmati[4] + kb*hmati[5] + kc*hmati[6]);
        double gz = tpi*(ka*hmati[7] + kb*hmati[8] + kc*hmati[9]);
        double g2 = gx*gx+gy*gy+gz*gz;
        if(g2<=ecut4){nPacked++;}
      }//endfor:kc
    }//endfor:kb
  }//endfor:ka

  //============================================================================
  // fill the k-vectors

  int *kx = (int *)cmalloc(nPacked*sizeof(int),"set_perd_corrs.C");
  int *ky = (int *)cmalloc(nPacked*sizeof(int),"set_perd_corrs.C");
  int *kz = (int *)cmalloc(nPacked*sizeof(int),"set_perd_corrs.C");

  int ic = 0;
  for(int ka=0;ka<=ka_max;ka++){
    int kb_strt = -kb_max;
    if(ka==0)kb_strt = 0;
    for(int kb=kb_strt;kb<=kb_max;kb++){
      int kc_strt = -kc_max;
      if(ka==0 && kb==0)kc_strt = 0;
      for(int kc=kc_strt;kc<=kc_max;kc++){
        double gx = tpi*(ka*hmati[1] + kb*hmati[2] + kc*hmati[3]);
        double gy = tpi*(ka*hmati[4] + kb*hmati[5] + kc*hmati[6]);
        double gz = tpi*(ka*hmati[7] + kb*hmati[8] + kc*hmati[9]);
        double g2 = gx*gx+gy*gy+gz*gz;
        if(g2<=ecut4){
          kx[ic]=ka;
          ky[ic]=kb;
          kz[ic]=kc;
          ic++;
        }/*endif*/
      }//endfor:kc
    }//endfor:kb
  }//endfor:ka

  //============================================================================
  // Set return values

  *kx_ret        = kx;
  *ky_ret        = ky;
  *kz_ret        = kz;
  *nPacked_ret   = nPacked;

  //============================================================================
}//end routine
//============================================================================

#endif //GLENN_PERIODIC_CORRECTION

