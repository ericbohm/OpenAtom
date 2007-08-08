/** \file para_grp_parse.C
 *
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "charm++.h"
#include "ckcomplex.h"

#include "para_grp_parse.h"
#include "../../src_piny_physics_v1.0/include/class_defs/piny_constants.h"

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Decompose the lines
//==========================================================================

void ParaGrpParse::get_chareG_line_prms(int nktot, int nchareG,int nline,
     int *npts,int *istrt_lgrp,int *iend_lgrp,int *npts_lgrp,int *nline_lgrp,
					bool statespace)

//==========================================================================
   {//begin routine
//==========================================================================
// Find the best load balancing factor

  PRINT_LINE_STAR;
  if(statespace){
      PRINTF("Statically load balancing pts and lines in state G-space : \n");
  }else{
      PRINTF("Statically load balancing pts and lines in rho G-space : \n");
  }//endif
  PRINTF("  There are %d pts %d lines and %d chunks\n",nktot,nline,nchareG);
  PRINT_LINE_DASH;printf("\n");

  if(nchareG>nline){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Dude, too much juice on the chunks. Chill on gExpandFact\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

  double dev_min  = 0.0;
  int nbal_min    = 0;
  int ibal_min    = 0;
  int ifirst      = 0;
  int mmm         = (nktot/nchareG);
  for(int ibal=0;ibal<=mmm-1;ibal++){
    int nbal  = ibal;
//  for(int ibal=0;ibal<=128;ibal++){
//    int nbal  = ((ibal*nchareG)/32);
    int ntarg = (nktot/nchareG);
    if(ntarg > nbal){ntarg -= nbal;}
    int nmax  = 0;
    int nmin  = nktot;
    int nnow  = 0;
    int ic    = 0;
    for(int i=0;i<nline;i++){
      nnow += npts[i];
      if( (nnow>=ntarg) && (ic<(nchareG-1)) ){
        ic+=1;
        nmin = MIN(nmin,nnow);
        nmax = MAX(nmax,nnow);
        nnow = 0;
      }//endif
    }//endfor
    nmin = MIN(nmin,nnow);
    nmax = MAX(nmax,nnow);
    double dev = 100.0*((double)(nmax-nmin))/((double)MAX(nmin,1));
    if(ic==nchareG-1){
     if(dev<dev_min || ifirst==0){
       ifirst   = 1;
       dev_min  = dev;
       nbal_min = nbal;
       ibal_min = ibal;
     }//endif
    }//endif
  }//endfor

  if(ifirst==0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Dude, too much juice on the chunks. Chill on gExpandFact\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

//==========================================================================
// Store the good decomposition

  int ntarg = (nktot/nchareG);
  if(ntarg > nbal_min){ntarg = ntarg-nbal_min;}

  int ic        = 0;
  npts_lgrp[0]  = 0;
  istrt_lgrp[0] = 0;
  nline_lgrp[0] = 0;
  for(int i=0;i<nline;i++){
    npts_lgrp[ic] += npts[i];
    nline_lgrp[ic]+= 1;
    if( (npts_lgrp[ic]>=ntarg) && (ic<(nchareG-1)) ){
       iend_lgrp[ic]  = i+1;
       ic+=1;
       npts_lgrp[ic]  = 0;
       nline_lgrp[ic] = 0;
       istrt_lgrp[ic] = i+1;
    }//endif
  }//endfor
  iend_lgrp[(nchareG-1)]  = nline;

//==========================================================================
// Output some statistics 

  int nmax      = 0;
  int nmin      = npts_lgrp[0];
  int nline_max = 0;
  int nline_min = nline_lgrp[0];
  for(int i=0;i<nchareG;i++){
    nmax = MAX(nmax,npts_lgrp[i]);
    nmin = MIN(nmin,npts_lgrp[i]);
    nline_max = MAX(nline_lgrp[i],nline_max);
    nline_min = MIN(nline_lgrp[i],nline_min);
  }//endfor
  double dev   = 100.0*((double)(nmax-nmin))/((double)MAX(nmin,1));
  double dev_l = 100.0*((double)(nline_max-nline_min))/((double)MAX(nline_min,1));

  if(nline_min==0 || nmin==0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Dude, too much juice on the chunks. Chill on gExpandFact\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif
  
  PRINTF("   Load balance statistics (%d %d)\n",ibal_min,nbal_min);
  PRINTF("     pt_max=%d : pt_min=%d : dev=%g\n",nmax,nmin,dev);
  PRINTF("     li_max=%d : li_min=%d : dev=%g\n",nline_max,nline_min,dev_l);

//==========================================================================
// Error check

  int ierr = 0;

  int nnn = npts_lgrp[0];
  int nc = nline_lgrp[0];
  if(istrt_lgrp[0]!=0) {ierr++;}
  for(int i=1;i<nchareG;i++){
    if(iend_lgrp[(i-1)]!=istrt_lgrp[i]){ierr++;CkPrintf("err.1\n");}
    nc += nline_lgrp[i];
    nnn += npts_lgrp[i];
  }//endfor
  if(iend_lgrp[(nchareG-1)]!=nline){ierr++;CkPrintf("err.2\n");}
  if(nc!=nline){ierr++;CkPrintf("err.3\n");}
  if(nnn!=nktot){ierr++;CkPrintf("err.4\n");}

  if(ierr!=0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Error in get_chareG_line_params\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

  printf("\n");
  PRINT_LINE_DASH;
  PRINTF("Static load balancing complete\n");
  PRINT_LINE_STAR;printf("\n");

//==========================================================================
   }
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

void ParaGrpParse::flip_data_set(int nktot, int *n_ret, int *kx, int *ky, int *kz,
                                complex *data,int iopt)

//==========================================================================
    {//begin routine 
//==========================================================================
// Count half plane kx=0 of piny data : check piny data

  int nplane0 = 0;
  for(int i=0;i<nktot;i++){
    if(kx[i]==0){nplane0++;}
  }//endfor
  int *kxt = (int *)malloc(nplane0*sizeof(int));
  int *kyt = (int *)malloc(nplane0*sizeof(int));
  int *kzt = (int *)malloc(nplane0*sizeof(int));

  complex *datat;
  if(iopt==1){datat = (complex *)malloc(nplane0*sizeof(complex));}

  for(int i=0;i<nplane0-1;i++){
    kxt[i]= kx[i];
    kyt[i]= ky[i];
    kzt[i]= kz[i];
    if(iopt==1){datat[i]= data[i];}
    if(kx[i]!=0){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Error while flipping piny dblpack data set\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif
  }//endfor
  kxt[(nplane0-1)]   = kx[(nktot-1)];
  kyt[(nplane0-1)]   = ky[(nktot-1)];
  kzt[(nplane0-1)]   = kz[(nktot-1)];
  if(iopt==1){datat[(nplane0-1)] = data[(nktot-1)];}

  if(kx[(nktot-1)]!=0 || ky[(nktot-1)] || kz[(nktot-1)]){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Error while flipping piny dblpack data set\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

//==========================================================================
// Expand the data set

  for(int i=nktot-2;i>=0;i--){
    kx[(i+nplane0)]   = kx[i];
    ky[(i+nplane0)]   = ky[i];
    kz[(i+nplane0)]   = kz[i];
    if(iopt==1){data[(i+nplane0)] = data[i];}
  }//endfor

//==========================================================================
// Create the bottom half of plane zero by symmetry : 

  int i1 = 0;
  for(int i=0;i<nplane0-1;i++){
    int ind = nplane0-i-2;
    kx[i]      =  kxt[ind];
    ky[i]      = -kyt[ind];
    kz[i]      = -kzt[ind];
    if(iopt==1){
      data[i].re =  datat[ind].re;
      data[i].im = -datat[ind].im;
    }//endif
    if(kx[i]!=0 || ky[i]<ky[i1]){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Error while flipping piny dblpack data set\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif
    i1 = i;
  }//endfor
  kx[(nplane0-1)]   = kxt[(nplane0-1)];
  ky[(nplane0-1)]   = kyt[(nplane0-1)];
  kz[(nplane0-1)]   = kzt[(nplane0-1)];
  if(iopt==1){data[(nplane0-1)] = datat[(nplane0-1)];}

  if(nktot>=nplane0+1){
    if(kx[nplane0]!= 0 || ky[nplane0]!=0 || kz[nplane0]!=1){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Error while flipping piny dblpack data set\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif
  }//endfor

#ifdef _GLENN_DBG_FLIP
  int nnn = MIN(nplane0+3,nktot);
  for(int i=0;i<nnn;i++){
    PRINTF(" %d : %d %d %d \n",i,kx[i],ky[i],kz[i]);
  }
#endif

//==========================================================================
// Exit

  (*n_ret) = (nktot+nplane0-1);

  free(kxt);
  free(kyt);
  free(kzt);

  if(iopt==1){free(datat);}

//==========================================================================
   }//end routine
//==========================================================================


