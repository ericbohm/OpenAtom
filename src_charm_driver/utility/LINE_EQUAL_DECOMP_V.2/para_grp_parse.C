#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "para_grp_parse.h"

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Decompose the lines
//==========================================================================

void ParaGrpParse::get_plane_line_prms(int nktot, int nplane,int nline,int *npts,
                int *kx_line, int *ky_line,
		int *istrt_lgrp,int *iend_lgrp,int *npts_lgrp,
 	        int *nline_lgrp,
                int *kx_str_lgrp,int *kx_end_lgrp,
                int *ky_str_lgrp,int *ky_end_lgrp)

//==========================================================================
   {
//==========================================================================
// Find the best load balancing factor

  double dev_min  = 0.0;
  int nbal_min    = 0;
  int ibal_min    = 0;
  for(int ibal=1;ibal<=32;ibal++){
    int nbal  = ((ibal*nplane)/32);
    int ntarg = (nktot/nplane);
    if(ntarg > nbal){ntarg -= nbal;}
    int nmax  = 0;
    int nmin  = nktot;
    int nnow  = 0;
    int ic    = 0;
    for(int i=0;i<nline;i++){
      nnow += npts[i];
      if( (nnow>=ntarg) && (ic<(nplane-1)) ){
        ic+=1;
        nmin = MIN(nmin,nnow);
        nmax = MAX(nmax,nnow);
        nnow = 0;
      }//endif
    }//endfor
    nmin = MIN(nmin,nnow);
    nmax = MAX(nmax,nnow);
    double dev = 100.0*((double)(nmax-nmin))/((double)MAX(nmin,1));
    if(ic==nplane-1){
     if(dev<dev_min || ibal==1){
       dev_min  = dev;
       nbal_min = nbal;
       ibal_min = ibal;
     }//endif
    }//endif
  }//endfor

//==========================================================================
// Store the good decomposition

  int ntarg = (nktot/nplane);
  if(ntarg > nbal_min){ntarg = ntarg-nbal_min;}

  int ic        = 0;
  npts_lgrp[0]  = 0;
  istrt_lgrp[0] = 0;
  nline_lgrp[0] = 0;
  for(int i=0;i<nline;i++){
    npts_lgrp[ic] += npts[i];
    nline_lgrp[ic]+= 1;
    if( (npts_lgrp[ic]>=ntarg) && (ic<(nplane-1)) ){
       iend_lgrp[ic]  = i+1;
       ic+=1;
       npts_lgrp[ic]  = 0;
       nline_lgrp[ic] = 0;
       istrt_lgrp[ic] = i+1;
    }//endif
  }//endfor
  iend_lgrp[(nplane-1)]  = nline;

  for(int i=0;i<nplane;i++){
    kx_str_lgrp[i] = kx_line[istrt_lgrp[i]];
    kx_end_lgrp[i] = kx_line[(iend_lgrp[i]-1)];
    ky_str_lgrp[i] = ky_line[istrt_lgrp[i]];
    ky_end_lgrp[i] = ky_line[(iend_lgrp[i]-1)];
  }//endfor

//==========================================================================
// Output some statistics 

  int nmax      = 0;
  int nmin      = npts_lgrp[0];
  int nline_max = 0;
  for(int i=0;i<nplane;i++){
    nmax = MAX(nmax,npts_lgrp[i]);
    nmin = MIN(nmin,npts_lgrp[i]);
    nline_max = MAX(nline_lgrp[i],nline_max);
  }//endfor
  double dev = 100.0*((double)(nmax-nmin))/((double)MAX(nmin,1));
  printf("max=%d : min=%d : dev=%g : %d %d\n",nmax,nmin,dev,nbal_min,ibal_min);

//==========================================================================
// Error check

  int ierr = 0;

  int nnn = npts_lgrp[0];
  int nc = nline_lgrp[0];
  if(istrt_lgrp[0]!=0) {ierr++;}
  for(int i=1;i<nplane;i++){
    if(iend_lgrp[(i-1)]!=istrt_lgrp[i]){ierr++;}
    nc += nline_lgrp[i];
    nnn += npts_lgrp[i];
  }//endfor
  if(iend_lgrp[(nplane-1)]!=nline){ierr++;}
  if(nc!=nline){ierr++;}
  if(nnn!=nktot){ierr++;}

  if(ierr!=0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Error in get_plane_line_params\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

//==========================================================================
   }
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

void ParaGrpParse::flip_data_set(int nktot, int *n_ret, int *kx, int *ky, int *kz
 //                                 ,complex *data
                                 )

//==========================================================================
    {//begin routine 
//==========================================================================
// Count half plane kx=0 of piny data : check piny data

  PRINTF("Flip the Complex data!!\n");

  int nplane0 = 0;
  for(int i=0;i<nktot;i++){
    if(kx[i]==0){nplane0++;}
  }//endfor
  int *kxt = (int *)malloc(nplane0*sizeof(int));
  int *kyt = (int *)malloc(nplane0*sizeof(int));
  int *kzt = (int *)malloc(nplane0*sizeof(int));
//  complex *datat = (complex *)malloc(nplane0*sizeof(complex));

  for(int i=0;i<nplane0-1;i++){
    kxt[i]= kx[i];
    kyt[i]= ky[i];
    kzt[i]= kz[i];
//    datat[i]= data[i];
    if(kx[i]!=0){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Error while flipping piny dblpack data set\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif
  }//endfor
  kxt[(nplane0-1)] = kx[(nktot-1)];
  kyt[(nplane0-1)] = ky[(nktot-1)];
  kzt[(nplane0-1)] = kz[(nktot-1)];
//    datat[(nplane0-1)]= data[(nktot-1)];
  if(kx[(nktot-1)]!=0 || ky[(nktot-1)] || kz[(nktot-1)]){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Error while flipping piny dblpack data set\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

//==========================================================================
// Expand the data set

  for(int i=nktot-2;i>=0;i--){
    kx[(i+nplane0)] = kx[i];
    ky[(i+nplane0)] = ky[i];
    kz[(i+nplane0)] = kz[i];
    //data[(i+nplane0)] =  data[i];
  }//endfor

//==========================================================================
// Create the bottom half of plane zero by symmetry : 

  int i1 = 0;
  for(int i=0;i<nplane0-1;i++){
    kx[i] =  kxt[nplane0-i-2];
    ky[i] = -kyt[nplane0-i-2];
    kz[i] = -kzt[nplane0-i-2];
//   data[i].re =  datat[i].re;
//   data[i].im = -datat[i].im;
    if(kx[i]!=0 || ky[i]<ky[i1]){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Error while flipping piny dblpack data set\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif
    i1 = i;
  }//endfor
  kx[(nplane0-1)] = kxt[(nplane0-1)];
  ky[(nplane0-1)] = kyt[(nplane0-1)];
  kz[(nplane0-1)] = kzt[(nplane0-1)];
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

//==========================================================================
   }//end routine
//==========================================================================
