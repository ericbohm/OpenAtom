/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** \file test_para_grp_parse.C
 *
 */
//==============================================================================



#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "para_grp_parse.h"

void create_kvec(int **,int **,int **, int ,int *);

void get_glenn_prms(int ,int *,int *,int *,int *,int **,int **,
                    int **,int **,int **);

void test_flip(int,int, int *,int *,int*);

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
   int main(){
//==========================================================================

   int nktot;
   int nplane;
   int *kx,*ky,*kz;

//==========================================================================
// create the model data set

   PRINTF("\n");
   PRINTF("====================================================\n");
   PRINTF("          Input data                                \n");
   PRINTF("----------------------------------------------------\n");
   PRINTF("Enter nplane : ");
   scanf("%d",&nplane);

   int ntemp;
   create_kvec(&kx,&ky,&kz,nplane,&ntemp);
   ParaGrpParse::flip_data_set(ntemp,&nktot,kx,ky,kz);
   test_flip(nktot,nplane,kx,ky,kz);

   PRINTF("      nktot  : %d\n",nktot);

   PRINTF("====================================================\n");

//==========================================================================
// Parse the data out by plane index

   PRINTF("\n");
   PRINTF("====================================================\n");
   PRINTF("          Line decomposition                        \n");
   PRINTF("----------------------------------------------------\n");

   int nline;
   int *istrt;
   int *iend;
   int *npts;
   int *kx_line;
   int *ky_line;

   get_glenn_prms(nktot,kx,ky,&nplane,&nline,&istrt,&iend,&npts,&kx_line,&ky_line);
   PRINTF("nplane %d nline %d\n",nplane,nline);

   int *istrt_lgrp   = (int *)malloc(nplane*sizeof(int));
   int *iend_lgrp    = (int *)malloc(nplane*sizeof(int));
   int *npts_lgrp    = (int *)malloc(nplane*sizeof(int));
   int *nline_lgrp   = (int *)malloc(nplane*sizeof(int));
   int *kx_str_lgrp  = (int *)malloc(nplane*sizeof(int));
   int *kx_end_lgrp  = (int *)malloc(nplane*sizeof(int));
   int *ky_str_lgrp  = (int *)malloc(nplane*sizeof(int));
   int *ky_end_lgrp  = (int *)malloc(nplane*sizeof(int));

   ParaGrpParse::get_plane_line_prms(nktot,nplane,nline,npts,kx_line, ky_line,
                    istrt_lgrp,iend_lgrp,npts_lgrp,nline_lgrp,
                    kx_str_lgrp,kx_end_lgrp,ky_str_lgrp,ky_end_lgrp);

   for(int i=0;i<nplane;i++){
     PRINTF("i=%d : n=%d istrt=%d iend=%d nline=%d ",
             i,npts_lgrp[i],istrt_lgrp[i],iend_lgrp[i],nline_lgrp[i]);
     PRINTF("kx_str %d kx_end %d ky_str %d ky_end %d\n",
             kx_str_lgrp[i],kx_end_lgrp[i],ky_str_lgrp[i],ky_end_lgrp[i]);  
   }//endfor

   PRINTF("====================================================\n");

   return 1;
//-------------------------------------------------------------------------
   }//end routine 
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Create model data points in 3D with spherical truncation
//==========================================================================

void create_kvec(int **kx_out,int **ky_out,int **kz_out,
                 int nplane,int *nktot_out)

//==========================================================================
      {//begin routine
//==========================================================================

  int nktot = 0;
  double aplane1 = (double)(nplane-1);
  for(int ikx=0;ikx<=nplane;ikx++){
    int kymin = -nplane;
    if(ikx==0)kymin=0;
    for(int iky=kymin;iky<=nplane;iky++){
      int kzmin = -nplane;
      if(ikx==0 && iky==0)kzmin=1;
      for(int ikz=kzmin;ikz<=nplane;ikz++){
        double aka = (double)ikx;
        double akb = (double)iky;
        double akc = (double)ikz;
        double g   = sqrt(aka*aka + akb*akb + akc*akc);
        if(g<=aplane1){nktot++;}
      }//endfor
    }//endfor
  }//endfor
  nktot++;

  int *kx = (int *)malloc(2*nktot*sizeof(int));
  int *ky = (int *)malloc(2*nktot*sizeof(int));
  int *kz = (int *)malloc(2*nktot*sizeof(int));

  int ic = 0;
  for(int ikx=0;ikx<=nplane;ikx++){
    int kymin = -nplane;
    if(ikx==0)kymin=0;
    for(int iky=kymin;iky<=nplane;iky++){
      int kzmin = -nplane;
      if(ikx==0 && iky==0)kzmin=1;
      for(int ikz=kzmin;ikz<=nplane;ikz++){
        double aka = (double)ikx;
        double akb = (double)iky;
        double akc = (double)ikz;
        double g   = sqrt(aka*aka + akb*akb + akc*akc);
        if(g<=aplane1){kx[ic]=ikx; ky[ic]=iky; kz[ic]=ikz; ic++;}
      }//endfor
    }//endfor
  }//endfor
  kx[ic]=0;ky[ic]=0;kz[ic]=0;
  ic++;

  *nktot_out = nktot;
  *kx_out    = kx;
  *ky_out    = ky;
  *kz_out    = kz;

//-------------------------------------------------------------------------
   }//end routine 
//==========================================================================



//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Decompose the nktot data points into chunks based on plane index
//==========================================================================

void get_glenn_prms(int nktot,int *kx,int *ky,int *nplane_ret,
                    int *nline_ret,int **istrt_ret,int **iend_ret,
                    int **npts_ret,int **kx_line_ret,int **ky_line_ret)

//==========================================================================
   { //begin routine
//==========================================================================
// count the planes
  
   int nplane = 1;
   int nline  = 1;
   for(int i=1;i<nktot;i++){
     if(kx[i]!=kx[(i-1)]){nplane++;}
     if(ky[i]!=ky[(i-1)] || kx[i]!=kx[(i-1)]){nline++;}
   }

   int *iend    = (int *)malloc(nline*sizeof(int));
   int *istrt   = (int *)malloc(nline*sizeof(int));
   int *npts    = (int *)malloc(nline*sizeof(int));
   int *kx_line = (int *)malloc(nline*sizeof(int));
   int *ky_line = (int *)malloc(nline*sizeof(int));

   int ic     = 0;
   npts[0]    = 1;
   istrt[0]   = 0;
   iend[0]    = 1;
   kx_line[0] = kx[0];
   ky_line[0] = ky[0];
   for(int i=1;i<nktot;i++){
    if(kx[i]!=kx[(i-1)] || ky[i]!=ky[(i-1)]){
       iend[ic]  = i;
       ic       += 1;
       npts[ic]  = 0;
       istrt[ic] = i;
       kx_line[ic] = kx[i];
       ky_line[ic] = ky[i];
    }//endif
    npts[ic]++;
  }//endfor
  iend[(nline-1)]=nktot;

  *nplane_ret  = nplane;
  *nline_ret   = nline;
  *istrt_ret   = istrt;
  *iend_ret    = iend;
  *npts_ret    = npts;
  *kx_line_ret = kx_line;
  *ky_line_ret = ky_line;

//-------------------------------------------------------------------------
   }//end routine 
//==========================================================================



//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Create model data points in 3D with spherical truncation
//==========================================================================

void test_flip(int nktot, int nplane, int *kx,int *ky,int *kz)

//==========================================================================
      {//begin routine
//==========================================================================

  double aplane1 = (double)(nplane-1);

  int ic = 0;
  for(int ikx=0;ikx<=nplane;ikx++){
    int kymin = -nplane;
    for(int iky=kymin;iky<=nplane;iky++){
      int kzmin = -nplane;
      for(int ikz=kzmin;ikz<=nplane;ikz++){
        double aka = (double)ikx;
        double akb = (double)iky;
        double akc = (double)ikz;
        double g   = sqrt(aka*aka + akb*akb + akc*akc);
        if(g<=aplane1){
          if(ic>=nktot){
            printf("flip error 1\n");exit(1);
	  }//endif
          if(kx[ic]!=ikx || ky[ic]!=iky || kz[ic]!=ikz){
            printf("flip error 2\n");exit(1);
	  }//endif
          ic++;
        }//endif
      }//endfor
    }//endfor
  }//endfor
  if(ic!=nktot){printf("flip error 3\n");exit(1);}

//==========================================================================
      }//end routine
//==========================================================================
