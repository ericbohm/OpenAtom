//==========================================================================
/** \file para_grp_parse.C
 *
 */
//==========================================================================
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "para_grp_parse.h"

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Decompose the nktot data points into chunks based on plane index
//==========================================================================

void ParaGrpParse::get_plane_decomp_prms(
                int index,int nktot,int nplane, int *kz,
                int *n_ret,int *istrt_ret,int *iend_ret
               )

//==========================================================================
   { //begin routine
//==========================================================================
// Error check input

  if(nktot<0 || index< 0 || nplane < 0 || index >= nplane){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("Input error to get_plane_decomp_parmsn");
     PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
  }//endif

  for(int i=0;i<nktot;i++){
    if(kz[i]<0 || kz[i]>=nplane){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Error in get_plane_decomp_parmsn");
      PRINTF("Decomposing spherical cutoff only.\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif
  }//endfor

  for(int i=1;i<nktot;i++){
    if(kz[i]<kz[(i-1)]){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Error in get_plane_decomp_parmsn");
      PRINTF("The g-space pts must lie in ascending order\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif
  }//endfor

//==========================================================================
// Compute decomposition

  int n     = 0;
  int istrt = -1;
  int iend  = 0;
  for(int i=0;i<nktot;i++){
    if(kz[i]==index){
       n++;
       if(istrt<0){istrt=i;}
       iend=i+1;
    }//endif
  }//endfor
  if(n==0){istrt=nktot+1; iend=nktot;}

//==========================================================================
// set return values

  *n_ret     = n;
  *istrt_ret = istrt;
  *iend_ret  = iend;

//-------------------------------------------------------------------------
   }//end routine 
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Find the data pts of the plane decomp inside the equal decomp
//    The number of planes in this portion of the equal decomp
//    The index of each plane in this portion.
//    The start and end of each plane in this portion.
//    The number of pts of each plane in this portion.
//==========================================================================

void ParaGrpParse::get_equal_to_plane_info(
                    int nsize, int nplane,int *kz,int *nplane_pln_ret,
                    int **iplane_pln_ret,int **nplane_pln_pts_ret,
                    int **iplane_pln_strt_ret,int **iplane_pln_iend_ret
                   )

//==========================================================================
   {//begin routine
//==========================================================================
// check for input errors

  for(int i=0;i<nsize;i++){
    if(kz[i]<0 || kz[i]>=nplane){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Error in get_equal_to_plane_info\n");
      PRINTF("Decomposing spherical cutoff only.\n");    
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif
  }//endfor

  for(int i=1;i<nsize;i++){
    if(kz[i]<kz[(i-1)]){
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("Error in get_equal_to_plane_info\n");
      PRINTF("The g-space pts must lie in ascending order\n");
      PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif
  }//endfor

//==========================================================================
// Get the mappling of equal decomposition to plane decomposition

  int nplane_pln = 1;
  int ikz        = kz[0];
  for(int i=0;i<nsize;i++){if(ikz!=kz[i]){nplane_pln++; ikz=kz[i];}}
  int *iplane_pln      = (int *)malloc(nplane_pln*sizeof(int));
  int *nplane_pln_pts  = (int *)malloc(nplane_pln*sizeof(int));
  int *iplane_pln_strt = (int *)malloc(nplane_pln*sizeof(int));
  int *iplane_pln_iend = (int *)malloc(nplane_pln*sizeof(int));

  int ic              = 0;
  iplane_pln[ic]      = kz[0];
  iplane_pln_strt[ic] = 0;
  for(int i=0;i<nsize;i++){
    if(iplane_pln[ic]!=kz[i]){
      nplane_pln_pts[ic]  = i-iplane_pln_strt[ic];
      iplane_pln_iend[ic] = i;
      ic++;
      iplane_pln[ic]      = kz[i];
      iplane_pln_strt[ic] = i;
    }//endif
  }//endfor
  iplane_pln_iend[ic] = nsize;
  nplane_pln_pts[ic]  = nsize-iplane_pln_strt[ic];

//==========================================================================
// Check the output

  int ierr = 0;

  int nc   = nplane_pln_pts[0];;
  if(nplane_pln > nplane){ierr++;}
  if(iplane_pln_strt[0]!=0) {ierr++;}
  if(iplane_pln[0]>nplane < iplane_pln[0]<0){ierr++;}
  for(int i=1;i<nplane_pln;i++){
    if(iplane_pln[i]<=iplane_pln[(i-1)]){ierr++;}
    if(iplane_pln[i]>nplane || iplane_pln[i]< 0){ierr++;}
    if(iplane_pln_iend[(i-1)]!=iplane_pln_strt[i]){ierr++;}
    nc += nplane_pln_pts[i];
  }
  if(iplane_pln_iend[(nplane_pln-1)]!=nsize){ierr++;}
  if(nc!=nsize){ierr++;}

  if(ierr!=0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Error in get_equal_to_plane_info\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }

//==========================================================================
// Assign return values


  *nplane_pln_ret      = nplane_pln;
  *iplane_pln_ret      = iplane_pln;
  *nplane_pln_pts_ret  = nplane_pln_pts;
  *iplane_pln_strt_ret = iplane_pln_strt;
  *iplane_pln_iend_ret = iplane_pln_iend;

//-------------------------------------------------------------------------
   }//end routine 
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Decompose the ntot data points into nprt chunks
//==========================================================================

void ParaGrpParse::get_equal_decomp_prms(
                           int ntot, int nprt, int index, 
                           int *n_ret, int *istrt_ret, int *iend_ret
                           )

//==========================================================================
   {// begin routine
//==========================================================================
// Error check input

  if(ntot<0 || nprt < 0 || index < 0 || index >= nprt){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("Input error in get_equal_decomp_prms\n");
     PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
  }//endif

//==========================================================================
// Compute decomposition

  int n      = (ntot/nprt);   //nprt is number of partitions
  int m      = (ntot % nprt);
  int istrt  = n*index;
  if(index>=m){istrt += m;}
  if(index< m){istrt += index;}
  if(index< m){n++;}
  int iend   = n+istrt;
  if(iend==istrt){istrt++;}// if nprt>ntot && index>=ntot then istrt=iend 
                           // For clarity increment istrt so istrt>iend

  *n_ret     = n;
  *istrt_ret = istrt;
  *iend_ret  = iend;

//==========================================================================
   }//end routine
//==========================================================================




//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// Find the data pts of the equal decomp inside the plane decomp
//==========================================================================

void ParaGrpParse::get_plane_to_equal_info(
                        int n,int istrt,int iend, int nplane, int nktot,
                        int *nplane_eql_ret, int **iplane_eql_ret, 
                        int **nplane_eql_pts_ret, int **iplane_eql_strt_ret,
                        int **iplane_eql_iend_ret
                      )

//==========================================================================
   { //begin routine
//==========================================================================
// Error check input

  if(n<0 || istrt < 0 || iend > nktot || n != (iend-istrt) || nplane < 0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Error in get_plane_to_equal_info");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
  }//endif

//==========================================================================
// Compute mapping from plane decomp to equal decomp

  int nplane_eql = 0;
  for(int i=0;i<nplane;i++){
    int neql,istrt_eql,iend_eql;
    get_equal_decomp_prms(nktot,nplane,i,&neql,&istrt_eql,&iend_eql);
    if(iend_eql > istrt && istrt_eql < iend){nplane_eql++;}
    if(istrt_eql >=iend){break;}
  }//endfor

  int *iplane_eql      = (int *)malloc(nplane_eql*sizeof(int));
  int *nplane_eql_pts  = (int *)malloc(nplane_eql*sizeof(int));
  int *iplane_eql_strt = (int *)malloc(nplane_eql*sizeof(int));
  int *iplane_eql_iend = (int *)malloc(nplane_eql*sizeof(int));

  int ic=0;
  for(int i=0;i<nplane;i++){
    int neql,istrt_eql,iend_eql;
    get_equal_decomp_prms(nktot,nplane,i,&neql,&istrt_eql,&iend_eql);
    if(iend_eql > istrt && istrt_eql < iend){
      iplane_eql[ic]      = i;
      iplane_eql_strt[ic] = MAX(istrt,istrt_eql)-istrt;
      iplane_eql_iend[ic] = MIN(iend,iend_eql)  -istrt;
      nplane_eql_pts[ic]  = iplane_eql_iend[ic]-iplane_eql_strt[ic];
      ic++;
    }//endif
    if(istrt_eql >=iend){break;}
  }//endfor

//==========================================================================
// Error check output

  int ierr = 0;

  int nc   = nplane_eql_pts[0];;
  if(nplane_eql > nplane){ierr++;}
  if(iplane_eql_strt[0]!=0) {ierr++;}
  if(iplane_eql[0]>nplane < iplane_eql[0]<0){ierr++;}
  for(int i=1;i<nplane_eql;i++){
    if(iplane_eql[i]<=iplane_eql[(i-1)]){ierr++;}
    if(iplane_eql[i]>nplane < iplane_eql[i]<0){ierr++;}
    if(iplane_eql_iend[(i-1)]!=iplane_eql_strt[i]){ierr++;}
    nc += nplane_eql_pts[i];
  }
  if(iplane_eql_iend[(nplane_eql-1)]!=n){ierr++;}
  if(nc!=n){ierr++;}

  if(ierr!=0){
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("Error in get_plane_to_equal_info\n");
    PRINTF("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

//==========================================================================
// Assign the return values

  *nplane_eql_ret      = nplane_eql;
  *iplane_eql_ret      = iplane_eql;
  *nplane_eql_pts_ret  = nplane_eql_pts;
  *iplane_eql_strt_ret = iplane_eql_strt;
  *iplane_eql_iend_ret = iplane_eql_iend;

//-------------------------------------------------------------------------
   }//end routine 
//==========================================================================
