//========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pup.h"
#include "mathlib.h"

//========================================================================

#ifdef _IBM_ESSL_
extern "C" {void dcft(int *,double *, int *,int *,double *y,int *,int *,
                      int *,int *,int *,double *,
		      double *, int *,double *, int *);}
extern "C" {void dcrft(int *,double *,int *,double *,int *,
                       int *,int *,int *,double *,
    	               double *, int *,double *,int *);}
extern "C" {void drcft(int *,double *,int *,double *,int *,
                       int *,int *,int *,double *,
 	               double *, int *,double *,int *);}
#endif

//========================================================================


//========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================
void dcftWrap(int *init,complex *x,int *istride,int *iskip,
                        complex *y,int *ostride,int *oskip,
              int *nfft,int *num,int *isign,double *scale,
              double *work1, int *nwork1,double *work2, int *nwork2){
//========================================================================

#ifdef _IBM_ESSL_
    double *xx    = reinterpret_cast<double *> (x);
    double *yy    = reinterpret_cast<double *> (y);
    dcft(init,xx,istride,iskip,yy,ostride,oskip,
	 nfft,num,isign,scale,
         work1,nwork1,work2,nwork2);
#else
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("IBM ESSL is not available\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    exit(1);
#endif

//------------------------------------------------------------------------
  }//end routine
//========================================================================


//========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================
void drcftWrap(int *init,double *x,int *istride,complex *y,int *ostride,
               int *nfft,int *num,int *isign,double *scale,
    	       double *work1, int *nwork1,double *work2,int *nwork2){
//========================================================================

#ifdef _IBM_ESSL_
    double *yy    = reinterpret_cast<double *> (y);
    drcft(init,x,istride,yy,ostride,
          nfft,num,isign,scale,
          work1,nwork1,work2,nwork2);
#else
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("IBM ESSL is not available\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    exit(1);
#endif

//------------------------------------------------------------------------
  }//end routine
//========================================================================


//========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//========================================================================
void dcrftWrap(int *init,complex *x,int *istride,double *y,int *ostride,
               int *nfft,int *num,int *isign,double *scale,
	       double *work1, int *nwork1,double *work2,int *nwork2){
//========================================================================

#ifdef _IBM_ESSL_
    double *xx    = reinterpret_cast<double *> (x);
    dcrft(init,xx,istride,y,ostride,
          nfft,num,isign,scale,
          work1,nwork1,work2,nwork2);
#else
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("IBM ESSL is not available\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    exit(1);
#endif

//------------------------------------------------------------------------
  }//end routine
//========================================================================
