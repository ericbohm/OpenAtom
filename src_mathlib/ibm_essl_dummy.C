#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pup.h"
#include "mathlib.h"

void dcft(int *init,complex *x,int *istride,int *iskip,complex *y,int *ostride,int *oskip,
          int *nfft,int *num,int *isign,double *scale,
          double *work1, int *nwork1,double *work2, int *nwork2){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("IBM ESSL is not available\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    exit(1);
}//end routine

void drcft(int *init,double *x,int *istride,complex *y,int *ostride,
           int *nfft,int *num,int *isign,double *scale,
	   double *work1, int *nwork1,double *work2,int *nwork2){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("IBM ESSL is not available\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    exit(1);
}//end routine

void dcrft(int *init,complex *x,int *istride,double *y,int *ostride,
           int *nfft,int *num,int *isign,double *scale,
	   double *work1, int *nwork1,double *work2,int *nwork2){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("IBM ESSL is not available\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    exit(1);
}//end routine
