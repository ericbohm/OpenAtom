#ifndef _MathlibInc_
#define _MathlibInc_

#ifdef FORTRANUNDERSCORE
#define DURAND durand_
#define DGEFA dgefa_
#define DGESL dgesl_
#define DGEMM dgemm_
#define DGEMV dgemv_
#define GENMATMUL genmatmul_
#define RS rs_
#define ZGEMV zgemv_
#define DSPEV dspev_
#ifdef _CP_USE_BLAS_GATHER_
#define ZGTHR zgthr_
#define ZSCTR zsctr_
#endif
#endif

#ifdef FORTRANUNDERSCORE_OFF
#define DURAND durand
#define DGEFA dgefa
#define DGESL dgesl
#define DGEMM dgemm
#define DGEMV dgemv
#define GENMATMUL genmatmul
#define RS rs
#define DSPEV dspev
#define ZGEMV zgemv
#ifdef _CP_USE_BLAS_GATHER_
#define ZGTHR zgthr
#define ZSCTR zsctr
#endif
#endif

#ifdef FORTRANUNDERSCORE_DBLE
#define DURAND _durand_
#define DGEFA _dgefa_
#define DGESL _dgesl_
#define DGEMM _dgemm_
#define DGEMV _dgemv_
#define GENMATMUL _genmatmul_
#define RS _rs_
#define DSPEV _dspev_
#define ZGEMV _zgemv_
#ifdef _CP_USE_BLAS_GATHER_
#define ZGTHR _zgthr_
#define ZSCTR _zsctr_
#endif
#endif

#ifdef _CRAY
#define DURAND DURAND
#define DGEFA DGEFA
#define DGESL DGESL
#define DGEMM DGEMM
#define DGEMV DGEMV
#define GENMATMUL GENMATMUL
#define RS RS
#define DSPEV DSPEV
#define ZGEMV ZGEMV
#ifdef _CP_USE_BLAS_GATHER_
#define ZGTHR ZGTHR
#define ZSCTR ZSCTR
#endif
#endif

#include "ckcomplex.h"

//========================================================================
// IBM ESSL FFT STUFF 
void dcftWrap(int *,complex *,int *,int *,complex *,int *,int *, int *,int *,int *,double *,
          double *, int *,double *,int *);
void dcrftWrap(int *,complex *,int *,double *,int *,int *,int *,int *,double *,
	   double *, int *,double *,int *);
void drcftWrap(int *,double *,int *,complex *,int *,int *,int *,int *,double *,
	   double *, int *,double *,int *);
//========================================================================

void lst_sort_clean(int , int *, int *);
void sort_commence(int , int *,int *);
void sort_commence_piny(int , int *,int *);

double altRandom(long *);

extern "C" {void DURAND(double *,int *, double *,int *);}

extern "C" {void DGEFA(double *, int *, int *, int *, int * );}

extern "C" {void DGESL(double *, int *, int *, int *, double *, int * );}

extern "C" {void DGEMM (char *, char *, int *, int *, int *,double *,double *,
                        int *, double *, int *, double *, double *, int * );}

extern "C" void DGEMV(char *,int *,int *,double *,double *,int *,
                      double *,int *,double *,double *,int *);

extern "C" {void  GENMATMUL(double *,int *,int *,double *,int *,int *,
                             double *,int *,int *,int *,int *,
			     double *,double *);}

extern "C" {void RS(int *,int *,double [],double [],int *,double [],
                    double [],double [],int *);}

extern "C" void ZGEMV(char *,int *,int *,complex *,complex *,int *,
                      complex *,int *,complex *,complex *,int *);

extern "C" {void  GENMATMUL(double *,int *,int *,double *,int *,int *,
                             double *,int *,int *,int *,int *,
			     double *,double *);}

extern "C" {void DSPEV(char *,char *,int *,double *,double *,double *,int *,
                       double *,int *);}
#ifdef _CP_USE_BLAS_GATHER_
extern "C" void ZGTHR(int *nz, complex *y, complex*x, int *indx);
extern "C" void ZSCTR(int *nz, const complex *x, int *indx,complex*y);
#endif


#endif
