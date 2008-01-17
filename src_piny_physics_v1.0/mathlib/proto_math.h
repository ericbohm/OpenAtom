#include "../../../src_mathlib/mathlib.h"

#ifdef FORTRANUNDERSCORE
#define DCFFTI_GENERIC dcffti_generic_
#define DCFFTF_GENERIC dcfftf_generic_
#define DCFFTB_GENERIC dcfftb_generic_
#endif

#ifdef _CRAY
#define DCFFTI_GENERIC DCFFTI_GENERIC
#define DCFFTF_GENERIC DCFFTF_GENERIC
#define DCFFTB_GENERIC DCFFTB_GENERIC
#endif

#ifdef FORTRANUNDERSCORE_OFF
#define DCFFTI_GENERIC dcffti_generic
#define DCFFTF_GENERIC dcfftf_generic
#define DCFFTB_GENERIC dcfftb_generic
#endif

#ifdef FORTRANUNDERSCORE_DBLE
#define DCFFTI_GENERIC _dcffti_generic_
#define DCFFTF_GENERIC _dcfftf_generic_
#define DCFFTB_GENERIC _dcfftb_generic_
#endif


/*----------------------------------------------------------------------*/

extern "C" {void DCFFTI_GENERIC(int *,double *);}

extern "C" {void DCFFTF_GENERIC(int *,double *,double *);}

extern "C" {void DCFFTB_GENERIC(int *,double *,double *);}

/*----------------------------------------------------------------------*/

void gethinv(double *, double *, double *, int );

double getdeth(double *);

double ddot1(int ,double *,int,double *,int);

double dsum1(int ,double *,int);

double gerf(double);

double gerfc(double);

double surf_corr(double);

double dsurf_corr(double);

double d2surf_corr(double);

void matmul_2(double *, double *, double *, int );

void matmul_2s(double *, double *, int );

void matmul_3(double *, double *);

void matmul_tt(double *, double *, double *, int );

void matmul_t(double *, double *, double *, int );

void matmul_t2(double *, double *, double *, int );

void diag33(double *, double *, double *, double *, double *);

void cputime(double *);

/*----------------------------------------------------------------------*/

void gaussran(int, long *, long *, double *, double *);

double ran_essl(double *);

double altRandom(long *);

/*----------------------------------------------------------------------*/

void sort_commence(int , int [],int []);

/*----------------------------------------------------------------------*/
