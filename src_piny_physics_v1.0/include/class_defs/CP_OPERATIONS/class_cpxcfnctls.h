//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                          CP exchange correlation (of whatever/however)
//==========================================================================

#ifndef _CPXCFNCTS_
#define _CPXCFNCTS_

class CPXCFNCTS{

 //---------------------------------------------------------------------------
 public:

 //---------------------------------------------------------------------------
 //con-destruct:

   CPXCFNCTS(){};
  ~CPXCFNCTS(){};

 //---------------------------------------------------------------------------
 // functions 

static void CP_exc_calc(
          const int , const int , const int , const int ,
          double *,double *,double *,double *);

static void CP_div_rho_gspace_calc( 
          complex* ,const int *, const int *, const int *, 
          int size, complex* , complex* , complex* );

 static void CP_getGGAFunctional(
          const int , const int ,const int ,const int ,
          double *,double *, double *, double *, 
          double *,int ,double *);

static void CP_white_byrd_gspace_calc(
          complex *, complex *, complex *, int *, int *, int *, 
          const int , const int , const int , const int , complex *);

static void becke_gcx_lda(double ,double , double *,double *,double *,double );

static void lyp_gcc(double ,double ,double *,double *,double *);

static void CP_fetch_hmati(double **, double *);

//---------------------------------------------------------------------------
  }; //CPXCFNCTS
//==========================================================================

#endif

