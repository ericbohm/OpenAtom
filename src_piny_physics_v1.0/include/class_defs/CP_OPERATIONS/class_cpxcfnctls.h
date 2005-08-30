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
          complex *,complex *,double *,double *);

static void CP_div_rho_gspace_calc( 
          complex* ,const int *, const int *, const int *, 
          int size, complex* , complex* , complex* );

 static void CP_getGGAFunctional(
          const int , const int ,const int ,const int ,
          complex *,complex *, complex *, complex *, 
          complex *,int ,double *);

static void CP_white_byrd_gspace_calc(
          complex *, complex *, complex *, int *, int *, int *, 
          const int , const int , const int , const int , complex *);

static void becke_gcx_lda(double ,double , double *,double *,double *,double );

//---------------------------------------------------------------------------
  }; //CPXCFNCTS
//==========================================================================

#endif

