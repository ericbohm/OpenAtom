//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                          CP Integration                                    
//==========================================================================

#ifndef _CPORTHOG_
#define _CPORTHOG_

class CPORTHOG{

 //---------------------------------------------------------------------------
 public:

 //---------------------------------------------------------------------------
 //con-destruct:

   CPORTHOG(){};
  ~CPORTHOG(){};

 //---------------------------------------------------------------------------
 // functions 

static void CP_orthoTransform(double *, const double *, int );
static void get_iter_Tmat(const double *,double *,int );
static void get_diag_Tmat(const double *,double *,int );
static void get_unit_Tmat(double *,int );

//---------------------------------------------------------------------------
   }; //CPORTHOG
//==========================================================================

#endif
