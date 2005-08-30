#include "../../../../include/Atoms.h"

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                          CP Local energy (however you calculate it) 
//==========================================================================

#ifndef _CPLOCAL_
#define _CPLOCAL_

class CPLOCAL{

 //---------------------------------------------------------------------------
 public:

 //---------------------------------------------------------------------------
 //con-destruct:
  CPLOCAL(){};
 ~CPLOCAL(){};

 //---------------------------------------------------------------------------
 // functions


static void CP_hart_eext_calc(
           const int , complex *,const int , Atom *,
           complex *, double *,double *,double *,
           const int *, const int *, const int *,int );

static void CP_get_vpsnow(int *,int ,
                          double ,double ,double ,
                          double *,double *,double *,double *,
                          double *,int *,int ,int ,double);

//---------------------------------------------------------------------------
   };  //CPLOCAL
//==========================================================================

#endif
