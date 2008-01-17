#include "../../../../include/Atoms.h"

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                          CP Integration                                    
//==========================================================================

#ifndef _CPRSPACEION_
#define _CPRSPACEION_

class CPRSPACEION{

 //---------------------------------------------------------------------------
 public:

 //---------------------------------------------------------------------------
 //con-destruct:

   CPRSPACEION(){};
  ~CPRSPACEION(){};

 //---------------------------------------------------------------------------
 // Functions

static void CP_getionforce(const int ,FastAtoms *,int ,int ,double *, double *, double *,double *);

static void atm_recip_corr(int ,double *, double *, double *,
                           double *, double *, double *,double *,
                           double, double *,double *,int ,
  			   double *,int , int );


//---------------------------------------------------------------------------
   }; //CPRSPACEION
//==========================================================================

#endif
