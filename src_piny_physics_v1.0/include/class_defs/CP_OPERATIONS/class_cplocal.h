#include "../../../../include/Atoms.h"

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                          CP Local energy (however you calculate it) 
//==========================================================================

#ifndef _CPLOCAL_
#define _CPLOCAL_

#include "../../../../include/Atoms.h"
#include "../../../../include/eesDataClass.h"

class CPLOCAL{

 //---------------------------------------------------------------------------
 public:

 //---------------------------------------------------------------------------
 //con-destruct:
  CPLOCAL(){};
 ~CPLOCAL(){};

 //---------------------------------------------------------------------------
 // functions
  static void CP_hart_eext_calc(int , complex *,int , Atom *,complex *, 
                          double *,double *,double *,int *, int *, int *,int );

  static void CP_get_vpsnow(int *,int ,double ,double ,double ,
                            double *,double *,double *,double *,
                            double *,int *,int ,int ,double);
 //---------------------------------------------------------------------------

  static void getEesPrms(int *, int *, int *, int *, int *);

  static void eesSetEesWghtGgrp(int , int *, int *, int *,double *, double *, 
                                int ,int ,int ,int );

  static void eesAtmBsplineRgrp(Atom *, int *, RHORHARTDATA *);

  static void eesPackGridRchare(int , int , double *, int ,int **, double **, int *);

  static void eesHartEextGchare(int , int , complex *, complex *, complex *, complex *,
                                double *, double *,double *,double *,
                                int *,int *,int *,int);

  static void eesEwaldGchare(int ,complex *,double *, double *,double *,
                             int *, int *, int *,int);

  static void eesAtmForceRchare(int , Atom *,int ,int **, double **, double **, 
                                double **, int *, double *, int ,int );

//---------------------------------------------------------------------------
   };  //CPLOCAL
//==========================================================================

#endif
