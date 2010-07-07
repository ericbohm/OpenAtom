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
  static void CP_hart_eext_calc(int , complex *,int , FastAtoms *,complex *, 
                          double *,double *,double *,int *, int *, int *, double *,int,int nfreq_cmi_update=4);

  static void CP_get_vpsnow(int *,int ,double ,double ,double ,
                            double *,double *,double *,double *,
                            double *,int *,int ,int ,double);
 //---------------------------------------------------------------------------

  static void getEesPrms(int *, int *, int *, int *, int *);

  static void eesSetEesWghtGgrp(int , int *, int *, int *,double *, double *, 
                                int ,int ,int ,int );

  static void eesAtmBsplineRgrp(FastAtoms *, int *, RHORHARTDATA **);

  static void eesPackGridRchare(int , int , double *, int ,int , int ***, double ***, int *, 
                                int **, int );

  static void eesHartEextGchare(int , int , complex *, complex *, complex *, complex *,
                                double *, double *,double *,double *,
                                int *,int *,int *,double *,int,int nfreq_cmi_update = 100);

  static void eesEwaldGchare(int ,complex *,double *, double *,double *,
                             int *, int *, int *, double *,int,int nfreq_cmi_update = 100);

  static void eesAtmForceRchare(int , FastAtoms *,int , 
                                int ***, double ***, double ***, double ***, 
                                int *, int **,  double *, int ,int ,int );
//---------------------------------------------------------------------------
   };  //CPLOCAL
//==========================================================================

#endif
