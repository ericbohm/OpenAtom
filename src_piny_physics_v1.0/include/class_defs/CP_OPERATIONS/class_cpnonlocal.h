//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                          CP Non Local energies (however calculated)
//==========================================================================

#ifndef _CPNONLOCAL_
#define _CPNONLOCAL_

#include "../../../../include/Atoms.h"
#include "../../../../include/eesDataClass.h"

class CPNONLOCAL{

 //---------------------------------------------------------------------------
 public:

 //---------------------------------------------------------------------------
 //con-destruct:
   CPNONLOCAL(){};
  ~CPNONLOCAL(){};

 //---------------------------------------------------------------------------
 // functions

static void CP_eke_calc(int, int ,complex *,complex *,
                        int *, int *, int *,double *,int);

static void CP_enl_matrix_calc(int , complex *, int *, int *, int *, 
                               complex *, complex *,complex *,complex *,
                               complex *, complex *,complex *, complex *, 
                               int , int , int , int );

static void CP_enl_force_calc(complex* , int , int *, int *, int *, 
                          complex *,complex *, int , int ,int , int);

static void CP_enl_atm_forc_calc(int , int , Atom *,
                          complex *,complex *,complex *,complex *,double *,int );

static void CP_calc_Struct_Fact(int , int *, int *, int *, 
                                complex *, complex *, complex *, complex *,
                                Atom *,  int , int ,int );

static void get_rad_proj(int ,int ,double ,double ,double ,int , int , 
                         double *,double *,double *,double *,double *,double *);

static void get_grp_params(int , int , int ,int *, int *, int *,int *);

//---------------------------------------------------------------------------

static void getEesPrms(int *, int *, int *,int *, int *);

static void eesSetEesWghtGgrp(int , int *, int *, int *,double *, double *, 
                              int ,int ,int ,int );

static void eesSplProjectorGgrp(int , int *, int *, int *,double *,int *);

static void eesAtmBsplineRgrp(Atom *, int *, RPPDATA *);

static void eesProjGchare(int , complex *,int *,int *, int *,int , int , int ,
                          double *, double *, double *, double *,complex *,int *,double *);

static void eesYlmOnD(int ,int ,int ,int *,int *,int *,double *,double *, double *, 
                      double *,double *);

static void eesZmatRchare(double *, int ,double *,int **, double **,int *,int );

static void eesEnergyAtmForcRchare(int , double *, double *, int **,double **,
                                   double **,double **,double *, int *, int ,Atom *);

static void eesPsiForcRchare(int , double *,int **,double **,double *, int *,int );

static void eesPsiForcGspace(int , int , int ,complex *,complex *, double *,double *);

//---------------------------------------------------------------------------
   }; //CPNONLOCAL
//==========================================================================

#endif 
