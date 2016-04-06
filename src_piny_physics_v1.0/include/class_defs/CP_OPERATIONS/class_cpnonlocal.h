//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                          CP Non Local energies (however calculated)
//==========================================================================

#ifndef _CPNONLOCAL_
#define _CPNONLOCAL_

#include "../../../../include/Atoms.h"
#include "../../../../include/RunDescriptor.h"
#include "../../../../include/eesDataClass.h"
#include "../../../include/class_defs/Interface_ctrl.h"
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
        int *, int *, int *,double **,double *,int,int,int,int nfreq_cmi_update = 400);

    static void CP_enl_matrix_calc(int , complex *, int *, int *, int *, 
        complex *, complex *,complex *,complex *,
        complex *, complex *,complex *, complex *, 
        int , int , int , int );

    static void CP_enl_force_calc(complex* , int , int *, int *, int *, 
        complex *,complex *, int , int ,int , int);

    static void CP_enl_atm_forc_calc(int , int , FastAtoms *,
        complex *,complex *,complex *,complex *,double *,int, int);

    static void CP_calc_Struct_Fact(int , int *, int *, int *, 
        complex *, complex *, complex *, complex *,
        FastAtoms *,  int , int ,int, PSSCRATCH * );

    static void get_rad_proj(int ,int ,double ,double ,double ,int , int , 
        double *,double *,double *,double *,double *,double *);

    static void get_grp_params(int , int , int ,int *, int *, int *,int *);

    //---------------------------------------------------------------------------

    static void getEesPrms(int *, int *, int *,int *, int *);

    static void eesSetEesWghtGgrp(int , int *, int *, int *,double *, double *, 
        int ,int ,int ,int );

    static void eesSplProjectorGgrp(int ,int *,int *,int *,double **,int **);

    static void eesAtmBsplineRgrp(FastAtoms *, int *, RPPDATA **, PSSCRATCH *);

    static void eesProjGchare(int , complex *,int *,int *, int *,int , int , int ,
        double *, double *, double *, double *,complex *,int **,
			      double **,int, int, int, int nfreq_cmi_update, int);

    static void eesYlmOnD(int ,int ,int ,int *,int *,int *,double *,double *, double *, 
        double *,double *,int, int nfreq_cmi_update);

    static void eesZmatRchare(double *, int ,double *,int **, double **,int *,int ,int);

    static void eesZmatRchareC(complex *, int , complex *, 
        int **, double **, int *, int ,int);

    static void eesEnergyAtmForcRchare(int , double *, double *, 
        int **,double **,double **,double **,double **,
        double *, double *, int *, int *,int **,
        int , int , int, FastAtoms *);

    static void eesEnergyAtmForcRchareC(int, double *, complex *, 
        int **, double **, double **, double **,double **,
        complex *, complex *, int *, int *, int **, 
        int, int, int, FastAtoms *);

    static void eesPsiForcGspace(int, int , int ,int, complex *,complex *, double *,double *,
				 int *,int *, int *,int , int, int, int nfreq_cmi_update, int);

    static void genericSetKvector(int ,int *,int *,int *,
        int ,RunDescriptor *,GCHAREPKG *,int ,int ,int ,int,
        double **, double** );

    //---------------------------------------------------------------------------
}; //CPNONLOCAL
//==========================================================================

#endif 
