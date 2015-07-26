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
#include <vector>
#include "fft_types.h"
class PSSCRATCH;

class CPLOCAL{

  //---------------------------------------------------------------------------
  public:

    //---------------------------------------------------------------------------
    //con-destruct:
    CPLOCAL(){};
    ~CPLOCAL(){};

    //---------------------------------------------------------------------------
    // functions
    static void CP_hart_eext_calc(int , std::vector< gridPoint > & , complex *,
        int, int , FastAtoms *, complex *, double *, double *, double *, int,
        PSSCRATCH *, int );

    static void CP_get_vpsnow(int *,int ,double ,double ,double ,
        double *,double *,double *,double *,
        double *,int *,int ,int ,double);
    //---------------------------------------------------------------------------

    static void getEesPrms(int *, int *, int *, int *, int *);

    static void eesSetEesWghtGgrp(int , std::vector< gridPoint > &, double *,
        double *, int, int, int, int);

    static void eesAtmBsplineRgrp(FastAtoms *, std::vector<RHORHARTDATA>
        &, PSSCRATCH *);

    static void eesPackGridRchare(int, int, double*, int, int, int, int, int,
        int ***, double ***, int **, int **);

    static void eesHartEextGchare(int , int , complex *, complex *, complex *,
        complex *, double *, double *,double *,double *,
        std::vector< gridPoint > &, double *, int,
        int nfreq_cmi_update = 100);

    static void eesEwaldGchare(int ,complex *,double *, double *,double *,
        std::vector< gridPoint > &, int, double *,int,int nfreq_cmi_update = 100);

    static void eesAtmForceRchare(int , FastAtoms *,int ,
        int ***, double ***, double ***, double ***,
        int **, int **,  double *, int, int, int, int, int );
    //---------------------------------------------------------------------------
};  //CPLOCAL
//==========================================================================

#endif
