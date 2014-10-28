
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                          CP Integration                                    
//==========================================================================

#ifndef _CPINTEGRATE_
#define _CPINTEGRATE_

class CPINTEGRATE{

  public:
    //==========================================================================
    //con-destruct:

    CPINTEGRATE(){};
    ~CPINTEGRATE(){};

    //==========================================================================
    // Functions

    //-------------------------------------------------------------------------
    // cp_integrate.C
    static void CP_integrate(int , int ,int , complex *,complex *,complex *,double *,
        int *,int *, int *, int , int , int ,
        double ***,double ***,double***,double ***,
        double *,double *, double *, double *,double ,
        double , double *,int ,int ,int,
        double *,double*,double,double,double*,int *,int *,int,int nfreq_cmi_update=400);

    static void CP_create_mass(int ,int *, int *, int *, double *,int);

    static void CheckCoefGradMag(double);

    //-------------------------------------------------------------------------
    // cp_min_STD.C
    static void CP_integrate_min_STD(int , int ,complex *,complex *,
        int *,int *, int *,double *,int nfreq_cmi_update);

    //-------------------------------------------------------------------------
    // cp_min_CG.C
    static void CP_integrate_min_CG(int ,int ,complex *,complex *,complex *,
        int *,int *, int *,double *, double, int nfreq_cmi_update);

    static void CP_fovlap_calc(int, int, complex *,double *);

    //-------------------------------------------------------------------------
    // cp_velSampl.C
    static void CPSmplVel(int ,double *,complex *,int ,int ,int, 
        double *,double ***, double ***, double ***,double *, 
        double ,int , int ,int ,int , double ,double);

    static void sampl1DVelOneT(int , double* ,double* ,double );

    static void cpSamplNHC(int ,int ,int,double*** ,double ***, double ***,double *, double *, 
        double  );

    static void CPrandomPsi(int , int , complex *);
    //-------------------------------------------------------------------------
    // cp_dynamics.C
    static void CP_integrate_dyn(int , int ,int ,complex *,complex *,complex *,double *, 
        int *, int *,int *,int , int ,int,
        double ***,double ***,double ***,double ***,
        double *,double *, double *, double *, double ,
        double *,int ,int ,int , double *, double *,
        double ,double , double *,int *,int *,int );

    static void cp_evolve_vel(int , complex *, complex *,double *,int , int , int,
        double ***,double ***,double ***,double ***,
        double *,double *, double *, double *, double ,
        int ,int ,int , int ,int ,double ,double,double *,
        int *, int * ,int);

    static void get_fictKE(int ,complex *, double *,int , int ,int,
        double ***,double ***,double ***,double *,
        int ,int ,int , double *,double *,double *);

    //-------------------------------------------------------------------------
    // cp_isokin.C

    static void fetchNHCsize(int *,int *,int *);

    static void initCPNHC(int ,int , int, int, int,double *,
        double *,double *,double *, double *,
        double *,double *, double *,double *,int *,int *);

    static void cp_isoNHC_update(int ,complex *,double *,int ,int ,int ,
        double ***,double ***, double ***,double ***,
        double *,double *, double *, double *, 
        double, double , double , double* ,int *,int *,int);

    static void set_yosh(int ,double ,double *,double *,double *,double *);

    //---------------------------------------------------------------------------
}; //CPINTEGRATE
//==========================================================================

#endif
