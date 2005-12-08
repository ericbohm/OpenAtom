
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                          CP Integration                                    
//==========================================================================

#ifndef _CPINTEGRATE_
#define _CPINTEGRATE_

class CPINTEGRATE{

 //---------------------------------------------------------------------------
 public:

 //---------------------------------------------------------------------------
 //con-destruct:

   CPINTEGRATE(){};
  ~CPINTEGRATE(){};

 //---------------------------------------------------------------------------
 // Functions


  static void CP_integrate(int , int ,int , complex *,complex *,complex *,double *,
                           int *,int *, int *, int , int , double **,double **,
                           double *,double ,double , double ,double *,int ,int ,int,
                           double *,double,double,double);

  static void CP_create_mass(int ,int *, int *, int *, double *,int);

  static void CheckCoefGradMag(double);

  static void CP_fovlap_calc(int, int, complex *,double *);

  static void CP_integrate_dyn(int , int ,int ,complex *,complex *,complex *,double *, 
                               int *, int *,int *,int , int , double **,double **,
                               double *,double , double ,  double *,int ,int ,int,
                               double *,double,double,double);

  static void initCPNHC(int ,int ,int , int *, int *, double *,	double *,double *,
                        double*,double*,double*,int,int);

  static void CPSmplVel(int,double *,complex *,int ,int,double ,double **,
                        double ,int, int, int,int,double,double,double);

  static void sampl1DVelOneT(int , double* ,double* ,double );

  static void cpSamplNHC(int , int ,double** ,double , double );

  static void cp_evolve_vel(int , complex *, complex *, double *,int , int , 
                            double **,double **,double *,double ,double ,
     		            int ,int ,int ,int ,int,double,double,double);

  static void CP_integrate_min_STD(int , int ,complex *,complex *,
                                   int *,int *, int *,double *);

  static void CP_integrate_min_CG(int ,int ,complex *,complex *,complex *,
                                  int *,int *, int *,double *, double);

  static void get_fictKE(int ,complex *, double *,int , int ,
			 double **,double ,int ,int ,int,double *,
                         double*,double);

  static void cp_isoNHC_update(int ,complex *,double *,
                      int ,int ,double *,double **,double **,double ,double ,
                      double,double,double);

  static void cp_isoVel_update(int,complex *,complex*,double *,int,double **,
                               double ,double ,double ,double ,double );

  inline static void evolve_v_vNHC(int ,double *,double *,double **,double ,double ,double ,
                            double ,double );

  inline static void get_forc_NHC(int ,int ,double **,double **,double ,double );

  inline static void evolve_vNHCM(int ,int,double **,double **,double );

  inline static void evolve_vNHC(int ,int ,double **,double **,double ,double );

  inline static void evolve_pNHC(int ,int ,double *,double **,double ,double ,double );

  static void set_yosh(int ,double ,double *,double *,double *,double *);

//---------------------------------------------------------------------------
   }; //CPINTEGRATE
//==========================================================================

#endif
