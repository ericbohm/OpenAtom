
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

  static void CP_integrate(const int , const int ,int,complex *,complex *,complex *,
                           const int *,const int *, const int *,
                           double *,double ,double *);

  static void CP_integrate_min_STD(const int , const int ,complex *,complex *,
                           const int *,const int *, const int *,const double *);

  static void CP_integrate_min_CG(const int ,const int ,complex *,complex *,complex *,
                           const int *,const int *, const int *,const double *,
                           double);

  static void CP_create_mass(const int ,const int *, const int *, const int *,
                             double *,const int);

  static void CheckCoefGradMag(const double);

  static void CP_fovlap_calc(const int, const int, complex *,double *);

  static void CP_integrate_dyn(const int ,const int , int, complex *,complex *,
                               complex *,const int *, const int *, const int *,
                               double *, double *);

  static void CP_integrate_half_vel(int ,int ,complex *,complex *,complex *,
                                    int *,int *,int *,double *);

  static void initCPNHC(int,int ,int *,double *,double *,double *);
  static void CPSmplVel(int,double *,complex *,int ,double *,double *,
                        double ,int);
  static void sampl1DVelOneT(int , double* ,double* ,double );


//---------------------------------------------------------------------------
   }; //CPINTEGRATE
//==========================================================================

#endif
