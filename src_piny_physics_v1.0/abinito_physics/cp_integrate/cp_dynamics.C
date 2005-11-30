#include "standard_include.h"
#include "ckcomplex.h"

#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_mdintegrate.h"

#include "../class_defs/CP_OPERATIONS/class_cpintegrate.h"
#include "../../mathlib/proto_math.h"

#define ISOKIN_OPT 0

//============================================================================
// Dynamics controller
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPINTEGRATE::CP_integrate_dyn(int ncoef, int istate,int iteration,
              complex *forces,complex *vpsi,complex *psi,double *cmass, 
              int *k_x, int *k_y,int *k_z,int len_nhc, int num_nhc,
              double **fNHC,double **vNHC,double *xNHC,double mNHC,double kTCP,
              double *fictEke,int nkx0_red,int nkx0_uni,int nkx0_zero)
//============================================================================
    { // Begin Function 
//----------------------------------------------------------------------------
// Local Variables

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP *cp = CP::get();
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

   double dt   = gentimeinfo->dt;
   double dt2  = dt*0.5;
   int istrt   = nkx0_red;
   int iend    = nkx0_uni+nkx0_red;
   int applyorder;

//============================================================================
// Update the velocities from time, t-dt/2, to time, t (if not the intial step).
// Compute the Fictious Kinetic Energy at time, t (e.g. t=0 on initial step).

   if(iteration>1){
     applyorder = 2;
     cp_evolve_vel(ncoef,forces,vpsi,cmass,len_nhc,num_nhc,fNHC,vNHC,
                   xNHC,mNHC,kTCP,nkx0_red,nkx0_uni,nkx0_zero,applyorder,iteration);
   }//endif

   get_fictKE(ncoef,vpsi,cmass,len_nhc,num_nhc,vNHC,mNHC,
              nkx0_red,nkx0_uni,nkx0_zero,fictEke);

//============================================================================
// Update the velocities from time, t, to time, t+dt/2. 

   applyorder = 1;
   cp_evolve_vel(ncoef,forces,vpsi,cmass,len_nhc,num_nhc,fNHC,vNHC,
                 xNHC,mNHC,kTCP,nkx0_red,nkx0_uni,nkx0_zero,applyorder,iteration);

//============================================================================
// Update the positions to the next step (t+dt)

   for(int i=istrt;i<ncoef;i++){
     psi[i] += vpsi[i]*dt;
   }//endfor

//============================================================================
  } // End function 
//============================================================================


//============================================================================
// Fancy Dan velocity integration under CP
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPINTEGRATE::cp_evolve_vel(int ncoef, complex *forces, complex *vpsi,
                    double *cmass,int len_nhc, int num_nhc, double **fNHC,
                    double **vNHC,double *xNHC,double mNHC,double kTCP,
                    int nkx0_red,int nkx0_uni,int nkx0_zero,
                    int applyorder,int iteration)
//============================================================================
   { // Begin Function
//----------------------------------------------------------------------------
// Local Variables and error checking

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP *cp = CP::get();
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

   double dt      = gentimeinfo->dt;
   double dt2     = dt*0.5;
   int istrt      = nkx0_red;
   int iend       = nkx0_uni+nkx0_red;
   int igo        = 1;
   if(applyorder==2 && iteration==1){igo=0;}

//============================================================================
// Isokin-NHC apply first

#ifdef NEXT_STEP
   if(applyorder==1 && ISOKIN_OPT==1){
     cp_NHC_apply();
   }//endif
#endif

//============================================================================
// Evolve velocity with forces

   if(ISOKIN_OPT==0 && igo==1){
     for(int i=istrt;i<ncoef;i++){
       vpsi[i] += forces[i]*(dt2/cmass[i]);     
     }//endfor
   }//endif

#ifdef NEXT_STEP
   if(ISOKIN_OPT==1 && igo==1){
     cp_isokin_update();
   }//endif
#endif

//============================================================================
// Isokin-NHC apply second

#ifdef NEXT_STEP
   if(applyorder==2 && igo==1){
     cp_NHC_apply();
   }//endif
#endif

//============================================================================
  } // End function
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPINTEGRATE::get_fictKE(int n,complex *v, double *m,int num, int len,
          double **vNHC,double mNHC,int nkx0_red,int nkx0_uni,int nkx0_zero,
          double *EKin_ret)
//============================================================================
  {// Begin Function 
//============================================================================

  int istrt0     = nkx0_red;
  int istrt      = nkx0_red+nkx0_zero;
  int iend       = nkx0_red+nkx0_uni;

  double EKin = 0.0;
  if(ISOKIN_OPT==1){
    for(int i=0;i<num;i++){
    for(int j=0;j<len;j++){
      EKin += mNHC*vNHC[i][j]*vNHC[i][j];
    }}//endfor
  }//endif
  for(int i=istrt0;i<istrt;i++){EKin += v[i].getMagSqr()*m[i];}
  for(int i=istrt;i<iend;i++)  {EKin += v[i].getMagSqr()*(2.0*m[i]);}
  for(int i=iend;i<n;i++)      {EKin += v[i].getMagSqr()*m[i];}
  EKin       *=0.5;

  (*EKin_ret) = EKin;

//============================================================================
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPINTEGRATE::initCPNHC(int npts,int maxLen,int maxNum, int *len_nhc_ret,
                            int *num_nhc_ret, double *kTCP_ret,
                            double *tau_ret,double *mNHC_ret)
//============================================================================
  {// Begin Function 
//============================================================================

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP *cp = CP::get();
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

  int len_nhc   = cptherm_info->len_c_nhc;
  int num_nhc   = 3;
  double Text   = cpopts->te_ext;
  double tau_in = cpopts->cp_tau_nhc;

  double kT     = Text/BOLTZ;
  double tau    = tau_in/TIME_CONV;

  if(len_nhc>maxLen){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Maximum NHC length is %d < %d\n",maxLen,len_nhc);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
  }//endif
  if(num_nhc>maxNum){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Maximum NHC number is %d < %d\n",maxNum,num_nhc);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
  }//endif

//============================================================================

  (*len_nhc_ret) = len_nhc;
  (*num_nhc_ret) = num_nhc;
  (*kTCP_ret)    = kT;
  (*tau_ret)     = tau;
  
  double pre     = tau*tau*kT;
  (*mNHC_ret)    = pre;

//============================================================================
   }// End function
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPINTEGRATE::CPSmplVel(int n,double *m,complex *v,int len_nhc,int num_nhc,
                            double mNHC,double **vNHC,double kT,int istart_typ_cp,
                            int nkx0_red,int nkx0_uni,int nkx0_zero)
//============================================================================
    {//Begin Function
//============================================================================
// Sample the velocities 

  if(istart_typ_cp<3){
     double *vtmp = new double [n];
       sampl1DVelOneT(n,vtmp,m,kT);
       for(int i=0;i<n;i++){v[i].re = vtmp[i];}
       sampl1DVelOneT(n,vtmp,m,kT);
       for(int i=0;i<n;i++){v[i].im = vtmp[i];}
     delete [] vtmp;
  }//endif
  cpSamplNHC(len_nhc,num_nhc,vNHC,mNHC,kT);

//============================================================================
// Scale the velocities
  
  int istrt0 = nkx0_red;
  int istrt  = nkx0_red+nkx0_zero;
  int iend   = nkx0_red+nkx0_uni;

  double temp = 0.0;
  for(int i=0;i<num_nhc;i++){temp += mNHC*vNHC[i][0]*vNHC[i][0];}
  for(int i=istrt0;i<istrt;i++){temp += v[i].getMagSqr()*m[i];}
  for(int i=istrt;i<iend;i++){temp += v[i].getMagSqr()*(2.0*m[i]);}
  for(int i=iend;i<n;i++){temp += v[i].getMagSqr()*m[i];}

  int nfree    = 2*(n-istrt) + num_nhc;
  double scale = sqrt( (((double)(nfree))*kT/temp) );

  for(int i=0;i<num_nhc;i++){vNHC[i][0] *= scale;}
  for(int i=0;i<n;i++){v[i].re *= scale;  v[i].im *= scale;}

//============================================================================
    }// End function
//============================================================================


//==================================================================== 
// Velocities  
//==================================================================== 
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==================================================================== 
void CPINTEGRATE::sampl1DVelOneT(int n, double* v,double* mass,
                                   double  kT)
//=================================================================== 
   {//begin routine 
//===================================================================

  MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();                  
#include "../class_defs/allclass_strip_mdintegrate.h"
   double *qseed  =  &(mdvel_samp->qseed);
   int    *iseed  =  &(mdvel_samp->iseed);
   int    *iseed2 =  &(mdvel_samp->iseed2);

//===================================================================
// Sample unit gaussian

   double *temp  = new double[n];
   double *ptemp = temp-1;

   gaussran(n,iseed,iseed2,qseed,ptemp);
   for(int i=0;i<n;i++){v[i] = temp[i]; }

   delete [] temp;

//===================================================================
// Apply the width

   for(int i=0;i<n;i++){
     double width = sqrt(kT/mass[i]); //kT has boltz
     v[i] *= width;
   }//endfor 

//------------------------------------------------------------------
 } //end routine 
//=================================================================== 

//==================================================================== 
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==================================================================== 
void CPINTEGRATE::cpSamplNHC(int len,int num,double** v,double mass,
                                   double  kT)
//=================================================================== 
   {//begin routine 
//===================================================================

  MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();                  
#include "../class_defs/allclass_strip_mdintegrate.h"
   double *qseed  =  &(mdvel_samp->qseed);
   int    *iseed  =  &(mdvel_samp->iseed);
   int    *iseed2 =  &(mdvel_samp->iseed2);

//===================================================================
// Sample unit gaussian

   int n = num*len;
   double *temp  = new double[n];
   double *ptemp = temp-1;
   gaussran(n,iseed,iseed2,qseed,ptemp);
   int iii=0;
   for(int i=0;i<num;i++){
   for(int j=0;j<len;j++){
     v[i][j]=temp[iii];
     iii++;
   }}//endfor

   delete []temp;

//===================================================================
// Add the width

   double width = sqrt(kT/mass); //kT has boltz
   for(int i=0;i<num;i++){
   for(int j=0;j<len;j++){
     v[i][j] *= width;
   }}

//------------------------------------------------------------------
 } //end routine 
//=================================================================== 

