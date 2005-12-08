#include "standard_include.h"
#include "ckcomplex.h"

#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_mdintegrate.h"

#include "../class_defs/CP_OPERATIONS/class_cpintegrate.h"
#include "../../mathlib/proto_math.h"

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPINTEGRATE::initCPNHC(int npts,int maxLen,int maxNum, int *len_nhc_ret,
                        int *num_nhc_ret, double *kTCP_ret,
                        double *tau_ret,double *mNHC_ret,double *degfree_ret,
                        double *degfreeNHC_ret,double *gammaNHC_ret,int ncoef_true,
                        int ncoef_zero)
//============================================================================
  {// Begin Function 
//============================================================================

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP *cp = CP::get();
#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

  int len_nhc   = cptherm_info->len_c_nhc;
  int num_nhc   = 1;
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
  if(len_nhc<2){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Minimum NHC len_nhc is 2 > %d\n",maxNum,num_nhc);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
  }//endif

  double degfree    = (double)(2*ncoef_true-ncoef_zero+num_nhc-1);
  double degfreeNHC = (double)( (len_nhc-1)*num_nhc );
  double gammaNHC   = degfree/(degfree+1.0);

//============================================================================

  (*len_nhc_ret) = len_nhc;
  (*num_nhc_ret) = num_nhc;
  (*kTCP_ret)    = kT;
  (*tau_ret)     = tau;
  
  double pre     = tau*tau*kT;
  (*mNHC_ret)    = pre;

  (*degfree_ret)    = degfree;
  (*degfreeNHC_ret) = degfreeNHC;
  (*gammaNHC_ret)   = gammaNHC;

//============================================================================
   }// End function
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CPINTEGRATE::CPSmplVel(int n,double *m,complex *v,int len_nhc,int num_nhc,
                            double mNHC,double **vNHC,double kT,int istart_typ_cp,
                            int nkx0_red,int nkx0_uni,int nkx0_zero,
                            double degfree,double degfreeNHC,double gammaNHC)
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

  for(int i=istrt0;i<istrt;i++){v[i].im=0.0;} // zero imaginary part of g=0

  double temp = 0.0;
  for(int i=0;i<num_nhc;i++)   {temp += gammaNHC*mNHC*vNHC[i][0]*vNHC[i][0];}
  for(int i=istrt0;i<istrt;i++){temp += v[i].getMagSqr()*m[i];}
  for(int i=istrt;i<iend;i++)  {temp += v[i].getMagSqr()*(2.0*m[i]);}
  for(int i=iend;i<n;i++)      {temp += v[i].getMagSqr()*m[i];}

  double scale = sqrt( (degfree*kT)/temp );

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
