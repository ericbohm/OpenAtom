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
void CPINTEGRATE::CPSmplVel(int n,double *m,complex *v,int len_nhc,int num_nhc,int nck_nhc,
                            double *mNHC,double ***vNHC,double ***xNHC, double ***xNHCP,
                            double *a2NHC, double kT,int istart_typ_cp,
                            int nkx0_red,int nkx0_uni,int nkx0_zero,
                            double degfree,double degfreeNHC)
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
  cpSamplNHC(len_nhc,num_nhc,nck_nhc,vNHC,xNHC,xNHCP,mNHC,a2NHC,kT);

//============================================================================
// Scale the velocities
  
  if(istart_typ_cp<3){
    int istrt0 = nkx0_red;
    int istrt  = nkx0_red+nkx0_zero;
    int iend   = nkx0_red+nkx0_uni;

    for(int i=istrt0;i<istrt;i++){v[i].im=0.0;} // zero imaginary part of g=0

    double temp = 0.0;
    for(int i=istrt0;i<istrt;i++){temp += v[i].getMagSqr()*m[i];}
    for(int i=istrt;i<iend;i++)  {temp += v[i].getMagSqr()*(2.0*m[i]);}
    for(int i=iend;i<n;i++)      {temp += v[i].getMagSqr()*m[i];}

    double scale = sqrt( (degfree*kT)/temp );
    for(int i=0;i<n;i++){v[i].re *= scale;  v[i].im *= scale;}
  }//endif

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
   long   *iseed  =  &(mdvel_samp->iseed);
   long   *iseed2 =  &(mdvel_samp->iseed2);

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
void CPINTEGRATE::cpSamplNHC(int len,int num,int nck,
                             double*** v,double ***x,double ***xp,
                             double *mass, double *a2,  double  kT)
//=================================================================== 
   {//begin routine 
//===================================================================

  MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();                  
#include "../class_defs/allclass_strip_mdintegrate.h"
   double *qseed  =  &(mdvel_samp->qseed);
   long   *iseed  =  &(mdvel_samp->iseed);
   long   *iseed2 =  &(mdvel_samp->iseed2);

//===================================================================
// Sample unit gaussian

   int n = num*len*nck;
   double *temp  = new double[n];
   double *ptemp = temp-1;

   gaussran(n,iseed,iseed2,qseed,ptemp);
   int iii=0;
   for(int k=0;k<nck;k++){
   for(int i=0;i<num;i++){
   for(int j=0;j<len;j++){
     v[k][i][j]=temp[iii];
     iii++;
   }}}//endfor

   gaussran(n,iseed,iseed2,qseed,ptemp);

   iii=0;
   for(int k=0;k<nck;k++){   
   for(int i=0;i<num;i++){
   for(int j=0;j<len;j++){
     xp[k][i][j]=temp[iii];
     iii++;
   }}}//endfor

   delete []temp;

//===================================================================
// Add the width

   for(int k=0;k<nck;k++){
   for(int j=0;j<len;j++){
     double width = sqrt(kT/mass[j]); //kT has boltz
     for(int i=0;i<num;i++){
       v[k][i][j]  *= width;
       xp[k][i][j] *= sqrt(a2[i]);
     }//endfor
   }}//endfor

//===================================================================
// Apply isokinetic constraint to each dimer

   if(len!=2){
     PRINTF("@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@\n");
     PRINTF("Deprecated len_nhc_cp %d \n",len);
     PRINTF("len_nhc_cp must be equal to 2 now\n");
     PRINTF("@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@\n");
     EXIT(1);
   }//endif

   for(int k=0;k<nck;k++){
   for(int i=0;i<num;i++){
     double ekin = mass[0]*v[k][i][0]*v[k][i][0]
                 + mass[1]*v[k][i][1]*v[k][i][1];
     double sc   = sqrt(kT/ekin);
     x[k][i][0] = 0.0;
     v[k][i][0] *= sc;
     v[k][i][1] *= sc;
   }}//endfor

//------------------------------------------------------------------
  } //end routine 
//=================================================================== 
