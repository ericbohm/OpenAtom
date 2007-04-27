//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

#include "standard_include.h"
#include "ckcomplex.h"

#include "../class_defs/CP_OPERATIONS/class_cporthog.h"
#include "../proto_defs/proto_math.h"

#define _USE_EISPACK_
//#define _USE_LAPACK_

//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// Control creation of inverse square root of S :  T = S^{-1/2}
//============================================================================

void CPORTHOG::CP_orthoTransform(double *T, const double *S, int nstate) 

//============================================================================
   {/* Begin function */
//----------------------------------------------------------------------------

  int imeth = 1;
  if(nstate>64){imeth=2;}
#ifdef _CP_DEBUG_NO_DIAG_
  imeth = 0;
#endif

  switch(imeth){
    case 0: CPORTHOG::get_unit_Tmat(T,nstate);   //wrong but very fast
            break;
    case 1: CPORTHOG::get_diag_Tmat(S,T,nstate); //correct but slow
            break;
    case 2: CPORTHOG::get_iter_Tmat(S,T,nstate); //fast but finite accuracy
            break;
  }//endswitch

//----------------------------------------------------------------------------
  }/* End function */
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// Set Tmax to the Unit matrix : remove cputime overhead of diag to test
//                               parallel performance
//============================================================================
void CPORTHOG::get_unit_Tmat(double *Tunit,int nstate){
   int nstate_sq = nstate*nstate;
   memset(Tunit,0,nstate_sq*sizeof(double));
   for(int i=0;i<nstate;i++){int ind = i+i*nstate;Tunit[ind] = 1.0;}
}
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// Diagonalize S and use Eigenvalues and Eigenvectors to get S^{-1/2}=T
//============================================================================

void CPORTHOG::get_diag_Tmat(const double *S,double *T,int nstate)

//============================================================================
  {//begin routine
//============================================================================
// I) Get some scratch

   double cpu1,cpu2;
   cputime(&cpu1);

   int nstate_sq     = nstate*nstate;
   double *umat      = new double[nstate_sq];
   double *scr_mat1  = new double[nstate_sq];
   double *scr_mat2  = new double[nstate_sq];
   double *s_eigs    = new double[nstate];
   double *scr1      = new double[3*nstate];        
   double *scr2      = new double[3*nstate];

//==========================================================================
// II. Diagonalize S using rs_ FORTRAN diagonalization routine

  int ifound = 0;
  int ierr   = 0;

  //----------------------------------------------------------------------
  // Use LAPACK : Captain Jack is Happy.
#ifdef _USE_LAPACK_
   ifound ++;
   for(int i = 1; i <= nstate; i++){
   for(int j = 1; j <= i; j++){
     int ind  = (i-1) + (j-1)*nstate;
     int ind2 = (i-1) + (j-1)*(2*nstate-j)/2;
     scr_mat1[ind2] = S[ind];
   }}//endfor
   char Vstuff    ='V';
   char Lstuff    ='L';
   DSPEV(&Vstuff,&Lstuff,&nstate,scr_mat1,s_eigs,umat,&nstate,scr1,&ierr);
#endif

  //----------------------------------------------------------------------
  // Use EISPACK : Captain Jack is Unhappy.
#ifdef _USE_EISPACK_
   ifound ++;
   int info = 1;
   for(int i = 0; i < nstate_sq; i++){ scr_mat1[i] = S[i];}
   RS(&nstate,&nstate,scr_mat1,s_eigs,&info,umat,scr1,scr2,&ierr);
#endif

  //----------------------------------------------------------------------
  // Error checking
   if(ifound!=1 || ierr != 0){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("Error trying to diagonalize S : %d %d\n",ifound,ierr);
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
   }//endif

//==========================================================================
// III. Compute inverse square root of eigenvalues:  Occupation numbers 
//      are HACKED!!!!!

  //----------------------------------------------------------------------
  // A) Construct diagonal matrix of eigenvalues : sqrt(2/lamba)

   for(int i = 0; i < nstate; i++){s_eigs[i] = sqrt(2.0/s_eigs[i]);}
   memset(scr_mat1,0,sizeof(double)*nstate_sq);
   for(int i = 0; i < nstate; i++){
     int ind = i*nstate+i;
     scr_mat1[ind]=s_eigs[i];
   }/* endfor */

  //------------------------------------------------------------------------
  // B) Transform matrix back to original representation using eigenvectors
  //    to create S^{-1/2}

   double alpha = 1.0; double beta = 0.0;
   int itransp  = 0;   int inorm   = 1;

   GENMATMUL(scr_mat1,&nstate,&inorm,umat,&nstate,&itransp,scr_mat2,
             &nstate,&nstate,&nstate,&nstate,&alpha,&beta);
   GENMATMUL(umat,&nstate,&inorm,scr_mat2,&nstate,&inorm,T,
             &nstate,&nstate,&nstate,&nstate,&alpha,&beta);

//============================================================================
// IV) Free allocated temporary memory

   delete [] umat;
   delete [] scr_mat1;
   delete [] scr_mat2;
   delete [] s_eigs;
   delete [] scr1;
   delete [] scr2;

   cputime(&cpu2);
   PRINTF("cpu time diag : %g\n",cpu2-cpu1);

//============================================================================
  } /* End function */
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// Schulz iteration for inverse sqrt root : quadratic convergence!
//============================================================================

void CPORTHOG::get_iter_Tmat(const double *S,double *Titer,int nstate)

//============================================================================
  {//begin routine
//============================================================================
// I) Get some scratch

   double cpu1,cpu2; 
   cputime(&cpu1);

   int nstate_sq     = nstate*nstate;
   double *scr_mat1  = new double[nstate_sq];
   double *scr_mat2  = new double[nstate_sq];
   double *scr_mat3  = new double[nstate_sq];

//============================================================================
// II) Schulz iteration

  //--------------------------------------------------------------------
  // A) Initialize with hacked occupation numbers of 2

   for(int i=0;i<nstate_sq;i++){scr_mat1[i] = S[i]/2.0;}
   memset(Titer,0,nstate_sq*sizeof(double));
   for(int i=0;i<nstate;i++){int ind = i+i*nstate;Titer[ind] = 1.0;}

  //--------------------------------------------------------------------
  // B) Iterate

   int iter        = 0;
   double tol_now  = 1.0;
   while (tol_now > 1.0e-15 && iter<10){

     iter++;
     //--------------------------------
     // scr_mat2 =  3*I - Titer*scr_mat1 
     int itransp  = 0;    int inorm    = 1;
     double alpha = -1.0; double beta  = 1.0;
     memset(scr_mat2,0,nstate_sq*sizeof(double));
     for(int i=0;i<nstate;i++){int ind = i+i*nstate;scr_mat2[ind]=3.0;}
     GENMATMUL(Titer,&nstate,&inorm,scr_mat1,&nstate,&itransp,scr_mat2,
               &nstate,&nstate,&nstate,&nstate,&alpha,&beta);
     //--------------------------------
     // scr_mat1 = 0.5*scr_mat1*scr_mat2 = 0.5*scr_mat3*scr_mat2
     alpha = 0.5;  beta  = 0.0;
     CmiMemcpy(scr_mat3,scr_mat1,nstate_sq*sizeof(double));
     GENMATMUL(scr_mat3,&nstate,&inorm,scr_mat2,&nstate,&itransp,scr_mat1,
               &nstate,&nstate,&nstate,&nstate,&alpha,&beta);
     //--------------------------------
     // Titer = 0.5*scr_mat2*Titer = 0.5*scr_mat2*scr_mat3
     alpha = 0.5;  beta  = 0.0;
     CmiMemcpy(scr_mat3,Titer,nstate_sq*sizeof(double));
     GENMATMUL(scr_mat2,&nstate,&inorm,scr_mat3,&nstate,&itransp,Titer,
               &nstate,&nstate,&nstate,&nstate,&alpha,&beta);
     //--------------------------------
     // tolerence check
     tol_now = 0.0;
     for(int i=0;i<nstate_sq;i++){
       double tmp=scr_mat3[i]-Titer[i];
       tol_now += tmp*tmp;
     }//endfor
     tol_now /= ((double)nstate_sq);
     tol_now = sqrt(tol_now);
     PRINTF("iter %d : tol %g\n",iter,tol_now);
     
   }//endwhile

   if(tol_now>1.0e-15){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("Iterative computation of S^{-1/2} failed in 10 itertions\n");
     PRINTF("The tolerence is now %g\n",tol_now);
     PRINTF("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     EXIT(1);
   }//endif

/*==========================================================================*/
// III) Clean up

   delete [] scr_mat1;
   delete [] scr_mat2;
   delete [] scr_mat3;

   cputime(&cpu2);
   PRINTF("cpu time iter : %g\n",cpu2-cpu1);

/*==========================================================================*/
  }//end routine
/*==========================================================================*/
