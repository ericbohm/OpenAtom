// member functions for GPP calculations
#include <cstdlib>
#include <iostream>
#include "include/ckcomplex.h"
#include "class_defs/gpp.h"
#include "my_fftw.h"
#include "fft_routines.h"
#include "class_defs/sysinfo.h"

//================================================
//        GPPOPTS class member functions               
//================================================
// constructor
GPPOPTS::GPPOPTS(SYSINFO& _sys){

  // get system information
  sys = &_sys;
  
  // read input files
  char fileName[1000];

  sprintf(fileName, "gpp.in");

  readInputFile(fileName);

  // read rho data
  readRho(rhoFile);

  // do fft to get G space density
  fft_R_to_G();
  
#ifdef DEBUG
  printf("rhog(0) data: %lg  %lg\n", rhoData[0].re, rhoData[0].im);
#endif  
}





// read input files
void GPPOPTS::readInputFile(char *fromFile){

  FILE *fp = fopen(fromFile,"r");

  fscanf(fp, "%d", &gpp_is_on);
  fscanf(fp, "%d", &qespresso);
  fscanf(fp, "%s", &rhoFile);
  fscanf(fp, "%d", &numq);
  fscanf(fp, "%d", &numEnergy);
  Energy = new double [numEnergy];
  for(int i=0; i<numEnergy; i++){
    fscanf(fp,"%lg",&Energy[i]);
  }
  fclose(fp);  
}





// read density data from the file (rho.dat as default)
// density is saved in real-space
void GPPOPTS::readRho(char *fromFile){
    
  FILE *fp = fopen(fromFile,"r");
  // read number of grids 
  fscanf(fp, "%d %d %d", &nr[0], &nr[1], &nr[2]);
  int ndata = nr[0]*nr[1]*nr[2];
  // malloc
  ga = new int [ndata];
  gb = new int [ndata];
  gc = new int [ndata];
  rhoData = new complex [ndata];

  // scale rhoData if qespresso is true
  double scale;
  if(qespresso)
    scale = sys->vol/double(ndata); 
  
  else
    scale = 1.0; 

  // read data
  int counter=0;
  for(int i=0; i<nr[0]; i++){
    for(int j=0; j<nr[1]; j++){
      for(int k=0; k<nr[2]; k++){
      	fscanf(fp, "%lg", &(rhoData[counter].re));
      	if( qespresso ){ rhoData[counter].re *= scale; }
      	counter += 1;
      }
    }
  }
#ifdef DEBUG
  double dsum = 0;
  for(int i=0; i<ndata; i++){
    dsum += rhoData[i].re;
  }
  printf(" sum of rho : %lg\n", dsum);
#endif
}





void GPPOPTS::fft_R_to_G(){

  // create 3D FFTW plan
  int forward = -1;
  setup_fftw_3d(nr,forward);

  // in_pointer has been declared in my_fftw.h 
  put_into_fftbox(nr, rhoData, in_pointer);
  
  // call fftw
  do_fftw();

  // out_pointer has been declared in my_fftw.h
  int ndata = nr[0]*nr[1]*nr[2];
  double scale  = 1.0;  // the scale actually doesn't matter to the final results because we want rho(g-g')/rho(0)
  fftbox_to_array(ndata, out_pointer, rhoData, scale);

  // get g index
  fftidx_to_gidx(ga, gb, gc, nr);

}








//===============================================
//        GPPDATA class member functions          
//===============================================

// S(i,j) = sqrt(Vc_i) * { S(i,j) - delta_i,j } * sqrt(Vc_j)
void GPPDATA::calc_Smtrx(){
  
  // set S matrix
  for(int i=0; i < S->nrow; i++){
    for(int j=0; j < S->ncol; j++){
      double Vcsq_i = sqrt( vc->data[i] );
      double Vcsq_j = sqrt( vc->data[j] );
      if(i==j){ S->get(i,j) = S->get(i,j) - double(1); }
      S->get(i,j) = Vcsq_i * S->get(i,j) * Vcsq_j;
    }// end column
  }// end row
  
}



#include "include/mylapack.h"
#include "util.h"

void GPPDATA::eigenDecomp(){

#ifdef USE_LAPACK
  char jobz = 'V';  // compute eigenvalues and eigenvectors
  char uplo = 'L';  // Lower triangle
  int n = ng;       // number of g vectors (=ncol=nrow)
  complex* WORK;
  int lwork = -1;
  double* RWORK;
  int info;
  int lwmax = n*250;

#ifdef DEBUG
  CMATRIX input = CMATRIX(ng, ng);
  for(int i=0; i<ng*ng; i++){
    input.m[i] = S->m[i];
  }
  mymessage("S matrix is safely saved");
#endif

  // m'alloc
  WORK = new complex [lwmax];
  RWORK = new double [3*n-2];

  // since LAPACK is column major, we need transpose
  S->transpose(); 
  
  // 1. inquire optimal size of lwork
  myHEEV(&jobz, &uplo, &n, S->m, &n, eigval, WORK, &lwork, RWORK, &info);
  // find minimum value
  if ( lwmax >= int( WORK[0].re ) ){
    lwork = int( WORK[0].re );
  }
  else{
    lwork = lwmax;
  }
  // 2. Do eigendecomposition
  
  // S will contain eigenvectors!
  myHEEV(&jobz, &uplo, &n, S->m, &n, eigval, WORK, &lwork, RWORK, &info);

  if ( info == 0 ){
    mymessage("eigendecomposition done successfully");
  }
  if ( info > 0 ){
    mymessage("eigendecomposition failed to converge");
  }
  if ( info < 0 ){
    char *msg;
    sprintf(msg, "%d th element has illigal value", -info);
    mymessage(msg);
  }

#ifdef DEBUG
  /*
  mymessage("printintg eigenvalues");
  for(int i=0; i<n; i++){
    printf("%lg\n",eigval[i]);
  }
  */
#endif


  
#else

  Die("Eigendecomposition requires LAPACK routine. Set USE_LAPACK flag on");

#endif
  // transpose again... 
  S->transpose();
  // so now each column of S matrix is an eigenvector!
  
#ifdef DEBUG
  /*
  // let's print some values   !! CHECK IF IT REALLY GIVES THE ORIGINAL S MATRIX !!  ->  yes it does!!!
  printf("Now we calculate S matrix from eigenvalues and eigenvectors and original S\n");
  int nprint=5;
  for(int g1=0; g1<nprint; g1++){
    for(int g2=0; g2<nprint; g2++){
      complex Sggp=0;
      for(int alpha=0; alpha<ng; alpha++){
	Sggp += S->get(g1,alpha) * eigval[alpha] * S->get(g2,alpha).conj();
      }
      printf("(%d,%d) %lg %lg   vs.   %lg %lg\n",g1,g2,Sggp.re, Sggp.im, input.get(g1,g2).re, input.get(g1,g2).im);
    }
  }
  */
#endif

}





#include "mtrxop_3x3.h"
void GPPDATA::calc_Omsq(){

  // calculate plasma frequency square
  // \omega_{pl}^2 = 4*\pi*\rho(0)   rho(0) = number of electrons/volume
  double Wpl2;
  Wpl2 = double(4) * M_PI * opt->rhoData[0].re/opt->sys->vol;

#ifdef DEBUG
  mymessage("printing Omega_pl^2");
  printf("%lg\n", Wpl2);
#endif
  
  double factor;
  factor = -Wpl2 * double(4) * M_PI / (opt->sys->vol * opt->sys->nkpt);

  // change q vector to cartesian coordiates
  double qcart[3];
  for(int i=0; i<3; i++){
    if(qIndex == 0){
      qcart[i] = vc->sys.b1[i]*(vc->sys.qvec[qIndex][0] + vc->sys.shift[0])
	+ vc->sys.b2[i]*(vc->sys.qvec[qIndex][1]+vc->sys.shift[1])
	+ vc->sys.b3[i]*(vc->sys.qvec[qIndex][2]+vc->sys.shift[2]);
    }
    else{
      qcart[i] = vc->sys.b1[i]*vc->sys.qvec[qIndex][0] + vc->sys.b2[i]*vc->sys.qvec[qIndex][1] + vc->sys.b3[i]*vc->sys.qvec[qIndex][2];
    }
    qcart[i] *= double(2) * M_PI / opt->sys->alat;
  }
  
#ifdef DEBUG
  printf("At this %d index, q vector is %lg %lg %lg\n", qIndex, vc->sys.qvec[qIndex][0], vc->sys.qvec[qIndex][1], vc->sys.qvec[qIndex][2]);
  printf("and their cartesian coordinates: %lg %lg %lg\n",qcart[0], qcart[1], qcart[2]);
#endif

  CMATRIX Mggp(ng, ng); //save (q+G)*(q+G')/(|q+G|^2*|q+G'|^2)*rho(G-G')/rho(0)

  for(int ig1=0; ig1<ng; ig1++){
    for(int ig2=0; ig2<ng; ig2++){
      
      // g vector at ig1, ig2
      double g1cryst[3], g2cryst[3];

      // g vector in cartesian coordinates at ig1, ig2
      double g1[3], g2[3];
      g1cryst[0] = double( vc->ga[ig1] );
      g1cryst[1] = double( vc->gb[ig1] );
      g1cryst[2] = double( vc->gc[ig1] );
      g2cryst[0] = double( vc->ga[ig2] );
      g2cryst[1] = double( vc->gb[ig2] );
      g2cryst[2] = double( vc->gc[ig2] );

      // change g1 and g2 to cartesian coordiates
      for(int i=0; i<3; i++){
	g1[i] = vc->sys.b1[i]*g1cryst[0] + vc->sys.b2[i]*g1cryst[1] + vc->sys.b3[i]*g1cryst[2];
	g2[i] = vc->sys.b1[i]*g2cryst[0] + vc->sys.b2[i]*g2cryst[1] + vc->sys.b3[i]*g2cryst[2];
	g1[i] *= 2*M_PI/vc->sys.alat;
	g2[i] *= 2*M_PI/vc->sys.alat;
      }
      
      // calculate q+G, q+G'
      for(int i=0; i<3; i++){
	g1[i] += qcart[i];
	g2[i] += qcart[i];
      }

      // find rho(G-G')
      int gdiff[3];
      gdiff[0] = vc->ga[ig1] - vc->ga[ig2];
      gdiff[1] = vc->gb[ig1] - vc->gb[ig2];
      gdiff[2] = vc->gc[ig1] - vc->gc[ig2];

      // index for (G-G') vector
      int gdiffIndex;
      bool gdiffTrue;
      // find the gdiffIndex
      find_rho_gdiff(gdiff, gdiffIndex, gdiffTrue);

      // calculate M_GG' matrix element
      if( gdiffTrue ){
	Mggp.get(ig1,ig2).re = dot_product(g1,g2) / (dot_product(g1,g1)*dot_product(g2,g2)) * opt->rhoData[gdiffIndex].re / opt->rhoData[0].re;
      }
      else{
	Mggp.get(ig1,ig2).re = double(0);
      }
      
    }//end ig2
  }//end ig1

  // do matrix multiplication
  // omega^2 = V^-1 * Mggp * V
  // eigenvectors are stored in S (columns are eigenvector)

  // call LAPACK

  char transformT = 'T'; // transpose
  char transform = 'N';  // do nothing
  double Lalpha = double(1.0);  // scale A*B by this factor
  double Lbeta = double(0.0); // scale initial value of C by this factor
  CMATRIX M1(ng, ng);    // temporary space

  // 1. M1 = Mggp * V
  S->transpose();
  myGEMM( &transformT, &transform, &ng, &ng, &ng, &Lalpha, Mggp.m, &ng, S->m, &ng, &Lbeta, M1.m, &ng);

  // 2. Mggp = V^-1 * M1
  S->ctranspose();
  myGEMM( &transform, &transform, &ng, &ng, &ng, &Lalpha, S->m, &ng, M1.m, &ng, &Lbeta, Mggp.m, &ng);

  for(int i=0; i<ng; i++){
    omsq[i] = Mggp.get(i,i).re * factor / eigval[i];

    // print omega^2
#ifdef DEBUG
    /*
    if(i==0){mymessage("Printing Omsq value\n");}
    printf("omsq: %lg\n",omsq[i]);
    */
#endif
  }

  // NOTE:
  // now S contains eigenvectors (but conjugated) in columns
    
    
}





void GPPDATA::find_rho_gdiff(int gdiff[3], int &gindex, bool &gdiffTrue){

  int ndata = opt->nr[0]*opt->nr[1]*opt->nr[2];
  gdiffTrue = false;
  
  for(int i=0; i<ndata; i++){
    //we need to compare gdiff with density g index
    if( gdiff[0] == opt->ga[i] && gdiff[1] == opt->gb[i] && gdiff[2] == opt->gc[i] ){
      gindex = i;
      gdiffTrue = true;
    }
  }
}




void GPPDATA::calc_Sw(int iw){

  complex E; // evaluation frequency
  const double Hartree = 27.211385;
  E.re = opt->Energy[iw]/Hartree;  // input value is eV
  E.im = 0.2/Hartree; // broadening factor
  // TODO  broadening factor should be an input value. Change it soon.

  // loop over row (G)
  for(int r=0; r<ng; r++){
    // loop over column (G')
    for(int c=0; c<ng; c++){
      // loop over eigenvalues
      for(int a=0; a<ng; a++){
	// 1/(omsa-E) = tmp
	complex tmp;
	tmp.re = omsq[a] - (E.re)*(E.re) + (E.im)*(E.im);
	tmp.re = tmp.re / ( (omsq[a]-E.re*E.re+E.im*E.im)*(omsq[a]-E.re*E.re+E.im*E.im)+4*E.re*E.re*E.im*E.im );
	tmp.im = 2*E.re*E.im / ( (omsq[a]-E.re*E.re+E.im*E.im)*(omsq[a]-E.re*E.re+E.im*E.im)+4*E.re*E.re*E.im*E.im );
	Sw[iw]->get(r,c) += S->get(r,a) * eigval[a] * omsq[a] * tmp * S->get(c,a).conj();
	/*
#ifdef DEBUG
	// printing out values
	if(r==0 && c==0){
	  if (a==0 | a==1 | a==2 | a==4 | a==ng-1){
	    double swreal = Sw[iw]->get(r,c).re*(1/vc->data[0])+1;
	    double swimag = Sw[iw]->get(r,c).im*(1/vc->data[0]);
	    printf("%.12f %.12f\n", swreal, swimag);
	  }
	}
#endif	
	*/
      }
    }
  }
  /*
#ifdef DEBUG
  // let's print some values
  for(int r=0; r<3; r++){
    for(int c=0; c<3; c++){
      double Vr = 1 / vc->data[r];
      double Vc = 1 / vc->data[c];
      Vr = sqrt( Vr );
      Vc = sqrt( Vc );
      Sw[iw]->get(r,c) *= Vr*Vc;
      if(r==c){ Sw[iw]->get(r,c) += 1;}
      printf("S(w) at g=%d, g'=%d : %lg, %lg\n", r, c, Sw[iw]->get(r,c).re, Sw[iw]->get(r,c).im);
    }
  }
#endif
  */	

}
