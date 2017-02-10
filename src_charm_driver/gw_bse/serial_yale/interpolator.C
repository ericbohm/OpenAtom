// member functions for P interpolation
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "include/ckcomplex.h"
#include "class_defs/interpolator.h"
#include "util.h"


//=================================================
//    INTERPOLATOR class member function
//=================================================

void INTERPOLATOR::make_zlist(){

  double margin = 0.00000001; // in Hartree Unit, (QE's energy unit is in Hartree)
  EoccMin = EoccMin - abs(EoccMin*margin);
  EoccMax = EoccMax + abs(EoccMax*margin);

  /* for case 1 to 3
  // allocate z array: It contains nzpt+1 data point!
  z = new double [nzpt+1];
  z[0] = EoccMin;
  z[nzpt] = EoccMax;
  */

  /* 
  // CASE 1: EVEN sampling on x-axis (not so good)

  // calculate interval
  double delta;
  delta = ( EoccMax - EoccMin ) / double(nzpt);

  for(int i=1; i<nzpt; i++){
    z[i] = z[0] + double(i) * delta;
  }
  */
  
  /* 
  // CASE 2: uneven sampling on x-axis, but even in y-axis
  // based on 1/(EunoccMin-z) function
  double vspace;

  vspace = ( 1/(EunoccMin-z[nzpt]) - 1/(EunoccMin-z[0]) )/nzpt;
  printf("vertical spacing: %lg\n",vspace);

  
  for(int i=1; i<nzpt; i++){
    z[i] = EunoccMin - 1/( 1/(EunoccMin-z[i-1]) + vspace );
  }
  */
  
  /*
  // CASE 3: uneven sampling based on maximum error analysis
  // we need some new parameters
  double wdth;
  wdth = (EoccMax - EoccMin)/Egap;
  double maxerr; // maximum error -> depends on nzpt;

  maxerr = pow ( (1-sqrt(1/(1+wdth)))/nzpt, 2 ); 

  for(int i=1; i<nzpt; i++){
    double ztmp = pow( 1-sqrt(maxerr)*(nzpt-i), 2 );
    ztmp = 1/ztmp;
    z[i] = EoccMax - (ztmp - 1)*Egap;
  }
  */
  // add 1 to nzpt
  nzpt += 1;

  
  // CASE 4: uneven sampling with weight of the spacing (exponential grid)
  // now we want to find nzpt based on the % error analysis
  double ptol; // % tolerance for sum_c 1/(Ec-Ev)
  ptol = 100;
  int satisfied; // whether we achieved the tolerance or not
  satisfied = 0; // initial value = false
  double wdth;
  wdth = (EoccMax - EoccMin);

#ifdef DEBUG
  printf("percentage error tolerance for sum_c 1/(Ec-Ev)\n");
#endif

  // starting at nzpt = 2 if linear interpolation
  // starting at nzpt = 3 for quadratic interpolation

  // outer loop for znpt
  // z[0] = 1, z[N]= 1 + wdth/Egap
  // total number of points = N+1 = nzpt

#ifdef DEBUG
  printf("we are doing linear interpolation\n");
#endif

  for(int i=0; i<10000; i++){
    satisfied = 0; // initially, tolerance is NOT satisfied
    nzpt = 2 + i; // for linear interpolation

    // setting z[] array
    int N; // temporary number for number of spacing
    N = nzpt-1; 
    z = new double [nzpt];
    double A; 
    A = log(1+wdth/Egap)/N;
    for(int iz=0;iz<N+1;iz++){
      z[iz] = exp(A*(N-iz));
    }
    // transform z[] to real energy value
    for(int iz=0;iz<nzpt;iz++){
      z[iz] = EoccMax - (z[iz]-1)*Egap;
    }
    // set the first and last value again
    z[0] = EoccMin;
    z[N] = EoccMax;

    // now we need to test if this nzpt is enough to get desirable error
    double* exact; // exact sum (i.e., sum_c 1/(Ec-Ev) at this Ev)
    double* errEvi; // percentage error at this Ev 
    double* zdata; // sum_c 1/(Ec-Ez) at this Ez
    
    errEvi = new double [nocc];


    // a. calculate exact value
    exact = new double [nocc];
    for(int iv=0;iv<nocc;iv++){
      exact[iv]=0; // initiallization
      // loop over conduction bands
      for(int ic=0;ic<nunocc;ic++){
	double Ec = psiR->eig[nocc+ic];
	double Ev = psiRq->eig[iv];
	exact[iv] += 1/(Ec-Ev);
      }
    }

    // b. evalute sum_c 1/(Ec-Ez) for interpolation
    zdata = new double [nzpt];
    for(int iz=0;iz<nzpt;iz++){
      zdata[iz] = 0;
      for(int ic=0;ic<nunocc;ic++){
	double Ec = psiR->eig[nocc+ic];
	double Ez = z[iz];
	zdata[iz] += 1/(Ec-Ez);
      }
    }

    // c. let's interpolate
    for(int iv=0; iv<nocc; iv++){
      double Ev = psiRq->eig[iv];
      for(int iz=0; iz<nzpt; iz++){
	if( Ev == z[iz] ){
	  errEvi[iv] = ( zdata[iz] - exact[iv] ) / exact[iv] * 100;
	  break;
	}
	if( z[iz] < Ev && Ev < z[iz+1] ){
	  double x0,x1,y0,y1;
	  x0 = z[iz];
	  x1 = z[iz+1];
	  y0 = zdata[iz];
	  y1 = zdata[iz+1];
	  double intpval; // interpolated value at this Ev
	  intpval = (y1-y0)/(x1-x0)*(Ev-x0) + y0;
	  errEvi[iv] = abs( intpval - exact[iv] ) / exact[iv] *100;
	  break;
	}
      }
      if( errEvi[iv] < ptol){
	satisfied = 1;
      }
      else{
	satisfied = 0;
	break; // break this loop and try next nzpt
      }
      
    }
    if ( satisfied == 1 ){
#ifdef DEBUG
      printf("We got nzpt! nzpt for %lg %% is %d\n",ptol,nzpt);
      // let's print out the largest error with this number of points
      double errmax;
      for(int iv=0;iv<nocc;iv++){
	if( errmax < errEvi[iv] ){
	  errmax = errEvi[iv];
	}
      }
      printf("The maximum %% error for sum_c 1/(Ec-Ev) is %lg\n",errmax);
#endif      
      break;
    }

  }// end loop for searching nzpt which gives the desired tolerance



#ifdef DEBUG
  //print_values();
  printf("EoccMin: %lg  EoccMax: %lg\n",EoccMin, EoccMax);
  printf("energy values to be evaluated:\n");
  for(int i=0;i<nzpt;i++){
    printf("z[%d]: %lg\n",i,z[i]);
  }
#endif


}

void INTERPOLATOR::initialize(int izpt){
  for(int ir=0; ir<nr; ir++){
    colG[izpt][ir] = 0;
  }
}

void INTERPOLATOR::G_calculator(int ic, int izpt){

  //initialize colG[izpt][:]
  initialize(izpt);

  for(int i=0; i<nunocc; i++){
    for(int ir=0; ir<nr; ir++){
      colG[izpt][ir] += psiR->coeff[i+nocc][ir].conj() * psiR->coeff[i+nocc][ic]
       	/ (psiR->eig[i+nocc] - z[izpt]) ;
    }
  }
  
}


void INTERPOLATOR::LinearInterpolation(double Eocc, complex* colG_Eocc){

  // output is saved to colG_Eocc
  for(int iz=0; iz<nzpt; iz++){
    if( Eocc == z[iz] ){
      for(int ir=0; ir<nr; ir++){
	colG_Eocc[ir] = colG[iz][ir];
      }
      break;
    }
    if( z[iz] < Eocc  &&  Eocc < z[iz+1] ){
      for(int ir=0; ir<nr; ir++){
	complex y0, y1;
	double x0, x1;
	y0 = colG[iz][ir];
	y1 = colG[iz+1][ir];
	x0 = z[iz];
	x1 = z[iz+1];
	colG_Eocc[ir] = (y1-y0)/(x1-x0)*(Eocc-x0) + y0;
      }
      break;
    }
  }
  
}



void INTERPOLATOR::QuadraticInterpolation(double Eocc, complex* colG_Eocc){

  // output is saved to colG_Eocc
  if(nzpt<3){
    Die("Quadratic Interpolation requires more than 3 points.\n Increase the number of points.");
  }
  for(int iz=0; iz<nzpt; iz++){
    if( Eocc == z[iz] ){
      for(int ir=0; ir<nr; ir++){
	colG_Eocc[ir] = colG[iz][ir];
      }
      break;
    }
    if( z[iz] < Eocc  &&  Eocc < z[iz+1] ){
      // check all cases
      double x0, x1, x2;
      complex y0, y1, y2;
      for(int ir=0; ir<nr; ir++){
	if( iz < 2 ){
	  x0 = z[0];
	  x1 = z[1];
	  x2 = z[2];
	  y0 = colG[0][ir];
	  y1 = colG[1][ir];
	  y2 = colG[2][ir];
	}
	else if( iz > nzpt-3 ){
	  x0 = z[nzpt-3];
	  x1 = z[nzpt-2];
	  x2 = z[nzpt-1];
	  y0 = colG[nzpt-3][ir];
	  y1 = colG[nzpt-2][ir];
	  y2 = colG[nzpt-1][ir];
	}
	else{
	  double test = (z[iz] + z[iz+1]) / 2;
	  if ( Eocc < test ){
	    x0 = z[iz-1];
	    x1 = z[iz];
	    x2 = z[iz+1];
	    y0 = colG[iz-1][ir];
	    y1 = colG[iz][ir];
	    y2 = colG[iz+1][ir+1];
	  }
	  else{
	    x0 = z[iz];
	    x1 = z[iz+1];
	    x2 = z[iz+2];
	    y0 = colG[iz][ir];
	    y1 = colG[iz+1][ir];
	    y2 = colG[iz+2][ir];
	  }
	}
	// finally, interpolate!
	colG_Eocc[ir] = (Eocc-x1)*(Eocc-x2)/((x0-x1)*(x0-x2))*y0 +
	  (Eocc-x0)*(Eocc-x2)/((x1-x0)*(x1-x2))*y1 +
	  (Eocc-x0)*(Eocc-x1)/((x2-x0)*(x2-x1))*y2;
      }//end ir loop
      break;
    }//end if
  }//end iz loop
}	



void INTERPOLATOR::print_G(int col){

  printf("***");
  printf("Printing out column of G at %d th column\n",col);

  char fname[100];
  sprintf(fname, "G_column_%d",col);
  FILE* fp = fopen(fname,"w");

  for(int ir=0; ir<nr; ir++){
    // printing real values
    for(int iz=0; iz<nzpt; iz++){
      fprintf(fp,"%lg  ", colG[iz][ir].re);
    }
    fprintf(fp,"\n");
    // printing imaginary values
    for(int iz=0; iz<nzpt; iz++){
      fprintf(fp,"%lg  ", colG[iz][ir].im);
    }
    fprintf(fp,"\n");

  }
  fclose(fp);

}



void INTERPOLATOR::print_values(){
  mymessage("printing interpolator class values");
  printf("q index: %d\n", thisqIndex);
  printf("number of z point: %d\n", nzpt);
  printf("EoccMax/EoccMin:  %lg,  %lg\n",EoccMax, EoccMin);
  printf("EunoccMin:  %lg\n",EunoccMin);
  printf("number of occupied/unoccupied states: %d, %d\n", nocc, nunocc);
  printf("number of real-space grid for states: %d\n", nr);
}






void INTERPOLATOR::G_eachband(int r1, int r2){
#ifdef VERBOSE
  printf("calculate G for each conduction band at (%d, %d)\n", r1,r2);
#endif
  char fname[100];
  sprintf(fname,"./Gdata_q1_k1/row%d_col%d.dat",r1,r2);

  FILE* fp = fopen(fname,"w");

  fprintf(fp,"#z(Hartree)");
  for(int ic=0; ic<nunocc; ic++){
    fprintf(fp,"%d    ",ic);
  }
  fprintf(fp,"\n");
  
  for(int iz=0; iz<nzpt; iz++){
    fprintf(fp,"%lg  ",z[iz]);
    for(int ic=0; ic<nunocc; ic++){
      complex data;
      data = psiR->coeff[ic+nocc][r1].conj() * psiR->coeff[ic+nocc][r2]
       	/ (psiR->eig[ic+nocc] - z[iz]) ;
      fprintf(fp,"%lg  ", data.re);
    }
    fprintf(fp,"\n");
  }

  fclose(fp);
  
  
}
