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


  // allocate z array: It contains nzpt+1 data point!
  z = new double [nzpt+1];
  z[0] = EoccMin;
  z[nzpt] = EoccMax;


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
  
  
  // CASE 3: uneven sampling based on maximum error analysis
  // we need some new parameters
  double wdth;
  wdth = (EoccMax - EoccMin)/Egap;
  double maxerr; // maximum error -> depends on nzpt;

  maxerr = pow ( (1-sqrt(1/(1+wdth)))/nzpt, 2 ); 
  printf("EunoccMin %lg, wdth %lg\n",EunoccMin, wdth);
  printf("maximum error will be %lg. The error is estimated for 1/z function.\n",maxerr);

  for(int i=1; i<nzpt; i++){
    double ztmp = pow( 1-sqrt(maxerr)*(nzpt-i), 2 );
    ztmp = 1/ztmp;
    z[i] = EoccMax - (ztmp - 1)*Egap;
  }
  
  // add 1 to nzpt
  nzpt += 1;

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
