//===================================================
//  Check alegbra for finite difference scheme
//===================================================
//***************************************************
//===================================================
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main();
void taylor3(double *,double ,double ,double ,double ,double );

//===================================================
//***************************************************
//===================================================
int main(){
//===================================================
  double x,v,a,j,dt;
  double x1,x2,x3;
  double dtm1,dtm2,dtm3;
  double vp;
  double c0,c1,c2,c3;
  double s0,s1,s2,s3,s4;
//===================================================
// Tell the user what's going down

  printf("\n");
  printf("================================================\n");
  printf("Checking 1st and 2nd order finite difference   \n");
  printf("methods to extract the velocity at time, t, from\n");
  printf("positions at times t, t-dt, t-2dt.\n");
  printf("================================================\n");
  printf("\n");

//===================================================
// position, velocity, acceleration, jerk, time step

  x    = 0.145714;
  v    = 0.156721;
  a    = 0.174155;
  j    = 0.185672;
  dt   = 0.671457;

  dtm1 = -dt;
  dtm2 =  2.0*dtm1;
  dtm3 =  3.0*dtm1;

//========================================================
// Create positions 1,2,3 time steps back by 3rd order Taylor Series

  taylor3(&x1,dtm1,x,v,a,j);
  taylor3(&x2,dtm2,x,v,a,j);
  taylor3(&x3,dtm3,x,v,a,j);

//===================================================
// Compute velocity by finite difference to 2nd order

  // 2nd order coefficients
  c0 = 11.0/6.0;
  c1 = -3.0;
  c2 =  1.5;
  c3 = -1.0/3.0;
  
  // check cancelation of each power in dt
  s0 =       c3+    c2+c1+c0;  // = 0  x
  s1 =  -3.0*c3-2.0*c2-c1;     // = 1  v
  s2 =  +9.0*c3+4.0*c2+c1;     // = 0  a
  s3 = -27.0*c3-8.0*c2-c1;     // = 0  j
  s4 =  81.0*c3+16.0*c2+c1;    //!= 0  error term

  // velocity computation
  vp = (c3*x3 + c2*x2 +c1*x1 + c0*x)/dt;

  printf("================================================\n");
  printf("   2nd order: c={%g,%g,%g,%g}\n",c0,c1,c2,c3);
  printf("------------------------------------------------\n");
  printf("2nd order coef chk: %.10g %.10g %.10g %.10g %.10g\n",s0,s1,s2,s3,s4/24.0);
  printf("2nd order vel  chk: %.10g %.10g\n",v,vp);
  printf("================================================\n");
  printf("\n");

//===============================================================
// Create positions 1,2 time steps back by 2nd order Taylor Series

  j = 0;  // zero 3rd order term and use 3rd order routine
  taylor3(&x1,dtm1,x,v,a,j);
  taylor3(&x2,dtm2,x,v,a,j);

//===================================================
// Compute velocity by finite difference to 1st order

  // 1st order coefficients
  c0 =  1.5;
  c1 = -2.0;
  c2 =  0.5;
  
  // check cancelation of each power in dt
  s0 =      c2+c1+c0;       // = 0  x
  s1 = -2.0*c2-c1;          // = 1  v
  s2 =  4.0*c2+c1;          // = 0  a
  s3 = -8.0*c2-c1;          //!= 0  error term

  // velocity computation
  vp = (c2*x2 +c1*x1 + c0*x)/dt;

  printf("================================================\n");
  printf("   1st order: c={%g,%g,%g}\n",c0,c1,c2);
  printf("------------------------------------------------\n");
  printf("1st order coef chk: %.10g %.10g %.10g %.10g\n",s0,s1,s2,s3/6.0);
  printf("1st order vel  chk: %.10g %.10g\n",v,vp);
  printf("================================================\n");

//--------------------------------------------------
  return 1;
//===================================================
}/*end main*/
//===================================================


//===================================================
//***************************************************
//===================================================
void taylor3(double *xn,double dt,double x,double v,double a,double j){
//===================================================
  double twfac=2.0; double thfac=2.0*3.0;

  xn[0] = x + v*dt + a*dt*dt/twfac + j*dt*dt*dt/thfac;

//===================================================
  }/*end taylor3*/
//===================================================
