//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
class PIBeadAtoms{
  public:
    ~PIBeadAtoms(){}
    PIBeadAtoms(int );

    int numBeads;
    void compute_PIMD_Fu();
    void compute_PIMD_u();
    void compute_PIMD_x();
    void output_PIMD_u();
    void output_PIMD_x();
    void zero_PIMD_u();
    void zero_PIMD_x();
    void zero_PIMD_fu();
    void energy_PIMD_u();
    void energy_PIMD_x();
    void modelpot_PIMD_x(double *);
    void checkUforce();

    double *x,*y,*z;
    double *xu,*yu,*zu;
    double *fx,*fy,*fz;
    double *fxu,*fyu,*fzu;
    double *rat1,*rat2,*veig;
};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
int main()
  //============================================================================
{//begin routine
  //============================================================================

  int numBeads = 10;
  PIBeadAtoms  pibeads(numBeads);

  //============================================================================
  // Check x to u and u to x

  pibeads.output_PIMD_x();
  pibeads.energy_PIMD_x();

  //  pibeads.zero_PIMD_u();   // not necessary but useful for debugging
  pibeads.compute_PIMD_u();
  pibeads.output_PIMD_u();
  pibeads.energy_PIMD_u();

  //  pibeads.zero_PIMD_x();  // not necessary but useful for debugging
  pibeads.compute_PIMD_x();
  pibeads.output_PIMD_x();
  pibeads.energy_PIMD_x();

  //============================================================================
  // Check the u forces

  //  pibeads.zero_PIMD_fu();  // not necessary but useful for debugging
  pibeads.compute_PIMD_Fu();
  pibeads.checkUforce();

  //============================================================================

  return 1;
  exit(0);

  //============================================================================
}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
PIBeadAtoms::PIBeadAtoms(int numBeads_dum)
  //============================================================================
{// begin routine
  //============================================================================
  //  printf("{%d}[%d] PIBeadAtoms::PIBeadAtoms numbeads %d\n",thisInstance.proxyOffset,thisIndex, numBeads);

  numBeads = numBeads_dum;

  x= new double[numBeads];
  y= new double[numBeads];
  z= new double[numBeads];
  xu= new double[numBeads];
  yu= new double[numBeads];
  zu= new double[numBeads];
  fx= new double[numBeads];
  fy= new double[numBeads];
  fz= new double[numBeads];
  fxu= new double[numBeads];
  fyu= new double[numBeads];
  fzu= new double[numBeads];

  rat1 = new double[numBeads];
  rat2 = new double[numBeads];
  veig = new double[numBeads];

  //============================================================================
  // Initialize magic ratios

  for(int ip=2;ip<=numBeads;ip++){
    int ip1 = ip-1;
    rat1[ip1] = ((double)(ip1))/((double)(ip));
    rat2[ip1] = 1.0/((double)(ip));
    veig[ip1] = ((double)(ip))/((double)(ip1));
  }//endfor

  //============================================================================

  long int seedval = 91435715;
  srand48(seedval);
  for(int i =0;i<numBeads;i++){
    x[i] = drand48() - 0.5;
    y[i] = drand48() - 0.5;
    z[i] = drand48() - 0.5;
    xu[i] = 0.0;
    yu[i] = 0.0;
    zu[i] = 0.0;
    fx[i] = -x[i];
    fy[i] = -y[i];
    fz[i] = -z[i];
  }//endfor

  //============================================================================
}//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::compute_PIMD_Fu()
  //============================================================================
{// begin routine
  //============================================================================

  if(numBeads<=0){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("The number of beads must be >=0\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    exit(0);
  }//endif

  //============================================================================

  if(numBeads>1){
    fxu[0]  = 0.0;  fyu[0]  = 0.0;  fzu[0]  = 0.0;
    for(int ip=0;ip<numBeads;ip++){
      fxu[0] += fx[ip]; fyu[0] += fy[ip]; fzu[0] += fz[ip];
    }/*endfor*/
    fxu[1]  = fx[1];  fyu[1]  = fy[1];  fzu[1] = fz[1];
    for(int ip=1;ip<numBeads-1;ip++){
      int ip1 = ip+1;
      fxu[ip1] = fx[ip1] + rat1[ip]*fxu[ip];
      fyu[ip1] = fy[ip1] + rat1[ip]*fyu[ip];
      fzu[ip1] = fz[ip1] + rat1[ip]*fzu[ip];
    }/*endfor*/
  }else{
    fxu[0]  = fx[0];  fyu[0]  = fy[0];  fzu[0] = fz[0];
  }/*endif*/

  //============================================================================
}//end routine
//============================================================================




//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::compute_PIMD_u()
  //============================================================================
{// begin routine
  //============================================================================

  if(numBeads<=0){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("The number of beads must be >=0\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    exit(0);
  }//endif

  //============================================================================

  if(numBeads>1){
    xu[0]   = x[0];  yu[0]   = y[0];  zu[0]   = z[0];
    for(int ip=1;ip<numBeads-1;ip++){
      int ip1 = ip+1;
      double xstar = rat1[ip]*x[ip1] + rat2[ip]*x[0];
      double ystar = rat1[ip]*y[ip1] + rat2[ip]*y[0];
      double zstar = rat1[ip]*z[ip1] + rat2[ip]*z[0];
      xu[ip] =  x[ip] - xstar;
      yu[ip] =  y[ip] - ystar;
      zu[ip] =  z[ip] - zstar;
    }/*endfor:ip*/
    int ip = numBeads-1;
    xu[ip] =  x[ip] - x[0];
    yu[ip] =  y[ip] - y[0];
    zu[ip] =  z[ip] - z[0];
  }else{
    xu[0]  = x[0];  yu[0]  = y[0];  zu[0] = z[0];
  }//endif

  //============================================================================
}//end routine
//============================================================================




//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::compute_PIMD_x()
  //============================================================================
{// begin routine
  //============================================================================

  if(numBeads<=0){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("The number of beads must be >=0\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    exit(0);
  }//endif

  //============================================================================


  if(numBeads>1){
    x[0] = xu[0];  y[0] = yu[0];  z[0] = zu[0];
    int ip = numBeads-1;
    x[ip] = xu[ip] + x[0]; y[ip] = yu[ip] + y[0]; z[ip] = zu[ip] + z[0];
    for(int ip=numBeads-2;ip>=1;ip--){
      int ip1 = ip+1;
      double xadd = rat1[ip]*x[ip1]+xu[0]*rat2[ip];
      double yadd = rat1[ip]*y[ip1]+yu[0]*rat2[ip];
      double zadd = rat1[ip]*z[ip1]+zu[0]*rat2[ip];
      x[ip] = xu[ip] + xadd;
      y[ip] = yu[ip] + yadd;
      z[ip] = zu[ip] + zadd;
    }/*endfor*/
  }else{
    x[0]  = xu[0];  y[0]  = yu[0];  z[0] = zu[0];
  }//endif

  //============================================================================
}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::output_PIMD_x()
  //============================================================================
{// begin routine
  //============================================================================

  printf("\n=================================\n");
  for(int i =0;i<numBeads;i++){
    printf("x[%d]: %g %g %g\n",i,x[i],y[i],z[i]);
  }//endfor
  printf("=================================\n");

  //============================================================================
}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::output_PIMD_u()
  //============================================================================
{// begin routine
  //============================================================================

  printf("\n=================================\n");
  for(int i =0;i<numBeads;i++){
    printf("u[%d]: %g %g %g\n",i,xu[i],yu[i],zu[i]);
  }//endfor
  printf("=================================\n");

  //============================================================================
}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::zero_PIMD_u()
  //============================================================================
{// begin routine
  //============================================================================

  for(int i =0;i<numBeads;i++){
    xu[i] = 0.0; yu[i] = 0.0; zu[i] = 0.0;;
  }//endfor

  //============================================================================
}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::zero_PIMD_fu()
  //============================================================================
{// begin routine
  //============================================================================

  for(int i =0;i<numBeads;i++){
    fxu[i] = 0.0; fyu[i] = 0.0; fzu[i] = 0.0;;
  }//endfor

  //============================================================================
}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::zero_PIMD_x()
  //============================================================================
{// begin routine
  //============================================================================

  for(int i =0;i<numBeads;i++){
    x[i] = 0.0; y[i] = 0.0; z[i] = 0.0;;
  }//endfor

  //============================================================================
}//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::energy_PIMD_x()
  //============================================================================
{// begin routine
  //============================================================================
  double ep;

  int ip  = 0;
  int ip1 = numBeads-1;
  ep = ( (x[ip]-x[ip1])*(x[ip]-x[ip1])
      +(y[ip]-y[ip1])*(y[ip]-y[ip1])
      +(z[ip]-z[ip1])*(z[ip]-z[ip1]) 
      );
  for(int ip=1;ip<numBeads;ip++){
    int ip1 = ip-1;
    ep += ( (x[ip]-x[ip1])*(x[ip]-x[ip1])
        +(y[ip]-y[ip1])*(y[ip]-y[ip1])
        +(z[ip]-z[ip1])*(z[ip]-z[ip1]) 
        );
  }//endfor
  printf("\n=================================\n");
  printf("x energy %.10g\n",ep);
  printf("=================================\n");

  //============================================================================
}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::energy_PIMD_u()
  //============================================================================
{// begin routine
  //============================================================================

  double ep = 0.0;
  for(int ip =1;ip<numBeads;ip++){
    double fk = veig[ip];
    ep += fk*(xu[ip]*xu[ip]
        +yu[ip]*yu[ip]
        +zu[ip]*zu[ip]
        );
  }//endfor
  printf("\n=================================\n");
  printf("u energy %.10g\n",ep);
  printf("=================================\n");

  //============================================================================
}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::checkUforce()
  //============================================================================
{// begin routine
  //============================================================================

  printf("\n=================================\n");

  double potpx,potmx;
  double potpy,potmy;
  double potpz,potmz;
  double delta = 1.e-05;
  for(int ip =0;ip<numBeads;ip++){

    xu[ip] += delta;
    compute_PIMD_x();
    modelpot_PIMD_x(&potpx);
    xu[ip] -= 2.0*delta;
    compute_PIMD_x();
    modelpot_PIMD_x(&potmx);
    xu[ip] += delta;

    yu[ip] += delta;
    compute_PIMD_x();
    modelpot_PIMD_x(&potpy);
    yu[ip] -= 2.0*delta;
    compute_PIMD_x();
    modelpot_PIMD_x(&potmy);
    yu[ip] += delta;

    zu[ip] += delta;
    compute_PIMD_x();
    modelpot_PIMD_x(&potpz);
    zu[ip] -= 2.0*delta;
    compute_PIMD_x();
    modelpot_PIMD_x(&potmz);
    zu[ip] += delta;
    printf("fu[%d]: %.10g %.10g : %.10g %.10g : %.10g %.10g\n",
        ip,fxu[ip],0.5*(potmx-potpx)/delta,
        fyu[ip],0.5*(potmy-potpy)/delta,
        fzu[ip],0.5*(potmz-potpz)/delta);
  }//endfor

  printf("=================================\n");

  //============================================================================
}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::modelpot_PIMD_x(double *pot_out)
  //============================================================================
{// begin routine
  //============================================================================

  double pot = 0.0;
  double fk = 0.5;
  for(int ip =0;ip<numBeads;ip++){
    pot += fk*(x[ip]*x[ip]
        +y[ip]*y[ip]
        +z[ip]*z[ip]
        );
  }//endfor

  pot_out[0] = pot;

  //============================================================================
}//end routine
//============================================================================
