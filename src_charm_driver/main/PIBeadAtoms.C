/*
** PIBeadAtoms.C
** 
** Made by (Eric Bohm)
** Login   <bohm@alacrity>
** 
** Started on  Tue Mar  2 10:15:33 2010 Eric Bohm
** Last update Sun May 12 01:17:25 2002 Speed Blue

*/


#include "AtomsCompute.h" 
#include "PIBeadAtoms.h"

extern CkVec <CProxy_AtomsCompute>                   UatomsComputeProxy;

// NOTE: thisIndex == atomIndex

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
PIBeadAtoms::PIBeadAtoms(UberCollection _thisInstance, int _numBeads, int _natm) :thisInstance(_thisInstance), 
										  numBeads(_numBeads), natm (_natm)
//============================================================================
    {// begin routine
//============================================================================
//  CkPrintf("{%d}[%d] PIBeadAtoms::PIBeadAtoms numbeads %d\n",thisInstance.proxyOffset,thisIndex, numBeads);

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

  rat1[0] = 0.0;
  rat2[0] = 0.0;
  veig[0] = 0.0;
  for(int ip=2;ip<=numBeads;ip++){
   int ip1 = ip-1;
   rat1[ip1] = ((double)(ip1))/((double)(ip));
   rat2[ip1] = 1.0/((double)(ip));
   veig[ip1] = ((double)(ip))/((double)(ip1));
  }//endfor

//============================================================================
// Initialize communication counters

  acceptCount_Fx=0;
  acceptCount_u=0;
  acceptCount_x=0;
  
//============================================================================
  }//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::accept_PIMD_Fx(AtomXYZMsg *msg)
//============================================================================
    {// begin routine
  fx[msg->index]=msg->x[thisIndex];
  fy[msg->index]=msg->y[thisIndex];
  fz[msg->index]=msg->z[thisIndex];
  delete msg;
  //  atomdest[PIBeadIndex]=atomdest;
  acceptCount_Fx++;
  //  CkPrintf("{%d}[%d] PIBeadAtoms::accept_PIMD_Fx (%d of %d)\n",thisInstance.proxyOffset,thisIndex, acceptCount_Fx, numBeads);

  // 
  if(acceptCount_Fx==numBeads){
      compute_PIMD_Fu();
      acceptCount_Fx=0;
      UberCollection instance=thisInstance;
      for(int bead=0;bead<numBeads; bead++){
	  instance.idxU.x=bead;
	  int proxyOffset=instance.setPO();
	  //	  CkPrintf("{%d}[%d] PIBeadAtoms::accept_PIMD_Fx sending atomsGrp{%d}.accept_PIMD_fu\n",thisInstance.proxyOffset,thisIndex, proxyOffset);
	  //	  UatomsGrpProxy[proxyOffset][atomdest].accept_PIMD_fu(fxu[bead],
	  //	  fyu[bead], fzu[bead], thisIndex);
	  // this atom index has to send the Fu to everyone
	  UatomsComputeProxy[proxyOffset].accept_PIMD_Fu(fxu[bead], fyu[bead], fzu[bead], thisIndex);
      }//endfor
    }//endif
//============================================================================
    }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::accept_PIMD_Fx_and_x(AtomXYZMsg *msg)
//============================================================================
  {// begin routine
//============================================================================

  fx[msg->index]=msg->x[thisIndex];
  fy[msg->index]=msg->y[thisIndex];
  fz[msg->index]=msg->z[thisIndex];
  x[msg->index]=msg->x[(thisIndex+natm)];
  y[msg->index]=msg->y[(thisIndex+natm)];
  z[msg->index]=msg->z[(thisIndex+natm)];
  delete msg;

//============================================================================
//  atomdest[PIBeadIndex]=atomdest;

  acceptCount_Fx++;
  //  CkPrintf("{%d}[%d] PIBeadAtoms::accept_PIMD_Fx (%d of %d)\n",thisInstance.proxyOffset,thisIndex, acceptCount_Fx, numBeads);

  if(acceptCount_Fx==numBeads){
      compute_PIMD_Fu();
      compute_PIMD_u();
      acceptCount_Fx=0;
      UberCollection instance=thisInstance;
      for(int bead=0;bead<numBeads; bead++){
	  instance.idxU.x=bead;
	  int proxyOffset=instance.setPO();
	  //	  CkPrintf("{%d}[%d] PIBeadAtoms::accept_PIMD_Fx sending atomsGrp{%d}.accept_PIMD_fu\n",thisInstance.proxyOffset,thisIndex, proxyOffset);
	  //	  UatomsGrpProxy[proxyOffset][atomdest].accept_PIMD_fu(fxu[bead],
	  //	  fyu[bead], fzu[bead], thisIndex);
	  // this atom index has to send the Fu to everyone
	  UatomsComputeProxy[proxyOffset].accept_PIMD_Fu_and_u(fxu[bead], fyu[bead], fzu[bead],xu[bead], yu[bead], zu[bead], thisIndex);
      }//endfor
    }//endif
//============================================================================
    }//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::accept_PIMD_u(double _xu, double _yu, double _zu, int PIBeadIndex)
//============================================================================
    {// begin routine
//============================================================================

  xu[PIBeadIndex]=_xu;
  yu[PIBeadIndex]=_yu;
  zu[PIBeadIndex]=_zu;
  acceptCount_u++;
  //  CkPrintf("{%d}[%d] PIBeadAtoms::accept_PIMD_u (%d of %d)\n",thisInstance.proxyOffset,thisIndex, acceptCount_u, numBeads );
  if(acceptCount_u==numBeads){
      compute_PIMD_x();
      acceptCount_u=0;
      UberCollection instance=thisInstance;
      for(int bead=0;bead<numBeads; bead++){
	    instance.idxU.x=bead;
	    int proxyOffset=instance.setPO();
	    int atomdest=0;
	    //	    UatomsGrpProxy[proxyOffset][atomdest].accept_PIMD_x(x[bead],y[bead], z[bead], thisIndex);
	    UatomsComputeProxy[proxyOffset].accept_PIMD_x(x[bead], y[bead], z[bead], thisIndex);
	}//endfor
    }//endif

//============================================================================
  }//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::accept_PIMD_x(double _x, double _y, double _z, int PIBeadIndex )
//============================================================================
    {// begin routine
//============================================================================

  x[PIBeadIndex]=_x;
  y[PIBeadIndex]=_y;
  z[PIBeadIndex]=_z;
  acceptCount_x++;
  //  CkPrintf("{%d}[%d] PIBeadAtoms::accept_PIMD_x (%d of %d) \n",thisInstance.proxyOffset,thisIndex, acceptCount_x, numBeads);
  if(acceptCount_x==numBeads){
      compute_PIMD_u();
      acceptCount_x=0;
      UberCollection instance=thisInstance;
      for(int bead=0;bead<numBeads; bead++){
	    instance.idxU.x=bead;
	    int proxyOffset=instance.setPO();
	    int atomdest=0;
	    //	    UatomsGrpProxy[proxyOffset][atomdest].accept_PIMD_u(xu[bead],yu[bead],zu[bead],  thisIndex);
	    UatomsComputeProxy[proxyOffset].accept_PIMD_u(xu[bead],yu[bead],zu[bead], thisIndex);
	}//endfor
    }//endif

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
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("The number of beads must be >=0\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
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
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("The number of beads must be >=0\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
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
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("The number of beads must be >=0\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
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
// useful for debugging
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::output_PIMD_x()
//============================================================================
  {// begin routine
//============================================================================

  CkPrintf("\n=================================\n");
  for(int i =0;i<numBeads;i++){
    CkPrintf("x[%d]: %g %g %g\n",i,x[i],y[i],z[i]);
  }//endfor
  CkPrintf("=================================\n\n");

//============================================================================
  }//end routine
//============================================================================


//============================================================================
// useful for debugging
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::output_PIMD_u()
//============================================================================
  {// begin routine
//============================================================================

  CkPrintf("\n=================================\n");
  for(int i =0;i<numBeads;i++){
    CkPrintf("u[%d]: %g %g %g\n",i,xu[i],yu[i],zu[i]);
  }//endfor
  CkPrintf("=================================\n\n");

//============================================================================
  }//end routine
//============================================================================


//============================================================================
// useful for debugging
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
// useful for debugging
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
// useful for debugging
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
// useful for debugging
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void PIBeadAtoms::energy_PIMD_x()
//============================================================================
  {// begin routine
//============================================================================

  int ip  = 0;
  int ip1 = numBeads-1;
  double ep;
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
  CkPrintf("\n=================================\n");
  CkPrintf("x energy %.10g\n",ep);
  CkPrintf("=================================\n");

//============================================================================
  }//end routine
//============================================================================


//============================================================================
// useful for debugging
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
  CkPrintf("\n=================================\n");
  CkPrintf("u energy %.10g\n",ep);
  CkPrintf("=================================\n");

//============================================================================
  }//end routine
//============================================================================



//============================================================================
// useful for debugging
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
// useful for debugging
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

#include "PIBeadAtoms.def.h"
