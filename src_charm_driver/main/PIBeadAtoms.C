/*
** PIBeadAtoms.C
** 
** Made by (Eric Bohm)
** Login   <bohm@alacrity>
** 
** Started on  Tue Mar  2 10:15:33 2010 Eric Bohm
** Last update Sun May 12 01:17:25 2002 Speed Blue

*/

#include "PIBeadAtoms.h"
#include "groups.h" 

extern CkVec <CProxy_AtomsGrp>                   UatomsGrpProxy;

// NOTE: thisIndex == atomIndex

PIBeadAtoms::PIBeadAtoms(UberCollection _thisInstance, int _numBeads) :thisInstance(_thisInstance), numBeads(_numBeads)
{
  //  CkPrintf("{%d}[%d] PIBeadAtoms::PIBeadAtoms \n",thisInstance.proxyOffset,thisIndex);
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
  acceptCount_Fx=0;
  acceptCount_u=0;
  acceptCount_x=0;
}

void PIBeadAtoms::accept_PIMD_Fx(double _Fx, double _Fy, double _Fz, int PIBeadIndex){
  fx[PIBeadIndex]=_Fx;
  fy[PIBeadIndex]=_Fy;
  fz[PIBeadIndex]=_Fz;
  acceptCount_Fx++;
  //  CkPrintf("{%d}[%d] PIBeadAtoms::accept_PIMD_Fx \n",thisInstance.proxyOffset,thisIndex);
  if(acceptCount_Fx==numBeads)
    {
      compute_PIMD_Fu();
      acceptCount_Fx=0;
      UberCollection instance=thisInstance;
      for(int bead=0;bead<numBeads; bead++)
      {
	  instance.idxU.x=bead;
	  int proxyOffset=instance.setPO();
	  int atomdest=0;
	  UatomsGrpProxy[proxyOffset][atomdest].accept_PIMD_fu(fxu[bead], fyu[bead], fzu[bead], thisIndex);
      }
    }
}

void PIBeadAtoms::accept_PIMD_u(double _xu, double _yu, double _zu, int PIBeadIndex)
{
  //  CkPrintf("{%d}[%d] PIBeadAtoms::accept_PIMD_u \n",thisInstance.proxyOffset,thisIndex);

  xu[PIBeadIndex]=_xu;
  yu[PIBeadIndex]=_yu;
  zu[PIBeadIndex]=_zu;
  acceptCount_u++;

  if(acceptCount_u==numBeads)
    {
      compute_PIMD_x();
      acceptCount_u=0;
      UberCollection instance=thisInstance;
      for(int bead=0;bead<numBeads; bead++)
	{
	    instance.idxU.x=bead;
	    int proxyOffset=instance.setPO();
	    int atomdest=0;
	    UatomsGrpProxy[proxyOffset][atomdest].accept_PIMD_x(x[bead], y[bead], z[bead], thisIndex);
	}
    }
}
void PIBeadAtoms::accept_PIMD_x(double _x, double _y, double _z, int PIBeadIndex )
{
  //  CkPrintf("{%d}[%d] PIBeadAtoms::accept_PIMD_x \n",thisInstance.proxyOffset,thisIndex);
  x[PIBeadIndex]=_x;
  y[PIBeadIndex]=_y;
  z[PIBeadIndex]=_z;
  acceptCount_x++;
  if(acceptCount_x==numBeads)
    {
      compute_PIMD_u();
      acceptCount_x=0;
      UberCollection instance=thisInstance;
      for(int bead=0;bead<numBeads; bead++)
	{
	    instance.idxU.x=bead;
	    int proxyOffset=instance.setPO();
	    int atomdest=0;
	    UatomsGrpProxy[proxyOffset][atomdest].accept_PIMD_u(xu[bead],yu[bead],zu[bead], thisIndex);
	}

    }
}
#include "PIBeadAtoms.def.h"
