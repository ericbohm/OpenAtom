/*
** PIBeadAtoms.C
** 
** Made by (Eric Bohm)
** Login   <bohm@alacrity>
** 
** Started on  Tue Mar  2 10:15:33 2010 Eric Bohm
** Last update Sun May 12 01:17:25 2002 Speed Blue

*/


#include "groups.h" 
#include "PIBeadAtoms.h"

extern CkVec <CProxy_AtomsGrp>                   UatomsGrpProxy;

// NOTE: thisIndex == atomIndex

PIBeadAtoms::PIBeadAtoms(UberCollection _thisInstance, int _numBeads) :thisInstance(_thisInstance), numBeads(_numBeads)
{
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
  acceptCount_Fx=0;
  acceptCount_u=0;
  acceptCount_x=0;
  
}

void PIBeadAtoms::accept_PIMD_Fx(AtomXYZMsg *msg){
  fx[msg->index]=msg->x[thisIndex];
  fy[msg->index]=msg->y[thisIndex];
  fz[msg->index]=msg->z[thisIndex];
  delete msg;
  //  atomdest[PIBeadIndex]=atomdest;
  acceptCount_Fx++;
  //  CkPrintf("{%d}[%d] PIBeadAtoms::accept_PIMD_Fx (%d of %d)\n",thisInstance.proxyOffset,thisIndex, acceptCount_Fx, numBeads);

  // 
  if(acceptCount_Fx==numBeads)
    {
      compute_PIMD_Fu();
      acceptCount_Fx=0;
      UberCollection instance=thisInstance;
      for(int bead=0;bead<numBeads; bead++)
      {
	  instance.idxU.x=bead;
	  int proxyOffset=instance.setPO();
	  //	  CkPrintf("{%d}[%d] PIBeadAtoms::accept_PIMD_Fx sending atomsGrp{%d}.accept_PIMD_fu\n",thisInstance.proxyOffset,thisIndex, proxyOffset);
	  //	  UatomsGrpProxy[proxyOffset][atomdest].accept_PIMD_fu(fxu[bead],
	  //	  fyu[bead], fzu[bead], thisIndex);
	  // this atom index has to send the Fu to everyone
	  UatomsGrpProxy[proxyOffset].accept_PIMD_Fu(fxu[bead], fyu[bead], fzu[bead], thisIndex);
      }
    }
}

void PIBeadAtoms::accept_PIMD_u(double _xu, double _yu, double _zu, int PIBeadIndex)
{


  xu[PIBeadIndex]=_xu;
  yu[PIBeadIndex]=_yu;
  zu[PIBeadIndex]=_zu;
  acceptCount_u++;
  //  CkPrintf("{%d}[%d] PIBeadAtoms::accept_PIMD_u (%d of %d)\n",thisInstance.proxyOffset,thisIndex, acceptCount_u, numBeads );
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
	    //	    UatomsGrpProxy[proxyOffset][atomdest].accept_PIMD_x(x[bead],y[bead], z[bead], thisIndex);
	    UatomsGrpProxy[proxyOffset].accept_PIMD_x(x[bead], y[bead], z[bead], thisIndex);
	}
    }
}
void PIBeadAtoms::accept_PIMD_x(double _x, double _y, double _z, int PIBeadIndex )
{

  x[PIBeadIndex]=_x;
  y[PIBeadIndex]=_y;
  z[PIBeadIndex]=_z;
  acceptCount_x++;
  //  CkPrintf("{%d}[%d] PIBeadAtoms::accept_PIMD_x (%d of %d) \n",thisInstance.proxyOffset,thisIndex, acceptCount_x, numBeads);
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
	    //	    UatomsGrpProxy[proxyOffset][atomdest].accept_PIMD_u(xu[bead],yu[bead],zu[bead],  thisIndex);
	    UatomsGrpProxy[proxyOffset].accept_PIMD_u(xu[bead],yu[bead],zu[bead], thisIndex);
	}

    }
}
#include "PIBeadAtoms.def.h"