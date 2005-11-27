//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file CP_Rho_GHartExc.C
 *
 *  This is a description of the "life" of a CP_Rho_GHartExc  object
 * 
 *  At the start of the program, the constructor CP_Rho_GHartExc is
 *  called.  We create our own rho_gs slab because we need the
 *  kvectors and the fft code.
 *
 *  Each iteration, the CP_Rho_GpacePlane object sends us the same rho
 *  data it will use in its own FFTs.  We then call the very intensive
 *  HartExtVksG function on the data.  We contribute the energy to the
 *  reduction, fft the vks, and ship the vks to rhoReal.
 *
 *  Then we're done until the next iteration.
 * 
 *  There is no RthThread control loop here because there is no
 *  meaningful flow of entry methods.  We get a message, we calculate,
 *  we send, we're done.  This object exists solely as a way to
 *  parallelize the uncomfortably long HartExtVks computation.
 */ 
//============================================================================

#include "charm++.h"
#include <iostream.h>
#include <fstream.h>
#include <math.h>

#include "../../include/debug_flags.h"
#include "util.h"
#include "groups.h"
#include "cpaimd.h"
#include "sim_subroutines.h"
#include "CP_State_Plane.h"

#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cplocal.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h"

//============================================================================
extern int sizeX;
extern Config config;
extern CProxy_CP_Rho_RealSpacePlane rhoRealProxy;
extern CProxy_CPcharmParaInfoGrp scProxy;
extern CProxy_AtomsGrp atomsGrpProxy;
extern ComlibInstanceHandle commGHartInstance;

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
/**
 *  This object just gets a rho message, computes GHartExt, and sends
 *  vks.  
 */
//============================================================================
CP_Rho_GHartExt::CP_Rho_GHartExt(size2d sizeYZ)
//============================================================================
{//begin routine
//============================================================================

  iopt=0;
  CkAssert(sizeX>0); //check for startup wackiness
  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;      
  CkVec <RunDescriptor> *sortedRunDescriptors = sim->RhosortedRunDescriptors;
  rho_gs.sizeX    = sizeX;
  rho_gs.sizeY    = sizeYZ[0];
  rho_gs.sizeZ    = sizeYZ[1];
  rho_gs.xdim     = rho_gs.sizeX;
  rho_gs.ydim     = rho_gs.sizeY;
  rho_gs.zdim     = 1;
  int x = thisIndex.x;
  rho_gs.numRuns  = sortedRunDescriptors[x].size();
  CkAssert(rho_gs.numRuns>0);

  rho_gs.numLines = sortedRunDescriptors[x].size()/2;
  rho_gs.numFull  = (rho_gs.numLines)*rho_gs.sizeZ;
  rho_gs.size     = rho_gs.numFull;
  rho_gs.runs     = new RunDescriptor[rho_gs.numRuns];
  rho_gs.numPoints = 0;
  for (int r = 0; r < rho_gs.numRuns; r++) {
    rho_gs.numPoints += sortedRunDescriptors[x][r].length;
    rho_gs.runs[r]    = sortedRunDescriptors[x][r];
  }//endfor
  int nPacked;
  rho_gs.setKVectors(&nPacked);
  CkAssert(nPacked==rho_gs.numPoints);
  rho_gs.packedVks = (complex *)fftw_malloc(nPacked*sizeof(complex));
  int numFull=rho_gs.numFull;
  rho_gs.Vks       = (complex *)fftw_malloc(numFull*sizeof(complex));
  rho_gs.packedRho = (complex *)fftw_malloc(nPacked*sizeof(complex));
  // we'd rather use a pointer to the inbound message
  //  rho_gs.packedRho = NULL;
  // we don't use these gs_slab objects
  rho_gs.divRhoX   = NULL;
  rho_gs.divRhoY   = NULL;
  rho_gs.divRhoZ   = NULL;


  setMigratable(false);

  rhoRealProxy_com = rhoRealProxy;
  if(config.useCommlib){
      ComlibAssociateProxy(&commGHartInstance, rhoRealProxy_com);
  }
}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_Rho_GHartExt::~CP_Rho_GHartExt()
{
}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GHartExt::acceptData(RhoGHartMsg *msg){
//============================================================================
// compute hart+Ext and vks(g) using rho(g) 
  
  rho_gs.packedRho = msg->data;  // why bother copying this readonly data?
  HartExtVksG();
  rho_gs.packedRho=NULL;  
  delete msg;  

}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GHartExt::HartExtVksG() { 
//============================================================================
// get the variables
   AtomsGrp *ag = atomsGrpProxy.ckLocalBranch(); // find me the local copy
   int natm          = ag->natm;
   Atom *atoms       = ag->atoms;
   double ehart_ret=0.0;
   double eext_ret=0.0;
   double ewd_ret=0.0;
   int numPoints= rho_gs.numPoints;
   int numLines= rho_gs.numLines;
   int numFull= rho_gs.numFull;
   complex *rho=rho_gs.packedRho;
   complex *vks=rho_gs.packedVks;
   complex *Vks=rho_gs.Vks;
   bzero(vks, numPoints*sizeof(complex));
   int *k_x=rho_gs.k_x;
   int *k_y=rho_gs.k_y;
   int *k_z=rho_gs.k_z;

//============================================================================
// compute vks(g) from hart eext and reduce eext and ehart

#ifndef CMK_OPTIMIZE
   double  StartTime=CmiWallTimer();
#endif    

   CPLOCAL::CP_hart_eext_calc(numPoints, rho, natm, atoms, vks,
                              &ehart_ret,&eext_ret,&ewd_ret,k_x,k_y,k_z,
                              thisIndex.x);
#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(HartExcVksG_, StartTime, CmiWallTimer());    
#endif
   double e[3];
   e[0] = ehart_ret;
   e[1] = eext_ret;
   e[2] = ewd_ret;
   contribute(3 * sizeof(double),e,CkReduction::sum_double);
#ifdef _CP_DEBUG_RHOG_VKSA_
     char myFileName[100];
     sprintf(myFileName, "Vks_Gspace_%d%d.out", thisIndex.x,thisIndex.y);
     FILE *fp = fopen(myFileName,"w");
       for (int i = 0; i < numPoints; i++){ 
              fprintf(fp," %d %d %d : %g %g\n",
                 k_x[i],k_y[i],k_z[i],
                 vks[i].re,vks[i].im);
       }//endfor
     fclose(fp);
#endif

//============================================================================
// partly fft vks  fft_gz(vks)

     bzero(Vks, numFull*sizeof(complex));
     rho_gs.expandRhoGSpace(Vks, vks);
     int ioptvks=4;
     rho_gs.doFwFFTGtoR(ioptvks, thisIndex.x); 

#ifdef _CP_DEBUG_RHOG_VKSA_
     sprintf(myFileName, "Vks_GspaceAFFT_%d%d.out", thisIndex.x,thisIndex.y);
     fp = fopen(myFileName,"w");
       for (int i = 0; i < numFull; i++){ 
              fprintf(fp," %g %g\n",
                 Vks[i].re, Vks[i].im);
       }//endfor
     fclose(fp);
#endif

//============================================================================
// send vks_hart_ext back to rho_real where fft(gx,gy) will be performed

     sendVks();

//---------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GHartExt::sendVks() { 
//============================================================================

#ifdef _CP_DEBUG_RHOG_VERBOSE_
  CkPrintf("Communicating data from RhoGHart to RhoR : %d %d\n",
	   thisIndex.x,thisIndex.y);
#endif

//============================================================================
// Do a Comlib Dance

  if (config.useCommlib){
      commGHartInstance.beginIteration();
  }
  int numLines=rho_gs.numLines;
  
//============================================================================

  int sizeZ=rho_gs.sizeZ;
  for(int z=0; z < sizeZ; z++) {

    RhoRSFFTMsg *msg = new (numLines,8*sizeof(int)) RhoRSFFTMsg;
    msg->size        = numLines;     // number of z-lines in this batch
    msg->senderIndex = thisIndex.x;  // line batch index
    msg->iopt        = iopt;         // iopt always 0 for us
    // FIGURE OUT WHAT iopt should be
    if(config.prioFFTMsg){
       CkSetQueueing(msg, CK_QUEUEING_IFIFO);
       *(int*)CkPriorityPtr(msg) = config.rhorpriority + thisIndex.x + thisIndex.y;
    }//endif

    // beam out all points with same z to chare array index z
    complex *data = msg->data;
    complex *Vks = rho_gs.Vks;
    for (int i=0,j=z; i<numLines; i++,j+=sizeZ){data[i] = Vks[j];}
    rhoRealProxy_com(z,0).acceptHartVks(msg);

  }//endfor

//============================================================================
// Complete the commlib dance
    
  if (config.useCommlib){
      commGHartInstance.endIteration();
  }

//---------------------------------------------------------------------------
  }// end routine
//============================================================================


