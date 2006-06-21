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
#include "cpaimd.h"
#include "groups.h"
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
extern CProxy_CP_Rho_GHartExt rhoGHartExtProxy;
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

  CkAssert(sizeX>0); //check for startup wackiness
  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;      
  CkVec <RunDescriptor> *sortedRunDescriptors = sim->RhosortedRunDescriptors;

  iopt            = 0;
  iteration       = 0;
  rhoGHelpers     = config.rhoGHelpers;
  rho_gs.sizeX    = sizeX;
  rho_gs.sizeY    = sizeYZ[0];
  rho_gs.sizeZ    = sizeYZ[1];
  rho_gs.xdim     = rho_gs.sizeX;
  rho_gs.ydim     = rho_gs.sizeY;
  rho_gs.zdim     = 1;

//==================================================================================
// Decomposition rhoG lines into slices of size rhoGHelper

  ind_x            = thisIndex.x;
  ind_xdiv         = (ind_x/rhoGHelpers);
  ind_xrem         = (ind_x%rhoGHelpers);
  int numLines_tot = sortedRunDescriptors[ind_xdiv].size()/2;

  getSplitDecomp(&istrt_lines,&iend_lines,&numLines,
                 numLines_tot,rhoGHelpers,ind_xrem);

//==================================================================================
// Carve out your rundescriptor, make the k-vectors , malloc the memory

  rho_gs.numLines  = numLines;
  rho_gs.numRuns   = (numLines*2);
  rho_gs.numFull   = (numLines*rho_gs.sizeZ);
  rho_gs.size      = rho_gs.numFull;
  rho_gs.runs      = new RunDescriptor[(rho_gs.numRuns)];
  rho_gs.numPoints = 0;
  for (int r = (2*istrt_lines),s=0; r < (2*iend_lines); r++,s++) {
    rho_gs.numPoints += sortedRunDescriptors[ind_xdiv][r].length;
    rho_gs.runs[s]    = sortedRunDescriptors[ind_xdiv][r];
  }//endfor

  int nPacked;
  rho_gs.setKVectors(&nPacked);
  rho_gs.nPacked=nPacked;
  int numFull=rho_gs.numFull;

  CkAssert(nPacked==rho_gs.numPoints);


  rho_gs.packedVks = (complex *)fftw_malloc(nPacked*sizeof(complex));
  rho_gs.Vks       = (complex *)fftw_malloc(numFull*sizeof(complex));
  rho_gs.packedRho = (complex *)fftw_malloc(nPacked*sizeof(complex));
  rho_gs.divRhoX   = NULL;
  rho_gs.divRhoY   = NULL;
  rho_gs.divRhoZ   = NULL;
  rho_gs.Rho       = NULL;

//==================================================================================
// Set some proxies, set the migratable flag

  setMigratable(false);

  rhoRealProxy_com = rhoRealProxy;
  if(config.useGHartInsRhoRP){
     ComlibAssociateProxy(&commGHartInstance,rhoRealProxy_com);
  }//endif
  usesAtSync = CmiTrue;
  if(config.lbdensity){
    setMigratable(true);
  }else{
    setMigratable(false);
  }//endif


//---------------------------------------------------------------------------
  }//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_Rho_GHartExt::~CP_Rho_GHartExt(){
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GHartExt::pup(PUP::er &p){
  ArrayElement2D::pup(p);
  rho_gs.pup(p);
  p|iopt;
  p|iteration;
  p|ind_x;
  p|ind_xdiv;
  p|ind_xrem;
  p|rhoGHelpers; 
  p|istrt_lines; 
  p|iend_lines;
  p|numLines;
  p|rhoRealProxy_com;
  //  if(config.useCommlib)
  //    ComlibResetProxy(&rhoRealProxy_com);

//---------------------------------------------------------------------------
}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GHartExt::acceptData(RhoGHartMsg *msg){
//============================================================================
// Copy out the data and wait until we can use it
  
  int ncoef  = rho_gs.numPoints;
  CkAssert(ncoef==msg->size);
  memcpy(rho_gs.packedRho,msg->data,sizeof(complex)*ncoef);
  delete msg;  

// Check the flow of control to see if we can use the data.
   if(atomsGrpProxy.ckLocalBranch()->iteration == iteration){
      HartExtVksG();
   }else{
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Flow of Control Error in HartExtVks : atoms slow\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
   }//endif
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GHartExt::HartExtVksG() { 
//============================================================================
// Get the variables

   AtomsGrp *ag     = atomsGrpProxy.ckLocalBranch(); // find me the local copy
   int natm         = ag->natm;
   Atom *atoms      = ag->atoms;
   double ehart_ret = 0.0;
   double eext_ret  = 0.0;
   double ewd_ret   = 0.0;
   int numPoints    = rho_gs.numPoints;
   int numLines     = rho_gs.numLines;
   int numFull      = rho_gs.numFull;
   complex *rho     = rho_gs.packedRho;
   complex *vks     = rho_gs.packedVks;
   complex *VksExpd = rho_gs.Vks;
   int *k_x         = rho_gs.k_x;
   int *k_y         = rho_gs.k_y;
   int *k_z         = rho_gs.k_z;

//============================================================================
// compute vks(g) from hart eext and reduce eext and ehart

#ifndef CMK_OPTIMIZE
   double  StartTime=CmiWallTimer();
#endif    

   bzero(vks,numPoints*sizeof(complex));
   CPLOCAL::CP_hart_eext_calc(numPoints, rho, natm, atoms, vks,
                              &ehart_ret,&eext_ret,&ewd_ret,k_x,k_y,k_z,
                              thisIndex.x);
   iteration ++;
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

     bzero(VksExpd, numFull*sizeof(complex));
     rho_gs.expandRhoGSpace(VksExpd, vks);
     int ioptvks=4;
     rho_gs.doFwFFTGtoR(ioptvks, thisIndex.x); 

#ifdef _CP_DEBUG_RHOG_VKSA_
     sprintf(myFileName, "Vks_GspaceAFFT_%d%d.out", thisIndex.x,thisIndex.y);
     fp = fopen(myFileName,"w");
       for (int i = 0; i < numFull; i++){ 
              fprintf(fp," %g %g\n",
                 VksExpd[i].re, VksExpd[i].im);
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

  if (config.useGHartInsRhoRP){
      commGHartInstance.beginIteration();
  }
  
//============================================================================

  int sizeZ=rho_gs.sizeZ;
  for(int z=0; z < sizeZ; z++) {

    RhoHartRSFFTMsg *msg = new (numLines,8*sizeof(int)) RhoHartRSFFTMsg;
    msg->size           = numLines;   // number of z-lines in this batch
    msg->senderBigIndex = ind_xdiv;   // big line batch index
    msg->senderStrtLine = istrt_lines;// where my lines start in big batch
    msg->iopt           = iopt;       // iopt always 0 for us
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
    
  if (config.useGHartInsRhoRP){
      commGHartInstance.endIteration();
  }

//---------------------------------------------------------------------------
  }// end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void getSplitDecomp(int *istrt_ret,int *iend_ret,int *n_ret,
                    int ntot, int ndiv,int idiv) 
//============================================================================
  {//begin routine
//============================================================================

   if(idiv>=ndiv || ntot< ndiv){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Incorrect input to RhoGHart collection creator.\n");
     CkPrintf("idiv %d ndiv %d, ntot %d.\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
   }//endif

   int n     = (ntot/ndiv);
   int r     = (ntot%ndiv);

   int istrt = n*idiv;
   if(idiv>=r){istrt += r;}
   if(idiv<r) {istrt += idiv;}
   if(idiv<r) {n++;}
   int iend  = n+istrt;

   if(n==0){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("No lines in a RhoGHart collection!!\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
   }//endif
  
   (*n_ret)     = n;
   (*istrt_ret) = istrt;
   (*iend_ret)  = iend;

//---------------------------------------------------------------------------
  }//end routine
//============================================================================
