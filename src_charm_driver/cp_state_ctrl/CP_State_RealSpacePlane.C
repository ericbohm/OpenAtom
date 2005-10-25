//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file CP_State_RealSpacePlane.C
 * This is a description of the "life" of a CP_State_RealSpacePlane object.
 *
 * At the start of the program, the constructor CP_State_RealSpacePlane() is called.
 * The CP_State_GSpacePlane objects send data to CP_State_RealSpacePlane 
 * after doing FFT in the y direction. The data is sent through the doFFT() method. 
 * In this method FFTs in the z and x directions are performed. After this the 
 * squared magnitudes of the psi_r values are sent to RealSpaceDensity. The 
 * CP_State_RealSpacePlane object is idle until further messages are sent to it.
 *
 * The idle period of the CP_State_RealSpacePlane is ended by a message from 
 * RealSpaceDensity objects - doProduct(). In this method the Slab of data
 * in CP_State_RealSpacePlane, psi_r, is multiplied with the data sent from 
 * RealSpaceDensity. Then, inverse ffts in z and x directions are performed.
 * Then the slab of data is send to CP_State_GSpacePlanes so that the inverse fft
 * in the x direction can be performed.
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

//============================================================================

extern CProxy_CP_State_GSpacePlane gSpacePlaneProxy;
extern CProxy_CP_Rho_RealSpacePlane rhoRealProxy;
extern CProxy_CP_State_RealSpacePlane realSpacePlaneProxy;
extern CProxy_CPcharmParaInfoGrp scProxy;
extern CkReduction::reducerType complexVectorAdderType;
extern CkGroupID mCastGrpId;
extern Config config;
extern int sizeX;
extern ComlibInstanceHandle ssInstance;
extern ComlibInstanceHandle mssInstance;
extern CProxy_CP_State_ParticlePlane particlePlaneProxy;
extern CProxy_FFTcache fftCacheProxy;

//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
RTH_Routine_locals(CP_State_RealSpacePlane,run)
RTH_Routine_code(CP_State_RealSpacePlane,run) {
//============================================================================

  while(1) {
  
    RTH_Suspend(); 
    c->doFFT();
    RTH_Suspend(); 
    c->doProduct();

  } //end while not done

//--------------------------------------------------------------------------
   } RTH_Routine_end(CP_State_RealSpacePlane,run)
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealSpacePlane::run () {
  run_thread = RTH_Runtime_create(RTH_Routine_lookup(CP_State_RealSpacePlane,run),this);
  RTH_Runtime_resume(run_thread);
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_State_RealSpacePlane::CP_State_RealSpacePlane(size2d size, int gSpaceUnits, 
                                                 int realSpaceUnits) {
//============================================================================
    count = 0;
    initRealStateSlab(&rs, size, gSpaceUnits, realSpaceUnits, thisIndex.x, thisIndex.y);
    flagsRecd = false;
    sendFFTDataSize = 0;
    setMigratable(false);
    gproxy = gSpacePlaneProxy;
    if (config.useCommlib){
      ComlibAssociateProxy(&mssInstance,gproxy);
    }//endif

    run();
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealSpacePlane::setNumPlanesToExpect(int num)
{
    rs.numPlanesToExpect = num;
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void
reduce(void *data)
{
    CkIndex2D idx = *(CkIndex2D *) data;
    realSpacePlaneProxy(idx.x, idx.y).ckLocal()->doReduction();
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void CP_State_RealSpacePlane::doFFT(RSFFTMsg *msg) {

//============================================================================

    int size               = msg->size; 
    int Index              = msg->senderIndex;
    complex *partiallyFFTd = msg->data;
    int nchareG            = scProxy.ckLocalBranch()->cpcharmParaInfo->nchareG;
    int sizeY              = scProxy.ckLocalBranch()->cpcharmParaInfo->sizeY;
    int sizeZ              = scProxy.ckLocalBranch()->cpcharmParaInfo->sizeZ;
    int **tranUnpack       = scProxy.ckLocalBranch()->cpcharmParaInfo->index_tran_upack;
    int *nline_per_chareG  = scProxy.ckLocalBranch()->cpcharmParaInfo->nlines_per_chareG;

    int planeSize          = rs.size;

    count++;
    if (count > nchareG) {
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Mismatch in allowed gspace chare arrays : %d %d %d %d\n",
                count,nchareG,thisIndex.x,thisIndex.y);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

//============================================================================

    if(config.conserveMemory && count==1){rs.allocate();}
    complex *planeArr = rs.planeArr;
    if(count==1){memset(planeArr,0,planeSize*sizeof(complex));} 

// You have received packed data (x,y) from processor sendIndex
// Every real space chare receives the same x,y indicies.
// For double pack, x=0,1,2,3,4 ...  y= {-K ... K}
// The x increase with processor number. The y are split.
// The rundescriptor contains all we need to unpack the data.
// For doublepack : nffty*run[i][j].x + run[i][j].y
// we store this stuff in the convenient package
// Pictorially a half cylinder is sent which is unpacked into
// a half cube for easy FFTing. Y is the inner index.

    if(size!=nline_per_chareG[Index]){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Dude, %d != %d for chare %d %d\n",size,nline_per_chareG[Index],
                   thisIndex.y,Index);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    for(int i=0;i< size;i++){planeArr[tranUnpack[Index][i]] = partiallyFFTd[i];}

    delete msg;

    if (count == nchareG) {count=0;RTH_Runtime_resume(run_thread);}

//============================================================================
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealSpacePlane::doFFT()
{
        flagsRecd = true;
        count = 0;
#ifdef _CP_RS_VERBOSE_
        ckout << "Real Space " << thisIndex.x << " " << thisIndex.y << " doing FFT" << endl;
#endif
        doReduction();
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealSpacePlane::doReduction(){
//============================================================================
// double *data = rs.doRealFwFFT();  // 2-D forward FFT, in x*Z dimension,
// data will contain the squared magnitudes of the fftd data.    
// magnitude is computed and stored in *data
//============================================================================
// Perform the FFT

#ifndef CMK_OPTIMIZE    
      double StartTime=CmiWallTimer();
#endif

      double *data = fftCacheProxy.ckLocalBranch()->doRealFwFFT(rs.planeArr);

#ifndef CMK_OPTIMIZE
      traceUserBracketEvent(doRealFwFFT_, StartTime, CmiWallTimer());
#endif    

//============================================================================
// Perform the Reduction to get the density

//    CkPrintf("In RealSpacePlane[%d %d] doReduction %d\n", thisIndex.x, thisIndex.y,CmiMemoryUsage());

    if (config.useGReduction) {

      // we should be able to cache this setup stuff and only do it
      // once per iteration
        CkMulticastMgr *mcastGrp = 
            CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch(); 
        CkCallback cb(CkIndex_CP_Rho_RealSpacePlane::acceptDensity(0),
                      CkArrayIndex1D(thisIndex.y),
                      rhoRealProxy.ckGetArrayID());

#ifndef CMK_OPTIMIZE    
      StartTime=CmiWallTimer();
#endif
        mcastGrp->contribute(rs.planeSize[0] * sizeX * 
                             sizeof(double), data, CkReduction::sum_double, 
                             cookie, cb);
#ifndef CMK_OPTIMIZE
      traceUserBracketEvent(DoFFTContribute_, StartTime, CmiWallTimer());
      //delete [] data;
#endif    
    }else {
        rhoRealProxy[thisIndex.y].acceptDensity(rs.planeSize[0] * sizeX
                                                , data, thisIndex.y);
    }//endif

    delete [] data;

//============================================================================
    }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * In this function we multiply (element-to-element) psi_r with vks
 * We follow this up with finishing a part of the inverse fft (i.e., 
 * ifft in the x-direction. Once that is don, send the "half-baked"
 * data to the CP_State_GSpacePlanes
 */
//============================================================================
void CP_State_RealSpacePlane::doProduct(ProductMsg *msg) {
    /* for debugging*/

    if (msg->datalen != rs.planeSize[0] * sizeX)
        ckout << msg->datalen << endl;
	
    size=msg->datalen;
    vks = msg->data;

    RSDummyResume *pmsg= new (8*sizeof(int)) RSDummyResume;
    CkSetQueueing(pmsg, CK_QUEUEING_IFIFO);
    *(int*)CkPriorityPtr(pmsg) = config.rsifftpriority;
    thisProxy(thisIndex.x,thisIndex.y).resumeProduct(pmsg);
}
//============================================================================

void CP_State_RealSpacePlane::ResumeFromSync(){
    if(config.useCommlib)
	ComlibResetProxy(&gproxy);
}

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealSpacePlane::resumeProduct(RSDummyResume *msg){
    delete msg;
    RTH_Runtime_resume(run_thread);
    }


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealSpacePlane::doProduct(int Size, const double *Vks)
{
  CkAbort("Depricated doProduct \n");
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/* doing the product */
//============================================================================
void CP_State_RealSpacePlane::doProduct() {
//============================================================================
// check size

//#ifdef _CP_RS_VERBOSE_
//    CkPrintf("In RealSpacePlane[%d %d] doProduct %d\n", thisIndex.x, thisIndex.y,CmiMemoryUsage());
//#endif

//===================================================================
// debugging

#ifdef _CP_DEBUG_VKS_RSPACE_
    if(thisIndex.x==0 && thisIndex.y == 0){
      FILE *fp = fopen("vks_real_y0_state0.out","w");
      for(int i=0;i<size;i++){
        fprintf(fp,"%g\n",vks[i]);
      }//endfor
      fclose(fp);
    }
#endif    

//===================================================================
// Log the 2D-FFT call and do it : in place

#ifndef CMK_OPTIMIZE
    double StartTime=CmiWallTimer();
#endif

    fftCacheProxy.ckLocalBranch()->doRealBwFFT(vks,rs.planeArr,
                                       thisIndex.x,thisIndex.y);
    complex *vks_on_state = rs.planeArr;

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(doRealBwFFT_, StartTime, CmiWallTimer());
#endif

//===================================================================
// Perform the transpose and then the blast off the final 1D-FFT

    int nchareG    = scProxy.ckLocalBranch()->cpcharmParaInfo->nchareG;
    int sizeY      = scProxy.ckLocalBranch()->cpcharmParaInfo->sizeY;
    int sizeZ      = scProxy.ckLocalBranch()->cpcharmParaInfo->sizeZ;
    int planeSize  = sizeY*sizeX;
    int **tranpack = scProxy.ckLocalBranch()->cpcharmParaInfo->index_tran_upack;
    int *nlines_per_chareG = scProxy.ckLocalBranch()->cpcharmParaInfo->nlines_per_chareG;
    if (config.useCommlib){mssInstance.beginIteration();}
    for (int ic = 0; ic < nchareG; ic ++) { // chare arrays to which we will send
      int sendFFTDataSize = nlines_per_chareG[ic];
      GSIFFTMsg *msg = new (sendFFTDataSize, 8 * sizeof(int)) GSIFFTMsg; 
      msg->size      = sendFFTDataSize;
      msg->offset    = thisIndex.y;    // z-index
      complex *data  = msg->data;
      if(config.prioFFTMsg){
	  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	  *(int*)CkPriorityPtr(msg) = config.gsifftpriority+thisIndex.x*rs.planeSize[0];
      }//endif
      for(int i=0;i<sendFFTDataSize;i++){data[i] = vks_on_state[tranpack[ic][i]];}
      gproxy(thisIndex.x, ic).doIFFT(msg); // send the message
    }//end for : chare sending
    if (config.useCommlib){mssInstance.endIteration();}
    if(config.conserveMemory){
      rs.destroy();
    }else{
      rs.zeroOutPlanes();
    }//endif

//============================================================================
   }//end routine : CP_State_RealSpacePlane::doProduct()
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/* Setting up the multicast trees for Gengbin's library */
//============================================================================
void 
CP_State_RealSpacePlane::init(ProductMsg *msg)
{
    int i; contribute(sizeof(int), &i, CkReduction::sum_int);
    CkGetSectionInfo(cookie, msg);
    delete msg;
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealSpacePlane::printData() {
    char str[20];
    sprintf(str, "RSP%dx%d.out", thisIndex.x, thisIndex.y);
    FILE *outfile = fopen(str, "w");
    int ioff = 0;
    for (int i = 0; i < rs.planeSize[0]; i++) {
      for (int j = 0; j < rs.planeSize[1]; j++){
         fprintf(outfile, "%10.5lf", rs.planeArr[(j+ioff)].getMagSqr());
         fprintf(outfile, "\n");
      }
      ioff += rs.planeSize[0];
    }
    fclose(outfile);
}
//============================================================================
