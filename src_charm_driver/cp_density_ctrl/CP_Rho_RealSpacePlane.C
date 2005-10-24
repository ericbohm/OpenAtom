//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//  This is the description of the "life" of a CP_Rho_RealSpacePlane object.
// 
//  At the start of the program, the constructor CP_Rho_RealSpacePlane() is called.
// 
//  The CP_State_RealSpacePlanes send the squared magnitudes of the psi_r values 
//  using the acceptData() method. The squared magnitudes are summed across states.
// A copy of this data is made, inverse fft is done in the z and x directions 
// and sent to rhoGDensity. The other copy is processed using CP_exc_calc. 
// Then the CP_Rho_RealSpacePlane object waits for a reply from the RhoGDensity 
// object.
//
// The reply from RhoGDensity comes in the form of the method 
// acceptDensityForSumming(). The data obtained from this reply is taken and
// forward fft in z and x directions is performed. The resulting data is 
// summed with the result of CP_exc_calc. The sum is sent to the 
// CP_State_RealSpacePlane objects.
//
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
#include <fftlib.h>
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h"

//============================================================================

extern CProxy_CP_State_RealSpacePlane realSpacePlaneProxy;
extern CProxy_CP_Rho_RealSpacePlane rhoRealProxy;
extern CProxy_CP_Rho_GSpacePlane rhoGProxy;
extern CProxy_CPcharmParaInfoGrp scProxy;
extern ComlibInstanceHandle commRealInstance;
extern ComlibInstanceHandle mcastInstance;
extern CkGroupID mCastGrpId;
extern Config config;
extern int nstates;


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
RTH_Routine_locals(CP_Rho_RealSpacePlane,run)
//---------------------------------------------------------------------------
RTH_Routine_code(CP_Rho_RealSpacePlane,run) {
//---------------------------------------------------------------------------
  while(1) {
  

    RTH_Suspend(); 
    //ckout<<"US1 RHO_REAL_SPACE:"<<c->thisIndex << " acceptDensity"<<endl;
    c->acceptDensity();
     
    RTH_Suspend(); 
    //ckout<<"RHO_REAL_SPACE:"<<c->thisIndex << " doneFFT0"<<endl;
    c->doneFFT();
	
    RTH_Suspend(); 
    //ckout<<"RHO_REAL_SPACE:"<<c->thisIndex << " doneFFT1"<<endl;
    c->doneFFT();

    RTH_Suspend(); 
    ////ckout<<"RHO_REAL_SPACE:"<<c->thisIndex << " doneFFT2"<<endl;
    c->doneFFT();
    
    RTH_Suspend(); 
    //ckout<<"RHO_REAL_SPACE: "<<c->thisIndex << " doneFFT3"<<endl;
    c->doneFFT();
    if(!(c->gotAllRhoEnergy&&(c->doneDoingFFT)))// PRE: doneDoingFFT==TRUE
    {
      //      ckout<<"RHO_REAL_SPACE: "<<c->thisIndex << " pending rho "<< c->gotAllRhoEnergy << "fft " << c->doneDoingFFT <<endl;
	RTH_Suspend(); 
    }
    //    ckout<<"US4 RHO_REAL_SPACE: "<<c->thisIndex << " acceptEnergyForSumming"<<endl;
    c->acceptEnergyForSumming();

  } //end while not done
//--------------------------------------------------------------------------
   } RTH_Routine_end(CP_Rho_RealSpacePlane,run)
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::run () {
  run_thread = RTH_Runtime_create(RTH_Routine_lookup(CP_Rho_RealSpacePlane,run),this);
  RTH_Runtime_resume(run_thread);
}
//============================================================================

//! return tru if input is power of 2
bool is_pow2(int input)
{
    int y=0;
    for(int x=0;x<32;x++)
    {
	y= 1<<x;
	if(y==input)
	    return true;
    }
    return false;
}



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// This class (array) accepts the real space densities from all the states,
// adds them up and remembers them for the next stage in computation
//
//============================================================================
CP_Rho_RealSpacePlane::CP_Rho_RealSpacePlane(int xdim, size2d yzdim, 
				     int numRealSpace, int numRhoG, bool _useCommlib, 
				     ComlibInstanceHandle _fftcommInstance) 
//============================================================================
   {//begin routine
//============================================================================  

    initRhoRealSlab(&rho_rs, xdim, yzdim[0], yzdim[1], numRealSpace, 
                    numRhoG, thisIndex);
    count = 0;
    countFFTdata = 0;

    // volumeFactor = box volume / fft points
    //volumeFactor=0.1013576503098923;//
    CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
    volumeFactor = ((double)(sim->vol)) / config.numFFTPoints;

    probScale=(1.0/(xdim * yzdim[0] * yzdim[1] * volumeFactor));
    numMcastSent = 0;

    CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();

    // make sections in the realSpacePlane array. These will be used when 
    // computing real-space densities and multicasting v_ks values
    CkArrayIndexMax *elems = new CkArrayIndexMax[nstates];
    int j;
    // section i has al the portions with all 
    CkArrayIndex2D idx(0, thisIndex);
    if(is_pow2(nstates))
    {
	for (j = 0; j < nstates; j++) {
	    idx.index[0] = j^(thisIndex%nstates);
	    elems[j] = idx;
	}
    }
    else
    {
	for (j = 0; j < nstates; j++) {
	    idx.index[0] = (j+thisIndex)%nstates;
	    elems[j] = idx;
	}
    }

    realSpaceSectionProxy = CProxySection_CP_State_RealSpacePlane::
        ckNew(realSpacePlaneProxy.ckGetArrayID(), elems, nstates);

    realSpaceSectionCProxy = CProxySection_CP_State_RealSpacePlane::
        ckNew(realSpacePlaneProxy.ckGetArrayID(), elems, nstates);
    
    //if (config.useGMulticast) {
    realSpaceSectionProxy.ckDelegate
        (CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch());
    //}        
    mcastGrp->setSection(realSpaceSectionProxy);
    
    ProductMsg *dummyProductMessage = new (0) ProductMsg;    
    // inform density element of this section proxy.
    realSpaceSectionProxy.init(dummyProductMessage);

    //if (config.useCommlibMulticast) {
    ComlibAssociateProxy(&mcastInstance,realSpaceSectionCProxy);
//    ComlibInitSectionID(realSpaceSectionCProxy.ckGetSectionID());
    //} 
    rhoGProxy_com = rhoGProxy;
    if (config.useCommlib) {
	ComlibAssociateProxy(&commRealInstance,rhoGProxy_com);          
    }//endif

    delete [] elems;
        
    doneGradientCorrection = doneRhoGStuff = false;
    vectorFFTCount = 0;
    
    /* Gradient Correction stuff */
    CkAssert(config.rhoGPPC == 1); // read up fftlib docs
    int srcDim[2], dstDim[2];
    // TODO: read the fftlib docs and check if the following 4 assignments are
    // correct. For now, this does not matter, since all 4 numbers are equal
    srcDim[0] = rho_rs.sizeY; // source is RhoG
    srcDim[1] = rho_rs.sizeX;
    dstDim[0] = rho_rs.sizeZ; // destination is CP_Rho_RealSpacePlane
    dstDim[1] = rho_rs.sizeX;
    rhoIRX = (complex *) fftw_malloc(rho_rs.size*sizeof(complex));
    NormalFFTinfo fftinfos0(srcDim, 
                                    dstDim, false, rhoIRX, 


				    COMPLEX_TO_COMPLEX, 1, 1);
    rhoIRY = (complex *) fftw_malloc(rho_rs.size*sizeof(complex));
    NormalFFTinfo fftinfos1(srcDim, 
                                    dstDim, false, rhoIRY, 
				    COMPLEX_TO_COMPLEX, 1, 1);
    rhoIRZ = (complex *) fftw_malloc(rho_rs.size*sizeof(complex));
    NormalFFTinfo fftinfos2(srcDim, 
                                    dstDim, false, rhoIRZ, 
				    COMPLEX_TO_COMPLEX, 1, 1);
    gradientCorrection = (complex *) fftw_malloc(rho_rs.size*sizeof(complex));
    NormalFFTinfo fftinfos3(srcDim, 
				    dstDim, false, gradientCorrection, 
				    COMPLEX_TO_COMPLEX, 1, 1);

    //    fftuseCommlib = _useCommlib;
    //    fftcommInstance = _fftcommInstance;

    setup(fftinfos0, rhoGProxy, rhoRealProxy, _useCommlib, _fftcommInstance);
    setup(fftinfos1, rhoGProxy, rhoRealProxy, _useCommlib, _fftcommInstance);
    setup(fftinfos2, rhoGProxy, rhoRealProxy, _useCommlib, _fftcommInstance);
    setup(fftinfos3, rhoGProxy, rhoRealProxy, _useCommlib, _fftcommInstance);

    // this is a destination array, so
    // setup fwd1DPlan and bwd2DPlan
    // TODO: verify
    fwd1DPlan = fftw_create_plan(rho_rs.sizeX, FFTW_FORWARD, 
                                 FFTW_MEASURE|FFTW_IN_PLACE|FFTW_USE_WISDOM);
    bwd2DPlan = fftw2d_create_plan(rho_rs.sizeZ, rho_rs.sizeX, FFTW_BACKWARD, 
                                   FFTW_MEASURE|FFTW_IN_PLACE|FFTW_USE_WISDOM);
    doneDoingFFT=false;
    gotAllRhoEnergy=false;
    setMigratable(false);
    //ckout<<"starting run"<<endl;
    run();

//============================================================================
   }//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_Rho_RealSpacePlane::~CP_Rho_RealSpacePlane()
{
    fftw_free(rhoIRX);
    fftw_free(rhoIRY);
    fftw_free(rhoIRZ);
    fftw_free(gradientCorrection);
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::doneFFT(int Id)
{
  id=Id;
  if(Id==3)
    {
      //ckout<<"  setting doneDoingFFT"<<endl;
      doneDoingFFT=true;
      //      CkPrintf("[%d] has completed FFT\n",thisIndex);
      //doneFFT();
    }
  RTH_Runtime_resume(run_thread);

}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::doneFFT(){
    int i;
    if (id <= 2) {
        vectorFFTCount++;
        if (vectorFFTCount == 3) {
            vectorFFTCount = 0;
	    rho_rs.exc_gga_ret = 0.0;
            for (i = 0; i < rho_rs.size; i++){
                gradientCorrection[i].re = 0.0;
                gradientCorrection[i].im = 0.0;
	    }/*endfor*/
#define GGA_ON
#ifdef GGA_ON
            // compute the functional
            int npts         = rho_rs.size;
            int nf1          = rho_rs.sizeX;
            int nf2          = rho_rs.sizeY;
            int nf3          = rho_rs.sizeZ;
            complex *density = rho_rs.density;
            double *exc_gga_ret = &(rho_rs.exc_gga_ret);

            CPXCFNCTS::CP_getGGAFunctional(npts,nf1,nf2,nf3,density,
             rhoIRX,rhoIRY,rhoIRZ,gradientCorrection,thisIndex,exc_gga_ret);
	    // add (dfxc_dn_x + dfxc_dn_c) to vks
            for (i = 0; i < rho_rs.size; i++){ 
                rho_rs.Vks[i].re += gradientCorrection[i].re;
                rho_rs.Vks[i].im  = 0.0;
	    }/*endfor*/
#endif
#ifdef CMK_VERSION_BLUEGENE
     CmiNetworkProgress();
#endif
            doIFFT(0,0);
#ifdef CMK_VERSION_BLUEGENE
     CmiNetworkProgress();
#endif
            doIFFT(1,1);
#ifdef CMK_VERSION_BLUEGENE
     CmiNetworkProgress();
#endif
            doIFFT(2,2);
        }
    }
    else {
        // this is White-Byrd correction factor that
        // is added to the KS potential in real space.
#ifdef GGA_ON
        for (i = 0; i < rho_rs.size; i++){
	  rho_rs.Vks[i].re -= gradientCorrection[i].re;
            rho_rs.Vks[i].im  = 0.0;
        }/*endfor*/
#else
        rho_rs.exc_gga_ret = 0.0;
        for (i = 0; i < rho_rs.size; i++){
            gradientCorrection[i].re = 0.0;
            gradientCorrection[i].im = 0.0;
        }/*endfor*/
#endif

        doneGradientCorrection = true;

	double exc[2];
	exc[0]=rho_rs.exc_ret;
	exc[1]=rho_rs.exc_gga_ret;

	contribute(2*sizeof(double),exc, CkReduction::sum_double);

        if (doneRhoGStuff) 	    
            doMulticast();	
    }
//============================================================================
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// Here the density from all the states is added up. The data from all the
// states is received via an array section reduction.
//
//============================================================================
void CP_Rho_RealSpacePlane::acceptDensity(CkReductionMsg *msg) {

    CkAssert(msg->getSize() == rho_rs.sizeZ * rho_rs.sizeX * sizeof(double));
    double *realValues = (double *) msg->getData(); 

#ifdef _CP_DEBUG_RHO_
    if(thisIndex==0){
      char myFileName[MAX_CHAR_ARRAY_LENGTH];
      if(config.doublePack){
	sprintf(myFileName, "dpRho_Real_%d_%d.out", thisIndex, CkMyPe() );
      }else{
	sprintf(myFileName, "Rho_Real_%d_%d.out", thisIndex, CkMyPe() );
      }//endif
      FILE *fp = fopen(myFileName,"a");
        for (int i = 0; i <rho_rs.sizeZ*rho_rs.sizeX; i++){
          fprintf(fp,"%g\n",realValues[i]*probScale);
        }//endfor
      fclose(fp);
    }//endif
#endif
    for(int i = 0; i < rho_rs.sizeZ*rho_rs.sizeX; i++){
      rho_rs.doFFTonThis[i] = complex(realValues[i] * probScale, 0);	
    }//endfor
    delete msg;

//    CkPrintf("In RhoRealSpacePlane[%d] accept density %d %d\n",thisIndex, CkMyPe(),CmiMemoryUsage());
    //energyComputation();
    //acceptDensity();
    RTH_Runtime_resume(run_thread);

//============================================================================
  }//end routine
//============================================================================

void CP_Rho_RealSpacePlane::ResumeFromSync()
{

    if(config.useCommlib)
    {
	ComlibResetSectionProxy(&realSpaceSectionCProxy);
	ComlibResetProxy(&rhoGProxy_com);
    }//endif


}
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::acceptDensity() {
	energyComputation();
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// 
// Make a copy of densities and compute one part of the energies using
// CP_exc_calc
//
//============================================================================
void CP_Rho_RealSpacePlane::energyComputation()
//============================================================================
    {//begin routine 
//============================================================================

// CkPrintf("In RhoRealSpacePlane[%d] energyComp %d %d\n",thisIndex, CkMyPe(),CmiMemoryUsage());

    CmiMemcpy(rho_rs.density, rho_rs.doFFTonThis, 
           rho_rs.sizeX * rho_rs.sizeZ * sizeof(complex));
    // doing the rho_r subroutine here:
    
    /*
     * Conflict with rho_rs.doBwFFT()
     *
     */
    complex *density = rho_rs.density;
    complex *Vks     = rho_rs.Vks;
    int size         = rho_rs.size;
    int nf1          = rho_rs.sizeX;
    int nf2          = rho_rs.sizeY;
    int nf3          = rho_rs.sizeZ;
    double *exc_ret  = &(rho_rs.exc_ret);
    double *muxc_ret = &(rho_rs.muxc_ret);

    // perform exchange correlation computation
    // doesn't depend on gradient of the density in realspace
    CPXCFNCTS::CP_exc_calc(size,nf1,nf2,nf3,density,Vks,exc_ret,muxc_ret);
    // get Rho(G) to compute Hartree and E_EXT as well as div(rho)

#ifdef CMK_VERSION_BLUEGENE
     CmiNetworkProgressAfter(1);
#endif


    startTranspose();
    //contribute(sizeof(int), &count, CkReduction::sum_int);
//============================================================================
    }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// do FFT on rho_rs.doFFTonThis
//
//============================================================================
void CP_Rho_RealSpacePlane::startTranspose()
//============================================================================
   {//begin routine
//============================================================================
// CkPrintf("In RhoRealSpacePlane[%d] startTranspose %d %d\n",thisIndex, CkMyPe(),CmiMemoryUsage());
    int planeSize = rho_rs.sizeX * rho_rs.sizeZ;
    int dataSize = rho_rs.sizeX;
    int maxk = rho_rs.sizeZ;
    int offset = thisIndex % maxk;
#ifndef CMK_OPTIMIZE
    double StartTime=CmiWallTimer();
#endif
    // do z-x plane fft in-place, on data rho_rs.doFFTonThis
    complex *planePtr = rho_rs.doFFTonThis;

#ifndef CMK_OPTIMIZE
    StartTime=CmiWallTimer();
#endif

    fftwnd_one(rho_rs.fft2dBwPlan, (fftw_complex *)planePtr, NULL);

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(RhoRtoGxzFFT_, StartTime, CmiWallTimer());    
#endif

        
    // send chunk to RhoGDensity
    int l, k;
    if (config.useCommlib) {
        commRealInstance.beginIteration();
    }//endif

//    CkPrintf("In RhoRealSpacePlane[%d] startTranspose %d %d %d\n",thisIndex,
//                   rho_rs.sizeZ,rho_rs.sizeX,CmiMemoryUsage());
    for (l = 0; l < rho_rs.sizeZ; l++){
      complex *dataPtr = planePtr + l * rho_rs.sizeX;
      rhoGProxy_com[l].acceptData(dataSize, dataPtr, thisIndex, l, 0);
    }//endfor

    if (config.useCommlib){commRealInstance.endIteration();}

//============================================================================
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::acceptDensity(int size, double *arr, int ycoord) {
    /* for debugging: should never be true*/ 

   if (count >= nstates) {
        ckerr << "Phase Three " << thisIndex << " got one plane too many" << endl;
        return;
    }
    
    int z, x;
    // loop to copy received data: realspace density
    for(z = 0; z < rho_rs.sizeZ; z++){
      for(x = 0; x < rho_rs.sizeX; x++){
        rho_rs.doFFTonThis[z * rho_rs.sizeX + x].re +=
                          probScale * arr[z * rho_rs.sizeX + x];
      }//endfor
    }//endfor
    
    count++;
    if (count == nstates) {
        count = 0;
        energyComputation();
    }//endif

//============================================================================
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::doMulticast(){
//============================================================================
    doneGradientCorrection = doneRhoGStuff = false;
    // Broadcast the data we have to all the states
    if (config.useGMulticast || config.useCommlibMulticast) {
        int s;
        int dataSize = rho_rs.sizeX * rho_rs.sizeZ;
        ProductMsg *msg = new (dataSize, 0) ProductMsg;
        //CmiMemcpy(msg->data, rho_rs.Vks, sizeof(complex) * dataSize);

        for(int count = 0; count < dataSize; count ++)
            msg->data[count] =  rho_rs.Vks[count].re;

        msg->idx = thisIndex;
        msg->datalen = dataSize;
        msg->hops = 0;

        if (config.useCommlibMulticast) 
            mcastInstance.beginIteration();
        
        if(config.useCommlibMulticast)
            realSpaceSectionCProxy.doProduct(msg);
        else
            realSpaceSectionProxy.doProduct(msg);

        if (config.useCommlibMulticast)
            mcastInstance.endIteration();
    }              
    else
        CkAbort("No multicast strategy\n");
    
    //re-initialize the array into which we will sum the squared densities.
    //rho_rs.vks does not need to be re-initialized as its contents will be
    //over-written by those of rho_rs.doFFTonThis. 
    // This is to zero out the sums for next interation
#ifdef _CP_DEBUG_VKS_RSPACE_
    if(thisIndex==0){
      FILE *fp = fopen("vks_real_y0.out","w");
      for(int i=0;i<rho_rs.size;i++){
        fprintf(fp,"%g\n",rho_rs.Vks[i].re);
      }//endfor
      fclose(fp);
    }
#endif
    bzero(rho_rs.doFFTonThis,sizeof(complex)*rho_rs.size);
    bzero(rho_rs.Vks,sizeof(complex)*rho_rs.size);
//============================================================================
   }//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::resumeThread(PPDummyMsg *dmsg) {

  delete dmsg;
  //ckout<<"-->testing doneDoing "<<doneDoingFFT<<endl;
  
  CkAbort("how the hell did we get here?");
  if(doneDoingFFT){
    //ckout<<"-->Resuming Thread1"<<endl;
    RTH_Runtime_resume(run_thread);
  }else{
    PPDummyMsg *msg = new(8*sizeof(int)) PPDummyMsg;
    CkSetQueueing(msg, CK_QUEUEING_IFIFO);
    *(int*)CkPriorityPtr(msg) = config.rhorpriority;
    rhoRealProxy(thisIndex).resumeThread(msg);
  }//endif

}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// In this method the energies from rho-g are received
//
//============================================================================
void CP_Rho_RealSpacePlane::acceptEnergyForSumming(int size, complex *densities,
                                                   int zindex) 
//============================================================================
   {//begin routine
//============================================================================

    int planeSize = rho_rs.sizeX * rho_rs.sizeZ;

    CmiMemcpy(rho_rs.doFFTonThis + zindex*rho_rs.sizeX,densities,
		 rho_rs.sizeX * sizeof(complex));

    countFFTdata++;
    if (countFFTdata == rho_rs.sizeZ) {
	gotAllRhoEnergy=true;
	//	CkPrintf("[%d] has all [%d] rho energy inputs\n",thisIndex,countFFTdata);
	//	  acceptEnergyForSumming();
	if(doneDoingFFT)
	  RTH_Runtime_resume(run_thread);
    }else{// do nothing 
/*	PPDummyMsg *msg = new(8*sizeof(int)) PPDummyMsg;
*	CkSetQueueing(msg, CK_QUEUEING_IFIFO);
*	*(int*)CkPriorityPtr(msg) = config.rhorpriority;
*	rhoRealProxy(thisIndex).resumeThread (msg);
*/
    }//endif

//============================================================================    
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::acceptEnergyForSumming() {
	doneDoingFFT = false;
	gotAllRhoEnergy=false;
        countFFTdata = 0;
        // do the forward plane ffts here
        rho_rs.doFwFFT();
        
        // sum energies
        int i;
        for (i = 0; i < rho_rs.size; i++)
	  rho_rs.Vks[i] += rho_rs.doFFTonThis[i];   
        
		doneRhoGStuff = true;
		if (doneGradientCorrection) {
		  doMulticast();
		}
                
}
//============================================================================





