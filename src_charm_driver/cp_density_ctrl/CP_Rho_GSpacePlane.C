//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
//  This is a description of the "life" of a CP_Rho_GSpacePlane  object
// 
//  At the start of the program, the constructor CP_Rho_GSpacePlane is called.
//  The RealSpaceDensity objects send data to CP_Rho_GSpacePlane using the 
//  acceptData() method. Inverse ffts in z and x directions are performed 
//  before the data is received, so here inverse fft in y direction is 
//  performed. This data is processed using the CP_hart_eext_calc. Then forward
//  fft in the y direction is performed and data is send back to 
//  RealSpaceDensity objects.
// 
//  The CP_Rho_GSpacePlaneHelper objects essentially help to split the work involved
//  in g-space density computation. They receive their share of the work
//  through the method recvCP_Rho_GSpacePlanePart() and send the processed 
//  data back to CP_Rho_GSpacePlane objects using the recvProcessedPart() method.
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
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cplocal.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h"

//============================================================================

extern Config config;
extern CProxy_CP_Rho_RealSpacePlane rhoRealProxy;
extern int nstates;
extern CProxy_CPcharmParaInfoGrp scProxy;
extern CProxy_AtomsGrp atomsGrpProxy;
extern ComlibInstanceHandle commInstance;
extern CProxy_ComlibManager mgrProxy;
extern CProxy_CP_Rho_GSpacePlane rhoGProxy;
extern CProxy_CP_Rho_GSpacePlaneHelper rhoGHelperProxy;

//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_Rho_GSpacePlane::CP_Rho_GSpacePlane(int xdim, size2d sizeYZ, 
	       int numRealSpace, int numRealSpaceDensity, bool _useCommlib,
	       ComlibInstanceHandle _fftcommInstance) 
//============================================================================
   {//begin routine
//============================================================================

    CkAssert(numRealSpaceDensity == 1);
    initRhoGSlab(&rho_gs, xdim, sizeYZ[0], sizeYZ[1], numRealSpace, 
                 numRealSpaceDensity, thisIndex);

    count = 0;
    helperWidth = sizeYZ[1]/config.rhoGHelpers;

    if( (sizeYZ[1] % config.rhoGHelpers) !=0 ){
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("Helper size must be a mod of %d.\n",sizeYZ[1]);
       CkPrintf("Please fix your cpaimd_config.\n");
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       CkExit();
    }//endif

    helperCount = 0;
    CkAssert(config.rhoGPPC == 1); // read up fftlib docs
    int srcDim[2], dstDim[2];

    // TODO: read the fftlib docs and check if the following 4 assignments are
    // correct. For now, this does not matter, since all 4 numbers are equal
    srcDim[0] = rho_gs.sizeY; // source is CP_Rho_GSpacePlane
    srcDim[1] = rho_gs.sizeX;
    dstDim[0] = rho_gs.sizeZ; // dst is CP_Rho_RealSpacePlane
    dstDim[1] = rho_gs.sizeX;

    rhoIGX = new complex[rho_gs.size];
    NormalFFTinfo fftinfos0(srcDim, dstDim, 
                                    true, rhoIGX, COMPLEX_TO_COMPLEX, 1, 1);

    rhoIGY = new complex[rho_gs.size];
    NormalFFTinfo fftinfos1(srcDim, dstDim, 
                                    true, rhoIGY, COMPLEX_TO_COMPLEX, 1, 1);

    rhoIGZ = new complex[rho_gs.size];
    NormalFFTinfo fftinfos2(srcDim, dstDim, 
                                    true, rhoIGZ, COMPLEX_TO_COMPLEX, 1, 1);
    gradientCorrection = new complex[rho_gs.size];
    NormalFFTinfo fftinfos3(srcDim, 
                                    dstDim, true, gradientCorrection,
				    COMPLEX_TO_COMPLEX, 1, 1);

    setup(fftinfos0, rhoGProxy, rhoRealProxy, _useCommlib, _fftcommInstance);
    setup(fftinfos1, rhoGProxy, rhoRealProxy, _useCommlib, _fftcommInstance);
    setup(fftinfos2, rhoGProxy, rhoRealProxy, _useCommlib, _fftcommInstance);
    setup(fftinfos3, rhoGProxy, rhoRealProxy, _useCommlib, _fftcommInstance);

    
    // this is a src array, so
    // setup bwd1DPlan and fwd2DPlan
    // TODO: verify
    bwd1DPlan = fftw_create_plan(rho_gs.sizeY, FFTW_BACKWARD, 
                                 FFTW_MEASURE|FFTW_IN_PLACE|FFTW_USE_WISDOM);

    fwd2DPlan = fftw2d_create_plan(rho_gs.sizeY, rho_gs.sizeX, FFTW_FORWARD, 
                                   FFTW_MEASURE|FFTW_IN_PLACE|FFTW_USE_WISDOM);
    vectorIFFTCount = 0;
    
    k_x = new int[rho_gs.size];		
    k_y = new int[rho_gs.size];		
    k_z = new int[rho_gs.size];		

    computeK(rho_gs.size,rho_gs.startx, rho_gs.starty, rho_gs.startz,
             rho_gs.xdim, rho_gs.ydim, rho_gs.zdim);
    
    setMigratable(false);

    //    fftuseCommlib = _useCommlib;
    //    fftcommInstance = _fftcommInstance;

//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_Rho_GSpacePlane::~CP_Rho_GSpacePlane()
{
    delete [] rhoIGX;
    delete [] rhoIGY;
    delete [] rhoIGZ;
    delete [] k_x;
    delete [] k_y;
    delete [] k_z;

    delete [] gradientCorrection;
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GSpacePlane::computeK(const int size,
                    const int startx, const int starty, const int startz,
                    const int xdim, const int ydim, const int zdim) 
//============================================================================
   {//begin routine
//============================================================================

  int ksize = rho_gs.xdim;

  memset(k_x, 0, size * sizeof(int));
  memset(k_y, 0, size * sizeof(int));
  memset(k_z, 0, size * sizeof(int));

  int y, z, x, index = 0;
  int xcoord, ycoord, zcoord;
  for (y = 0, ycoord = starty; y < ydim; y++, ycoord++)
    for (z = 0, zcoord = startz; z < zdim; z++, zcoord++)
      for (x = 0, xcoord = startx; x < xdim; x++, xcoord++) {

        k_x[index] = xcoord;         
        if(k_x[index] > ksize/2)
          k_x[index] -= ksize;

        k_y[index] = zcoord;
        if(k_y[index] > ksize/2)
          k_y[index] -= ksize;

        k_z[index] = ycoord;
        if(k_z[index] > ksize/2)
          k_z[index] -= ksize;

        index++;
      }
//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * ycoord  : this will identify the chunk which is sending the data
 * planeNum: this identifies the plane a portion of which we have got
 * zcoord  : this identifies the portion of the plane that we hava got
 */
//============================================================================
void CP_Rho_GSpacePlane::acceptData(int size, complex *densities, int yindex, 
                      int zcoord, int planeNum)
//============================================================================
{
	//cout<<"-->accepting data"<<endl;
	// copy the data at the right place
    int yoffset = yindex + planeNum;
    memcpy(rho_gs.chunk + yoffset * size, densities, size * sizeof(complex));  //<-

    count++;
    acceptData();
} 
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GSpacePlane::acceptData() { 
//============================================================================

    if (count == rho_gs.sizeY) {
        count = 0;

        // do the line inv-fft, in y-D, in-place on rho_rs.chunk 
#ifndef CMK_OPTIMIZE
	double StartTime=CmiWallTimer();
#endif

	rho_gs.doBwFFT(thisIndex);

#ifndef CMK_OPTIMIZE
	traceUserBracketEvent(RhoRtoGyFFT_, StartTime, CmiWallTimer());    
#endif


        CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
        double volumeFactor = ((double)(sim->vol)) / config.numFFTPoints;
	for(int i=0; i<rho_gs.size; i++) {
	  rho_gs.chunk[i] *= complex(volumeFactor,0);
	}
#ifdef _CP_DEBUG_VKS_GSPACE_
        if(thisIndex==0){
          char myFileName[100];
   	  sprintf(myFileName, "Rho_Gspace_%d.out", thisIndex);
 	  FILE *fp = fopen(myFileName,"w");
  	    for (int i = 0; i < rho_gs.size; i++){ 
              fprintf(fp," %d %d %d : %g %g\n",k_x[i],k_y[i],k_z[i],
                          rho_gs.chunk[i].re,rho_gs.chunk[i].im);
 	    }//endfor
          fclose(fp);
	}//endif
#endif
        // Make a copy for gradient correction
        memcpy(gradientCorrection, rho_gs.chunk, rho_gs.size * sizeof(complex));

        CPXCFNCTS::CP_div_rho_gspace_calc(gradientCorrection,k_x, k_y, k_z, 
                                          rho_gs.size,rhoIGX,rhoIGY,rhoIGZ);
#ifndef CMK_OPTIMIZE
	StartTime=CmiWallTimer();
#endif

        doFFT(0,0); //rhoIGX

#ifndef CMK_OPTIMIZE
	traceUserBracketEvent(RhoDivRhoXFFT_, StartTime, CmiWallTimer());    
#endif

#ifndef CMK_OPTIMIZE
	StartTime=CmiWallTimer();
#endif

        doFFT(1,1); //rhoIGY

#ifndef CMK_OPTIMIZE
	traceUserBracketEvent(RhoDivRhoYFFT_, StartTime, CmiWallTimer());    
#endif

#ifndef CMK_OPTIMIZE
	StartTime=CmiWallTimer();
#endif

        doFFT(2,2); //rhoIGZ

#ifndef CMK_OPTIMIZE
	traceUserBracketEvent(RhoDivRhoZFFT_, StartTime, CmiWallTimer());    
#endif

        if( (rho_gs.sizeY % helperWidth) !=0 ){
          CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
          CkPrintf("GHelper size must be a mod of %d.\n",rho_gs.sizeY);
          CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
          CkExit();
        }//endif

        int y;
        for (y = 0; y < rho_gs.sizeY; y += helperWidth){
	  rhoGHelperProxy(thisIndex, y).recvRhoGPart(helperWidth*rho_gs.sizeX,
                                             (rho_gs.chunk+y*rho_gs.sizeX));
	}//endfor

    }//endif : routine received all parts expected

//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// For the parallel fft library
//============================================================================
void CP_Rho_GSpacePlane::doneIFFT(int id)
{
    vectorIFFTCount++;
    doneIFFT();
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GSpacePlane::doneIFFT()
{
    if (vectorIFFTCount == 3) {
        // all 3 IFFTs are done, use the 3 vectors to get new data
        vectorIFFTCount = 0;
        CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;      
        int nf1 = sim->sizeX;
        int nf2 = sim->sizeY;
        int nf3 = sim->sizeZ;
        int npts = rho_gs.size;
        CPXCFNCTS::CP_white_byrd_gspace_calc(rhoIGX,rhoIGY,rhoIGZ, 
              k_x,k_y,k_z,npts,nf1,nf2,nf3,gradientCorrection);
        doFFT(3, 3);
    }
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GSpacePlane::recvProcessedPart(int size, complex *points, int pos)
{
    memcpy(rho_gs.chunk + pos*rho_gs.sizeX, points, size*sizeof(complex)); 
    helperCount++;
    recvProcessedPart();
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GSpacePlane::recvProcessedPart(){
//============================================================================

    if (helperCount == config.rhoGHelpers) {
        helperCount = 0;
#ifndef CMK_OPTIMIZE
	double StartTime=CmiWallTimer();
#endif
        rho_gs.doFwFFT();  

#ifndef CMK_OPTIMIZE
	traceUserBracketEvent(VksofGFFT_, StartTime, CmiWallTimer());    
#endif

        // now send data back to rhoReal
        int dataSize = rho_gs.sizeX;

        int maxk = rho_gs.sizeY;
        int offset = thisIndex % maxk;

        CProxy_CP_Rho_RealSpacePlane rhoRealProxy_com = rhoRealProxy;;
        if (config.useCommlib) {
            commInstance.beginIteration();
            ComlibDelegateProxy(&rhoRealProxy_com);
        }
        
        for (int y = 0; y < rho_gs.sizeY; y ++){
            rhoRealProxy_com[y].acceptEnergyForSumming
                (dataSize, rho_gs.chunk + y*rho_gs.sizeX, thisIndex);
	}//endfor
        
        if (config.useCommlib)
            commInstance.endIteration();
    }// helpr = config.rhoGhelpers

//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_Rho_GSpacePlaneHelper::CP_Rho_GSpacePlaneHelper(int sizeX, size2d sizeYZ, 
                                                   int pos)
//============================================================================
   {//begin routine
//============================================================================
    helperWidth = sizeYZ[0]/config.rhoGHelpers;

    if( (sizeYZ[0] % config.rhoGHelpers) !=0 ){
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("Helper size must be a mod of %d.\n",sizeYZ[0]);
       CkPrintf("Please fix your cpaimd_config.\n");
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       CkExit();
    }//endif

    helperPos = pos; //starting y coord in YX plane
    rho_gs.size = helperWidth * sizeX;
    rho_gs.chunk = new complex[rho_gs.size];
    temp = new complex[rho_gs.size];
    k_x = new int[rho_gs.size];
    k_y = new int[rho_gs.size];
    k_z = new int[rho_gs.size];
    rho_gs.sizeX = sizeX;
    rho_gs.sizeY = sizeYZ[0];
    rho_gs.sizeZ = sizeYZ[1];
    rho_gs.xdim  = sizeX;

    int ksize = rho_gs.xdim;
    int y, x, i = 0;
    for (y = 0; y < helperWidth; y++)
        for (x = 0; x < sizeX; x++){
	    k_x[i] = x;

	    if(k_x[i] > ksize/2)
	      k_x[i] -= ksize;
	    
            k_z[i] = y + helperPos;

	    if(k_z[i] > ksize/2)
	      k_z[i] -= ksize;
	    
            k_y[i] = thisIndex.x;

	    if(k_y[i] > ksize/2)
	      k_y[i] -= ksize;

	    
            i++;
        }
    
    setMigratable(false);
//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_Rho_GSpacePlaneHelper::~CP_Rho_GSpacePlaneHelper()
//============================================================================
{
    delete [] k_x;	
    delete [] k_y;	
    delete [] k_z;	
    delete [] rho_gs.chunk;
    delete [] temp;
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GSpacePlaneHelper::recvRhoGPart(int size, complex *points)
//============================================================================
   {//begin routine
//============================================================================
   CkAssert(size == rho_gs.size);

   AtomsGrp *ag = atomsGrpProxy.ckLocalBranch(); // find me the local copy

   memcpy(rho_gs.chunk, points, size * sizeof(complex));
   complex *chunk    = rho_gs.chunk;
   int nsize         = (int) (rho_gs.size);
   int natm          = ag->natm;
   Atom *atoms       = ag->atoms;
   double *ehart_ret = &(rho_gs.ehart_ret);
   double *eext_ret  = &(rho_gs.eext_ret);
   double *ewd_ret   = &(rho_gs.ewd_ret);

#ifdef RHO_DBLE_PACK
   int mydoublePack = config.doublePack;
#else
   int mydoublePack = 0;  
#endif
   CPLOCAL::CP_hart_eext_calc(nsize,chunk,natm,atoms,temp,
                              ehart_ret,eext_ret,ewd_ret,k_x,k_y,k_z,
                              mydoublePack);
    // temp is Vks from EHartee and Eext in G Space !!!


    //ckout<<"Contributing"<<thisIndex.x<<endl;

   double e[3];
   e[0] = rho_gs.ehart_ret;
   e[1] = rho_gs.eext_ret;
   e[2] = rho_gs.ewd_ret;

   contribute(3 * sizeof(double), e, CkReduction::sum_double);
   //DEBUG HACK
   //bzero(temp, size * sizeof(complex));
   rhoGProxy[thisIndex.x].recvProcessedPart(size, temp, helperPos);  
    // fwd FFT is done to produce Vks in Real space
    // The result of first 1D FFT is store in rho_gs->chunk
    // ????
//---------------------------------------------------------------------------
   }//end routine
//============================================================================



