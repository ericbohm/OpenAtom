//=========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=========================================================================
/** \file CP_State_ParticlePlane.C
 * Life-cycle of a CP_State_ParticlePlane:
 *
 * The particle-plane array is a shadow of the g-space planes. This means
 * that particlePlaneProxy(s, p) exists on the same processor as
 * gSpacePlaneProxy(s, p). Also,  both chare arrays have the same number of
 * objects.
 *
 * The life-cycle of this object involves computation of a portion of
 * the Z-matrix, reduction of the Z matrix over g-space and computing
 * forces using the particle data.  ComputeZ is triggered by the
 * arrival of the structure factor for each atom group at the local
 * sfcache.
 *
 * The computation of the Z-matrix is done in the method computeZ().
 * Immediately after this, the matrix is reduced using the reduceZ() method.
 * Through the reduceZ() method, a portion of the reduced matrix is 
 * available at the zeroth plane of each state.
 * 
 * The reduced matrix is sent from the particle-plane 0 to all the particle
 * planes in the state using the getForces() method, which computes
 * forces and adds these forces to the forces compputed using electron
 * density in the g-space planes. 
 *
 * Important information:
 * 1. The dimensionality of the Z-matrix depends upon how may orbital levels
 * are considered (l = 0,1...). The current code assumes that only l = 0
 * is allowed. This has to be changed.
 * 
 * 2. The special point gx=gy=gz=0 is not dealt with.
 */ 
//=========================================================================

#include "charm++.h"
#include "util.h"
#include "cpaimd.h"
#include "groups.h"
#include "fftCacheSlab.h"
#include "eesCache.h"
#include "CP_State_Plane.h"
#include "StructFactorCache.h"
#include "ckmulticast.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cplocal.h"

//=========================================================================
extern CProxy_main                       mainProxy;
extern CProxy_CP_State_GSpacePlane       gSpacePlaneProxy;
extern CProxy_AtomsGrp                   atomsGrpProxy;
extern CProxy_CPcharmParaInfoGrp         scProxy;
extern CProxy_CP_State_ParticlePlane     particlePlaneProxy;
extern CProxy_CP_State_RealParticlePlane realParticlePlaneProxy;
extern CProxy_StructFactCache            sfCacheProxy;
extern CProxy_EnergyGroup                egroupProxy; //energy group proxy
extern CProxy_eesCache                   eesCacheProxy;
extern CProxy_FFTcache                   fftCacheProxy;

extern CkGroupID            mCastGrpId;
extern ComlibInstanceHandle gssPInstance;

extern int    nstates;
extern int    nchareG;
extern Config config;

//#define _CP_DEBUG_STATE_GPP_VERBOSE_

//=========================================================================



//============================================================================
// Reduction client for N^3 nonlocal energy computation
//============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void printEnl(void *param, void *msg){

  CkReductionMsg *m=(CkReductionMsg *)msg; //unpack
  double d = ((double *)m->getData())[0];
  delete m;

  CkPrintf("ENL         = %5.8lf\n", d);   // tell the world

  gSpacePlaneProxy(0,0).computeEnergies(ENERGY_ENL, d);  //store it
}
//============================================================================



//============================================================================
// The glorious constructor
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_State_ParticlePlane::CP_State_ParticlePlane(
       int x, int y, int z, int xNL, int yNL, int zNL, 
       int _gSpacePPC, int numSfGrps_in,int natm_nl_in, int natm_nl_grp_max_in,
       int _nstates, int _nchareG, int _Gstates_per_pe, int _numNLiter,
       int _ees_nonlocal)
//============================================================================
  {//begin routine
//============================================================================

  sizeX                = x;
  sizeY                = y;
  sizeZ                = z;
  ngridaNL             = xNL; 
  ngridbNL             = yNL; 
  ngridcNL             = zNL; 
  ees_nonlocal         = _ees_nonlocal;
  nstates	       = _nstates;
  nchareG	       = _nchareG;
  Gstates_per_pe       = _Gstates_per_pe;
  numNLiter            = _numNLiter;
  myChareG             = thisIndex.y;
  gSpacePlanesPerChare = _gSpacePPC;

  iteration            = 0;
  iterNL               = 0;
  countNLIFFT          = 0;
  doneEnl              = 0;
  doneForces           = 0;
  enl                  = 0.0;
  energy_count         = 0;
  totalEnergy          = 0.0;
  sendDone             = 0;

  numSfGrps            = numSfGrps_in;
  natm_nl              = natm_nl_in;
  natm_nl_grp_max      = natm_nl_grp_max_in;
  count                = new int[numSfGrps];
  haveSFAtmGrp         = new int[numSfGrps];
  memset(haveSFAtmGrp,-1,numSfGrps*sizeof(int));
  bzero(count,numSfGrps*sizeof(int));

  zsize                = numSfGrps*natm_nl_grp_max;

  zmatrixSum    = NULL;
  zmatrix       = NULL;
  zmatrix_fx    = NULL; zmatrix_fy    = NULL; zmatrix_fz    = NULL;
  zmatrixSum_fx = NULL; zmatrixSum_fy = NULL; zmatrixSum_fz = NULL;

//============================================================================
// Comlib to talk to realPP : used when ees_enl_on==1

  realPP_proxy = realParticlePlaneProxy;
  if (config.useGssInsRealPP){
      ComlibAssociateProxy(&gssPInstance,realPP_proxy);
  }//endif

//============================================================================
//  Register with the SFCache

  if(ees_nonlocal==0){
#ifdef _CP_DEBUG_SF_CACHE_
  CkPrintf("PP [%d,%d] has %d numSfGrps\n",thisIndex.x, thisIndex.y, numSfGrps);
#endif
    StructFactCache *sfcache = sfCacheProxy.ckLocalBranch();
    for(int i=0;i<numSfGrps;i++){
#ifdef _CP_DEBUG_SF_CACHE_
      CkPrintf("PP [%d,%d] registers grp %i with SFC[%d]\n",
         thisIndex.x, thisIndex.y,i, CkMyPe());
#endif
      sfcache->registerPP(thisIndex.x, thisIndex.y,i);
    }//endfor
  }//endif

//============================================================================
// Create section proxy for ENL reduction ParticlePlane (any state#, reductionPlaneNum)
// This is used when ees_enl_on ==0 

  //--------------------------------------------------------------------------
  // Compute all reduction planes for all chares

#ifdef USE_TOPOMAP
  int *red_pl = new int[nstates];
  int l       = Gstates_per_pe;

  int pl = nstates / l;
  int pm = CkNumPes() / pl; if(pm==0){CkAbort("Choose a larger Gstates_per_pe\n");}
  int m  = (nchareG / pm);

  int planes_per_pe = m;
  for(int i=0; i<nstates;i++){
    red_pl[i]=((i % Gstates_per_pe)*planes_per_pe)%nchareG;
  }//endfor
  reductionPlaneNum    = red_pl[thisIndex.x];
#else
  reductionPlaneNum    = calcReductionPlaneNum(thisIndex.x);
#endif

  //--------------------------------------------------------------------------
  // If you are a reduction plane, set up your comm with your guys and the
  // the other reduction planes. Store your cookie to eat later.

  if(thisIndex.x==0 && thisIndex.y==reductionPlaneNum){

      CkArrayIndexMax *elems = new CkArrayIndexMax[nstates];

      CkArrayIndex2D idx(0, reductionPlaneNum);  // plane# = this plane#
      for (int j = 0; j < nstates; j++) {
	idx.index[0] = j;
#ifdef USE_TOPOMAP
	idx.index[1] = red_pl[j];
#else
	idx.index[1] = calcReductionPlaneNum(j);
#endif
	elems[j] = idx;
      }//endfor

      particlePlaneENLProxy = 
	CProxySection_CP_State_ParticlePlane::ckNew(particlePlaneProxy.ckGetArrayID(),
						    elems, nstates);
      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
      particlePlaneENLProxy.ckDelegate(mcastGrp);
      EnlCookieMsg *emsg= new EnlCookieMsg;
      mcastGrp->setSection(particlePlaneENLProxy);
      particlePlaneENLProxy.setEnlCookie(emsg);
      delete [] elems;

  }//endif
#ifdef USE_TOPOMAP
  delete [] red_pl;
#endif

//============================================================================
// No load balancing and no atsyncing either!

  setMigratable(false);
  usesAtSync           = CmiFalse;

//============================================================================
// report your status to main

  int constructed=1;
  contribute(sizeof(int), &constructed, CkReduction::sum_int, 
	     CkCallback(CkIndex_main::doneInit(NULL),mainProxy));

//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//  InitKVectors creates local k vectors for this chare and mallocs some memory.
/// It is invoked from CP_State_GSpacePlane::initGSpace()
//  It could be done in the constructor but this is really fine.
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::initKVectors(GStateSlab *gss){
//============================================================================

    numLines        = gss->numLines;
    numFullNL       = gss->numFullNL;
    ees_nonlocal    = gss->ees_nonlocal;
    gSpaceNumPoints = gss->numPoints;

    myForces =  (complex *)fftw_malloc(gSpaceNumPoints*sizeof(complex)); // forces on coefs
    memset(myForces, 0, gSpaceNumPoints * sizeof(complex));

    if(ees_nonlocal==1){
       dyp_re   = (double *)fftw_malloc(gSpaceNumPoints*sizeof(double));
       dyp_im   = (double *)fftw_malloc(gSpaceNumPoints*sizeof(double));
       projPsiG = (complex *)fftw_malloc(numFullNL *sizeof(complex)); // fft projector
       memset(projPsiG, 0, numFullNL * sizeof(complex));
    }//endif

    // This occurs AFTER GSP has registered
    registrationFlag = 0;
    if(ees_nonlocal==1){
       eesCache *eesData = eesCacheProxy.ckLocalBranch ();
       int ncoef = gSpaceNumPoints;
       int *k_x  = eesData->GspData[myChareG].ka;
       int *k_y  = eesData->GspData[myChareG].kb;
       int *k_z  = eesData->GspData[myChareG].kc;
       int mycoef= eesData->GspData[myChareG].ncoef;
       if(eesData->allowedGspChares[myChareG]==0 || mycoef != ncoef){
         CkPrintf("Plane %d of state %d toasy %d %d\n",myChareG,thisIndex.x,mycoef,ncoef);
         CkExit();
       }//endif
       eesData->registerCacheGPP(thisIndex.y,ncoef,k_x,k_y,k_z);
       int i=1;
       CkCallback cb(CkIndex_CP_State_ParticlePlane::registrationDone(NULL),
 		   particlePlaneProxy);
       contribute(sizeof(int),&i,CkReduction::sum_int,cb);
    }//endif

//--------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
// The destructor : Carefully about what it malloced and what is not
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_State_ParticlePlane::~CP_State_ParticlePlane(){

  fftw_free(myForces);

  if(ees_nonlocal==1){
     fftw_free(projPsiG); 
     fftw_free(dyp_re);
     fftw_free(dyp_im);
  }//endif

  if(ees_nonlocal==0){
    fftw_free(zmatrix); //FOO-BAR should i delete this?            
    fftw_free(zmatrix_fx); 
    fftw_free(zmatrix_fy); 
    fftw_free(zmatrix_fz); 
    delete [] count;
    delete [] haveSFAtmGrp;
  }//endif

}
//============================================================================


//============================================================================
// pup for migration
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::pup(PUP::er &p){
	ArrayElement2D::pup(p);
	p|registrationFlag;
        p|myChareG;
	p|iteration;
        p|iterNL;
        p|numNLiter;
        p|ees_nonlocal;
        p|ngridaNL;
        p|ngridbNL;
        p|ngridcNL;
	p|numFullNL;
	p|numLines;          // NL has same number as states but a copy is nice
	p|gSpaceNumPoints;   // NL has same number as states but a copy is nice
	p|doneGettingForces;
	p|sendDone;
	p|numSfGrps;
        p|natm_nl;
	p|natm_nl_grp_max;
        p|realPP_proxy;
        p|countNLIFFT;
        p|zsize;
	if (p.isUnpacking()) {
 	        myForces = (complex *)fftw_malloc(gSpaceNumPoints*sizeof(complex));
                if(ees_nonlocal==1){
                  dyp_re   = (double *)fftw_malloc(gSpaceNumPoints*sizeof(double));
                  dyp_im   = (double *)fftw_malloc(gSpaceNumPoints*sizeof(double));
                  projPsiG = (complex *)fftw_malloc(numFullNL*sizeof(complex)); // fft proj
		}//endif
                if(ees_nonlocal==0){
  		  zmatrixSum    = (complex *)fftw_malloc(zsize*sizeof(complex));
                  zmatrixSum_fx = (complex *)fftw_malloc(zsize*sizeof(complex));
                  zmatrixSum_fy = (complex *)fftw_malloc(zsize*sizeof(complex));
                  zmatrixSum_fz = (complex *)fftw_malloc(zsize*sizeof(complex));
		  zmatrix       = (complex *)fftw_malloc(zsize*sizeof(complex));
                  zmatrix_fx    = (complex *)fftw_malloc(zsize*sizeof(complex));
                  zmatrix_fy    = (complex *)fftw_malloc(zsize*sizeof(complex));
                  zmatrix_fz    = (complex *)fftw_malloc(zsize*sizeof(complex));
  		  count         = new int[numSfGrps];
		  haveSFAtmGrp  = new int[numSfGrps];
		}//endif
	}//endif
	p((char*)myForces,gSpaceNumPoints*sizeof(complex));
        if(ees_nonlocal==1){
  	  p((char*)projPsiG,numFullNL*sizeof(complex));
	  p(dyp_re,gSpaceNumPoints);
	  p(dyp_im,gSpaceNumPoints);
	}//endif
        if(ees_nonlocal==0){
  	  p(count,numSfGrps);
	  p(haveSFAtmGrp,numSfGrps);
          p((char*)zmatrixSum,zsize *sizeof(complex));
  	  p((char*)zmatrixSum_fx,zsize *sizeof(complex));
	  p((char*)zmatrixSum_fy,zsize *sizeof(complex));
	  p((char*)zmatrixSum_fz,zsize *sizeof(complex));
	  p((char*)zmatrix,zsize *sizeof(complex));
	  p((char*)zmatrix_fx,zsize *sizeof(complex));
	  p((char*)zmatrix_fy,zsize *sizeof(complex));
	  p((char*)zmatrix_fz,zsize *sizeof(complex));
	}//endif
	p|sizeX;
	p|sizeY;
	p|sizeZ;
	p|gSpacePlanesPerChare;
	p|energy_count;
	p|totalEnergy;
	p|enl;
	p|doneEnl;
	p|doneForces;
        p|enl_total;
        p|reductionPlaneNum;
	p|enlCookie;
	p|particlePlaneENLProxy;
	// "gspace" is not pup'ed since it is always assigned to
}
//============================================================================


//============================================================================
/* The N^3 routine computeZ is triggered by the arrival of the structure factor
    for each atom group in the local sfcache. */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::computeZ(PPDummyMsg *m){
//============================================================================    
#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
  if(thisIndex.x==0)
   CkPrintf("HI, I am gPP %d %d in computeZ \n",thisIndex.x,thisIndex.y);
#endif

    if(ees_nonlocal==1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Homes, ees nonlocal is on. You can't call computeZ\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    doneGettingForces = false;    //if you are here, you ain't done.
                                  //It is also reset in sendLambda 
                                  //which is when the results of the previous
                                  //iteration are no longer needed.

    CP_State_GSpacePlane *gsp = gSpacePlaneProxy(thisIndex.x, thisIndex.y).ckLocal();
    GStateSlab *gss           = &(gsp->gs);
    int atmIndex    = m->atmGrp;
    int sfindex     = m->sfindex;
    delete m;
   
//============================================================================    
// If you got here and you have the latest psi, correct iteration etc.
// This ensures you have reduced psi for dynamics which you really really need!
// OK, so the upshot is, get your part of the Zmatrix and then send to the reduction

    if(gsp->acceptedPsi && gsp->doneNewIter){
      // get ks
      eesCache *eesData        = eesCacheProxy.ckLocalBranch();
      int *k_x = eesData->GspData[myChareG].ka;
      int *k_y = eesData->GspData[myChareG].kb;
      int *k_z = eesData->GspData[myChareG].kc;

      // get sf
      StructFactCache *sfcache = sfCacheProxy.ckLocalBranch();
      complex *structureFactor,*structureFactor_fx,*structureFactor_fy,*structureFactor_fz;
      sfcache->getStructFact(thisIndex.y, atmIndex, &structureFactor, 
			     &structureFactor_fx, &structureFactor_fy, &structureFactor_fz);
      if(zmatrix==NULL){
	zmatrix    = (complex *)fftw_malloc(zsize*sizeof(complex));
	zmatrix_fx = (complex *)fftw_malloc(zsize*sizeof(complex));
	zmatrix_fy = (complex *)fftw_malloc(zsize*sizeof(complex));
	zmatrix_fz = (complex *)fftw_malloc(zsize*sizeof(complex));
	memset(zmatrix, 0, sizeof(complex)*zsize);
	memset(zmatrix_fx, 0, sizeof(complex)*zsize);
	memset(zmatrix_fy, 0, sizeof(complex)*zsize);
	memset(zmatrix_fz, 0, sizeof(complex)*zsize);
      }//endif

      // get zmat
      int mydoublePack = config.doublePack;
      int zoffset      = natm_nl_grp_max * atmIndex;
      CPNONLOCAL::CP_enl_matrix_calc(gSpaceNumPoints,gss->packedPlaneData, k_x,k_y,k_z, 
	     structureFactor,structureFactor_fx,structureFactor_fy,
	     structureFactor_fz,&zmatrix[zoffset],&zmatrix_fx[zoffset],&zmatrix_fy[zoffset],
	     &zmatrix_fz[zoffset],thisIndex.x,mydoublePack,numSfGrps,atmIndex);

      // reduce zmatrices over the planes of each state
      thisProxy(thisIndex.x, reductionPlaneNum).reduceZ(natm_nl_grp_max, atmIndex, 
	&zmatrix[zoffset],&zmatrix_fx[zoffset],&zmatrix_fy[zoffset],&zmatrix_fz[zoffset]);

      haveSFAtmGrp[atmIndex]=-1; //this one is done

//============================================================================    
//   If you to the `else', you don't have the latest psi (impossible now).
//   Set the flag so Psi can kick off computeZ
    }else{
	haveSFAtmGrp[atmIndex]=sfindex;
    }//endif

//============================================================================    
// Now, Hang out until reduced zmatrix is sent back.

//----------------------------------------------------------------------------
    }//end routine
//============================================================================


//============================================================================
// ReduceZ reduces Zmat over gspace for a particular state for N^3 method
//============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::reduceZ(int size, int atmIndex, complex *zmatrix_,
                   complex *zmatrix_fx_,complex *zmatrix_fy_,complex *zmatrix_fz_){
//============================================================================
// This method is invoked only CP_State_ParticlePlane(*,reductionPlaneNum)
//  CkPrintf("PP [%d %d] reduceZ\n",thisIndex.x, thisIndex.y);

#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
  if(thisIndex.x==0)
   CkPrintf("HI, I am gPP %d %d in redZ \n",thisIndex.x,thisIndex.y);
#endif

  if(ees_nonlocal==1){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Homes, ees nonlocal is on. You can't call reduceZ\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  count[atmIndex]++;
  int i,j;
  int zoffset= natm_nl_grp_max * atmIndex;
  if (zmatrixSum==NULL) {
      zmatrixSum = (complex *)fftw_malloc(zsize*sizeof(complex));
      memset(zmatrixSum, 0, zsize * sizeof(complex));
  }//endif

  if (zmatrixSum_fx==NULL) {
      zmatrixSum_fx = (complex *)fftw_malloc(zsize*sizeof(complex));
      memset(zmatrixSum_fx, 0, zsize * sizeof(complex));
  }//endif

  if (zmatrixSum_fy==NULL) {
      zmatrixSum_fy = (complex *)fftw_malloc(zsize*sizeof(complex));
      memset(zmatrixSum_fy, 0, zsize * sizeof(complex));
  }//endif

  if (zmatrixSum_fz==NULL) {
      zmatrixSum_fz = (complex *)fftw_malloc(zsize*sizeof(complex));
      memset(zmatrixSum_fz, 0, zsize * sizeof(complex));
  }//endif
  
  for (i = 0,j=zoffset; i < size; i++,j++){
    zmatrixSum[j]    += zmatrix_[i];
    zmatrixSum_fx[j] += zmatrix_fx_[i];
    zmatrixSum_fy[j] += zmatrix_fy_[i];
    zmatrixSum_fz[j] += zmatrix_fz_[i];
  }//endif

//============================================================================
// Now we have contributions from all the other gspace planes in this state. 
// Calculate the energies and atoms forces. Normalize ZmatrixSum

  if( ((count[atmIndex] == (sizeX/gSpacePlanesPerChare)/2 - 1)&&(!config.doublePack)) || 
      ((count[atmIndex] == nchareG)&&(config.doublePack)) ){ 

   //---------------------------------------------------------------------
   // accounting
    count[atmIndex] = 0;
    if(doneEnl==0){enl=0.0;}
    doneEnl++;

   //-----------------------------------------------------------------------
   // energy and atom forces
    AtomsGrp *ag         = atomsGrpProxy.ckLocalBranch();
    FastAtoms *fastAtoms = &(ag->fastAtoms);
    int mydoublePack = config.doublePack;
    double myenl     = 0.0;
    CPNONLOCAL::CP_enl_atm_forc_calc(numSfGrps,atmIndex,fastAtoms,&zmatrixSum[zoffset],
               &zmatrixSum_fx[zoffset],&zmatrixSum_fy[zoffset],&zmatrixSum_fz[zoffset],
               &myenl,mydoublePack);
    enl+=myenl;

#ifdef _CP_DEBUG_NLMAT_
    for (i = 0; i < zsize; i++){
      CkPrintf("Non-local matrix[%d][%d]={%g, %g}\n",i,thisIndex.x,
	       zmatrix_[i].re,zmatrix_[i].im );}
#endif

   //---------------------------------------------------------------
   // send the reduced zmatrix to the particleplanes of this state and 
   // compute the non-local forces on the coefs
    for (i = 0; i < nchareG; i ++){
      thisProxy(thisIndex.x, i).getForces(size,atmIndex, &zmatrixSum[zoffset]);
    }//endfor

   //---------------------------------------------------------------
   // If completely done with Z stuff send out our ENL
    if(thisIndex.y==reductionPlaneNum && doneEnl==numSfGrps){
       CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
       CkCallback cb=CkCallback(printEnl, NULL);
       mcastGrp->contribute(sizeof(double),(void*) &enl, 
			   CkReduction::sum_double, enlCookie, cb);
       doneEnl=0;
       enl=0.0;
       bzero(zmatrixSum, zsize * sizeof(complex));
       bzero(zmatrixSum_fx, zsize * sizeof(complex));
       bzero(zmatrixSum_fy, zsize * sizeof(complex));
       bzero(zmatrixSum_fz, zsize * sizeof(complex));
    }//endif

  }/*endif*/

//---------------------------------------------------------------------------
  }//end routine
//============================================================================



//============================================================================
// Compute psi forces for the N^3 method
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::getForces(int zmatSize, int atmIndex, 
                                       complex *zmatrixSum_loc)
//============================================================================
   {//begin routine
//============================================================================
//  Warninging and unpacking

#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
  if(thisIndex.x==0)
   CkPrintf("HI, I am gPP %d %d in getforces \n",thisIndex.x,thisIndex.y);
#endif
  if(ees_nonlocal==1){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Homes, ees nonlocal is on. You can't call getForces\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

  doneForces++;

//============================================================================
// Compute the non-local forces on the coefficients using PINY physics

  CP_State_GSpacePlane *gsp = gSpacePlaneProxy(thisIndex.x, thisIndex.y).ckLocal();
  GStateSlab *gss           = &(gsp->gs);
  eesCache *eesData         = eesCacheProxy.ckLocalBranch();
  StructFactCache *sfcache  = sfCacheProxy.ckLocalBranch();

  int mydoublePack = config.doublePack;
  int state_ind    = gss->istate_ind;                
  int *k_x         = eesData->GspData[myChareG].ka;
  int *k_y         = eesData->GspData[myChareG].kb;
  int *k_z         = eesData->GspData[myChareG].kc;
  complex *structureFactor,*structureFactor_fx,*structureFactor_fy,*structureFactor_fz;

  sfcache->getStructFact(thisIndex.y, atmIndex, &structureFactor, &structureFactor_fx, 
                         &structureFactor_fy, &structureFactor_fz);

  if (gSpaceNumPoints != 0){
    CPNONLOCAL::CP_enl_force_calc(zmatrixSum_loc,gSpaceNumPoints, k_x, k_y, k_z,
                                  structureFactor,myForces,state_ind,mydoublePack,
                                  numSfGrps, atmIndex);
  }//endif

//============================================================================
// Forces have been computed, let g-space plane know. This ends the N^3 method

  if(doneForces==numSfGrps){
     doneGettingForces = true;
     doneForces        = 0;
     gsp->acceptNLForces();
  }//endif : doneforces

//----------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
// Section reduction cookie N^3 method
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::setEnlCookie(EnlCookieMsg *m){

  CkGetSectionInfo(enlCookie,m);

#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
    CkPrintf("PP [%d %d] get enl cookie\n",thisIndex.x, thisIndex.y); 
#endif

}
//============================================================================


//============================================================================
/**
 * Spread the reduction plane numbers around to minimize map collisions for
 * the N^3 method
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
int CP_State_ParticlePlane::calcReductionPlaneNum(int state){
//============================================================================
  
  int nstatemax=nstates-1;
  int ncharemax=nchareG-1;
  int planeNum= state %ncharemax;
  if(planeNum<0){
     CkPrintf(" PP [%d %d] calc nstatemax %d ncharemax %d state %d planenum %d\n",
              thisIndex.x, thisIndex.y,nstatemax, ncharemax, state, planeNum); 
     CkExit();
  }//endif
  return planeNum;
}
//============================================================================


//============================================================================
// Resume from Sync which is a general routine
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::ResumeFromSync(){
  if(thisIndex.x==0 && thisIndex.y==reductionPlaneNum){
      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();       
      mcastGrp->resetSection(particlePlaneENLProxy);
  }//endif
}// end routine
//============================================================================



//==============================================================================
// The entry point to the Euler Exponential Spline non-local method.
// It is invoked serially by its friend, stateGspacePlane so its not an entry.
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_ParticlePlane::startNLEes(int iteration_in){
//==============================================================================
//  Increment counters, Check for inconsistancies

  iterNL++;  if(iterNL==1){iteration++;}

  doneGettingForces = false;   // If you are here, you are not done.
                               // The flag is also flipped in sendlambda
                               // because avoiding race conditions es muy importante

  if(ees_nonlocal==0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Yo dawg, ees nonlocal is off. You can't call startNLEes\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  if(iteration!=iteration_in){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Mismatch in time step between Gstate and Gpart %d %d for chare %d %d\n",
  	      iteration_in,iteration,thisIndex.y,thisIndex.x);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  if(iterNL>numNLiter){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Too many non-local iterations encountered %d %d for Gchare %d %d\n",
  	      iterNL,numNLiter,thisIndex.y,thisIndex.x);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  if(iterNL>1 && registrationFlag==0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Why am I not registered! PP chare %d %d\n",thisIndex.x,thisIndex.y);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
  }//endif

//==============================================================================
// Lets boogy : provided we have registered and we are not trying to do some
//              bogus iteration

#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
  if(thisIndex.x==0)
   CkPrintf("HI, I am gPP %d %d in startNLees : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif

  if(registrationFlag==1 && iteration<=config.maxIter){createNLEesFFTdata();}

  // Bogus iteration :
  //   Spin my wheels while atoms integrate, cleanExit and Energy reducts
  //   and ortho does who knows what.
  if(iteration>config.maxIter){
     CP_State_GSpacePlane *gsp = gSpacePlaneProxy(thisIndex.x,thisIndex.y).ckLocal();
     iterNL = 0;
     doneGettingForces = true; // False is flipped in cp_state_gspace_plane once
     gsp->acceptNLForcesEes(); // Let the lads know you are done
  }//endif

//-----------------------------------------------------------------------------
   }//end routine
//==============================================================================


//============================================================================
// In gspace, create the projector, projPsi
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::createNLEesFFTdata(){
//============================================================================
// Local pointers and warnings 

  if(ees_nonlocal==0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Yo dawg, ees nonlocal is off. You can't createNLEesFFTdata\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

   CP_State_GSpacePlane *gsp = gSpacePlaneProxy(thisIndex.x, thisIndex.y).ckLocal();
   GStateSlab *gss           = &(gsp->gs);
   eesCache *eesData         = eesCacheProxy.ckLocalBranch ();

//============================================================================
// Setup the projector and then launch the fft

   double *d_re   = eesData->GppData[myChareG].b_re;    // precomputed stuff
   double *d_im   = eesData->GppData[myChareG].b_im;
   int *ind_gspl  = eesData->GppData[myChareG].ind_gspl;
   double *h_gspl = eesData->GppData[myChareG].h_gspl;
   int *k_x       = eesData->GspData[myChareG].ka;
   int *k_y       = eesData->GspData[myChareG].kb;
   int *k_z       = eesData->GspData[myChareG].kc;

   int ncoef      = gSpaceNumPoints;
   complex *psi   = gss->packedPlaneData;
   int ihave_g0   = gss->ihave_g000;      // I have the special pt, g=0!
   int ind_g0     = gss->ind_g000;        // This is the index of   g=0 !

#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
   if(thisIndex.x==0)
    CkPrintf("HI, I am gPP %d %d in createNLEes : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif
    CPNONLOCAL::eesProjGchare(ncoef,psi,k_x,k_y,k_z,ihave_g0,ind_g0,iterNL,
                              d_re,d_im,dyp_re,dyp_im,projPsiG,ind_gspl,h_gspl,
                              thisIndex.x,thisIndex.y);

   FFTNLEesFwd();

//---------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
// FFT projPsi(gx,gy,gz) -> projPsi(gx,gy,z) 
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::FFTNLEesFwd(){

  if(ees_nonlocal==0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Yo dawg, ees nonlocal is off. You can't FFTNLEesFwd\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  // local pointers
  CP_State_GSpacePlane *gsp = gSpacePlaneProxy(thisIndex.x, thisIndex.y).ckLocal();
  GStateSlab *gss           = &(gsp->gs);
  FFTcache *fftcache        = fftCacheProxy.ckLocalBranch();

//============================================================================
// Do the FFT and then send it to r-space to complete the fft

#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
  if(thisIndex.x==0)
   CkPrintf("HI, I am gPP %d %d in FFTNLFwd : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif

  fftcache->doNlFftGtoR_Gchare(projPsiG,numFullNL,gSpaceNumPoints,numLines,
                               (gss->numRuns),(gss->runs),ngridcNL);
  sendToEesRPP();

//---------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
// send projPsi(gx,gy,z)  to Rspace chare where FFT is completed
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::sendToEesRPP(){
//============================================================================
// Do a Comlib Dance

  if(ees_nonlocal==0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Yo dawg, ees nonlocal is off. You can't sendtoEesRPP\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
  if(thisIndex.x==0)
   CkPrintf("HI, I am gPP %d %d in sendtoEesRPP : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif

  if (config.useGssInsRealPP){gssPInstance.beginIteration();}

//============================================================================
// Send your (x,y,z) to processors z.

   for(int z=0; z < ngridcNL; z++) {

     NLFFTMsg *msg    = new (numLines,8*sizeof(int)) NLFFTMsg;
     msg->size        = numLines;
     msg->senderIndex = thisIndex.y;  // planenumber
     msg->step        = iterNL;
    
     if(config.prioNLFFTMsg){
        CkSetQueueing(msg, CK_QUEUEING_IFIFO);
        *(int*)CkPriorityPtr(msg) = config.rsNLfftpriority + 
                                   thisIndex.x*ngridcNL + thisIndex.y;
     }//endif

     // beam out all points with same z to chare array index z
     complex *data = msg->data;
     for (int i=0,j=z; i<numLines; i++,j+=ngridcNL){data[i] = projPsiG[j];}
     realPP_proxy(thisIndex.x, z).recvFromEesGPP(msg);  // same state, realspace char[z]
#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif
   }//endfor

//============================================================================
// Turn off commlib
    
  if (config.useGssInsRealPP){gssPInstance.endIteration();}
  sendDone = 1;

//============================================================================
   }//end routine    
//============================================================================


//============================================================================
// Receive the projector back modified to generate psi forces
// A message cannot come back until I have sent!!
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::recvFromEesRPP(GSPPIFFTMsg *msg){
//============================================================================

  if(ees_nonlocal==0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Yo dawg, ees nonlocal is off. You can't recvFromEesRPP\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
    if(thisIndex.x==0)
     CkPrintf("HI, I am gPP %d %d in recvFromEesRPP : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif

//============================================================================
// unpack, check sizes, iteration numbers, and if you could be overwriting memory

  int size             = msg->size;
  int offset           = msg->offset;
  int iterNLnow        = msg->iterNL;
  complex *partlyIFFTd = msg->data;

  CkAssert(numLines == size);
  if(iterNLnow != iterNL || sendDone==0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Iteration mismatch in ParticlePlane %d %d %d\n",iterNL,iterNLnow,sendDone);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

//============================================================================
// Copy out the data and delete the message

  if(countNLIFFT==0) {memset(projPsiG, 0, sizeof(complex)*numFullNL);}

  // z=offset is inner index : collections of z-lines of constant (gx,gy)
  for(int i=0,j=offset; i< numLines; i++,j+=ngridcNL){projPsiG[j] = partlyIFFTd[i];}

  delete msg;
  // receive 1 message from each z-chare    

//============================================================================
// When everything is here, continue the computation

  countNLIFFT++;
  if (countNLIFFT == ngridcNL) {
    countNLIFFT = 0;
#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
    if(thisIndex.x==0)
     CkPrintf("HI, I am gPP %d %d in recvFromEesRPP : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif
    FFTNLEesBck();
    sendDone = 0;
  }//endif

//----------------------------------------------------------------------------
  }//end routine 
//============================================================================


//============================================================================
// Complete the FFT : projPsif(gx,gy,z) -> projPsif(gx,gy,gz)
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::FFTNLEesBck(){
//============================================================================

  if(ees_nonlocal==0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Yo dawg, ees nonlocal is off. You can't FFTEeesBck\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  // local pointers
  CP_State_GSpacePlane *gsp = gSpacePlaneProxy(thisIndex.x, thisIndex.y).ckLocal();
  GStateSlab *gss           = &(gsp->gs);
  FFTcache *fftcache        = fftCacheProxy.ckLocalBranch();

//============================================================================
// Do the FFT and then compute Psi forces

#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
  if(thisIndex.x==0)
   CkPrintf("HI, I am gPP %d %d in FFTNLeesBck : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif
  fftcache->doNlFftRtoG_Gchare(projPsiG,numFullNL,gSpaceNumPoints,numLines,
                               (gss->numRuns),(gss->runs),ngridcNL);
  computeNLEesForces();

//---------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
// Compute the Psi forces : If not done, start a new iteration
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::computeNLEesForces(){
//============================================================================

  CP_State_GSpacePlane *gsp = gSpacePlaneProxy(thisIndex.x, thisIndex.y).ckLocal();
  GStateSlab *gss           = &(gsp->gs);
  eesCache *eesData         = eesCacheProxy.ckLocalBranch ();

  int ncoef       = gSpaceNumPoints;
  int *k_x        = eesData->GspData[myChareG].ka;
  int *k_y        = eesData->GspData[myChareG].kb;
  int *k_z        = eesData->GspData[myChareG].kc;
  int ihave_g0    = gss->ihave_g000;
  int ind_g0      = gss->ind_g000;
  int nkx0        = gss->nkx0;
  complex *fPsiG  = myForces;

#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
  if(thisIndex.x==0)
   CkPrintf("HI, I am gPP %d %d in computeNLeesForc : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif

  CPNONLOCAL::eesPsiForcGspace(ncoef,ihave_g0,ind_g0,nkx0,projPsiG,fPsiG,dyp_re,dyp_im,
                               k_x,k_y,k_z,thisIndex.x,thisIndex.y,iterNL);

//============================================================================

  if(iterNL==numNLiter){
#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
    if(thisIndex.x==0)
     CkPrintf("HI, I am gPP %d %d done! : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif
    iterNL = 0;
    doneGettingForces = true; // False is flipped in cp_state_gspace_plane once
    gsp->acceptNLForcesEes(); // Let the lads know you are done
  }else{
    startNLEes(iteration);    // do another iteration
  }//endif

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//==========================================================================
// Make sure everyone is registered on the 1st time step
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
  void CP_State_ParticlePlane::registrationDone(CkReductionMsg *msg) {
//==========================================================================

  int sum = ((int *)msg->getData())[0];
  delete msg;

#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
  if(thisIndex.x==0)
    CkPrintf("HI, I am gPP %d %d in reg : %d\n",thisIndex.x,thisIndex.y,sum);
#endif

  registrationFlag = 1;
  if(iterNL==1){createNLEesFFTdata();}

  if(iterNL>1){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Homes, registeration must occur before the 1st time step in GPP\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

}
//==========================================================================
