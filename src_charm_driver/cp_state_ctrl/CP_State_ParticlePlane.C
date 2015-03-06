//=========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=========================================================================
/** \file CP_State_ParticlePlane.C
 * @defgroup Particle
 *    @{
 *
 * \brief Compute kinetic energy of the non-interacting electrons and non-local forces based on the particle view of the system, triggered by \ref GSpaceState and overlaps with those computations when possible.
 *
 * The kinetic energy of non-interacting electrons is expressed (in
 * the computer science parlance) as a ``point by point multiply'' and
 * reduction operation.
 * \f$$(\hbar^2/2m_e)\sum_{g_x,g_y,g_z\epsilon
 * |{\bf g}|<g_{cut}}\sum_s f_sg^2 |\Psi(s,g_x,g_y,g_z)|^2$\f$
 *
 * Life-cycle of a CP_State_ParticlePlane:
 *
 * The particle-plane array is a shadow of the g-space planes. This means
 * that UparticlePlaneProxy[thisInstance.proxyOffset](s, p) exists on the same processor as
 * UgSpacePlaneProxy[thisInstance.proxyOffset](s, p). Also,  both chare arrays have the same number of
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

#include "CP_State_ParticlePlane.h"
#include "CP_State_GSpacePlane.h"
#include "CP_State_Plane.h"
#include "fft_slab_ctrl/fftCacheSlab.h"
#include "structure_factor/StructFactorCache.h"
#include "main/cpaimd.h"
#include "main/AtomsCache.h"
#include "main/eesCache.h"
#include "utility/util.h"

#include "ckmulticast.h"
#include "charm++.h"
#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h"
#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cplocal.h"
extern CProxy_ENL_EKE_Collector                  ENLEKECollectorProxy;
//=========================================================================

extern CProxy_main                               mainProxy;
extern CProxy_InstanceController                 instControllerProxy;
extern CkVec <CProxy_CP_State_GSpacePlane>       UgSpacePlaneProxy;
extern CkVec <CProxy_GSpaceDriver>               UgSpaceDriverProxy;
extern CkVec <CProxy_AtomsCache>                   UatomsCacheProxy;
extern CkVec <CProxy_CP_State_ParticlePlane>     UparticlePlaneProxy;
extern CkVec <CProxy_CP_State_RealParticlePlane> UrealParticlePlaneProxy;
extern CkVec <CProxy_StructFactCache>            UsfCacheProxy;
extern CkVec <CProxy_eesCache>                   UeesCacheProxy;
extern CkVec <CProxy_FFTcache>                   UfftCacheProxy;

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
void CP_State_ParticlePlane::CP_State_ParticlePlane::printEnl(CkReductionMsg *msg){

  double d = ((double *)msg->getData())[0];
  delete msg;

  ENLEKECollectorProxy[thisInstance.idxU.z].acceptENL(d);
  //  CkPrintf("{%d} ENL         = %5.8lf\n", thisInstance.proxyOffset, d);   // tell the world
  //  CkAbort("fix CP_StateParticlePlane printEnl to be instance aware");
  UgSpacePlaneProxy[thisInstance.proxyOffset](0,0).computeEnergies(ENERGY_ENL, d);  //store it
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
    int _ees_nonlocal, UberCollection _instance) :
  sizeX(x), sizeY(y), sizeZ(z), ngridaNL(xNL), ngridbNL(yNL), ngridcNL(zNL),
  ees_nonlocal(_ees_nonlocal), nstates(_nstates), nchareG(_nchareG),
  Gstates_per_pe(_Gstates_per_pe), numNLiter(_numNLiter), 
  gSpacePlanesPerChare(_gSpacePPC), thisInstance(_instance)
  //============================================================================
{//begin routine
  //============================================================================

  myChareG             = thisIndex.y;
  ibead_ind            = thisInstance.idxU.x;
  kpoint_ind           = thisInstance.idxU.y;
  itemper_ind          = thisInstance.idxU.z;
  iteration            = 0;
  iterNL               = 0;
  countNLIFFT          = 0;
  doneEnl              = 0;
  doneForces           = 0;
  enl                  = 0.0;
  energy_count         = 0;
  totalEnergy          = 0.0;
  sendDone             = 0;
  istate_ind           = thisIndex.x;

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
  // No load balancing and no atsyncing either!

  setMigratable(false);
  usesAtSync           = false;

  //============================================================================
  // report your status to main

  int constructed=1;
  contribute(sizeof(int), &constructed, CkReduction::sum_int, 
      CkCallback(CkIndex_InstanceController::doneInit(NULL),CkArrayIndex1D(thisInstance.proxyOffset),instControllerProxy), thisInstance.proxyOffset);
#ifdef _CP_GS_DEBUG_COMPARE_VKS_
  savedprojpsiBf=NULL;
  savedprojpsiBfsend=NULL;
  savedprojpsiGBf=NULL;
#endif
  //---------------------------------------------------------------------------
}//end routine
//============================================================================


/** It is invoked from CP_State_GSpacePlane::initGSpace(). It also takes care of 
 * the array section creation which needs to wait for the doneInsertion phase
 * @todo: Made this public to remove friendship. Check if this can go back to being private.
 */
void CP_State_ParticlePlane::initKVectors()
{
  //============================================================================
  // Comlib to talk to realPP : used when ees_enl_on==1

  /// realPP proxies are created only if ees_nonlocal is turned on
  if(ees_nonlocal == 1)
  {
    realPP_proxy = UrealParticlePlaneProxy[thisInstance.proxyOffset];
#ifdef USE_COMLIB
    if (config.useGssInsRealPP)
      ComlibAssociateProxy(gssPInstance,realPP_proxy);
#endif
  }

  //============================================================================
  //  Register with the SFCache

  if(ees_nonlocal==0){
#ifdef _CP_DEBUG_SF_CACHE_
    CkPrintf("PP [%d,%d] has %d numSfGrps\n",thisIndex.x, thisIndex.y, numSfGrps);
#endif
    StructFactCache *sfcache = UsfCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
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



  int *red_pl = new int[nstates];
  int numProcs=CkNumPes();
  int *usedProc= new int[numProcs];
  memset(usedProc,0,sizeof(int)*numProcs);
  int charperpe=nstates/numProcs;
  if(nstates%numProcs!=0)  charperpe++;
  if(charperpe<1) charperpe=1;
  for(int state=0; state<nstates;state++){
    int plane=nchareG-1;
    while(plane>=0)
    {
      bool used=false;
      int thisstateplaneproc = GSImaptable[thisInstance.proxyOffset].get(state, plane)%numProcs;
      if(usedProc[thisstateplaneproc]>=charperpe);
      {
        used=true;
      }
      if(!used || plane==0)
      {
        red_pl[state]=plane;
        (usedProc[thisstateplaneproc])++;
        plane=-1;
      }
      plane--;
    }
  }
  for(int state=0; state<nstates;state++){
    int plane=0;
    while(plane<nchareG)
    {
      bool used=false;
      int thisstateplaneproc = GSImaptable[thisInstance.proxyOffset].get(state, plane)%CkNumPes();
      if(usedProc[thisstateplaneproc]>=charperpe);
      {
        used=true;
      }
      if(!used || (plane+1==nchareG))
      {
        usedProc[thisstateplaneproc]++;
        red_pl[state]=plane;
        plane=nchareG;
      }
      plane++;
    }
  }
  reductionPlaneNum = red_pl[thisIndex.x];
  delete [] usedProc;
  //--------------------------------------------------------------------------
  // If you are a reduction plane, set up your comm with your guys and the
  // the other reduction planes. Store your cookie to eat later.

  if(thisIndex.x==0 && thisIndex.y==reductionPlaneNum){

    CkArrayIndexMax *elems = new CkArrayIndexMax[nstates];

    CkArrayIndex2D idx(0, reductionPlaneNum);  // plane# = this plane#
    for (int j = 0; j < nstates; j++) {
      idx.index[0] = j;
      idx.index[1] = red_pl[j];
      // idx.index[1] = calcReductionPlaneNum(j); ifndef USE_TOPOMAP
      elems[j] = idx;
    }//endfor

    particlePlaneENLProxy = 
      CProxySection_CP_State_ParticlePlane::ckNew(UparticlePlaneProxy[thisInstance.proxyOffset].ckGetArrayID(),
          elems, nstates);
    CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
    particlePlaneENLProxy.ckDelegate(mcastGrp);
    EnlCookieMsg *emsg= new EnlCookieMsg;
    mcastGrp->setSection(particlePlaneENLProxy);
    particlePlaneENLProxy.setEnlCookie(emsg);
    delete [] elems;

  }//endif
  delete [] red_pl;

  GStateSlab *gss = &( UgSpacePlaneProxy[thisInstance.proxyOffset](thisIndex.x, thisIndex.y).ckLocal()->gs );
  numLines        = gss->numLines;
  numFullNL       = gss->numFullNL;
  ees_nonlocal    = gss->ees_nonlocal;
  gSpaceNumPoints = gss->numPoints;

  myForces =  (complex *)fftw_malloc(gSpaceNumPoints*sizeof(complex)); // forces on coefs
  bzero(myForces,gSpaceNumPoints * sizeof(complex));

  if(ees_nonlocal==1){
    dyp_re   = (double *)fftw_malloc(gSpaceNumPoints*sizeof(double));
    dyp_im   = (double *)fftw_malloc(gSpaceNumPoints*sizeof(double));
    projPsiG = (complex *)fftw_malloc(numFullNL *sizeof(complex)); // fft projector
    bzero(projPsiG, numFullNL * sizeof(complex));
  }//endif

  // This occurs AFTER GSP has registered
  registrationFlag = 0;
  if(ees_nonlocal==1){
    eesCache *eesData = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
    int ncoef = gSpaceNumPoints;
    int *k_x  = eesData->GspData[myChareG]->ka;
    int *k_y  = eesData->GspData[myChareG]->kb;
    int *k_z  = eesData->GspData[myChareG]->kc;
    int mycoef= eesData->GspData[myChareG]->ncoef;
    if(eesData->allowedGspChares[myChareG]==0 || mycoef != ncoef){
      CkPrintf("Plane %d of state %d toasy %d %d\n",myChareG,thisIndex.x,mycoef,ncoef);
      CkExit();
    }//endif
    eesData->registerCacheGPP(thisIndex.y,ncoef,k_x,k_y,k_z);
    int i=1;
    CkCallback cb(CkIndex_CP_State_ParticlePlane::registrationDone(NULL),UparticlePlaneProxy[thisInstance.proxyOffset]);
    contribute(sizeof(int),&i,CkReduction::sum_int,cb);
  }//endif
}//end routine


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
  p|istate_ind;
  p|ibead_ind; p|kpoint_ind; p|itemper_ind;
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
/**
 * Pushes all relevant computeZ() calls onto the runtime's queue. Computation for each atomindex  
 * is split into a computeZ call to increase the granularity of the computation and prevent it from 
 * blocking chares on the critical path as it would if it were a monolithic block of work.
 */ 
//============================================================================
void CP_State_ParticlePlane::launchComputeZs(){
  //============================================================================
#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
  CkPrintf("ParticlePlane[%d,%d] launching computeZs()\n",thisIndex.x,thisIndex.y);
#endif
  for(int i=0;i<config.numSfGrps;i++)
  {
    if(haveSFAtmGrp[i]>=0)
    {
      PPDummyMsg *pmsg = new (8*sizeof(int)) PPDummyMsg;
      pmsg->atmGrp  = i;
      pmsg->sfindex = haveSFAtmGrp[i];
      CkSetQueueing(pmsg, CK_QUEUEING_IFIFO);
      *(int*)CkPriorityPtr(pmsg) = config.sfpriority+i+config.numSfGrps;
      //lower than sf and sfcache
      thisProxy(thisIndex.x,thisIndex.y).computeZ(pmsg);
    }
  }
}//end routine
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

  CP_State_GSpacePlane *gsp = UgSpacePlaneProxy[thisInstance.proxyOffset](thisIndex.x, thisIndex.y).ckLocal();
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
    eesCache *eesData        = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
    int *k_x = eesData->GspData[myChareG]->ka;
    int *k_y = eesData->GspData[myChareG]->kb;
    int *k_z = eesData->GspData[myChareG]->kc;

    // get sf
    StructFactCache *sfcache = UsfCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
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

#if CMK_TRACE_ENABLED
    double  StartTime=CmiWallTimer();
#endif    

    CPNONLOCAL::CP_enl_matrix_calc(gSpaceNumPoints,gss->packedPlaneData, k_x,k_y,k_z, 
        structureFactor,structureFactor_fx,structureFactor_fy,
        structureFactor_fz,&zmatrix[zoffset],&zmatrix_fx[zoffset],&zmatrix_fy[zoffset],
        &zmatrix_fz[zoffset],thisIndex.x,mydoublePack,numSfGrps,atmIndex);
#if CMK_TRACE_ENABLED
    traceUserBracketEvent(enlMatrixCalc_, StartTime, CmiWallTimer());    
#endif

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
    AtomsCache *ag         = UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
    FastAtoms *fastAtoms = &(ag->fastAtoms);
    int mydoublePack = config.doublePack;
    double myenl     = 0.0;

#if CMK_TRACE_ENABLED
    double  StartTime=CmiWallTimer();
#endif    

    CPNONLOCAL::CP_enl_atm_forc_calc(numSfGrps,atmIndex,fastAtoms,&zmatrixSum[zoffset],
        &zmatrixSum_fx[zoffset],&zmatrixSum_fy[zoffset],&zmatrixSum_fz[zoffset],
        &myenl,mydoublePack,istate_ind);
#if CMK_TRACE_ENABLED
    traceUserBracketEvent(enlAtmForcCalc_, StartTime, CmiWallTimer());    
#endif

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

      CkCallback cb=CkCallback(CkIndex_CP_State_ParticlePlane::printEnl(NULL),CkArrayIndex2D(0,0), UparticlePlaneProxy[thisInstance.proxyOffset]);
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

  CP_State_GSpacePlane *gsp = UgSpacePlaneProxy[thisInstance.proxyOffset](thisIndex.x, thisIndex.y).ckLocal();
  GStateSlab *gss           = &(gsp->gs);
  eesCache *eesData         = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
  StructFactCache *sfcache  = UsfCacheProxy[thisInstance.proxyOffset].ckLocalBranch();

  int mydoublePack = config.doublePack;
  int state_ind    = gss->istate_ind;                
  int *k_x         = eesData->GspData[myChareG]->ka;
  int *k_y         = eesData->GspData[myChareG]->kb;
  int *k_z         = eesData->GspData[myChareG]->kc;
  complex *structureFactor,*structureFactor_fx,*structureFactor_fy,*structureFactor_fz;

  sfcache->getStructFact(thisIndex.y, atmIndex, &structureFactor, &structureFactor_fx, 
      &structureFactor_fy, &structureFactor_fz);

  if (gSpaceNumPoints != 0){

#if CMK_TRACE_ENABLED
    double  StartTime=CmiWallTimer();
#endif    

    CPNONLOCAL::CP_enl_force_calc(zmatrixSum_loc,gSpaceNumPoints, k_x, k_y, k_z,
        structureFactor,myForces,state_ind,mydoublePack,
        numSfGrps, atmIndex);
#if CMK_TRACE_ENABLED
    traceUserBracketEvent(enlForcCalc_, StartTime, CmiWallTimer());    
#endif

  }//endif

  //============================================================================
  // Forces have been computed, let g-space plane know. This ends the N^3 method

  if(doneForces==numSfGrps){
    doneForces        = 0;
    /// If the nonlocal force computations are barriered, contribute to the reduction
#ifdef BARRIER_CP_PARTICLEPLANE_NONLOCAL
    int nonLocDone = 1;
    contribute(sizeof(int),&nonLocDone,CkReduction::sum_int,
        CkCallback(CkIndex_GSpaceDriver::allDoneNLForces(NULL),UgSpaceDriverProxy[thisInstance.proxyOffset]));
    /// else, directly notify your driver
#else
    UgSpaceDriverProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).doneNLForces(); 
#endif
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
  //============================================================================
  // Do not delete msg. Its a nokeep.
  //============================================================================
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
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_ParticlePlane::lPrioStartNLEes(NLDummyMsg *msg){ 
  int iteration_in=msg->iteration;
  delete msg;
  startNLEes(iteration_in);
}
//==============================================================================


//==============================================================================
// The entry point to the Euler Exponential Spline non-local method.
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_ParticlePlane::startNLEes(int iteration_in){
  //==============================================================================
  //  Increment counters, Check for inconsistancies

  iterNL++;  if(iterNL==1){iteration++;}

  if(natm_nl==0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Duuuuude, natm_nl==0 don't start NL Ees\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

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
    ///TODO: Remove unused var gsp
    CP_State_GSpacePlane *gsp = UgSpacePlaneProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).ckLocal();
    iterNL = 0;
    /// If the nonlocal force computations are barriered, contribute to the reduction
#ifdef BARRIER_CP_PARTICLEPLANE_NONLOCAL
    int nonLocDone = 1;
    contribute(sizeof(int),&nonLocDone,CkReduction::sum_int,
        CkCallback(CkIndex_GSpaceDriver::allDoneNLForces(NULL),UgSpaceDriverProxy[thisInstance.proxyOffset]));
    /// else, directly notify your driver
#else
    UgSpaceDriverProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).doneNLForces(); 
#endif

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

  CP_State_GSpacePlane *gsp = UgSpacePlaneProxy[thisInstance.proxyOffset](thisIndex.x, thisIndex.y).ckLocal();
  GStateSlab *gss           = &(gsp->gs);
  eesCache *eesData         = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
  FFTcache *fftcache        = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();

  //============================================================================
  // Setup the projector and then launch the fft

  double *d_re   = eesData->GppData[myChareG]->b_re;    // precomputed stuff
  double *d_im   = eesData->GppData[myChareG]->b_im;
  int **ind_gspl = eesData->GppData[myChareG]->ind_gspl;
  double **h_gspl= eesData->GppData[myChareG]->h_gspl;
  int *k_x       = eesData->GspData[myChareG]->ka;
  int *k_y       = eesData->GspData[myChareG]->kb;
  int *k_z       = eesData->GspData[myChareG]->kc;

  int ncoef      = gSpaceNumPoints;
  complex *psi   = gss->packedPlaneData;
  int ihave_g0   = gss->ihave_g000;      // I have the special pt, g=0!
  int ind_g0     = gss->ind_g000;        // This is the index of   g=0 !

#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
  if(thisIndex.x==0)
    CkPrintf("HI, I am gPP %d %d in createNLEes : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif

  fftcache->getCacheMem("CP_State_ParticlePlane::createNLEesFFTdata");
  complex *projPsiGTmp = fftcache->tmpData; // store in tmp

#if CMK_TRACE_ENABLED
  double  StartTime=CmiWallTimer();
#endif    

  CPNONLOCAL::eesProjGchare(ncoef,psi,k_x,k_y,k_z,ihave_g0,ind_g0,iterNL,
      d_re,d_im,dyp_re,dyp_im,projPsiGTmp,ind_gspl,h_gspl,
      thisIndex.x,thisIndex.y,kpoint_ind,config.nfreq_cpnonlocal_eesfwd);
#if CMK_TRACE_ENABLED
  traceUserBracketEvent(eesProjG_, StartTime, CmiWallTimer());    
#endif

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
  CP_State_GSpacePlane *gsp = UgSpacePlaneProxy[thisInstance.proxyOffset](thisIndex.x, thisIndex.y).ckLocal();
  GStateSlab *gss           = &(gsp->gs);
  FFTcache *fftcache        = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
  eesCache *eesData         = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();

  //============================================================================
  // Do the FFT and then send it to r-space to complete the fft

#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
  if(thisIndex.x==0)
    CkPrintf("HI, I am gPP %d %d in FFTNLFwd : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif

  RunDescriptor *runs  = eesData->GspData[myChareG]->runs;
  complex *projPsiGTmp = fftcache->tmpData;
#if CMK_TRACE_ENABLED
  double StartTime=CmiWallTimer();
#endif

  fftcache->doNlFFTGtoR_Gchare(projPsiGTmp,projPsiG,numFullNL,gSpaceNumPoints,numLines,
      (gss->numRuns),runs,ngridcNL,myChareG);
  fftcache->freeCacheMem("CP_State_ParticlePlane::FFTNLEesFwd");

#if CMK_TRACE_ENABLED
  traceUserBracketEvent(doNlFFTGtoR_, StartTime, CmiWallTimer());    
#endif

#ifdef _CP_GS_DUMP_VKS_
  dumpMatrix("projPsiGb4send",(double *)projPsiG, 1, gSpaceNumPoints*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
#endif

#ifdef _CP_GS_DEBUG_COMPARE_VKS_
  if(savedprojpsiBfsend==NULL)
  { // load it
    savedprojpsiBfsend= new complex[gSpaceNumPoints];
    loadMatrix("projPsiGb4send",(double *)savedprojpsiBfsend, 1, gSpaceNumPoints*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
  }
  for(int i=0;i<gSpaceNumPoints;i++)
  {
    if(fabs(projPsiG[i].re-savedprojpsiBfsend[i].re)>0.0001)
    {
      fprintf(stderr, "PP [%d,%d] %d element projpsi  %.10g not %.10g\n",thisIndex.x, thisIndex.y,i, projPsiG[i].re, savedprojpsiBfsend[i].re);
    }
    CkAssert(fabs(projPsiG[i].re-savedprojpsiBfsend[i].re)<0.0001);
    CkAssert(fabs(projPsiG[i].im-savedprojpsiBfsend[i].im)<0.0001);
  }
#endif

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

#ifdef USE_COMLIB
#ifdef OLD_COMMLIB
  if (config.useGssInsRealPP){gssPInstance.beginIteration();}
#else
  if (config.useGssInsRealPP){ComlibBegin(realPP_proxy, iterNL);}
#endif
#endif

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
        thisIndex.x + thisIndex.y;
    }//endif

    // beam out all points with same z to chare array index z
    complex *data = msg->data;
    /// @todo: realPP_proxy is initialized in initKVectors only if nloc_on is true. Should we check the same here too?
    for (int i=0,j=z; i<numLines; i++,j+=ngridcNL){data[i] = projPsiG[j];}
    realPP_proxy(thisIndex.x, z).recvFromEesGPP(msg);  // same state, realspace char[z]
  }//endfor

  //============================================================================
  // Turn off commlib

#ifdef USE_COMLIB
#ifdef OLD_COMMLIB
  if (config.useGssInsRealPP){gssPInstance.endIteration();}
#else
  if (config.useGssInsRealPP){ComlibEnd(realPP_proxy, iterNL);}
#endif    
#endif
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

  // No need to zero because every value is set.
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
  CP_State_GSpacePlane *gsp = UgSpacePlaneProxy[thisInstance.proxyOffset](thisIndex.x, thisIndex.y).ckLocal();
  GStateSlab *gss           = &(gsp->gs);
  FFTcache *fftcache        = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
  eesCache *eesData         = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();

  //============================================================================
  // Do the FFT and then compute Psi forces

#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
  if(thisIndex.x==0)
    CkPrintf("HI, I am gPP %d %d in FFTNLeesBck : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif

  RunDescriptor *runs  = eesData->GspData[myChareG]->runs;
  fftcache->getCacheMem("CP_State_ParticlePlane::FFTNLEesBck");
  complex *projPsiGTmp = fftcache->tmpData;

#ifdef _CP_GS_DUMP_VKS_
  dumpMatrix("projPsiGb4",(double *)projPsiG, 1, numFullNL*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
#endif

#ifdef _CP_GS_DEBUG_COMPARE_VKS_
  if(savedprojpsiGBf==NULL)
  { // load it
    savedprojpsiGBf= new complex[numFullNL];
    loadMatrix("projPsiGb4",(double *)savedprojpsiGBf, 1, numFullNL*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
  }
  for(int i=0;i<numFullNL;i++)
  {
    if(fabs(projPsiG[i].re-savedprojpsiGBf[i].re)>0.0001)
    {
      fprintf(stderr, "PP [%d,%d] %d element projpsi  %.10g not %.10g\n",thisIndex.x, thisIndex.y,i, projPsiG[i].re, savedprojpsiGBf[i].re);
    }
    CkAssert(fabs(projPsiG[i].re-savedprojpsiGBf[i].re)<0.0001);
    CkAssert(fabs(projPsiG[i].im-savedprojpsiGBf[i].im)<0.0001);
  }
#endif

#if CMK_TRACE_ENABLED
  double  StartTime=CmiWallTimer();
#endif    

  fftcache->doNlFFTRtoG_Gchare(projPsiG,projPsiGTmp,numFullNL,gSpaceNumPoints,numLines,
      eesData->GspData[myChareG]->numRuns,runs,ngridcNL,myChareG);
#if CMK_TRACE_ENABLED
  traceUserBracketEvent(doNlFFTRtoG_, StartTime, CmiWallTimer());    
#endif

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

  CP_State_GSpacePlane *gsp = UgSpacePlaneProxy[thisInstance.proxyOffset](thisIndex.x, thisIndex.y).ckLocal();
  GStateSlab *gss           = &(gsp->gs);
  eesCache *eesData         = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
  FFTcache *fftcache        = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();

  int ncoef       = gSpaceNumPoints;
  int *k_x        = eesData->GspData[myChareG]->ka;
  int *k_y        = eesData->GspData[myChareG]->kb;
  int *k_z        = eesData->GspData[myChareG]->kc;
  int ihave_g0    = gss->ihave_g000;
  int ind_g0      = gss->ind_g000;
  int nkx0        = gss->nkx0;
  complex *fPsiG  = myForces; // these are zeroed in gspaceplane at the end of iteration.

#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
  if(thisIndex.x==0)
    CkPrintf("HI, I am gPP %d %d in computeNLeesForc : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif
  complex *projPsiGTmp = fftcache->tmpData;

#ifdef _CP_GS_DUMP_VKS_
  dumpMatrix("projPsiGtmpb4",(double *)projPsiGTmp, 1, ncoef*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
#endif

#ifdef _CP_GS_DEBUG_COMPARE_VKS_
  if(savedprojpsiBf==NULL)
  { // load it
    savedprojpsiBf= new complex[ncoef];
    loadMatrix("projPsiGtmpb4",(double *)savedprojpsiBf, 1, ncoef*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
  }
  for(int i=0;i<ncoef;i++)
  {
    if(fabs(projPsiGTmp[i].re-savedprojpsiBf[i].re)>0.0001)
    {
      fprintf(stderr, "PP [%d,%d] %d element projpsi  %.10g not %.10g\n",thisIndex.x, thisIndex.y,i, projPsiGTmp[i].re, savedprojpsiBf[i].re);
    }
    CkAssert(fabs(projPsiGTmp[i].re-savedprojpsiBf[i].re)<0.0001);
    CkAssert(fabs(projPsiGTmp[i].im-savedprojpsiBf[i].im)<0.0001);
  }
#endif

#if CMK_TRACE_ENABLED
  double  StartTime=CmiWallTimer();
#endif    

  CPNONLOCAL::eesPsiForcGspace(ncoef,ihave_g0,ind_g0,nkx0,projPsiGTmp,fPsiG,dyp_re,dyp_im,
      k_x,k_y,k_z,thisIndex.x,thisIndex.y,iterNL,config.nfreq_cpnonlocal_eesbk);
  fftcache->freeCacheMem("CP_State_ParticlePlane::computeNLEesForces");

#if CMK_TRACE_ENABLED
  traceUserBracketEvent(eesPsiForcGspace_, StartTime, CmiWallTimer());    
#endif    

  //============================================================================

  if(iterNL==numNLiter){
#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
    if(thisIndex.x==0)
      CkPrintf("HI, I am gPP %d %d done! : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif
    iterNL = 0;
    /// If the nonlocal force computations are barriered, contribute to the reduction
#ifdef BARRIER_CP_PARTICLEPLANE_NONLOCAL
    int nonLocDone = 1;
    contribute(sizeof(int),&nonLocDone,CkReduction::sum_int,
        CkCallback(CkIndex_GSpaceDriver::allDoneNLForces(NULL),UgSpaceDriverProxy[thisInstance.proxyOffset]));
    /// else, directly notify your driver
#else
    UgSpaceDriverProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).doneNLForces(); 
#endif
    // Gsp zeros myForces so it is ready next time
  }else{
    NLDummyMsg *msg = new(8*sizeof(int)) NLDummyMsg;
    CkSetQueueing(msg, CK_QUEUEING_IFIFO);
    *(int*)CkPriorityPtr(msg) = config.sfpriority;
    msg->iteration=iteration;
    thisProxy(thisIndex.x, thisIndex.y).lPrioStartNLEes(msg);
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
/*@}*/
#include "gParticlePlane.def.h"

