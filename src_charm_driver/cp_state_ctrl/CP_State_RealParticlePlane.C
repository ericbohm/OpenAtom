//=========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=========================================================================
/** \file CP_State_RealParticlePlane.C
 * Life-cycle of a CP_State_RealParticlePlane:
 *
 * Insert descriptive comment here please
 */ 
//=========================================================================

#include "utility/util.h"  ///@note: Putting this declaration further down seems to screw things. Whats the issue with fft malloc & free declarations?

#include "CP_State_ParticlePlane.h"
#include "CP_State_Plane.h"
#include "structure_factor/StructFactorCache.h"
#include "fft_slab_ctrl/fftCacheSlab.h"
#include "main/groups.h"
#include "main/eesCache.h"
#include "main/cpaimd.h"
#include "ckmulticast.h"
#include "charm++.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cplocal.h"

//=========================================================================
extern CProxy_InstanceController      instControllerProxy;
extern CProxy_main                       mainProxy;
extern CkVec <CProxy_CP_State_GSpacePlane>       UgSpacePlaneProxy;
extern CkVec <CProxy_AtomsGrp>                   UatomsGrpProxy;
extern CProxy_CPcharmParaInfoGrp         scProxy;
extern CkVec <CProxy_CP_State_ParticlePlane>     UparticlePlaneProxy;
extern CkVec <CProxy_StructFactCache>            UsfCacheProxy;
extern CkVec <CProxy_eesCache>                   UeesCacheProxy;
extern CkVec <CProxy_FFTcache>                   UfftCacheProxy;
extern CkVec <CProxy_CP_State_RealParticlePlane> UrealParticlePlaneProxy;
extern CkVec <MapType2>                          RPPImaptable;
extern CkGroupID            mCastGrpId;
extern ComlibInstanceHandle mssPInstance;
extern CkReduction::reducerType sumFastDoubleType;
extern int    nstates;
extern int    nchareG;
extern Config config;

//#define _CP_DEBUG_STATE_RPP_VERBOSE_

//=========================================================================


//============================================================================
// Energy reduction client for ees method!
//============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealParticlePlane::printEnlR(CkReductionMsg *m){

  //unpack
  double d = ((double *)m->getData())[0];
  itimeRed = m->getUserFlag();
  delete m;
  //output and save the data
  CkPrintf("{%d} ENL(EES)    = %5.8lf\n", thisInstance.proxyOffset,d);
  UgSpacePlaneProxy[thisInstance.proxyOffset](0,0).computeEnergies(ENERGY_ENL, d);  
}
//============================================================================


//============================================================================
// Energy reduction client for ees method!
//============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealParticlePlane::printEnlRSimp(double cp_enl_loc,int index,int itime_in){

  // upack
  countEnl++;
  if(countEnl==1){itimeRed=itime_in;}
  cp_enlTot += cp_enl_loc;

  if(itimeRed!=itime_in){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Yo dawg, ees nonlocal is racin' into the enl reduction\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
  CkPrintf("Enl contrib arrived from %d : %g %d %d %d\n",index,cp_enl_loc,countEnl,
            itime_in,itimeRed);
#endif

  // output and save the data
  if(countEnl==nstates){
    CkPrintf("{%d} ENL(EES)    = %5.8lf\n", thisInstance.proxyOffset, cp_enlTot);
    UgSpacePlaneProxy[thisInstance.proxyOffset](0,0).computeEnergies(ENERGY_ENL,cp_enlTot);  
    countEnl  = 0;
    cp_enlTot = 0.0;
  }//endif

}
//============================================================================



//============================================================================
// The constructor  : Called only once
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_State_RealParticlePlane::CP_State_RealParticlePlane(
                             int ngridA_in, int ngridB_in, int ngridC_in,
                             int numIterNL_in, int zmatSizeMax_in,
                             int Rstates_per_pe_in,int nChareG_in,
                             int ees_nonlocal_in, UberCollection _instance):
  thisInstance(_instance)
//============================================================================
  {//begin routine
//============================================================================
// The non-local guts variables

  ngridA         = ngridA_in;            // grid size along a
  ngridB         = ngridB_in;            // grid size along b
  ngridC         = ngridC_in;            // grid size along c

  nChareR        = ngridC;               // Real Space chares=# C-planes
  nChareG        = nChareG_in;           // G  Space chares
  Rstates_per_pe = Rstates_per_pe_in;    // Real Space topomap
  myPlane        = thisIndex.y;          // Real space plane number
  planeSize      = (ngridA+2)*ngridB;    // expanded plane size for FFTing
  planeSizeT     = ngridA*ngridB;        // true plane size 
  csize          = (ngridA+2)*ngridB/2;  // complex variable size
  CkAssert(csize*2==planeSize);
  numIterNl      = numIterNL_in;         // # of non-local iterations per time step
  zmatSizeMax    = zmatSizeMax_in;       // zmatrix size
  ees_nonlocal   = ees_nonlocal_in;

  rhoRTime         = 0;
  cp_enl           = 0.0;                // non-local energy
  cp_enlTot        = 0.0;
  count            = 0;                  // fft communication counter
  countEnl         = 0;                  // Energy counter
  iterNL           = 0;                  // Nl iteration counter
  itimeRed         = 0;                  // time of reduction
  itime            = 0;                  // time step;
  registrationFlag = 0;
  recvBlock        = 0;
  fftDataDone      = false;
  launchFFT        = false;
  countZ           = 0;      // Zmat communication counter



//============================================================================
// No migration: No atSync load-balancing act

  setMigratable(false);
  usesAtSync = CmiFalse;
#ifdef _CP_GS_DEBUG_COMPARE_VKS_
  savedprojpsiC=NULL;
  savedProjpsiCScr=NULL;
  savedProjpsiRScr=NULL;
  savedzmat=NULL;
  savedmn=NULL;
  saveddmn_x=NULL;
  saveddmn_y=NULL;
  saveddmn_z=NULL;
  savedigrid=NULL;
#endif
//----------------------------------------------------------------------------
  }//end routine
//============================================================================

//============================================================================
// Post construction initialization.
// you can't make sections until the multicast manager is done
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealParticlePlane::init(){

//============================================================================
// Malloc the projector memory, non-local matrix and register with your cache

  if(ees_nonlocal==1){
   //--------
   // malloc
    projPsiC    = (complex*) fftw_malloc(csize*sizeof(complex));
    projPsiR    = reinterpret_cast<double*> (projPsiC);
    zmat        = new double[zmatSizeMax];
    zmatScr     = new double[zmatSizeMax];
   //--------
   // Register
    eesCache *eesData  = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch (); 
    eesData->registerCacheRPP(thisIndex.y);
  }//endif


//============================================================================
// Choose reduction plane reasonably intelligently

  int *red_pl       = new int[nstates];
  //foreach state, try to find a new proc for each state's reduction plane
  int *usedProc= new int[CkNumPes()];
  memset(usedProc,0,sizeof(int)*CkNumPes());
  int charperpe=nstates/CkNumPes();
  if(nstates%CkNumPes()!=0)  charperpe++;
  if(charperpe<1) charperpe=1;
  for(int state=0; state<nstates;state++){
    int plane=0;
    while(plane<nChareR)
      {
        bool used=false;
        int thisstateplaneproc=RPPImaptable[thisInstance.proxyOffset].get(state,plane)%CkNumPes();
	if(usedProc[thisstateplaneproc]>charperpe);
	{
	  used=true;
	}

        if(!used || (plane+1==nChareR))
          {
	    usedProc[thisstateplaneproc]++;
            red_pl[state]=plane;
            plane=nChareR;
          }
        plane++;
      }
  }
  reductionPlaneNum = red_pl[thisIndex.x];
  delete [] usedProc;
	    /* old less reliable method
  int l             = Rstates_per_pe;
  int pl            = (nstates / l);
  int pm            = (CkNumPes() / pl);
  if(pm==0){CkAbort("Choose a larger Gstates_per_pe\n");}
  int planes_per_pe = (nChareR / pm);
  for(int i=0; i<nstates;i++){
    red_pl[i]= (  ((i % Rstates_per_pe)*planes_per_pe)% nChareR);
  }//endif
  reductionPlaneNum = red_pl[thisIndex.x];
	    */
  
  /*for(int i=0; i<nstates;i++){ ifndef USE_TOPOMAP
    red_pl[i] = calcReductionPlaneNum(thisIndex.x);
  }//endif
  reductionPlaneNum = calcReductionPlaneNum(thisIndex.x); */


//============================================================================
// Build section reductions

  //-----------------------------------------------------------
  // The plane section reduction to the reduction plane for zmat
  // Only the root of the reduction does setup dance.
  //-----------------------------------------------------------
  enlSectionComplete=false;
  planeRedSectionComplete=false;
  if(thisIndex.y==reductionPlaneNum){
    rPlaneSectProxy = 
       CProxySection_CP_State_RealParticlePlane::ckNew(thisProxy.ckGetArrayID(),
	  					       thisIndex.x,thisIndex.x,1,
						       0,nChareR-1,1);
    CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
    rPlaneSectProxy.ckDelegate(mcastGrp);
    mcastGrp->setSection(rPlaneSectProxy);
    EnlCookieMsg *emsg= new EnlCookieMsg;
    rPlaneSectProxy.setPlaneRedCookie(emsg);
  }//endif
  else
    { // if you aren't in the reduction plane you are setup
      planeRedSectionComplete=true;
    }

  //-----------------------------------------------------------
  // The reduction over states involving the reduction planes for cp_enl
  // Only the root of the reduction does setup dance, state=0
  //-----------------------------------------------------------

  if(thisIndex.y==reductionPlaneNum && thisIndex.x==0){
      CkArrayIndexMax *elems = new CkArrayIndexMax[nstates];
      CkArrayIndex2D idx(0, reductionPlaneNum);  // cheesy constructor
      for (int j = 0; j < nstates; j++) {
	idx.index[0] = j;  idx.index[1] = red_pl[j];
	elems[j] = idx;
      }//endfor
      rPlaneENLProxy = 
  	   CProxySection_CP_State_RealParticlePlane::ckNew(thisProxy.ckGetArrayID(),
  		  				           elems, nstates);
      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
      rPlaneENLProxy.ckDelegate(mcastGrp);
      mcastGrp->setSection(rPlaneENLProxy);
      EnlCookieMsg *emsg= new EnlCookieMsg;
      rPlaneENLProxy.setEnlCookie(emsg);
      delete [] elems;
  }//endif
  else
    { // if not part of enl, setup complete
      bool inEnl=false;
      for (int j = 0; j < nstates; j++) {
	if(thisIndex.x ==j && thisIndex.y== red_pl[j])
	  inEnl=true;
      }//endfor
      if(!inEnl)
	enlSectionComplete=true;
    }
  initComplete();
  delete [] red_pl;

//============================================================================
// Setup the comlib to talk to GPP

  gPP_proxy = UparticlePlaneProxy[thisInstance.proxyOffset];
#ifdef USE_COMLIB
  if (config.useMssInsGPP){
     ComlibAssociateProxy(mssPInstance,gPP_proxy);
  }//endif
#endif

}

//============================================================================
// All cookies initialized for enl section, only reduction root receives this
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealParticlePlane::enlSectDone(CkReductionMsg *m)
{
  //  CkPrintf("RPP[%d %d] gets enlSectDone\n",thisIndex.x, thisIndex.y);
  delete m;
  enlSectionComplete=true;
  initComplete();
}

//============================================================================
// All cookies initialized for zmat plane reduction
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealParticlePlane::planeRedSectDone(CkReductionMsg *m)
{
  //  CkPrintf("RPP[%d %d] gets rPlaneRedSectDone\n",thisIndex.x, thisIndex.y);
  delete m;
  planeRedSectionComplete=true;
  initComplete();
}

//============================================================================
// Initialization and registration done for this chare
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealParticlePlane::initComplete()
{
  // If the sections are done for enl and zmat, then we are also done
  // registering. The coalesced completion state is reported in to the
  // global initialization phase process.
  if(  enlSectionComplete && planeRedSectionComplete)
    {
      //      CkPrintf("RPP[%d %d] initComplete\n",thisIndex.x, thisIndex.y);
      int constructed=1;
      contribute(sizeof(int), &constructed, CkReduction::sum_int, 
		 CkCallback(CkIndex_InstanceController::doneInit(NULL),
			    CkArrayIndex1D(thisInstance.proxyOffset),
			    instControllerProxy), thisInstance.proxyOffset);
    }
}

//============================================================================
// The destructor : never called directly but I guess migration needs it
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_State_RealParticlePlane::~CP_State_RealParticlePlane(){

  if(ees_nonlocal==1){
    fftw_free(projPsiC);
    delete [] zmat;
    delete [] zmatScr;
  }//endif
  projPsiC    = NULL;
  projPsiR    = NULL;
  zmat        = NULL;
  zmatScr     = NULL;

}
//============================================================================


//============================================================================
// Recv FFT data/psi-projector from CP_StateParticlePlane
// kicking things off for this Chare.
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealParticlePlane::recvFromEesGPP(NLFFTMsg *msg){
//============================================================================
//  Make sure you belong here.

    if(ees_nonlocal==0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Yo dawg, ees nonlocal is off. RPP can't recvFromEesGPP\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

//============================================================================
//  Unpack the message, get some local variables and increment counters
#ifdef _NAN_CHECK_
  for(int i=0;i<msg->size ;i++)
    {
      CkAssert(finite(msg->data[i].re));
      CkAssert(finite(msg->data[i].im));
    }
#endif

    int size               = msg->size; 
    int iterNLNow          = msg->step; 
    int Index              = msg->senderIndex; // which g-space chare sent the data
    complex *partiallyFFTd = msg->data;

    CPcharmParaInfo *sim   = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
    int nchareG            = sim->nchareG;
    int **tranUnpack       = sim->index_tran_upackNL;
    int *nline_per_chareG  = sim->nlines_per_chareG;

    count++;
    // if we are starting a totally new time step, increment the time
    if(iterNL==0 && count==1){itime++; cp_enl=0.0;}

    // we doing a new iteration, refresh the memory, increment iterNL
    // In real space, we must zero out here because not every value
    // is initialized. e.g. some zero elements are not set.
    if(count==1){
      iterNL++; 
      bzero(projPsiR,planeSize*sizeof(double));
    }//endif

//============================================================================
// Consistency Checking

    if (count > nchareG) {
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf(
      "Mismatch in allowed # of NL-gspace chare arrays sending to NLR-chare: %d %d %d %d\n",
                count,nchareG,thisIndex.x,thisIndex.y);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    if(size!=nline_per_chareG[Index]){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Dude, size %d != %d for NL-Rchare %d %d from G-chare %d\n",
               size,nline_per_chareG[Index],thisIndex.y,Index,Index);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    if(iterNLNow!=iterNL || recvBlock==1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Dude, iter count %d != %d for NL-Rchare %d %d %d\n",iterNLNow,iterNL,
                   thisIndex.y,Index,recvBlock);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

//============================================================================
//   Copy out the data set and delete the message

// You have received packed data (x,y) from processor sendIndex
// Every real space chare receives the same x,y indicies.
// For double pack, x=0,1,2,3,4 ...  y= {-K ... K}
// The x increase with processor number. The y are split.
// The rundescriptor contains all we need to unpack the data.
// For doublepack : nffty*run[i][j].x + run[i][j].y
// we store this stuff in the convenient package
// Pictorially a half cylinder is sent which is unpacked into
// a half cube for easy FFTing. Y is the inner index.

    for(int i=0;i< size;i++){projPsiC[tranUnpack[Index][i]] = partiallyFFTd[i];}
    delete msg;

//============================================================================
// When you have collected the full data set, continue

    if(count == nchareG) {
      count     = 0;
      recvBlock = 1;
#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
      if(thisIndex.x==0)
       CkPrintf("HI, I am rPP %d %d in recvfromGpp : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif
      fftDataDone=true;
      if(launchFFT){
	thisProxy(thisIndex.x,thisIndex.y).FFTNLEesFwdR();
        if(rhoRTime!=itime){CkPrintf("Badddd launchFFT.1\n");CkExit();}
      }else{
        if(iterNL!=1){CkPrintf("Badddd launchFFT.2\n");CkExit();}
      }//endif
    }//endif

//----------------------------------------------------------------------------
  }// end routine
//============================================================================


//============================================================================
// Complete the Forward FFT of psi-projector sent from g-space chares
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealParticlePlane::FFTNLEesFwdR(){
  fftDataDone=false;
  int nplane_x = scProxy.ckLocalBranch()->cpcharmParaInfo->nplane_x;
#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
  if(thisIndex.x==0)
   CkPrintf("HI, I am rPP %d %d in FFTNL : %d %d %d %d\n",
     thisIndex.x,thisIndex.y,iterNL,ngridA,ngridB,nplane_x);
#endif
  // This is actually in place, projPsiR and projPsiC use the same location
#if CMK_TRACE_ENABLED
   double  StartTime=CmiWallTimer();
#endif    

  UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch()->doNlFFTGtoR_Rchare(projPsiC,projPsiR,
                                                    nplane_x,ngridA,ngridB,myPlane);
#if CMK_TRACE_ENABLED
   traceUserBracketEvent(doNlFFTGtoR_, StartTime, CmiWallTimer());    
#endif

#ifdef _CP_GS_DUMP_VKS_
    dumpMatrixDouble("projPsiC",(double *)projPsiC, 1, csize*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
#endif

#ifdef _CP_GS_DEBUG_COMPARE_VKS_
    // at this point we're checking that we did the fft correctly
    // which implicitly checks that we collated the input correctly

  if(savedprojpsiC==NULL)
    { // load it
      savedprojpsiC= new complex[csize];
      loadMatrixDouble("projPsiC",(double *)savedprojpsiC, 1, csize*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
    }
  for(int i=0;i<csize;i++)
    {
      if(fabs(projPsiC[i].re-savedprojpsiC[i].re)>0.0001)
	{
	  fprintf(stderr, "RPP [%d,%d] %d element projpsi  %.10g not %.10g\n",thisIndex.x, thisIndex.y,i, projPsiC[i].re, savedprojpsiC[i].re);
	}
      CkAssert(fabs(projPsiC[i].re-savedprojpsiC[i].re)<0.0001);
      CkAssert(fabs(projPsiC[i].im-savedprojpsiC[i].im)<0.0001);
    }
#endif

#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
  if(thisIndex.x==0)
   CkPrintf("HI, I am rPP %d %d in FFTNL.2 : %d %d %d %d\n",
     thisIndex.x,thisIndex.y,iterNL,ngridA,ngridB,nplane_x);
#endif

  if(registrationFlag==1){computeZmatEes();}

  if(registrationFlag==0 && iterNL!=1 || iterNL==0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Dude, iter count %d > 1 for NL-Rchare %d %d and no reg?\n",iterNL,
                   thisIndex.x,thisIndex.y);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
  }//endif

//----------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
// Compute the Zmat elements of this iteration : Spawn the section reduction
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealParticlePlane::computeZmatEes(){

   CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
   int iterNL1          = iterNL-1;           // silly C++ convention
   int *nmem_zmat       = sim->nmem_zmat;     // zmat size now
   int nZmat            = nmem_zmat[iterNL1]; 

//============================================================================ 
// Check out your B-splines from the cache and then compute Zmat

   eesCache *eesData  = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch (); 
   int *allowedRppChares = eesData->allowedRppChares;
   CkAssert(allowedRppChares[myPlane]==1);

   if(iterNL==1){eesData->queryCacheRPP(myPlane,itime,iterNL);}// query once a t-step

   int *plane_index = eesData->RppData[myPlane]->plane_index;
   int **igrid      = eesData->RppData[myPlane]->igrid;
   double **mn      = eesData->RppData[myPlane]->mn;

#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
   if(thisIndex.x==0)
    CkPrintf("HI, I am rPP %d %d in computeZmat : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif

#if CMK_TRACE_ENABLED
   double  StartTime=CmiWallTimer();
#endif    

   CPNONLOCAL::eesZmatRchare(projPsiR,iterNL,zmat,igrid,mn,
                             plane_index,thisIndex.x,myPlane);
#if CMK_TRACE_ENABLED
   traceUserBracketEvent(eesZmatR_, StartTime, CmiWallTimer());    
#endif

//============================================================================
// Launch section reduction : 

#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
   if(thisIndex.x==0)
    CkPrintf("HI, I am rPP %d %d in send to zmat-red : %d %d %d d\n",
         thisIndex.x,thisIndex.y,iterNL,reductionPlaneNum,nZmat,zmatSizeMax);
#endif

#define _FANCY_RED_METHOD_
#ifdef _FANCY_RED_METHOD_
   CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch(); 
   CkCallback cb(CkIndex_CP_State_RealParticlePlane::recvZMatEes(NULL),
                 CkArrayIndex2D(thisIndex.x,reductionPlaneNum),
                 UrealParticlePlaneProxy[thisInstance.proxyOffset].ckGetArrayID());
   mcastGrp->contribute((nZmat*sizeof(double)),zmat,sumFastDoubleType,
                         rPlaneRedCookie,cb, iterNL);
#else
   thisProxy(thisIndex.x,reductionPlaneNum).recvZMatEesSimp(nZmat,zmat,
                                               thisIndex.x,thisIndex.y,iterNL);
#endif


//----------------------------------------------------------------------------
    }//end routine
//============================================================================


//============================================================================
// A chare can be behind by 1 iteration only. This message can arrive before
// this chare has sent its zmat. This chare must therefore recv in zmatScr
//============================================================================
void CP_State_RealParticlePlane::recvZMatEesSimp(int size, double *_zmat,
                                        int state, int index,int iterNL_in){
//============================================================================

   CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
   int iterNL1          = iterNL_in-1;        // silly C++ convention
   int *nmem_zmat       = sim->nmem_zmat;     // zmat size now
   int nZmat            = nmem_zmat[iterNL1];

//============================================================================
// check for errors

   if(size !=  nZmat){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Dude, size != nzmat : %d : %d %d : %d %d\n",
                   iterNL,thisIndex.x,thisIndex.y,state,index);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
   }//endif
   if(iterNL!=iterNL_in && iterNL!=iterNL_in-1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Dude, wrong iterNls : %d %d : %d %d : %d %d\n",
                   iterNL,iterNL_in,thisIndex.x,thisIndex.y,state,index);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
   }//endif
   if(state != thisIndex.x){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Dude, you sent to the wrong state : %d %d  : %d %d : %d %d\n",
                   iterNL,iterNL_in,thisIndex.x,thisIndex.y,state,index);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
   }//endif

//============================================================================
// Recv the reduced zmatrix and chuck the message. Recv into scratch because
// it could be this chare has not yet contributed and you would overwrite
// critical data if you used zmat[]

    countZ++;
    if(countZ==1){bzero(zmatScr,nZmat*sizeof(double));}
    for(int i=0;i<nZmat;i++){zmatScr[i]+=_zmat[i];}

//============================================================================
// Bcast the puppy dog back out to your friends : Bow-Wow-Wow

    if(countZ==nChareR){
      if(iterNL != iterNL){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("Wrong iterNls strange racing: %d %d : %d %d : %d %d\n",
                   iterNL,iterNL_in,thisIndex.x,thisIndex.y,state,index);
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit();
      }//endif
#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
      CkPrintf("HI, I am rPP %d %d in recvZmatSimp blasting off to forces: %d %d\n",
          thisIndex.x,thisIndex.y,iterNL,index);
#endif
      countZ=0;
      for(int i=0;i<nChareR;i++){
	CompAtmForcMsg *rmsg = new (nZmat, 8*sizeof(int)) CompAtmForcMsg;
	if(config.prioNLFFTMsg){
	  CkSetQueueing(rmsg, CK_QUEUEING_IFIFO);
	  *(int*)CkPriorityPtr(rmsg) = config.gsNLfftpriority+thisIndex.x+iterNL_in;
	}//endif
	rmsg->init(nZmat,zmatScr,iterNL_in);
	thisProxy(thisIndex.x,i).computeAtmForcEes(rmsg);
      }//endfor
    }//endif

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
// Reduction client of zmat :  Here everyone must have the same iteration
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealParticlePlane::recvZMatEes(CkReductionMsg *msg){

   CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
   int iterNL_in        = msg->getUserFlag();
   //   int iterNL1          = iterNL_in-1;        // silly C++ convention
   int iterNL1          = iterNL-1;        // silly C++ convention
   int *nmem_zmat       = sim->nmem_zmat;     // zmat size now
   int nZmat            = nmem_zmat[iterNL1];

//============================================================================
// Recv the reduced zmatrix and chuck the message

#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
    if(thisIndex.x==0)
     CkPrintf("HI, I am rPP %d %d in recvZmat : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif
#ifdef _NAN_CHECK_
  for(int i=0;i<msg->getSize()/sizeof(double) ;i++)
    {
      CkAssert(finite(((double*) msg->getData())[i]));
    }
#endif

    CkAssert(msg->getSize() ==  nZmat* sizeof(double));
    double *realValues = (double *) msg->getData(); 
    for(int i=0;i<nZmat;i++){zmat[i]=realValues[i];}

    delete msg;

//============================================================================
// Bcast the puppy dog back out to your friends : Bow-Wow-Wow

#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
    if(thisIndex.x==0)
     CkPrintf("HI, I am rPP %d %d in recvZmat sending to compute: %d\n",
            thisIndex.x,thisIndex.y,iterNL);
#endif
   CompAtmForcMsg *rmsg=new (nZmat, 8*sizeof(int)) CompAtmForcMsg;
   if(config.prioNLFFTMsg){
	CkSetQueueing(rmsg, CK_QUEUEING_IFIFO);
	*(int*)CkPriorityPtr(rmsg) = config.gsNLfftpriority+thisIndex.x+iterNL_in;
   }//endif

   rmsg->init(nZmat,zmat,iterNL);
   CkAssert(rmsg->iterNL>0);
   rPlaneSectProxy.computeAtmForcEes(rmsg);

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
// Use Zmat and ProjPsi to get atmForces, Energy and psiforces
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealParticlePlane::computeAtmForcEes(CompAtmForcMsg *msg)
{
//============================================================================
// Do not delete msg. Its a nokeep.
//============================================================================
//============================================================================
// Unpack the message : Get some sizes

#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
   if(thisIndex.x==0)
    CkPrintf("HI, I am rPP %d %d in compteAtmforc : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif
   int nZmat_in=msg->nZmat;
   double *zmat_loc=msg->zmat;
   int iterNL_in=msg->iterNL;
   CkAssert(msg->iterNL>0);
   CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
   int iterNL1          = iterNL-1;          // silly C++ convention
   int *nmem_zmat       = sim->nmem_zmat;    // zmat memory size
   int nZmat            = nmem_zmat[iterNL1];
   //  int nZmat            = msg->nZmat;
   /*  In the reduction case each iteration arrives together
       So this coordination scheme is redundant */
   if(iterNL_in !=  iterNL || nZmat_in != nZmat){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("RPP [%d,%d] Dude, iteration mismatch : %d %d z %d %d\n",
	       thisIndex.x,thisIndex.y,
	       iterNL,iterNL_in,nZmat_in, nZmat );
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
   }//endif

   CmiMemcpy(zmat,zmat_loc,sizeof(double)*nZmat_in);

//============================================================================
// Check out your B-splines from the cache and then compute energy and forces

   eesCache *eesData   = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch (); 
   FFTcache *fftcache  = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();  

   int *allowedRppChares = eesData->allowedRppChares;
   CkAssert(allowedRppChares[myPlane]==1);

   int *plane_index    = eesData->RppData[myPlane]->plane_index;
   int **igrid         = eesData->RppData[myPlane]->igrid;
   int *nBreakJ        = eesData->RppData[myPlane]->nBreakJ;
   int **sBreakJ       = eesData->RppData[myPlane]->sBreakJ;
   double **dmn_x      = eesData->RppData[myPlane]->dmn_x;
   double **dmn_y      = eesData->RppData[myPlane]->dmn_y;
   double **dmn_z      = eesData->RppData[myPlane]->dmn_z;
   double **mn         = eesData->RppData[myPlane]->mn;
   double *projPsiRScr = fftcache->tmpDataR;
   fftcache->getCacheMem("CP_State_RealParticlePlane::computeAtmForcEes");

   AtomsGrp *ag         = UatomsGrpProxy[thisInstance.proxyOffset].ckLocalBranch();
   FastAtoms *fastAtoms = &(ag->fastAtoms);

#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
   if(thisIndex.x==0)
    CkPrintf("HI, I am rPP %d %d in compteAtmforc : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif

   // projPsiR     comes in with info for atoms
    // projPsiRSsr leaves with info for psiforces to (ngrid_a+2)*ngrid_b elements 
   int n_a, n_b, n_c, n_interp, nAtm;
   CPNONLOCAL::getEesPrms(&n_a,&n_b,&n_c,&n_interp,&nAtm);
   int n_interp2=n_interp*n_interp;
#ifdef _CP_GS_DUMP_VKS_
#ifdef _INSANEO_PARANOID_COMPARE_THAT_EATS_HUGE_MEMORY_
    dumpMatrix2DDouble("mn",mn, nAtm, n_interp2, thisIndex.y,thisIndex.x,thisIndex.x,iterNL,false);    
    dumpMatrix2DDouble("dmn_x",dmn_x, nAtm, n_interp2,thisIndex.y,thisIndex.x,thisIndex.x,iterNL,false);    
    dumpMatrix2DDouble("dmn_y",dmn_y, nAtm, n_interp2,thisIndex.y,thisIndex.x,thisIndex.x,iterNL,false);    

    dumpMatrix2DDouble("dmn_z",dmn_z, nAtm, n_interp2,thisIndex.y,thisIndex.x,thisIndex.x,iterNL,false);    
    dumpMatrix2DInt("igrid",igrid, nAtm, n_interp2,thisIndex.y,thisIndex.x,thisIndex.x,iterNL,false);    
#endif
#endif


#ifdef _CP_GS_DEBUG_COMPARE_VKS_
    // at this point we're checking that nobody poisoned projPsiR, igrid dmn_[xyz] mn 
    // projPsiR uses the same location as projpsiC, so we hold our nose
    // and cut and paste


  if(savedprojpsiC==NULL)
    { // load it
      savedprojpsiC= new complex[csize];
      loadMatrixDouble("projPsiC",(double *)savedprojpsiC, 1, csize*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
    }
#ifdef _INSANEO_PARANOID_COMPARE_THAT_EATS_HUGE_MEMORY_
  if(savedmn==NULL)
    { // load it
      savedigrid   = (int **)fftw_malloc(nAtm*sizeof(int*));

      savedmn    = (double **)fftw_malloc(nAtm*sizeof(double*));
      saveddmn_x = (double **)fftw_malloc(nAtm*sizeof(double*));
      saveddmn_y = (double **)fftw_malloc(nAtm*sizeof(double*));
      saveddmn_z = (double **)fftw_malloc(nAtm*sizeof(double*));
      // Argh! this crazy iterate from 1 stuff escaped from
      // the piny box of evil fortranisms
      for(int i=0;i<nAtm;i++){
	savedigrid[i]    = (int *)fftw_malloc(n_interp2*sizeof(int))-1;    
	double *tmp = (double *)fftw_malloc(4*n_interp2*sizeof(double));
	// by dint of more ugly trickery we store 4 array rows in 1
	// which provides a nice way for us to quietly violate logical
	// array boundaries without tripping over memory we don't own
	int ioff    = 0;
	savedmn[i]       = &tmp[ioff]-1;  ioff+=n_interp2;
	saveddmn_x[i]    = &tmp[ioff]-1;  ioff+=n_interp2;
	saveddmn_y[i]    = &tmp[ioff]-1;  ioff+=n_interp2;
	saveddmn_z[i]    = &tmp[ioff]-1;  
      }//endfor

      loadMatrix2DDouble("mn",savedmn, nAtm, n_interp2,thisIndex.y,thisIndex.x,thisIndex.x,iterNL,false);    
      loadMatrix2DDouble("dmn_x",saveddmn_x, nAtm, n_interp2,thisIndex.y,thisIndex.x,thisIndex.x,iterNL,false);    
      loadMatrix2DDouble("dmn_y",saveddmn_y, nAtm, n_interp2,thisIndex.y,thisIndex.x,thisIndex.x,iterNL,false);    
      loadMatrix2DDouble("dmn_z",saveddmn_z, nAtm, n_interp2,thisIndex.y,thisIndex.x,thisIndex.x,iterNL,false);    
      loadMatrix2DInt("igrid",savedigrid, nAtm, n_interp2,thisIndex.y,thisIndex.x,thisIndex.x,iterNL,false);    
    }

  for(int i=0;i<nAtm;i++)
    {
      if(plane_index[i]!=0){

	for(int j=1;j<=n_interp2;j++){
	  if(fabs(mn[i][j]-savedmn[i][j])>0.0001)
	    {
	      fprintf(stderr, "RPP [%d,%d] %d %d element mn  %.10g not %.10g\n",thisIndex.x, thisIndex.y,i,j, mn[i][j], savedmn[i][j]);
	    }
	  CkAssert(fabs(mn[i][j]-savedmn[i][j])<0.0001);
	  if(fabs(dmn_x[i][j]-saveddmn_x[i][j])>0.0001)
	    {
	      fprintf(stderr, "RPP [%d,%d] %d %d element dmn_x  %.10g not %.10g\n",thisIndex.x, thisIndex.y,i,j, dmn_x[i][j], saveddmn_x[i][j]);
	    }
	  CkAssert(fabs(dmn_x[i][j]-saveddmn_x[i][j])<0.0001);
	  if(fabs(dmn_y[i][j]-saveddmn_y[i][j])>0.0001)
	    {
	      fprintf(stderr, "RPP [%d,%d] %d %d element dmn_y  %.10g not %.10g\n",thisIndex.x, thisIndex.y,i,j, dmn_y[i][j], saveddmn_y[i][j]);
	    }
	  CkAssert(fabs(dmn_y[i][j]-saveddmn_y[i][j])<0.0001);
	  if(fabs(dmn_z[i][j]-saveddmn_z[i][j])>0.0001)
	    {
	      fprintf(stderr, "RPP [%d,%d] %d %d element dmn_z  %.10g not %.10g\n",thisIndex.x, thisIndex.y,i,j, dmn_z[i][j], saveddmn_z[i][j]);
	    }
	  CkAssert(fabs(dmn_z[i][j]-saveddmn_z[i][j])<0.0001);
	  if(igrid[i][j]!=savedigrid[i][j])
	    {
	      fprintf(stderr, "RPP [%d,%d] %d %d element igrid  %d not %d\n",thisIndex.x, thisIndex.y,i,j, igrid[i][j], savedigrid[i][j]);
	    }
	  CkAssert(igrid[i][j]==savedigrid[i][j]);
	}
      }
    }
#endif
#endif

#if CMK_TRACE_ENABLED
   double  StartTime=CmiWallTimer();
#endif    

   double cp_enl_now = 0.0;
   CPNONLOCAL::eesEnergyAtmForcRchare(iterNL,&cp_enl_now,zmat,igrid,mn,dmn_x,dmn_y,dmn_z,
 		      projPsiR,projPsiRScr,plane_index,nBreakJ,sBreakJ,
                      myPlane,thisIndex.x,fastAtoms);
   cp_enl += cp_enl_now;

#if CMK_TRACE_ENABLED
  traceUserBracketEvent(eesEnergyAtmForcR_, StartTime, CmiWallTimer());    
#endif

//============================================================================
// If we are done, send out the energy : HELP HELP Evil Section Multicast

   if(thisIndex.y==reductionPlaneNum && iterNL==numIterNl){
#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
     if(thisIndex.x==0)
       CkPrintf("HI, I am rPP %d %d sending Enl : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif

#ifdef _FANCY_RED_METHOD_
     CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
     //     CkCallback cb=CkCallback(printEnlR, NULL);
     CkCallback cb(CkIndex_CP_State_RealParticlePlane::printEnlR(NULL),
                 CkArrayIndex2D(0,reductionPlaneNum),
                 UrealParticlePlaneProxy[thisInstance.proxyOffset].ckGetArrayID());
     mcastGrp->contribute(sizeof(double),(void*) &cp_enl, 
		          CkReduction::sum_double,rEnlCookie, cb, itime);
#else
     thisProxy(0,0).printEnlRSimp(cp_enl,thisIndex.x,itime);
#endif
#ifdef CMK_BLUEGENEL
       CmiNetworkProgress();
#endif
   }//endif

#ifdef _CP_GS_DUMP_VKS_
    dumpMatrixDouble("zmat",(double *)zmat, 1, nZmat_in,thisIndex.y,thisIndex.x,thisIndex.x,iterNL,false);    
    dumpMatrixDouble("projPsiRScr",(double *)projPsiRScr, 1, planeSize,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
#endif

#ifdef _CP_GS_DEBUG_COMPARE_VKS_
    // check that we did the zmat gather/scatter correctly and 

  if(savedzmat==NULL)
    { // load it
      savedzmat= new double[nZmat_in];
      loadMatrixDouble("zmat",(double *)savedzmat, 1, nZmat_in,thisIndex.y,thisIndex.x,thisIndex.x,iterNL,false);    
    }
  for(int i=0;i<nZmat_in; i++)
    {
      if(fabs(zmat[i]-savedzmat[i])>0.0001)
	{
	  fprintf(stderr, "RPP [%d,%d] %d element projpsi  %.10g not %.10g\n",thisIndex.x, thisIndex.y,i, zmat[i], savedzmat[i]);
	}
      CkAssert(fabs(zmat[i]-savedzmat[i])<0.0001);

    }
  if(savedProjpsiRScr==NULL)
    { // load it
      savedProjpsiRScr= new double[planeSize];
      loadMatrixDouble("projPsiRScr",(double *)savedProjpsiRScr, 1, planeSize,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
    }
  for(int i=0;i<planeSize;i++)
    {
      if(fabs(projPsiRScr[i]-savedProjpsiRScr[i])>0.0001)
	{
	  fprintf(stderr, "RPP [%d,%d] %d element projpsi  %.10g not %.10g\n",thisIndex.x, thisIndex.y,i, projPsiRScr[i], savedProjpsiRScr[i]);
	  dumpMatrixDouble("badprojPsiRScr",(double *)projPsiRScr, 1, planeSize,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
	}
      CkAssert(fabs(projPsiRScr[i]-savedProjpsiRScr[i])<0.0001);

    }

#endif

   // zero the total enl energy if we are done.
   if(iterNL==numIterNl){
     cp_enl = 0.0;
     launchFFT=false;
   }

//============================================================================
// Time to make the Psiforces (donuts!)

  FFTNLEesBckR();
  // do not delete nokeep message
//----------------------------------------------------------------------------
}//end routine
//============================================================================

//============================================================================
// Control launch of FFT based on Rho having more of its act together
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealParticlePlane::launchFFTControl(int time_in){
  rhoRTime = time_in;
  launchFFT=true;
  if(fftDataDone){
    if(iterNL!=1){CkPrintf("Badddd launchFFT.3\n");CkExit();}
    if(rhoRTime!=itime){CkPrintf("Badddd launchFFT.4\n");CkExit();}
    thisProxy(thisIndex.x,thisIndex.y).FFTNLEesFwdR();
  }//endif
//----------------------------------------------------------------------------
}
//============================================================================


//============================================================================
// Do the FFT of projPsi(x,y,z) to get projPsi(gx,gy,z)
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealParticlePlane::FFTNLEesBckR(){
//============================================================================

  int nplane_x = scProxy.ckLocalBranch()->cpcharmParaInfo->nplane_x;

#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
  if(thisIndex.x==0)
   CkPrintf("HI, I am rPP %d %d in FFTNLbck : %d %d %d %d\n",
     thisIndex.x,thisIndex.y,iterNL,ngridA,ngridB,nplane_x);
#endif

  FFTcache *fftcache   = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();  
  double  *projPsiRScr = fftcache->tmpDataR;
  complex *projPsiCScr = fftcache->tmpData;


#ifdef _CP_GS_DEBUG_COMPARE_VKS_
  // double check that projPsiRScr is still ok
  if(savedProjpsiRScr==NULL)
    { // load it
      savedProjpsiRScr= new double[planeSize];
      loadMatrixDouble("projPsiRScr",(double *)savedProjpsiRScr, 1, planeSize,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
    }
  for(int i=0;i<planeSize;i++)
    {
      if(fabs(projPsiRScr[i]-savedProjpsiRScr[i])>0.0001)
	{
	  fprintf(stderr, "RPP [%d,%d] %d element projpsi  %.10g not %.10g\n",thisIndex.x, thisIndex.y,i, projPsiRScr[i], savedProjpsiRScr[i]);
	}
      CkAssert(fabs(projPsiRScr[i]-savedProjpsiRScr[i])<0.0001);
      CkAssert(fabs(projPsiRScr[i]-savedProjpsiRScr[i])<0.0001);
    }

#endif

#if CMK_TRACE_ENABLED
  double StartTime= CmiWallTimer();    
#endif
  /*rfftwnd_plan plan = fftcache->bwdXPlanNL.rfftwPlan;
  CkPrintf("outside rfft plan %p howmany %d in %p istride %d idist %d out %p ostride %d odist %d\n ",plan,  ngridB, projPsiRScr, 1, ngridA+2, NULL, 0 ,0 ) ;
  */
    fftcache->doNlFFTRtoG_Rchare(projPsiCScr,projPsiRScr,nplane_x,ngridA,ngridB,myPlane);
  /*

	   *out_in,  ostride_in, odist_in, split);

 rfftwnd_real_to_complex(
   	     plan,             // backward plan
	     ngridB,                      // these many 1D ffts
	     (fftw_real *)projPsiRScr,         // data set
             1,                          // stride
             (ngridA+2),                  // spacing between data sets
  	     NULL,0,0                   // input array is output array
           );     
  */

#ifdef _CP_GS_DUMP_VKS_
    dumpMatrixDouble("projPsiCScr",(double *)projPsiCScr, 1, (ngridA+2)*ngridB,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
#endif

#ifdef _CP_GS_DEBUG_COMPARE_VKS_
    // at this point we're checking that we did the fft correctly
    // which implicitly checks that we did the zmat gather/scatter correctly

  if(savedProjpsiCScr==NULL)
    { // load it
      savedProjpsiCScr= new complex[csize];
      loadMatrixDouble("projPsiCScr",(double *)savedProjpsiCScr, 1, csize*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
    }
  for(int i=0;i<csize;i++)
    {
      if(fabs(projPsiCScr[i].re-savedProjpsiCScr[i].re)>0.0001)
	{
	  fprintf(stderr, "RPP [%d,%d] %d element projpsi  %.10g not %.10g\n",thisIndex.x, thisIndex.y,i, projPsiCScr[i].re, savedProjpsiCScr[i].re);
	  fprintf(stderr, "HI, I am rPP %d %d in FFTNLbck : %d %d %d %d\n",
     thisIndex.x,thisIndex.y,iterNL,ngridA,ngridB,nplane_x);

	}
      CkAssert(fabs(projPsiCScr[i].re-savedProjpsiCScr[i].re)<0.0001);
      CkAssert(fabs(projPsiCScr[i].im-savedProjpsiCScr[i].im)<0.0001);
    }
#endif


#if CMK_TRACE_ENABLED
  traceUserBracketEvent(doNlFFTRtoG_, StartTime, CmiWallTimer());    
#endif

  sendToEesGPP();

//============================================================================
   }//end routine
//============================================================================


//============================================================================
// Send the PsiForce NL FFT back to GparticlePlane.
// Until ALL RealPP chares send forces, PP cannot start the next iteration
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealParticlePlane::sendToEesGPP(){
//============================================================================

  CPcharmParaInfo *sim   = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
  FFTcache *fftcache     = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();  
  int nchareG            = sim->nchareG;
  int **tranpack         = sim->index_tran_upackNL;
  int *nlines_per_chareG = sim->nlines_per_chareG;
  complex *projPsiCScr   = fftcache->tmpData;

//===================================================================
// Perform the transpose and then the blast off the final 1D-FFT

#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
  if(thisIndex.x==0)
   CkPrintf("HI, I am rPP %d %d in sendtoGPP : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif

#ifdef USE_COMLIB
#ifdef OLD_COMMLIB
  if(config.useMssInsGPP){mssPInstance.beginIteration();}
#else
//    if(config.useMssInsGPP){ComlibBegin(UparticlePlaneProxy[thisInstance.proxyOffset], iterNL);}
#endif
#endif


    for (int ic=0;ic<nchareG;ic++) { // chare arrays to which we will send
      int sendFFTDataSize = nlines_per_chareG[ic];
      GSPPIFFTMsg *msg    = new (sendFFTDataSize, 8 * sizeof(int)) GSPPIFFTMsg; 
      msg->iterNL         = iterNL;
      msg->size           = sendFFTDataSize;
      msg->offset         = thisIndex.y;    // c-plane-index
      complex *data       = msg->data;
      if(config.prioNLFFTMsg){
	  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	  *(int*)CkPriorityPtr(msg) = config.gsNLfftpriority+thisIndex.x;
      }//endif
      for(int i=0;i<sendFFTDataSize;i++){data[i] = projPsiCScr[tranpack[ic][i]];}
      UparticlePlaneProxy[thisInstance.proxyOffset](thisIndex.x, ic).recvFromEesRPP(msg); // send the message
#ifdef CMK_BLUEGENEL
       CmiNetworkProgress();
#endif
    }//end for : chare sending

#ifdef USE_COMLIB
#ifdef OLD_COMMLIB
  if(config.useMssInsGPP){mssPInstance.endIteration();}
#else
//    if(config.useMssInsGPP){ComlibEnd(UparticlePlaneProxy[thisInstance.proxyOffset], iterNL);}
#endif
#endif

//============================================================================
// If it looks like this is the end, reset my counters baby.
// If its not the end, GParticlePlane will invoke entry methods.

  recvBlock = 0; // I am done, I can recv messages now from PP
                 // If I recv with block==1, something is terribly wrong
  if(iterNL==numIterNl){
    iterNL = 0;
    cp_enl = 0;
    launchFFT=false;
  }//endif
   
  fftcache->freeCacheMem("CP_State_RealParticlePlane::sendToEesGPP");

//----------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
/**
 * spread the reduction plane numbers around to minimize map collisions
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
int CP_State_RealParticlePlane::calcReductionPlaneNum(int state){
//============================================================================
  
  int nstatemax=nstates-1;
  int ncharemax=nChareR-1;
  int planeNum= (state %ncharemax);
  if(planeNum<0){
      CkPrintf(" PP [%d %d] calc nstatemax %d ncharemax %d state %d planenum %d\n",
                 thisIndex.x, thisIndex.y,nstatemax, ncharemax, state, planeNum); 
      CkExit();
  }//endif
  return planeNum;

}//end routine
//============================================================================


//============================================================================
// Recv a dummy message, 1 integer,  and set the cookie monster
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealParticlePlane::setPlaneRedCookie(EnlCookieMsg *m){
//============================================================================
// Do not delete msg. Its a nokeep.
//============================================================================
#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
  if(thisIndex.x==0)
   CkPrintf("RPP[%d %d] gets rPlaneRedCookie\n",thisIndex.x, thisIndex.y);
#endif

  CkGetSectionInfo(rPlaneRedCookie,m);
  int cookiedone=1;
  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
   CkCallback cb(CkIndex_CP_State_RealParticlePlane::planeRedSectDone(NULL),
                 CkArrayIndex2D(thisIndex.x,reductionPlaneNum),
                 UrealParticlePlaneProxy[thisInstance.proxyOffset].ckGetArrayID());
   mcastGrp->contribute(sizeof(int),&cookiedone,CkReduction::sum_int,
     rPlaneRedCookie,cb);
     // if not root, we are set up
  if(thisIndex.y!=reductionPlaneNum)
    {
      planeRedSectionComplete=true;
      initComplete();
    }
}
//============================================================================


//============================================================================
// Recv a dummy message, 1 integer,  and set the cookie monster
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealParticlePlane::setEnlCookie(EnlCookieMsg *m){
//============================================================================
// Do not delete msg. Its a nokeep.
//============================================================================
#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
  if(thisIndex.x==0)
   CkPrintf("RPP[%d %d] gets rEnlCookie\n",thisIndex.x, thisIndex.y);
#endif
  CkGetSectionInfo(rEnlCookie,m);
  int cookiedone=1;
  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
  CkCallback cb(CkIndex_CP_State_RealParticlePlane::enlSectDone(NULL),
                 CkArrayIndex2D(0,reductionPlaneNum),
                 UrealParticlePlaneProxy[thisInstance.proxyOffset].ckGetArrayID());
     mcastGrp->contribute(sizeof(int), &cookiedone, 
		          CkReduction::sum_int,rEnlCookie, cb);
     // if not root, we are set up
  if(thisIndex.y!=reductionPlaneNum || thisIndex.x!=0){     
    enlSectionComplete=true;
    initComplete();
  }
}
//============================================================================

//==========================================================================
// Make sure everyone is registered on the 1st time step
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
  void CP_State_RealParticlePlane::registrationDone() {
//==========================================================================


#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
  if(thisIndex.x==0)
   CkPrintf("HI, I am realPPa %d %d \n",thisIndex.x,thisIndex.y);
#endif

  registrationFlag=1;
  if(iterNL==1){computeZmatEes();}

  if(iterNL>1){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Homes, registeration must occur before the 1st time step in RPP\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

}
//==========================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealParticlePlane::pup(PUP::er &p) {
//============================================================================

   p|rhoRTime;
   p|ees_nonlocal;
   p|nChareR; 
   p|nChareG;
   p|Rstates_per_pe;
   p|myPlane;
   p|registrationFlag;
   p|recvBlock;

   p|numIterNl;
   p|countEnl;
   p|count;
   p|iterNL;
   p|itime;
   p|itimeRed;
   p|countZ;
 
   p|ngridA;
   p|ngridB;
   p|ngridC;
   p|planeSize;
   p|planeSizeT;
   p|csize;
   p|zmatSizeMax;
   p|reductionPlaneNum;
   p|cp_enl; 
   p|cp_enlTot; 

   p|rPlaneSectProxy;
   p|rPlaneENLProxy;
   p|rPlaneRedCookie;
   p|rEnlCookie;
   p|gPP_proxy;

   if(ees_nonlocal==1){
     if (p.isUnpacking()) {
       projPsiC = (complex*) fftw_malloc(csize*sizeof(complex));
       projPsiR = reinterpret_cast<double*> (projPsiC);
       zmat     = new double[zmatSizeMax];
       zmatScr = new double[zmatSizeMax];
     }//endif
     p((char*)projPsiC,csize *sizeof(complex));
     p(zmat,zmatSizeMax);
     p(zmatScr,zmatSizeMax);
   }//endif

//---------------------------------------------------------------------------
  }//end routine
//============================================================================
