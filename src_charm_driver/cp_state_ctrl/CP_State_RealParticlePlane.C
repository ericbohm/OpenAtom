//=========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=========================================================================
/** \file CP_State_RealParticlePlane.C
 * Life-cycle of a CP_State_RealParticlePlane:
 *
 * Insert descriptive comment here please
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
extern CProxy_StructFactCache            sfCacheProxy;
extern CProxy_eesCache                   eesCacheProxy;
extern CProxy_FFTcache                   fftCacheProxy;
extern CProxy_CP_State_RealParticlePlane realParticlePlaneProxy;
extern CProxy_EnergyGroup                egroupProxy; //energy group proxy

extern CkGroupID            mCastGrpId;
extern ComlibInstanceHandle mssPInstance;

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
  delete m;

  //output and save the data
  CkPrintf("ENL(EES)    = %5.8lf\n", d);
  gSpacePlaneProxy(0,0).computeEnergies(ENERGY_ENL, d);  
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
    CkPrintf("ENL(EES)    = %5.8lf\n", cp_enlTot);
    gSpacePlaneProxy(0,0).computeEnergies(ENERGY_ENL,cp_enlTot);  
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
                             int ees_nonlocal_in)
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

  numIterNl      = numIterNL_in;         // # of non-local iterations per time step
  zmatSizeMax    = zmatSizeMax_in;       // zmatrix size
  ees_nonlocal   = ees_nonlocal_in;

  cp_enl           = 0.0;                // non-local energy
  cp_enlTot        = 0.0;
  count            = 0;                  // fft communication counter
  countEnl         = 0;                  // Energy counter
  iterNL           = 0;                  // Nl iteration counter
  itimeRed         = 0;                  // time of reduction
  itime            = 0;                  // time step;
  registrationFlag = 0;
  recvBlock        = 0;

  countZ = 0.0;      // Zmat communication counter

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
    eesCache *eesData  = eesCacheProxy.ckLocalBranch (); 
    eesData->registerCacheRPP(thisIndex.y);
   //--------
   // Tell your friends your are ready to boogy
    int i=1;
    CkCallback cb(CkIndex_CP_State_RealParticlePlane::registrationDone(NULL),
		  realParticlePlaneProxy);
    contribute(sizeof(int),&i,CkReduction::sum_int,cb);
  }//endif

//============================================================================
// Choose reduction plane reasonably intelligently

  int *red_pl       = new int[nstates];
#ifdef USE_TOPOMAP
  int l             = Rstates_per_pe;
  int pl            = (nstates / l);
  int pm            = (CkNumPes() / pl);
  if(pm==0){CkAbort("Choose a larger Gstates_per_pe\n");}
  int planes_per_pe = (nChareR / pm);
  for(int i=0; i<nstates;i++){
    red_pl[i]= (  ((i % Rstates_per_pe)*planes_per_pe)% nChareR);
  }//endif
  reductionPlaneNum = red_pl[thisIndex.x];
#else
  for(int i=0; i<nstates;i++){
    red_pl[i] = calcReductionPlaneNum(thisIndex.x);
  }//endif
  reductionPlaneNum = calcReductionPlaneNum(thisIndex.x);
#endif

//============================================================================
// Build section reductions

  //-----------------------------------------------------------
  // The plane section reduction to the reduction plane for zmat
  // Only the root of the reduction does setup dance.
  //-----------------------------------------------------------
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
  }//endif
  delete [] red_pl;

//============================================================================
// Setup the comlib to talk to GPP

  gPP_proxy = particlePlaneProxy;
  if (config.useMssInsGPP){
     ComlibAssociateProxy(&mssPInstance,gPP_proxy);
  }//endif

//============================================================================
// No migration: No atSync load-balancing act

  setMigratable(false);
  usesAtSync = CmiFalse;

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


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
      FFTNLEesFwdR();
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

  int nplane_x = scProxy.ckLocalBranch()->cpcharmParaInfo->nplane_x;
#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
  if(thisIndex.x==0)
   CkPrintf("HI, I am rPP %d %d in FFTNL : %d %d %d %d\n",
     thisIndex.x,thisIndex.y,iterNL,ngridA,ngridB,nplane_x);
#endif

  fftCacheProxy.ckLocalBranch()->doNlFFTGtoR_Rchare(projPsiC,projPsiR,
                                                    nplane_x,ngridA,ngridB);

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

   eesCache *eesData  = eesCacheProxy.ckLocalBranch (); 
   if(iterNL==1){eesData->queryCacheRPP(myPlane,itime,iterNL);}// query once a t-step

   int *plane_index = eesData->RppData[myPlane].plane_index;
   int **igrid      = eesData->RppData[myPlane].igrid;
   double **mn      = eesData->RppData[myPlane].mn;

#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
   if(thisIndex.x==0)
    CkPrintf("HI, I am rPP %d %d in computeZmat : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif

#ifndef CMK_OPTIMIZE
   double  StartTime=CmiWallTimer();
#endif    

   CPNONLOCAL::eesZmatRchare(projPsiR,iterNL,zmat,igrid,mn,
                             plane_index,thisIndex.x,myPlane);
#ifndef CMK_OPTIMIZE
   traceUserBracketEvent(eesZmatR_, StartTime, CmiWallTimer());    
#endif

//============================================================================
// Launch section reduction : 

#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
   if(thisIndex.x==0)
    CkPrintf("HI, I am rPP %d %d in send to zmat-red : %d %d %d d\n",
         thisIndex.x,thisIndex.y,iterNL,reductionPlaneNum,nZmat,zmatSizeMax);
#endif

#define _FANCY_RED_METHOD_OFF_
#ifdef _FANCY_RED_METHOD_
   CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch(); 
   CkCallback cb(CkIndex_CP_State_RealParticlePlane::recvZMatEes(NULL),
                 CkArrayIndex2D(thisIndex.x,reductionPlaneNum),
                 realParticlePlaneProxy.ckGetArrayID());
   mcastGrp->contribute((nZmat*sizeof(double)),zmat,CkReduction::sum_double,
                         rPlaneRedCookie,cb);
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
        thisProxy(thisIndex.x,i).computeAtmForcEes(nZmat,zmatScr,iterNL_in);
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
   int iterNL1          = iterNL-1;           // silly C++ convention
   int *nmem_zmat       = sim->nmem_zmat;     // zmat size now
   int nZmat            = nmem_zmat[iterNL1];

//============================================================================
// Recv the reduced zmatrix and chuck the message

#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
    if(thisIndex.x==0)
     CkPrintf("HI, I am rPP %d %d in recvZmat : %d\n",thisIndex.x,thisIndex.y,iterNL);
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

   rPlaneSectProxy.computeAtmForcEes(nZmat,zmat,iterNL);

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
// Use Zmat and ProjPsi to get atmForces, Energy and psiforces
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealParticlePlane::computeAtmForcEes(int nZmat_in, double *zmat_loc,
                                                   int iterNL_in){
//============================================================================
// Unpack the message : Get some sizes

#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
   if(thisIndex.x==0)
    CkPrintf("HI, I am rPP %d %d in compteAtmforc : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif

   CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
   int iterNL1          = iterNL-1;          // silly C++ convention
   int *nmem_zmat       = sim->nmem_zmat;    // zmat memory size
   int nZmat            = nmem_zmat[iterNL1];

   if(iterNL_in !=  iterNL || nZmat_in != nZmat){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Dude, iteration mismatch : %d %d : %d %d \n",iterNL,iterNL_in,
                   thisIndex.x,thisIndex.y);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
   }//endif

   CmiMemcpy(zmat,zmat_loc,sizeof(double)*nZmat);

//============================================================================
// Check out your B-splines from the cache and then compute energy and forces

   eesCache *eesData   = eesCacheProxy.ckLocalBranch (); 
   FFTcache *fftcache  = fftCacheProxy.ckLocalBranch();  
   int *plane_index    = eesData->RppData[myPlane].plane_index;
   int **igrid         = eesData->RppData[myPlane].igrid;
   double **dmn_x      = eesData->RppData[myPlane].dmn_x;
   double **dmn_y      = eesData->RppData[myPlane].dmn_y;
   double **dmn_z      = eesData->RppData[myPlane].dmn_z;
   double **mn         = eesData->RppData[myPlane].mn;
   double *projPsiRScr = fftcache->tmpDataR;

   AtomsGrp *ag         = atomsGrpProxy.ckLocalBranch();
   FastAtoms *fastAtoms = &(ag->fastAtoms);

#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
   if(thisIndex.x==0)
    CkPrintf("HI, I am rPP %d %d in compteAtmforc : %d\n",thisIndex.x,thisIndex.y,iterNL);
#endif

   // projPsiR     comes in with info for atoms
   // projPsiRSsr leaves with info for psiforces

#ifndef CMK_OPTIMIZE
   double  StartTime=CmiWallTimer();
#endif    

   CPNONLOCAL::eesEnergyAtmForcRchare(iterNL,&cp_enl,zmat,igrid,mn,dmn_x,dmn_y,dmn_z,
                       projPsiR,projPsiRScr,plane_index,myPlane,thisIndex.x,fastAtoms);
#ifndef CMK_OPTIMIZE
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
                 realParticlePlaneProxy.ckGetArrayID());
     mcastGrp->contribute(sizeof(double),(void*) &cp_enl, 
		          CkReduction::sum_double,rEnlCookie, cb);
#else
     thisProxy(0,0).printEnlRSimp(cp_enl,thisIndex.x,itime);
#endif
#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif
   }//endif

   // zero the total enl energy if we are done.
   if(iterNL==numIterNl){cp_enl = 0.0;}

//============================================================================
// Time to make the Psiforces (donuts!)

  FFTNLEesBckR();

//----------------------------------------------------------------------------
  }//end routine
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

  FFTcache *fftcache   = fftCacheProxy.ckLocalBranch();  
  double  *projPsiRScr = fftcache->tmpDataR;
  complex *projPsiCScr = fftcache->tmpData;
  fftcache->doNlFFTRtoG_Rchare(projPsiCScr,projPsiRScr,nplane_x,ngridA,ngridB);

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
  FFTcache *fftcache     = fftCacheProxy.ckLocalBranch();  
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

  if(config.useMssInsGPP){mssPInstance.beginIteration();}

    for (int ic=0;ic<nchareG;ic++) { // chare arrays to which we will send
      int sendFFTDataSize = nlines_per_chareG[ic];
      GSPPIFFTMsg *msg    = new (sendFFTDataSize, 8 * sizeof(int)) GSPPIFFTMsg; 
      msg->iterNL         = iterNL;
      msg->size           = sendFFTDataSize;
      msg->offset         = thisIndex.y;    // c-plane-index
      complex *data       = msg->data;
      if(config.prioNLFFTMsg){
	  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	  *(int*)CkPriorityPtr(msg) = config.gsNLfftpriority+thisIndex.x*planeSize;
      }//endif
      for(int i=0;i<sendFFTDataSize;i++){data[i] = projPsiCScr[tranpack[ic][i]];}
      particlePlaneProxy(thisIndex.x, ic).recvFromEesRPP(msg); // send the message
#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif
    }//end for : chare sending

  if(config.useMssInsGPP){mssPInstance.endIteration();}

//============================================================================
// If it looks like this is the end, reset my counters baby.
// If its not the end, GParticlePlane will invoke entry methods.

  recvBlock = 0; // I am done, I can recv messages now from PP
                 // If I recv with block==1, something is terribly wrong
  if(iterNL==numIterNl){
    iterNL = 0;
    cp_enl = 0;
  }//endif
   
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

#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
  if(thisIndex.x==0)
   CkPrintf("RPP[%d %d] gets rPlaneRedCookie\n",thisIndex.x, thisIndex.y);
#endif

  CkGetSectionInfo(rPlaneRedCookie,m);
}
//============================================================================


//============================================================================
// Recv a dummy message, 1 integer,  and set the cookie monster
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealParticlePlane::setEnlCookie(EnlCookieMsg *m){

#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
  if(thisIndex.x==0)
   CkPrintf("RPP[%d %d] gets rEnlCookie\n",thisIndex.x, thisIndex.y);
#endif
  CkGetSectionInfo(rEnlCookie,m);

}
//============================================================================

//==========================================================================
// Make sure everyone is registered on the 1st time step
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
  void CP_State_RealParticlePlane::registrationDone(CkReductionMsg *msg) {
//==========================================================================

  int sum = ((int *)msg->getData())[0];
  delete msg;

#ifdef _CP_DEBUG_STATE_RPP_VERBOSE_
  if(thisIndex.x==0)
   CkPrintf("HI, I am realPPa %d %d in reg : %d\n",thisIndex.x,thisIndex.y,sum);
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
