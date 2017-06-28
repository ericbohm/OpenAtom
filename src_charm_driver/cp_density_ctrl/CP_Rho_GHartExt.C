//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file CP_Rho_GHartExt.C
 *  This is a description of the "life" of a CP_Rho_GHartExc  object
 *
 *  At the start of the program, the constructor CP_Rho_GHartExt is
 *  called.  We create our own rho_gs slab because we need the
 *  kvectors and the fft code.
 *
 *  N^2 method:
 *   Each iteration, the CP_Rho_GpacePlane object sends us the same rho
 *   data it will use in its own FFTs.  We then call the very intensive
 *   HartExtVksG function on the data.  We contribute the energy to the
 *   reduction, fft the vks, and ship the vks to rhoReal. Done!
 *
 *
 *  N log N EES (Euler-Exponential Spline) method:
 *   Each iteration, the CP_Rho_GpacePlane object sends us the same rho
 *   data it will use in its own FFTs.  The CP_Rho_RHartExt object will
 *   send a partially FFTed EES approximated atom structure factor for an
 *   atom type. The atom types (water has two) are parallelized by nchareHartAtmT.
 *   Using an EES Atm SF and the density, eext energy, VKS and atom forces are computed.
 *   Atom forces need to go back to RHartEext to be completed and applied to that atoms.
 *   When all the atom types are down, the total atom SF is accumulated and used to
 *   compute the Ewald energy - hartree is computed at the same time. The vks is partly
 *   ffted and shipped back to rhoReal. The ewald is partly ffted and shipped back
 *   to CP_Rho_RHartExt to get the forces on the atoms. Done!
 *
 *  LSDA there are two contributions to the density rho=rhoUp+rhoDn
 *   Both have to arrive before you can proceed. 
 *   you also need vks to have been computed and sent back to down
 *   
 */
//============================================================================

#include <iostream>
#include <fstream>
#include <cmath>
using std::isnan;

#include "debug_flags.h"
#include "utility/util.h"
#include "main/cpaimd.h"
#include "main/AtomsCache.h"
#include "main/eesCache.h"
#include "cp_state_ctrl/CP_State_Plane.h"

#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cplocal.h"
#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h"
#include "src_piny_physics_v1.0/include/proto_defs/proto_cp_ewald_corrs.h"
#include "src_piny_physics_v1.0/include/class_defs/Interface_ctrl.h"
#include "main/PhysScratchCache.h"

#include "fft_charm.h"

extern CkVec <CProxy_CP_Rho_RealSpacePlane>     UrhoRealProxy;
extern CProxy_CPcharmParaInfoGrp                scProxy;
extern CkVec <CProxy_PhysScratchCache>     UpScratchProxy;
extern CkVec <CProxy_AtomsCache>                UatomsCacheProxy;
extern CkVec <CProxy_CP_Rho_GHartExt>           UrhoGHartExtProxy;
extern CkVec <CProxy_CP_Rho_RHartExt>           UrhoRHartExtProxy;
extern CkVec <CProxy_eesCache>                  UeesCacheProxy;
extern Config                                   config;
extern CPcharmParaInfo                          simReadOnly;
extern CkVec <CProxy_fft2d>                     Urho_fft_hartProxy;
extern CkVec < CkVec<CProxy_fft2d> >             Urho_fft_atmSFProxy;
extern CkVec <CkVec<CProxySection_CP_Rho_GHartExt> >        Ughart_sectionProxy;
extern CkGroupID                                mCastGrpId;

//#define _DEBUG_INT_TRANS_FWD_
//#define _CP_GHART_VERBOSE_ 1

/** @addtogroup Density
  @{
*/

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 *  This object just gets a rho message, computes GHartExt, and sends
 *  vks.
 */
//============================================================================
CP_Rho_GHartExt::CP_Rho_GHartExt( UberCollection _instance) :
  thisInstance(_instance)
{

  //============================================================================
  // Fetch and compute the useful variables : zero counters

  launchFlag      = 0;
  registrationFlag= 0;
  atmSFHere       = 0;
  iterAtmTyp      = 0;
  nsendAtmTyp     = 0;
  countVksTot     = 0;
  CountDebug      = 0;

  iteration       = 0;
  densityHere     = 0;
  ehart_ret       = 0.0;
  eext_ret        = 0.0;
  ewd_ret         = 0.0;
  countAtmSFtot   = 0;
  CountDebug      = 0;

  //RAZ: added spin variables:
  numAcceptDensity = 0;
  mySpinIndex      = thisInstance.idxU.s;
  cp_lsda          = simReadOnly.cp_lsda;

  //==================================================================================
  // AtmTyp parallelization
  int natmTypTot     = simReadOnly.natm_typ;
  int div        = (natmTypTot/config.nchareHartAtmT);
  int rem        = (natmTypTot % config.nchareHartAtmT);
  int max        = (thisIndex.y < rem ? thisIndex.y : rem);
  natmTyp        = (thisIndex.y < rem ? div+1 : div);
  atmTypoff      = div*thisIndex.y + max;

  if(simReadOnly.ees_eext_on == 0 && config.nchareHartAtmT > 1){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Parallel atom type not supported without ees Eext\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  Rho = Vks = NULL;
  atmSF = atmSFtot = NULL;
  VksRecv = atmSFtotRecv = NULL;
  atmSF = atmSFtot = atmSFtotRecv = NULL;
  VksRecv = NULL;

  //==================================================================================
  // Set some proxies, set the migratable flag

  setMigratable(false);

  usesAtSync = true;
  if(config.lbdensity){
    setMigratable(true);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Hold on there partner. Migratable true in CP_Rho_GHartExt?? \n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }else{
    setMigratable(false);
  }//endif

  //---------------------------------------------------------------------------
}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * Post constructor initialization - register in the ees cache if
 * necessary, contribute to cache reductions, set up some proxies, setmigratable,
 * periodic BC for wires,surfaces,clusters ...
 */
//============================================================================
void CP_Rho_GHartExt::init(){
  //==================================================================================
  // Register in the cache : contribute to a reduction to be sure everyone is done
  //query fft library to find my extents
  if(config.nchareHartAtmT == 1 || thisIndex.y == 1) {
    Charm_getOutputIndex(thisIndex.x, fft_hartoffset,
        Urho_fft_hartProxy[thisInstance.proxyOffset]);
    Charm_getOutputExtents(myGrid_start[MY_X], myGrid_end[MY_X],
                          myGrid_start[MY_Y], myGrid_end[MY_Y],
                          myGrid_start[MY_Z], myGrid_end[MY_Z],
                          Urho_fft_hartProxy[thisInstance.proxyOffset],
                          fft_hartoffset);

    Charm_getAllPoints(myPoints, Urho_fft_hartProxy[thisInstance.proxyOffset],
        fft_hartoffset);

    for(int i = 0; i < 3; i++) {
      myGrid_length[i] = myGrid_end[i] - myGrid_start[i];
    }

    numPoints = (*myPoints).size();
    myGrid_size = myGrid_length[MY_X] * myGrid_length[MY_Y];
  }

  if(simReadOnly.ees_eext_on) {
    if(thisIndex.y == 0) {
      Charm_getOutputIndex(thisIndex.x, fft_atmSFTotOffset,
          Urho_fft_atmSFProxy[config.nchareHartAtmT][thisInstance.proxyOffset]);
    }

    Charm_getOutputIndex(thisIndex.x, fft_atmSFOffset,
        Urho_fft_atmSFProxy[thisIndex.y][thisInstance.proxyOffset]);

    Charm_getOutputExtents(myGrid_start_ext[MY_X], myGrid_end_ext[MY_X],
        myGrid_start_ext[MY_Y], myGrid_end_ext[MY_Y],
        myGrid_start_ext[MY_Z], myGrid_end_ext[MY_Z],
        Urho_fft_atmSFProxy[thisIndex.y][thisInstance.proxyOffset],
        fft_atmSFOffset);

    Charm_getAllPoints(myPoints_ext,
        Urho_fft_atmSFProxy[thisIndex.y][thisInstance.proxyOffset],
        fft_atmSFOffset);

    for(int i = 0; i < 3; i++) {
      myGrid_length_ext[i] = myGrid_end_ext[i] - myGrid_start_ext[i];
    }

    numPoints_ext = (*myPoints_ext).size();
    myGrid_size_ext = myGrid_length_ext[MY_X] * myGrid_length_ext[MY_Y];

    if(config.nchareHartAtmT == 1 || thisIndex.y == 1) {
      if(numPoints_ext != numPoints) {
        CkAbort("numPoints do not match for density and ext\n");
      }
#if _CP_DEBUG_RHOG_VERBOSE_
      std::vector< gridPoint > & r_myPoints = (*myPoints);
      std::vector< gridPoint > & r_myPoints_ext = (*myPoints_ext);
      for(int i = 0; i < numPoints; i++) {
        if(r_myPoints_ext[i].d1 != r_myPoints[i].d1 ||
            r_myPoints_ext[i].d2 != r_myPoints[i].d2 ||
            r_myPoints_ext[i].d3 != r_myPoints[i].d3) {
          CkAbort("points do not match for density and ext\n");
        }
      }
#endif
    }
  } else {
    numPoints_ext = numPoints;
  }

  if(simReadOnly.ees_eext_on == 1) {
    numPoints = numPoints_ext;
    Rho  = (complex *) fftw_malloc(numPoints * sizeof(complex));
    Vks  = (complex *) fftw_malloc(numPoints * sizeof(complex));
    atmSF          = (complex *)fftw_malloc(myGrid_size_ext * sizeof(complex));
    Charm_setOutputMemory((void*)atmSF,
        Urho_fft_atmSFProxy[thisIndex.y][thisInstance.proxyOffset],
        fft_atmSFOffset);
    Charm_createOutputPlan(
        Urho_fft_atmSFProxy[thisIndex.y][thisInstance.proxyOffset],
        fft_atmSFOffset);
    atmSFtot       = (complex *)fftw_malloc(numPoints_ext * sizeof(complex));
    if(thisIndex.y == 0) {
      atmSFtotRecv = (complex *)fftw_malloc(myGrid_size_ext * sizeof(complex));
      Charm_setOutputMemory((void*)atmSFtotRecv,
          Urho_fft_atmSFProxy[config.nchareHartAtmT][thisInstance.proxyOffset],
          fft_atmSFTotOffset);
      Charm_createOutputPlan(
          Urho_fft_atmSFProxy[config.nchareHartAtmT][thisInstance.proxyOffset],
          fft_atmSFTotOffset);
    }//endif
    if(config.nchareHartAtmT == 1 || thisIndex.y == 1) {
      VksRecv        = (complex *)fftw_malloc(myGrid_size * sizeof(complex));
      Charm_setOutputMemory((void*)VksRecv,
          Urho_fft_hartProxy[thisInstance.proxyOffset],
          fft_hartoffset);
      Charm_createOutputPlan(Urho_fft_hartProxy[thisInstance.proxyOffset],
          fft_hartoffset);
    }//endif
  } else {
    /*TODO This Rho may be unnecessary - inplace copy to Vks? Atleast, we don't
     * need myGrid_size Rho.*/
    Rho  = (complex *) fftw_malloc(myGrid_size * sizeof(complex));
    Vks  = (complex *) fftw_malloc(myGrid_size * sizeof(complex));
    Charm_setOutputMemory((void*)Vks,
        Urho_fft_hartProxy[thisInstance.proxyOffset],
        fft_hartoffset);
    Charm_createOutputPlan(Urho_fft_hartProxy[thisInstance.proxyOffset],
        fft_hartoffset);
  }

#if 0 && GLENN_PERIODIC_CORRECTION
  perdCorr = (double *) fftw_malloc(numPoints * sizeof(double));
  setput_nd_eext_corrs(numPoints, myPoints, perdCorr);
#endif

  if(simReadOnly.ees_eext_on == 1) {
    if(config.nchareHartAtmT > 1) {
      if(thisIndex.y == 0) {
        CkMulticastMgr *mg = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
        mg->resetSection(Ughart_sectionProxy[thisIndex.y][thisInstance.proxyOffset]);
      }
    }

    eesCache *eesData  = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
    eesGHart_index = eesData->registerCacheGHart(thisIndex.x, myPoints_ext);
    CkCallback cb(CkIndex_CP_Rho_GHartExt::registrationDone(),
                  UrhoGHartExtProxy[thisInstance.proxyOffset]);
    contribute(cb);
  }//endif
  //---------------------------------------------------------------------------
}//end routine
//============================================================================

//============================================================================
// Destructor
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_Rho_GHartExt::~CP_Rho_GHartExt(){

  if(simReadOnly.ees_eext_on == 1){
    fftw_free(atmSF);
    fftw_free(atmSFtot);
    if(config.nchareHartAtmT == 1 || thisIndex.y == 1) { fftw_free(VksRecv); }
    if(thisIndex.y == 0) { fftw_free(atmSFtotRecv); }
  }
  fftw_free((void*)Vks);
  fftw_free((void*)Rho);
}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * The density arrives from RhoGspace - ONCE a time step (iteration).
 * Invoke compute of eext and/or Hartree energy is you have all the stuff you need.
 */
//============================================================================
void CP_Rho_GHartExt::acceptData(RhoGHartMsg *msg){
  //============================================================================
  // Check the flow of control to see if we can use the data.

#ifdef _CP_GHART_VERBOSE_
  CkPrintf("[%d] Ghart %d Here in acceptData numaccept density %d\n", thisInstance.proxyOffset, thisIndex.x, numAcceptDensity);
#endif

  int cp_min_opt = simReadOnly.cp_min_opt;


  //RAZ: this is only for spin up instance:
  if(mySpinIndex!=0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("No Ubers for CP_Rho_GHartExt!!  Good-bye.\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

 

#ifdef _NAN_CHECK_
  for(int i = 0; i < msg->size; i++){
    CkAssert(isnan(msg->data[i].re) == 0);
    CkAssert(isnan(msg->data[i].im) == 0);
  }
#endif

  //============================================================================
  // Copy out the data and flip arrival flag

  CkAssert(numPoints == msg->size);
  numAcceptDensity++;
  if(numAcceptDensity==1){ // do normal stuff
    iteration++;
    if(cp_min_opt == 0){
      if(UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch()->iteration != iteration-1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Flow of Control Error in GHartExtVks : atoms slow %d %d\n",
          UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch()->iteration,iteration);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
      }//endif
    }//endif
    if(simReadOnly.ees_eext_on == 1) {
      for(int p = 0; p < numPoints; p++) {
	Rho[p] = msg->data[p];
      }
    } else {
      for(int p = 0; p < numPoints; p++) {
	Rho[(*myPoints)[p].offset] = msg->data[p];
      }
    }
  }
  else // add it to the one from the other spin
    {
      if(simReadOnly.ees_eext_on == 1) {
	for(int p = 0; p < numPoints; p++) {
	  Rho[p] += msg->data[p];
	}
      } else {
	for(int p = 0; p < numPoints; p++) {
	  Rho[(*myPoints)[p].offset] += msg->data[p];
	}
      }
    }

  // spin case needs both up and dn before it computes
  if( (cp_lsda==1 && numAcceptDensity==2) || (cp_lsda==0) ){ 
    operateOnData();
    numAcceptDensity=0;
  }//endif

#ifdef _CP_DEBUG_RHOHART
  char myFileName[100];
  sprintf(myFileName, "Rho_GHart%d-%d.out", thisIndex.x,thisIndex.y);
  std::vector< gridPoint > & dpoints = (*myPoints);
  FILE *fp = fopen(myFileName,"w");
  for (int i = 0; i < numPoints; i++){
    fprintf(fp," %d %d %d : %g %g\n", dpoints[i].d3, dpoints[i].d2, dpoints[i].d1,
        msg->data[i].re, msg->data[i].im);
  }//endfor
  fclose(fp);
#endif

  delete msg;
//============================================================================
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * LSDA: The density arrives from acceptData ONCE a time step 
 */
//============================================================================
void CP_Rho_GHartExt::operateOnData(){
//============================================================================
// Check the flow of control to see if we can use the data.
#ifdef _CP_GHART_VERBOSE_
  CkPrintf("{%d} In CP_Rho_GHartExt[%d,%d] operateOnData, Memory %.2lf MB\n",
      thisInstance.proxyOffset, thisIndex.x, thisIndex.y, CmiMemoryUsage()/(1024.0 * 1024));
#endif

  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  int cp_min_opt = sim->cp_min_opt;

//============================================================================
// Flip arrival flag

  densityHere = 1;

  //============================================================================
  // If ees is off, go ahead with the N^2 method.
  // If ees is on, either go ahead or chill depending on whether atmSF is here

  if(simReadOnly.ees_eext_on == 0) {
    HartExtVksG();
  }
  if(simReadOnly.ees_eext_on == 1 && atmSFHere == 1){
    getHartEextEes();
  }
  //============================================================================
}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * Compute hartree eext and vks using the N^2 method
 */
//============================================================================
void CP_Rho_GHartExt::HartExtVksG() {
  //============================================================================
  // Get the variables

#ifdef _CP_GHART_VERBOSE_
  CkPrintf("[%d] Ghart %d Here in HartExtVksG \n", CkMyPe(), thisIndex.x);
#endif

  // find me the local copy
  AtomsCache *ag       = UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
  int natm             = ag->natm;
  FastAtoms *fastAtoms = &(ag->fastAtoms);

  ehart_ret     = 0.0;
  eext_ret      = 0.0;
  ewd_ret       = 0.0;

  //============================================================================
  // compute vks(g) from hart eext and reduce eext and ehart

#if CMK_TRACE_ENABLED
  double  StartTime=CmiWallTimer();
#endif

  memset(Vks, 0, myGrid_size*sizeof(fftw_complex));
  CPLOCAL::CP_hart_eext_calc(numPoints, (*myPoints), Rho, natm, simReadOnly.iperd,
      fastAtoms, Vks, &ehart_ret, &eext_ret, &ewd_ret, thisIndex.x,
      UpScratchProxy[thisInstance.proxyOffset].ckLocalBranch()->psscratch, config.nfreq_cplocal_hartext);

#if CMK_TRACE_ENABLED
  traceUserBracketEvent(HartExcVksG_, StartTime, CmiWallTimer());
#endif

#ifdef _CP_DEBUG_RHOG_VKSA_
  char myFileName[100];
  sprintf(myFileName, "Vks_Gspace_%d%d.out", thisIndex.x, thisIndex.y);
  FILE *fp = fopen(myFileName,"w");
  for (int i = 0; i < numPoints; i++) {
    gridPoint &g = (*myPoints)[i];
    fprintf(fp," %d %d %d : %g %g\n",
        g.d3, g.d2, g.d1, Vks[g.offset].re, Vks[g.offset].im);
  }//endfor
  fclose(fp);
#endif

#ifdef _CP_GHART_VERBOSE_
  CkPrintf("[%d] Ghart %d HartExtVksG completed \n", CkMyPe(), thisIndex.x);
#endif

  //============================================================================
  // Reduce the energies computed then FFT Vks into Real space (part way)

  double e[3];
  e[0] = ehart_ret;
  e[1] = eext_ret;
  e[2] = ewd_ret;
  contribute(3 * sizeof(double),e,CkReduction::sum_double);

  FFTVks();

  //---------------------------------------------------------------------------
}//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * Partly fft vks(gx,gy,gz) -> vks(gx,gy,z) then invoke transpose
 */
//============================================================================
void CP_Rho_GHartExt::FFTVks() {
  //============================================================================
  // Perform the FFT(gx,gy,gz) to FFT(gx,gy,z)

  Charm_doBackwardFFT(CkCallback(CkIndex_CP_Rho_RealSpacePlane::acceptHartVks(),
        UrhoRealProxy[thisInstance.proxyOffset]),
        Urho_fft_hartProxy[thisInstance.proxyOffset], fft_hartoffset);
  //============================================================================
}//end routine
//============================================================================


/*** FOLLOWING CALLS ARE ONLY INVOKED FOR EES method ***/
/*** FOLLOWING CALLS ARE ONLY INVOKED FOR EES method ***/
/*** FOLLOWING CALLS ARE ONLY INVOKED FOR EES method ***/
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/**
 * Make sure everyone is registered in the Cache on the 1st time step
 */
//==========================================================================
void CP_Rho_GHartExt::registrationDone() {
  //==========================================================================
#ifdef _CP_GHART_VERBOSE_
  CkPrintf("[%d] I am Ghart chare %d %d in reg : natmtyp %d\n", CkMyPe(),
      thisIndex.x, thisIndex.y, natmTyp);
#endif

  registrationFlag = 1;
  if(launchFlag == 1) {
    launchFlag = 0;
    FFTEesBck();
  }
}//end routine
//==========================================================================

//============================================================================
/**
 *  Auxilary functions to mark the ending of FFT and informing the section
 */
//============================================================================

void CP_Rho_GHartExt::doneAtmSF_FFT() {
  if(thisIndex.x != 0) {
    CkAbort("Rho_GHartExt::doneAtmSF_FFT called on non-zero indexed chare\n");
  }
  FFT_Done_Msg* msg = new FFT_Done_Msg;
  Ughart_sectionProxy[thisIndex.y][thisInstance.proxyOffset].doneAtmSF_Multicast(msg);
}

void CP_Rho_GHartExt::doneAtmSF_Multicast(FFT_Done_Msg *msg) {
  delete msg;
  recvAtmSFFromRhoRHart();
}

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 *  FFT data from RhoRhart : Euler Exponential spline based method
 */
//============================================================================
void CP_Rho_GHartExt::recvAtmSFFromRhoRHart(){
  //============================================================================
  launchFlag   = 1;

#ifdef _CP_RHART_VERBOSE_DUMP
  if(iterAtmTyp == 0) {
    std::vector< gridPoint > & r_myPoints_ext = (*myPoints_ext);
    for(int i = 0; i < numPoints; i++) {
      CkPrintf("%d %d %d %lf %lf\n", r_myPoints_ext[i].d3,
      r_myPoints_ext[i].d2, r_myPoints_ext[i].d1,
      atmSF[r_myPoints_ext[i].offset].re,
      atmSF[r_myPoints_ext[i].offset].im);
    }
  }
#endif

  if(registrationFlag == 1){
    launchFlag = 0;
    thisProxy(thisIndex.x, thisIndex.y).FFTEesBck();
  }
}//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** Finish FFting to G-space  ::
 *         2D)  atmSF(gx,gy,z) -> atmSF(gx,gy,gz)
 */
//============================================================================
void CP_Rho_GHartExt::FFTEesBck(){
  //============================================================================

  atmSFHere = 1;
  iterAtmTyp ++;

  if(iterAtmTyp > 1 && densityHere == 0) {
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("{%d} Yo dawg, the density has not arrived. You can't keep going. GhartEext\n", thisInstance.proxyOffset);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  if(densityHere == 1) {
    getHartEextEes();
  }
  //----------------------------------------------------------------------------
}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * compute HartreeEextEes
 */
//============================================================================
void CP_Rho_GHartExt::getHartEextEes(){
  //============================================================================
  // Output and Error check

#ifdef _CP_GHART_VERBOSE_
  CkPrintf("[%d] GHart %d %d starting getHartEextEes at %d on %d\n", CkMyPe(),
      thisIndex.x, thisIndex.y, iterAtmTyp);
#endif

  //============================================================================
  // Compute eext energy, hartree, total SF and VKS

  //----------------------------------------------------------
  // Local Variables
  eesCache *eesData  = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
  int iterAtmTypFull = iterAtmTyp + atmTypoff;

  double *b_re = eesData->RhoGHartData[eesGHart_index].b_re;
  double *b_im = eesData->RhoGHartData[eesGHart_index].b_im;

  //----------------------------------------------------------
  // Initialize
  if(iterAtmTyp == 1) {
    memset(Vks, 0, numPoints * sizeof(complex));  // no getting around these zeros
    memset(atmSFtot, 0, numPoints_ext * sizeof(complex));
    ehart_ret = 0.0;
    eext_ret  = 0.0;
    ewd_ret   = 0.0;
  }//endif

  //----------------------------------------------------------
  // Get the energy, vks, modifiy atmSF, contribute to total SF
#if CMK_TRACE_ENABLED
  double StartTime=CmiWallTimer();
#endif
  double *perdCorr = NULL; //no periodic corrections yet
  CPLOCAL::eesHartEextGchare(numPoints_ext, iterAtmTypFull, Rho, Vks, atmSF,
      atmSFtot, b_re, b_im, &ehart_ret, &eext_ret, (*myPoints_ext),
      perdCorr, thisIndex.x, config.nfreq_cplocal_eeshart);
#if CMK_TRACE_ENABLED
  traceUserBracketEvent(eesHartExcG_, StartTime, CmiWallTimer());
#endif

#ifdef _CP_GHART_VERBOSE_
  CkPrintf("[%d] Ghart %d %d Here in hartees: EHart %.10g Eext %.10g at %d\n",
      CkMyPe(),thisIndex.x,thisIndex.y,ehart_ret,eext_ret,iterAtmTyp);
#endif

  //============================================================================
  // If you have SFtot get the ewald energy  : A reduction is required when atmTyp parallel

  if(iterAtmTyp == natmTyp) {
    if(config.nchareHartAtmT == 1) {
#if CMK_TRACE_ENABLED
      double StartTime=CmiWallTimer();
#endif
      double *perdCorr = NULL;
      CPLOCAL::eesEwaldGchare(numPoints_ext, atmSFtot, b_re, b_im, &ewd_ret,
          (*myPoints_ext), 0, perdCorr, thisIndex.x, config.nfreq_cplocal_eesewald);
#if CMK_TRACE_ENABLED
      traceUserBracketEvent(eesEwaldG_, StartTime, CmiWallTimer());
#endif
#ifdef _CP_GHART_VERBOSE_
      CkPrintf("[%d] Ghart %d %d iter %d Here in hartees: EwaldRecip : %.10g\n",
          CkMyPe(), thisIndex.x,thisIndex.y, iterAtmTyp, ewd_ret);
#endif
    }//endif

  //============================================================================
  // Blast out the energies when done with atom type stuff : index=0 may have to wait

#ifdef _CP_GHART_VERBOSE_
    CkPrintf("[%d] Ghart %d reduces energies at %d on %d\n", CkMyPe(),
        iterAtmTyp);
#endif
    if(thisIndex.y != 0 || config.nchareHartAtmT == 1){
      double e[3];
      e[0] = ehart_ret;
      e[1] = eext_ret;
      e[2] = ewd_ret;
      contribute(3*sizeof(double),e,CkReduction::sum_double);
    }//endif

    //============================================================================
    // Perform the back FFT SF and SFtot to get atm forces and VKS to get e-forces

    //-----------------------------------------------------------------
    // Ewald needs contribs from all SF  : chare index 0 is large and in charge
    if(config.nchareHartAtmT == 1) {
      //first copy atmSFtot data to expanded space
      memset(atmSFtotRecv, 0, myGrid_size_ext * sizeof(complex));
      for(int p = 0; p < numPoints_ext; p++) {
        atmSFtotRecv[(*myPoints_ext)[p].offset] = atmSFtot[p];
      } //then fft to RHart
      FFTEesFwd(1);
    } else {
      RhoGHartMsg *msg = new (numPoints_ext) RhoGHartMsg;
      CmiMemcpy(msg->data, atmSFtot, numPoints_ext * sizeof(complex));
      msg->size = numPoints_ext;
      UrhoGHartExtProxy[thisInstance.proxyOffset](thisIndex.x, 0).acceptAtmSFTot(msg);
    }//endif

    //-----------------------------------------------------------------
    // Vks needs contribs from all SF : chare index 1 is large and in charge
    if(config.nchareHartAtmT == 1)  {
      //first copy vks data to expanded space
      memset(VksRecv, 0, myGrid_size * sizeof(complex));
      for(int p = 0; p < numPoints; p++) {
        VksRecv[(*myPoints)[p].offset] = Vks[p];
      } //then fft to RhoR
      FFTVks();
    }else{
      RhoGHartMsg *msg = new (numPoints) RhoGHartMsg;
      CmiMemcpy(msg->data, Vks, numPoints * sizeof(complex));
      msg->size = numPoints;
      UrhoGHartExtProxy[thisInstance.proxyOffset](thisIndex.x, 1).acceptVks(msg);
    }//endif
  }

#ifdef  _CP_DEBUG_RHOHART_ATMSF
  char myFileName[100];
  sprintf(myFileName, "AtmSF_GHart_%d_%d_%d.out", thisIndex.x, thisIndex.y,
      iterAtmTyp);
  FILE *fp = fopen(myFileName,"w");
  for (int i = 0; i < numPoints; i++) {
    gridPoint &g = (*myPoints_ext)[i];
    fprintf(fp," %d %d %d : %g %g\n",
        g.d3, g.d2, g.d1, atmSF[g.offset].re, atmSF[g.offset].im);
  }//endfor
  fclose(fp);
#endif

  //zero out non sphere part
  complex zero;
  zero.re = 0.0; zero.im = 0.0;
  std::vector< gridPoint > & points = (*myPoints_ext);
  int last_offset = -1;
  for(int p = 0; p < numPoints_ext; p++) {
    int offset = points[p].offset;
    if(offset != (last_offset + 1)) {
      for(int cur_off = last_offset + 1; cur_off < offset; cur_off++) {
        atmSF[cur_off] = zero;
      }
    }
    last_offset = offset;
  }
  for(int cur_off = last_offset + 1; cur_off < myGrid_size_ext; cur_off++) {
    atmSF[cur_off] = zero;
  }

  //-----------------------------------------------------------------
  // we always generate an SF : we have everthing for this atmtype
  //                            : flow of control says ``this guy goes last''
  //                            : This guy controls exit condition.
  FFTEesFwd(0);                // DON'T MOVE HIM
  //============================================================================
}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * Start FFTing back to R-space for atmSF
 */
//============================================================================
void CP_Rho_GHartExt::FFTEesFwd(int flag){
  //============================================================================
  if(flag == 0) {
    //fft call for atmSF
    if(config.nchareHartAtmT > 1) {
      //use section
      Charm_doBackwardFFT(
          CkCallback(CkIndex_CP_Rho_RHartExt::doneAtmSF_FFT(),
            UrhoRHartExtProxy[thisInstance.proxyOffset](0,0,thisIndex.y)),
          Urho_fft_atmSFProxy[thisIndex.y][thisInstance.proxyOffset],
          fft_atmSFOffset);
    } else {
      Charm_doBackwardFFT(
          CkCallback(CkIndex_CP_Rho_RHartExt::recvAtmForcFromRhoGHart(),
            UrhoRHartExtProxy[thisInstance.proxyOffset]),
          Urho_fft_atmSFProxy[thisIndex.y][thisInstance.proxyOffset],
          fft_atmSFOffset);
    }
  } else {
    //fft call for atmSFtot
    if(config.nchareHartAtmT > 1) {
    Charm_doBackwardFFT(
        CkCallback(CkIndex_CP_Rho_RHartExt::doneAtmSFTot_FFT(),
          UrhoRHartExtProxy[thisInstance.proxyOffset](0,0,0)),
        Urho_fft_atmSFProxy[config.nchareHartAtmT][thisInstance.proxyOffset],
        fft_atmSFTotOffset);
    } else {
      Charm_doBackwardFFT(
          CkCallback(CkIndex_CP_Rho_RHartExt::recvAtmForcTotFromRhoGHart(),
            UrhoRHartExtProxy[thisInstance.proxyOffset]),
          Urho_fft_atmSFProxy[config.nchareHartAtmT][thisInstance.proxyOffset],
          fft_atmSFTotOffset);
    }
  }

  int nsendExpect = natmTyp;
  if(thisIndex.y == 0) {
    nsendExpect++;
  } // chare 0 sends out SFtot, too.

  nsendAtmTyp++;
  if(nsendAtmTyp == nsendExpect){
    nsendAtmTyp = 0;
    iterAtmTyp  = 0;
    atmSFHere   = 0;
    densityHere = 0;
  }//endif

}//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * Collect the SF from all the atm type chares on chare 0 to compute Ewald sum
 */
//============================================================================
void CP_Rho_GHartExt::acceptAtmSFTot(RhoGHartMsg *msg){
  //============================================================================
  // Recv the contribs.
  CkAssert(thisIndex.y == 0);

  if(countAtmSFtot == 0) {
    memset(atmSFtotRecv, 0, myGrid_size_ext * sizeof(complex));
  }

  countAtmSFtot++;

  if(msg->size != numPoints_ext) {
    CkAbort("Contributed atmSFtot does not match in length\n");
  }

  for(int i = 0; i < msg->size; i++) {
    atmSFtotRecv[(*myPoints_ext)[i].offset] += msg->data[i];
  }

  //============================================================================
  // Once you have it all, compute the energy, contribute, fft back, send to Rhart

  if(countAtmSFtot == config.nchareHartAtmT){
    countAtmSFtot = 0;
    //---------------------------------------------------------
    // Compute ewald energy, modify the atmSFtot appropriately
    eesCache *eesData  = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
    double *b_re = eesData->RhoGHartData[eesGHart_index].b_re;
    double *b_im = eesData->RhoGHartData[eesGHart_index].b_im;
#if CMK_TRACE_ENABLED
    double StartTime=CmiWallTimer();
#endif
    double *perdCorr = NULL;
    CPLOCAL::eesEwaldGchare(numPoints_ext, atmSFtotRecv, b_re, b_im, &ewd_ret,
      (*myPoints_ext), 1, perdCorr, thisIndex.x,config.nfreq_cplocal_eesewald);
#if CMK_TRACE_ENABLED
    traceUserBracketEvent(eesEwaldG_, StartTime, CmiWallTimer());
#endif

    //---------------------------------------------------------
    // Contribute your energies now that you have them all
    double e[3];
    e[0] = ehart_ret;
    e[1] = eext_ret;
    e[2] = ewd_ret;
    contribute(3*sizeof(double),e,CkReduction::sum_double);

    //---------------------------------------------------------
    // FFT back, which generates a send back to RHart
    FFTEesFwd(1);
  }//endif
  delete msg;
  //============================================================================
}//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * Reduce VKS contribs from all the eext atm type chares onto chare 1
 * Invoke the FFT back to real space when it is all here
 */
//============================================================================
void CP_Rho_GHartExt::acceptVks(RhoGHartMsg *msg) {
  //============================================================================
  // Recv the contributions
  CkAssert(thisIndex.y == 1);

  if(countVksTot == 0) {
    memset(VksRecv, 0, myGrid_size * sizeof(complex));
  }

  countVksTot++;

  if(msg->size != numPoints) {
    CkAbort("Contributed vks does not match in length\n");
  }

  for(int i = 0; i < msg->size; i++) {
    VksRecv[(*myPoints)[i].offset] += msg->data[i];
  }

  //============================================================================
  // When all the guys have reported, you can do the fft and then send vks to RhoReal

  if(countVksTot == config.nchareHartAtmT) {
    countVksTot=0;
    FFTVks();
  }//endif
  delete msg;
  //============================================================================
}//end routine
//============================================================================

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GHartExt::pup(PUP::er &p){
  //============================================================================
  ArrayElement2D::pup(p);
  p|registrationFlag;
  p|launchFlag;
  p|atmTypoff;
  p|natmTyp;
  p|iterAtmTyp;
  p|nsendAtmTyp;
  p|atmSFHere;
  p|densityHere;
  p|ehart_ret;
  p|eext_ret;
  p|ewd_ret;
  p|countAtmSFtot;
  p|countVksTot;
  p|iteration;
  //TODO: this is still incomplete
  //---------------------------------------------------------------------------
}//endif
//============================================================================

//============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * Glenn's special exit for debugging - convenient to have around and not rewrite.
 */
//============================================================================
void CP_Rho_GHartExt::exitForDebugging(){
  //============================================================================
  CkPrintf("I am in the exitfordebuging rhoghartext puppy. Bye-bye\n");
  CkExit();
  //============================================================================
}//end routine
//============================================================================
/*@}*/
