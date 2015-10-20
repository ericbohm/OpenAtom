#include "charm++.h"
#include "utility/util.h"
#include "cpaimd.h"
#include "AtomsCache.h"
#include "AtomsCompute.h"
#include "InstanceController.h"

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/* ostensibly the InstanceController may need to know about everything */
extern CkVec < CkVec <int> > UplaneUsedByNLZ;
extern CProxy_TemperController temperControllerProxy;
extern CkVec <CProxy_CP_State_GSpacePlane>       UgSpacePlaneProxy;
extern CkVec <CProxy_GSpaceDriver>               UgSpaceDriverProxy;
extern CkVec <CProxy_CP_State_ParticlePlane>     UparticlePlaneProxy;
extern CkVec <CProxy_CP_State_RealParticlePlane> UrealParticlePlaneProxy;
extern CkVec <CProxy_CP_State_RealSpacePlane>    UrealSpacePlaneProxy;
extern CkVec <CProxy_CP_Rho_RealSpacePlane>      UrhoRealProxy;
extern CkVec <CProxy_CP_Rho_GSpacePlane>         UrhoGProxy;
extern CkVec <CProxy_CP_Rho_RHartExt>            UrhoRHartExtProxy;
extern CkVec <CProxy_CP_Rho_GHartExt>            UrhoGHartExtProxy;
extern CkVec <CProxy_AtomsCache>                 UatomsCacheProxy;
extern CkVec <CProxy_AtomsCompute>               UatomsComputeProxy;
extern CkVec <CProxy_EnergyGroup>                UegroupProxy;
extern CkVec <CProxy_FFTcache>                   UfftCacheProxy;
extern CkVec <CProxy_StructFactCache>            UsfCacheProxy;
extern CkVec <CProxy_StructureFactor>            UsfCompProxy;
extern CkVec <CProxy_eesCache>                   UeesCacheProxy;
extern CkVec <CProxy_OrthoHelper>                UorthoHelperProxy;
extern CProxy_TimeKeeper                         TimeKeeperProxy;
#ifdef CMK_BALANCED_INJECTION_API
extern CProxy_PlatformSpecific                   platformSpecificProxy;
extern uint16_t                                  origBIValue;
#endif
extern CkGroupID                                 mCastGrpId;
extern CkVec <UberCollection>                    UberAlles;
extern CProxy_ENL_EKE_Collector                  ENLEKECollectorProxy;
extern CPcharmParaInfo                          simReadOnly;
//============================================================================

InstanceController::InstanceController(int _fft_expected) {

  done_init=0;Timer=CmiWallTimer(); numKpointforces=0;
  fft_expected = _fft_expected;
  done_fft_creation = false;
  UberCollection instance=UberCollection(thisIndex);
  // 0th k, 0th spin makes this to lockdown everyone so the atoms
  // shared across all k and spin can start sanely
  if((config.UberMmax >1 || config.UberJmax>1) && instance.idxU.y==0 && instance.idxU.s==0)
  {
    // make section for k-points and spins
    CkVec <CkArrayIndex1D> elems;
    for(int kp =0; kp<config.UberJmax; kp++)
    {
      instance.idxU.y=kp;
      for(int spin =0; spin<config.UberMmax; spin++)
      {
        instance.idxU.s=spin;
        instance.setPO();
        elems.push_back(CkArrayIndex1D(instance.proxyOffset));
      }
    }

    CProxySection_InstanceController sectProxy=CProxySection_InstanceController::ckNew(thisProxy.ckGetArrayID(),elems.getVec(),elems.size());
    CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch(); 
    sectProxy.ckSectionDelegate(mcastGrp);
    ICCookieMsg *cookieme=new ICCookieMsg;
    sectProxy.initCookie(cookieme);
  }
}

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void InstanceController::init(){
  UberCollection instance=UberCollection(thisIndex);
  // 0th bead, 0th temper makes this to sync beads and tempers
  if((config.UberImax >1 || config.UberKmax>1) && instance.idxU.x==0 && instance.idxU.z==0)
  {
    // make section for beads and tempers
    int numDestinations=config.UberImax*config.UberKmax;
    CkArrayID *beadArrayIds= new CkArrayID[numDestinations];
    CkArrayIndex **elems  = new CkArrayIndex*[numDestinations];
    int *naelems = new int[numDestinations];
    for(int bead =0; bead<config.UberImax; bead++)
    {
      instance.idxU.x=bead;
      for(int temper =0; temper<config.UberKmax; temper++)
      {
        int index=bead*config.UberKmax+temper;
        elems[index]= new CkArrayIndex2D[1];
        naelems[index]=1;
        instance.idxU.z=temper;
        instance.setPO();
        CkPrintf("fmag sync section adding bead %d temper %d index %d proxyOffset %d\n",bead, temper, index, instance.proxyOffset);
        beadArrayIds[index]=UgSpacePlaneProxy[instance.proxyOffset].ckGetArrayID();
        elems[index][0]=CkArrayIndex2D(0,0);
      }
    }
    //finish setting this up      
    gTemperBeadProxy=CProxySection_CP_State_GSpacePlane(numDestinations, beadArrayIds, elems, naelems);
    CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch(); 
    gTemperBeadProxy.ckSectionDelegate(mcastGrp);
    ICCookieMsg *cookieme=new ICCookieMsg;
    CkCallback *cb = new CkCallback(CkIndex_InstanceController::fmagMinTest(NULL),CkArrayIndex1D(0),thisProxy);
    mcastGrp->setReductionClient(gTemperBeadProxy,cb);
    gTemperBeadProxy.initBeadCookie(cookieme);

  }
}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void InstanceController::fmagMinTest(CkReductionMsg *m){


  int result=( (int *) (m->getData()) )[0];
  // CkPrintf("[%d] fmagMinTest %d\n",thisIndex, result);
  ICCookieMsg *out=new ICCookieMsg;
  out->junk = result;
  delete m;
  gTemperBeadProxy.minimizeSync(out);

}


void InstanceController::instancesReady(CkReductionMsg *m){
  
  int result=( (int *) (m->getData()) )[0];
  delete m;
  CkPrintf("Starting all %d instances\n", result);
  for(int i=0;i<config.numInstances; i++)  UgSpaceDriverProxy[i].startControl();  
#ifdef CMK_BALANCED_INJECTION_API
  CkPrintf("Reseting Balanced Injection to %d \n", origBIValue);
  platformSpecificProxy.reset_BI();
#endif
}

void InstanceController::doneFFTCreation(idMsg *msg) {
  delete msg;
  fft_expected--;
  CkPrintf("Received FFT completion %d\n", fft_expected);
  if(fft_expected == 0) {
    done_fft_creation = true;
    if(done_init > 3) {
      initDensity();
    }
  }
}

void InstanceController::initDensity() {
  UrhoRealProxy[thisIndex].init();
  UrhoGProxy[thisIndex].init();
  UrhoGHartExtProxy[thisIndex].init();
  if(simReadOnly.ees_eext_on) { 
    UrhoRHartExtProxy[thisIndex].init();
  }
}

//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \addtogroup startup
# Parallel Startup Phases:

1.  Phase 0 is kicked off by reaching the end of main. This turns
execution completely over to the Charm++ scheduler, at which point it
will process the object constructor messages that we triggered in the
proxy creation in main.  Object construction will occur
in this order + all readonlies will be initialized.  + all groups
will be constructed + all arrays will be constructed

The ordering within those phases is non-deterministic, so we don't expect to
have control over the ordering of chare array construction.  The upshot of
this is that in order to safely make array sections we wait until the objects
are constructed and then call a second phase of initialization.  In practice
this means that arrays will contribute to reductions during construction and
the completion of those reductions will trigger a chain of section creation
which will eventually feed back into a reduction that reports to the global
startup phase ordering in InstanceController.  

2. Phase 2 and 3 are automatically triggered during the construction process.
These phases are ortho constructing proxies to sections of the
paircalculators.  In each case they construct a section and send a message on
that section to its elements to initialize a cookie.  Receipt of that cookie
increments a counter and when each PC element has received all the cookies it
expects, it contributes to a reduction which reports to
InstanceController::doneInit().  There is a phase for symmetric and asymmetric
calculator, they could complete in either order.

3.  Phase 4 triggers the post construction initialization of section proxies
and cache registrations in RhoReal RhoG RhoGHartExt.  The big ticket item here
is the sections of RealSpace made by RhoReal.  These operate in the previously
described fashion wherein you make a section, initialize the cookies with a
dummy message and report on completion via a reduction along the section. When
realspace has received as many cookies as there are rhoreal subplanes, it
contributes to a reduction reporting to InstanceController::doneInit.  

4.  Phase 5 is triggered by the completion of the RS sections. When EES is
enabled, phase 5 will launch the section construction and registration process
in RealParticlePlane.  The coalesced completion of eesCache, enlSection, and
planeRedSection initialization contributes to a single reduction reporting to
InstanceController::doneInit.  This phase always triggers the loading and
multicasting of the gspace state data from the statefiles. When all elements
of gspace are initialized with that data they contribute to a reduction which
reports to InstanceController::doneInit.  

5.  Phase 6 happens only if EES is enabled, it broadcasts registrationDone to
all RealParticlePlane elements.  

6.  Phase 7 (or 6 if no realparticleplane) means that all initialization is
complete and startup is effectively over.  Control is then turned over to the
gSpaceDriver::startControl.  Some chares will do a little local first
iteration initialization after this.  Semantically it should now be safe to
engage in any operation as the previous phases should have taken care of any
synchronized initialization issues.

7. Phase 8( or 7) if more than one instance is in play, then we barrier for
all intances to be ready.  If on Cray, we then reset the BALANCED_INJECTION to
a normal behavior value.

 */
/**@{*/
void InstanceController::doneInit(){
  CPcharmParaInfo *sim  = CPcharmParaInfo::get();
  CkPrintf("{%d} Done_init \n",thisIndex);
  // Also, when paircalc becomes completely instance unaware, it will fail. This single assert is not enough motivation
  // to provide instance info to the pc/ortho bubble. @todo: remove this assert
  int numPhases=5;
  if(CPcharmParaInfo::get()->ees_nloc_on==1)
    numPhases++;
  double newtime=CmiWallTimer();
  CkAssert(done_init<numPhases+1);

  if(done_init<numPhases){
    CkPrintf("{%d} Completed chare instantiation phase %d in %g\n",thisIndex,done_init+1,newtime-Timer);
  }
  if (done_init==3)
  { // kick off post constructor inits
    if(thisIndex==0) init();
    UberCollection thisInstance(thisIndex);
    if(thisInstance.idxU.y==0 && done_fft_creation) {
      initDensity();
    }
  }
  if (done_init == 4){
    // We do this after we know gsp, pp, rp, rpp exist
    //if(CPcharmParaInfo::get()->ees_nloc_on==1)
    //{UrealParticlePlaneProxy[thisIndex].init();}
    // kick off file reading in gspace
    CkPrintf("{%d} Initiating import of states\n",thisIndex);
    CkPrintf("{%d} IC uGSpacePlaneProxy[%d] is %d\n",thisIndex,thisIndex, CkGroupID(UgSpacePlaneProxy[thisIndex].ckGetArrayID()).idx);
    for(int s=0;s<simReadOnly.nstates;s++) {
      UgSpacePlaneProxy[thisIndex](s,UplaneUsedByNLZ[thisIndex][s]).readFile();
    } //endfor

  }//endif
  if (done_init == 5 && CPcharmParaInfo::get()->ees_nloc_on==1){
    CkPrintf("{%d} Completed chare data acquisition phase %d in %g\n",thisIndex, done_init+1,newtime-Timer);
    UrealParticlePlaneProxy[thisIndex].init();
  }

  if (done_init >= numPhases) {
    if (done_init == numPhases){ 
      //          PRINT_LINE_STAR;
      CkPrintf("{%d} Chare array launch and initialization complete       \n",thisIndex);
      if(CPcharmParaInfo::get()->cp_min_opt==1 && CPcharmParaInfo::get()->cp_bomd_opt==0){
        CkPrintf("{%d} Running Open Atom CP Minimization: \n",thisIndex);
      }
      if(CPcharmParaInfo::get()->cp_bomd_opt==1){
        CkPrintf("{%d} Running Open Atom CP BOMD Dynamics: \n", thisIndex);
      }
      if (CPcharmParaInfo::get()->cp_min_opt==0) {
        CkPrintf("{%d} Running Open Atom CP Dynamics: \n",thisIndex);
      }//endif
      //          PRINT_LINE_STAR; CkPrintf("\n");
      //          PRINT_LINE_STAR;

      if (CPcharmParaInfo::get()->ees_nloc_on==1){
        UrealParticlePlaneProxy[thisIndex].registrationDone();
      }

      if(config.numInstances>1) 
       {
         /* 
            barrier for all instances to be ready.  Not semantically
            required, but does make output less weird and avoid some startup
            related issues on some platforms.
         */
         int num=1;
         contribute(sizeof(int), &num, CkReduction::sum_int, 
                    CkCallback(CkIndex_InstanceController::instancesReady(NULL),CkArrayIndex1D(0),thisProxy), 0);
         
       }
      else
       {
         UgSpaceDriverProxy[thisIndex].startControl();
#ifdef CMK_BALANCED_INJECTION_API
         CkPrintf("Reseting Balanced Injection to %d\n", origBIValue);
         if(thisIndex==0)  platformSpecificProxy.reset_BI();
#endif
       }
    }//endif
  }
  Timer=newtime;
  ++done_init;
}
//============================================================================
/**@}*/

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void InstanceController::initCookie(ICCookieMsg *msg){
  CkGetSectionInfo(allKPcookie, msg);
  //    delete msg; nokeep
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void InstanceController::printEnergyHart(CkReductionMsg *msg){
  //============================================================================
  //  double ehart = 0, eext = 0.0, ewd = 0.0;
  void *data=msg->getData();
  double ehart = ((double *)data)[0];
  double eext = ((double *)data)[1];
  double ewd  = ((double *)data)[2];

  UgSpacePlaneProxy[thisIndex](0, 0).computeEnergies(ENERGY_EHART, ehart);
  UgSpacePlaneProxy[thisIndex](0, 0).computeEnergies(ENERGY_EEXT, eext);
  UgSpacePlaneProxy[thisIndex](0, 0).computeEnergies(ENERGY_EWD, ewd);
  delete msg;

  //============================================================================
}//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void InstanceController::printEnergyEexc(CkReductionMsg *msg){
  //============================================================================

  double eexc = 0;
  double egga = 0;
  void *data=msg->getData();
  eexc += ((double *)data)[0];
  egga += ((double *)data)[1];

  UgSpacePlaneProxy[thisIndex](0, 0).computeEnergies(ENERGY_EEXC, eexc);
  UgSpacePlaneProxy[thisIndex](0, 0).computeEnergies(ENERGY_EGGA, egga);

#define _GLENN_STUFF_OFF_
#ifdef _GLENN_STUFF_
  CkPrintf("exiting in printEnergyEexc\n");CkExit();
#endif
  delete msg;

  //============================================================================
}//end routine
//============================================================================

//============================================================================
// Print out Quantum KE and put all the energies into the message group 
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void InstanceController::printEnergyEke(CkReductionMsg *m){
  //============================================================================  
  double d = ((double *)m->getData())[0];
  delete m;
#ifdef _CP_DEBUG_SFNL_OFF_
  CkPrintf("EKE         = OFF FOR DEBUGGING\n");
#endif
  UberCollection thisInstance(thisIndex);
  if(config.UberKmax>1) // report to temper master
    ENLEKECollectorProxy[thisInstance.idxU.z].acceptEKE(d);
  UgSpacePlaneProxy[thisIndex](0,0).computeEnergies(ENERGY_EKE, d);
}



//============================================================================
// Print out Fict CP KE and send it to the energy group
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void InstanceController::printFictEke(CkReductionMsg *m){
  //============================================================================  
  CPcharmParaInfo *sim  = CPcharmParaInfo::get();
  double d0   = ((double *)m->getData())[0];
  double d1   = ((double *)m->getData())[1];
  double d2   = ((double *)m->getData())[2];
  double d3   = ((double *)m->getData())[3];
  double d4   = ((double *)m->getData())[4];
  delete m;

  if(CPcharmParaInfo::get()->cp_min_opt==0 && printToScreen){
    FILE *temperScreenFile = UatomsCacheProxy[thisIndex].ckLocalBranch()->temperScreenFile;
    int iteration= UatomsCacheProxy[thisIndex].ckLocalBranch()->iteration;

    fprintf(temperScreenFile,"InstanceController printing Fict stats computed in iteration %d\n",iteration);
    fprintf(temperScreenFile,"Iter [%d] Fict Temp   =  %.10g K\n", iteration, d0); // per g-chare temp
    fprintf(temperScreenFile,"Iter [%d] Fict Eke    =  %.10g K\n", iteration, d2); // total kinetic energy
    fprintf(temperScreenFile,"Iter [%d] Fict TempNHC=  %.10g K\n", iteration, d1); // per g-chare tempNHC
    fprintf(temperScreenFile,"Iter [%d] Fict EkeNHC =  %.10g K\n", iteration, d3); // total NHC kinetic energy
    fprintf(temperScreenFile,"Iter [%d] Fict PotNHC =  %.10g K\n", iteration, d4); // total potNHC
    fprintf(temperScreenFile,"Iter [%d] Fict EConv  =  %.10g K\n", iteration, d2+d3+d4);
  }//endif
  UgSpacePlaneProxy[thisIndex](0,0).computeEnergies(ENERGY_FICTEKE, d0);  

  //============================================================================
}//end routine
//============================================================================


//============================================================================
//  When ALL the cp forces are done, you can integrate the atoms.
//  Could be made slicker by using knowledge of which gpsace chares
//  are on each processor and coordinating that with atom launch.
//  cklocalbranch->count_psi_forc++; zeroed in atms grp constructor
//  and incremented could be incremented in launch atoms.
//  If(cklocalbranch->count_psi_forc==numgspacethisproc){
//    atomsgrpproxy(myid).startrealspaceforces();
//  }
//  That is, you only need to know that all the forces are complete
//  on each PROCECESSOR and you can go ahead to the atoms. The atoms
//  are a group which needs to contribute atom forces to a reduction.
//  If all the forces aren't in there, you can't start the reduction
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void InstanceController::allDoneCPForces(int tol_reached){
  //============================================================================
  // only the 0th instance of each k-point and spin is allowed to do this 

  CPcharmParaInfo *sim  = CPcharmParaInfo::get();
  printToScreen = (!sim->cp_bomd_opt || tol_reached);
  UberCollection thisInstance(thisIndex);
  if(config.UberJmax>1 || config.UberMmax >1 ) 
  {
    // contribute to the section which includes all k-points of this
    // bead
    UberCollection instance=thisInstance;
    instance.idxU.y=0;
    int offset=instance.calcPO();

    CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch(); 
    CkCallback cb(CkReductionTarget(InstanceController, allDoneCPForcesAllKPoint),CkArrayIndex1D(offset),thisProxy);
    mcastGrp->contribute(sizeof(int), &tol_reached, CkReduction::min_int, allKPcookie, cb);
  }
  else
  {
    FILE *temperScreenFile = UatomsCacheProxy[thisIndex].ckLocalBranch()->temperScreenFile;
    int iteration= UatomsCacheProxy[thisIndex].ckLocalBranch()->iteration;

    // The atoms cache iteration gets incremented at the end of iterations, but
    // the GSpace iteration is incremented at the beginning. Add 1 to make
    // output from the two consistent.
    if (printToScreen) {
      fprintf(temperScreenFile,"Iter [%d] allDoneCPForces bead %d\n",iteration+1,thisInstance.idxU.x);  
    }
    UatomsComputeProxy[thisIndex].startRealSpaceForces(tol_reached);
  }

  //============================================================================
}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void InstanceController::allDoneCPForcesAllKPoint(int tol_reached){
  UberCollection thisInstance(thisIndex);
  FILE *temperScreenFile = UatomsCacheProxy[thisIndex].ckLocalBranch()->temperScreenFile;
  int iteration= UatomsCacheProxy[thisIndex].ckLocalBranch()->iteration;

  // The atoms cache iteration gets incremented at the end of iterations, but
  // the GSpace iteration is incremented at the beginning. Add 1 to make
  // output from the two consistent.
  if (printToScreen) {
    fprintf(temperScreenFile,"Iter [%d] allDoneCPForces bead %d\n",iteration,thisInstance.idxU.x);  
  }
  UatomsComputeProxy[thisIndex].startRealSpaceForces(tol_reached);
}

//============================================================================

// When the simulation is done on each instance, make a clean exit  
// this gets called by each instance
void InstanceController::cleanExit(CkReductionMsg *m)
{
  FILE *temperScreenFile = UatomsCacheProxy[thisIndex].ckLocalBranch()->temperScreenFile;
  int iteration= UatomsCacheProxy[thisIndex].ckLocalBranch()->iteration;
  fprintf(temperScreenFile,"Iter [%d] called cleanExit\n",iteration);
  delete m;
  int exited=1;
  contribute(sizeof(int), &exited, CkReduction::sum_int, 
      CkCallback(CkIndex_InstanceController::cleanExitAll(NULL),CkArrayIndex1D(0),thisProxy), 0);
}


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void InstanceController::acceptNewTemperature(double temperature){
  // you are the instance controller for [temper,0,0.0] tell the rest of your
  // minions about this new temperature
  // NOTE: this should be done using a section
  //  CkPrintf("[%d] acceptNewTemperature\n",thisIndex);
  UberCollection anIndex(thisIndex);
  for(int integral=0; integral< config.UberImax; integral++)
  {
    anIndex.idxU.x=integral;
    for(int kpoint=0; kpoint< config.UberJmax; kpoint++)
    {
      anIndex.idxU.y=kpoint;
      for(int spin=0; spin< config.UberMmax; spin++) {
        anIndex.idxU.s=spin;
        anIndex.setPO();
        thisProxy[anIndex.proxyOffset].useNewTemperature(temperature);
      }
    }
  }

}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void InstanceController::useNewTemperature(double temperature){
  // broadcast the temp to your atoms and GSPs
  //  CkPrintf("[%d] useNewTemperature\n",thisIndex);
  atomsTempDone=false;
  gspTempDone=false;
  UatomsComputeProxy[thisIndex].acceptNewTemperature(temperature);
  UgSpacePlaneProxy[thisIndex].acceptNewTemperature(temperature);

}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void InstanceController::atomsDoneNewTemp(CkReductionMsg *m)
{
  atomsTempDone=true;
#define TEMPERBARRIER 1
#ifndef TEMPERBARRIER  
  if(gspTempDone)
    resumeFromTemper();
#else
  contribute(CkCallback(CkReductionTarget(TemperController,barrier),temperControllerProxy));
#endif
  delete m;
}

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void InstanceController::gspDoneNewTemp(CkReductionMsg *m)
{
  gspTempDone=true;
#ifndef TEMPERBARRIER  
  if(atomsTempDone)
    resumeFromTemper();
#else
  contribute(CkCallback(CkReductionTarget(TemperController,barrier),temperControllerProxy));
#endif
  delete m;
}

//in a nicer world this would be inlined
void InstanceController::resumeFromTemper()
{ 
  UegroupProxy[thisIndex].resumeFromTemper(); 
}

// When the simulation is done, make a clean exit  
// this gets called on the 0th element when everyone calls cleanExit
void InstanceController::cleanExitAll(CkReductionMsg *m)
{
  delete m;  
  CkPrintf("********************************************************************************\n");
  CkPrintf("\n"); CkPrintf("\n");
  CkPrintf("********************************************************************************\n");
  CkPrintf("         Open Atom Simulation Complete                \n");
  CkPrintf("********************************************************************************\n");
  CkExit();
}

#include "instanceController.def.h"

