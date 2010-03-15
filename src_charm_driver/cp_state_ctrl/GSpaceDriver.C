#include "GSpaceDriver.h"
#include "CP_State_GSpacePlane.h"
#include "CP_State_ParticlePlane.h"
#include "GSpaceRTH.h"
#include "main/TimeKeeper.h"

extern int nstates;
extern Config config;

extern CProxy_CPcharmParaInfoGrp 		scProxy;
extern CProxy_TimeKeeper 			TimeKeeperProxy;
extern CkVec <CProxy_CP_State_GSpacePlane>      UgSpacePlaneProxy;
extern CkVec <CProxy_CP_State_ParticlePlane> 	UparticlePlaneProxy;
extern CkVec <CProxy_StructureFactor> 		UsfCompProxy;
extern CkVec <CProxy_Ortho>                     UorthoProxy;


/// Constructor
GSpaceDriver::GSpaceDriver(const UberCollection _thisInstance): 
			thisInstance(_thisInstance),
			myGSpaceObj(0),
			paraInfo(0),
			isFirstStep(true),
			waitingForEnergy(false),
			waitingForAtoms(false),
			isAtomIntegrationDone(false),
			isEnergyReductionDone(false),
			areNLForcesDone(false),
			isPsiVupdateNeeded(false),
			controlThread(0)
{
	#ifdef DEBUG_CP_GSPACE_CREATION
		CkPrintf("GSpaceDriver[%d,%d] born\n",thisIndex.x,thisIndex.y);
	#endif
	/// Initialize flags and counters that record the control state
	paraInfo = scProxy.ckLocalBranch ()->cpcharmParaInfo;
	ees_nonlocal = paraInfo->ees_nloc_on;
        int natm_nl  = paraInfo->natm_nl;
        if(natm_nl==0){
	  areNLForcesDone=true;
        }//endif
}//end routine




///
GSpaceDriver::GSpaceDriver(CkMigrateMessage *msg): 
			myGSpaceObj(0),
			paraInfo(0),
			controlThread(0)
{}




/// PUP method
void GSpaceDriver::pup(PUP::er &p)
{
	p|isFirstStep;
	p|ees_nonlocal;
        p|areNLForcesDone;
	p|waitingForEnergy;
	p|waitingForAtoms;
	p|isAtomIntegrationDone;
	p|isEnergyReductionDone;
        p|isPsiVupdateNeeded;
	p|sfCompSectionProxy;
	if( p.isUnpacking() )
	{
		controlThread = RTH_Runtime_create(RTH_Routine_lookup(GSpaceDriver,driveGSpace),this);
		paraInfo = scProxy.ckLocalBranch ()->cpcharmParaInfo;
	}
	RTH_Runtime_pup(controlThread,p,this);
}




/** Creates a structure factor array section proxy that is used later in releaseSFComputeZ. 
 * @todo: When is the earliest that you can do this? Can we do this upon construction? Probably not, as we have no guarantee that
 * UsfCompProxy is up and ready.
 */ 
void GSpaceDriver::init()
{
	/// If EES nonlocals are off, each state=0 driver chare creates a structure factor array section
	if(thisIndex.x==0 && ees_nonlocal==0)
	{
		/// numDups must be less than the number of states because thats the maximum number of times you can duplicate a plane
		int numSfDups = config.numSfDups;
		if(numSfDups > nstates)
			numSfDups = nstates;
		/// Create a list of SF chares in this section
		CkVec <CkArrayIndex3D> sfelems;
		for(int dup=0; dup<numSfDups; dup++)
			for(int atm=0;atm<config.numSfGrps; atm++)
				sfelems.push_back(CkArrayIndex3D(atm,thisIndex.y,dup));
		/// Create the SF array section
		sfCompSectionProxy = CProxySection_StructureFactor::ckNew(	UsfCompProxy[thisInstance.proxyOffset].ckGetArrayID(),
																	(CkArrayIndexMax*)sfelems.getVec(), sfelems.size());
	} 
}




/** Called from InstanceController::doneInit which is a reduction client invoked by all (state,plane) 
 * objects after completion of CP_State_GspacePlane::initGSpace
 */
void GSpaceDriver::startControl()
{
    /// Get hold of the (local) GSpacePlane and ParticlePlane objects that we'll be working with
	myGSpaceObj = UgSpacePlaneProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).ckLocal();
	CkAssert(myGSpaceObj);
    myParticlePlaneObj = UparticlePlaneProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).ckLocal();
    CkAssert(myParticlePlaneObj);
    /// Do other initialization chores
	init();
    /// Create an RTH thread and give it control
	controlThread = RTH_Runtime_create(RTH_Routine_lookup(GSpaceDriver,driveGSpace),this);
	resumeControl();
}




// Forward declaration
void cleanExit(void *param, void *msg);

/// GSpace notifies me that its ready to exit by calling this method
void GSpaceDriver::readyToExit()
{
	/// Set a flag that indicates this chare is ready to exit
	myGSpaceObj->cleanExitCalled = 1;
	/// Check if this chare is not waiting for any steps in this iteration before hitting an exit reduction barrier
	if(isAtomIntegrationDone && isEnergyReductionDone)
	{
		#ifdef _CP_SUBSTEP_TIMING_
			#ifdef USE_HPM
				(TimeKeeperProxy.ckLocalBranch())->printHPM();
			#endif
		#endif
		int i=0;
		contribute(sizeof(int),&i,CkReduction::sum_int,CkCallback(cleanExit,NULL));
	}
}




/// GSpace notifies me that the energy reduction is done by calling this method
void GSpaceDriver::doneComputingEnergy(const int AtomsGrpIter) 				
{
	/// Ensure the iterations are synced
	if (myGSpaceObj->iteration != AtomsGrpIter)
        CkAbort("GSpaceDriver::doneComputingEnergy: GSpace iteration and Atoms group iteration number are not in sync. Aborting...");
	///
	isEnergyReductionDone = true;
	/// If GSpace has already called for an exit, check if we can exit again
	if (myGSpaceObj->cleanExitCalled==1)
		readyToExit();
	/// If we were waiting for the energy (and the atoms have been moved) resume the driver logic
	if (waitingForEnergy) 
	{
		waitingForEnergy = false; 
		if (!waitingForAtoms) resumeControl(); 
	} 
}





/** Probe for atom completion : Could have atom group invoke this function directly on all chares for which 
 * it is responsible. There is an atom group on each proc. Some number of driver chares reside on each proc. 
 * This information is known by main and could be put into atom constructor. The atom group
 * would then know which chares upon which to invoke this method. Careful, careful with migration with this 
 * alternative scheme. Enable config.localAtomBarrier to trigger that behavior
 */
void GSpaceDriver::doneMovingAtoms(const int AtomsGrpIter)
{
  /// Ensure the iterations are synced @todo: Should this be an if condition? It was when it lived in GSpace

  CkAssert(myGSpaceObj->iteration == AtomsGrpIter);

	///
	isAtomIntegrationDone = true;
	/// If GSpace has already called for an exit, check if we can exit again
	if (myGSpaceObj->cleanExitCalled==1)
		readyToExit();
	/// If we were waiting for the atom integration (and the energy computation was done) resume the driver logic
	if (waitingForAtoms) 
	{
		waitingForAtoms = false; 
		if (!waitingForEnergy) resumeControl(); 
	} 
}



/// GSpace notifies me when the non-local results have arrived
void GSpaceDriver::doneNLForces() 					
{
    areNLForcesDone = true;
    if (myGSpaceObj->doneDoingIFFT)
        resumeControl();
}




/// All GSpace objects have finished psi : For debugging only
void GSpaceDriver::allDonePsi(CkReductionMsg *msg)
{
    delete msg;
    resumeControl();
}




/// All GSpace objects have finished writing coefs : NECESSARY
void GSpaceDriver::allDoneWritingPsi(CkReductionMsg *msg)
{
    delete msg;
    resumeControl();
}




/// All GSpace objects have finished velocity rotation : For debugging only
void GSpaceDriver::allDonePsiV(CkReductionMsg *msg)
{
    #ifdef DEBUG_CP_GSPACE_PSIV
        CkPrintf("GSpaceDriver[%d,%d] allDonePsiV: PsiV update step complete in iteration %d. Barrier reduction reached.\n",thisIndex.x,thisIndex.y,myGSpaceObj->iteration);
    #endif
    delete msg;
    resumeControl();
}




/// All GSpace objects finished Inverse FFT : For debugging only
void GSpaceDriver::allDoneIFFT(CkReductionMsg *msg)
{
	delete msg;
    myGSpaceObj->doneDoingIFFT = true;
	resumeControl();
}




/// All ParticlePlane chares finished the nonlocal computations : For debugging only
void GSpaceDriver::allDoneNLForces(CkReductionMsg *msg)
{
    delete msg;
    areNLForcesDone = true;
    if (myGSpaceObj->doneDoingIFFT)
        resumeControl();
}




/// Ortho notifies us that GSpace needs a tolerance update (velocity rotation)
void GSpaceDriver::needUpdatedPsiV()
{
    isPsiVupdateNeeded = true;
    // Once all driver chares are notified, restart Ortho's backward path
    int foo=1;
    contribute(sizeof(int), &foo, CkReduction::min_int,CkCallback(CkIndex_Ortho::resumeV(NULL), UorthoProxy[thisInstance.proxyOffset]));
}




/// Trigger the nonlocal computations
void GSpaceDriver::startNonLocalEes(int iteration_loc)
{

    paraInfo = scProxy.ckLocalBranch ()->cpcharmParaInfo;
    int natm_nl  = paraInfo->natm_nl;

    if(iteration_loc!=myGSpaceObj->iteration)
        CkAbort("GSpaceDriver::startNonLocalEes - Iteration mismatch between GSpace and someone else who asked to launch NL computations\n");

    /// Set to false, just before I spawn the nonlocal work

    if(natm_nl!=0){
      areNLForcesDone = false;
      #define _NLEES_PRIO_START_
      #ifdef _NLEES_PRIO_START_OFF_
         myParticlePlaneObj->startNLEes(iteration);
      #else
         NLDummyMsg *msg = new(8*sizeof(int)) NLDummyMsg;
         CkSetQueueing(msg, CK_QUEUEING_IFIFO);
         *(int*)CkPriorityPtr(msg) = config.sfpriority;
         msg->iteration = myGSpaceObj->iteration;
         UparticlePlaneProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).lPrioStartNLEes(msg);
    #endif
    }//endif
}




#include "structure_factor/StructureFactorMessages.h"

void GSpaceDriver::releaseSFComputeZ()
{
	#ifdef _CP_DEBUG_SF_CACHE_
		CkPrintf("GSpaceDriver[%d,%d] Releasing SF computations\n",thisIndex.x,thisIndex.y);
	#endif
	#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
		CkPrintf("GSpaceDriver[%d,%d] Releasing SF computations\n",thisIndex.x,thisIndex.y);
	#endif
    
    /// Set to false, just before I spawn the nonlocal work
    areNLForcesDone = false;
    /// Tell the structure factor chares to do their work	
	if(thisIndex.x==0)
	{
		/// Multicast to all states of our plane and dups using the section proxy
		SFDummyMsg *msg = new(8*sizeof(int)) SFDummyMsg;
		CkSetQueueing(msg, CK_QUEUEING_IFIFO);
		*(int*)CkPriorityPtr(msg) = config.sfpriority;
		msg->iteration_src = myGSpaceObj->iteration;
		sfCompSectionProxy.computeSF(msg);
	}
	/// Call on corresponding particle plane chare to launch all the Z matrix computations
	myParticlePlaneObj->launchComputeZs();
}

#include "gSpaceDriver.def.h"

