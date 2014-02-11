#include "GSpaceDriver.h"
#include "CP_State_GSpacePlane.h"
#include "CP_State_ParticlePlane.h"
#include "main/TimeKeeper.h"

extern int nstates;
extern Config config;

extern CProxy_TimeKeeper 			TimeKeeperProxy;
extern CkVec <CProxy_CP_State_GSpacePlane>      UgSpacePlaneProxy;
extern CkVec <CProxy_CP_State_ParticlePlane> 	UparticlePlaneProxy;
extern CkVec <CProxy_StructureFactor> 			UsfCompProxy;
extern CProxy_InstanceController      instControllerProxy;

/** @addtogroup GSpaceState
    @{
*/
/// Constructor
GSpaceDriver::GSpaceDriver(const UberCollection _thisInstance): 
			thisInstance(_thisInstance),
			myGSpaceObj(0),
			isFirstStep(true),
			isPsiVupdateNeeded(false)
{
	#ifdef DEBUG_CP_GSPACE_CREATION
		CkPrintf("GSpaceDriver[%d,%d] born\n",thisIndex.x,thisIndex.y);
	#endif
	/// Initialize flags and counters that record the control
	/// state
        CPcharmParaInfo *sim  = CPcharmParaInfo::get();
		ees_nonlocal = sim->ees_nloc_on;
		natm_nl = sim->natm_nl;
		cp_min_opt = sim->cp_min_opt;
		gen_wave = sim->gen_wave;
}//end routine




///
GSpaceDriver::GSpaceDriver(CkMigrateMessage *msg): 
			myGSpaceObj(0)
{}




/// PUP method
void GSpaceDriver::pup(PUP::er &p)
{
	__sdag_pup(p);
	p|isFirstStep;
	p|ees_nonlocal;
	p|natm_nl;
	p|cp_min_opt;
	p|gen_wave;
	p|isPsiVupdateNeeded;
	p|sfCompSectionProxy;
}




/** Creates a structure factor array section proxy that is used later in releaseSFComputeZ. 
 * @todo: When is the earliest that you can do this? Can we do this upon construction? Probably not, as we have no guarantee that
 * UsfCompProxy is up and ready.
 */ 
void GSpaceDriver::init()
{
	/// If EES nonlocals are off, each state=0 driver chare
	/// creates a structure factor array section
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
	/// Call the SDAG method in charge of control flow
	driveGSpace();
}




/// GSpace notifies me that its ready to exit by calling this method
// TODO: This might not be needed anymore...GSpace can just set the flag.
// TODO: What if this msg isn't seen until the next iteration? Setting the flag in GSpace might be better?
void GSpaceDriver::readyToExit()
{
	/// Set a flag that indicates this chare is ready to exit
	myGSpaceObj->cleanExitCalled = 1;
}

/// Ortho notifies us that GSpace needs a tolerance update (velocity rotation)
void GSpaceDriver::needUpdatedPsiV()
{
    isPsiVupdateNeeded = true;
    // Once all driver chares are notified, restart Ortho's backward path
    int foo=1;
    contribute( sizeof(int), &foo, CkReduction::min_int, CkCallback(CkIndex_Ortho::resumeV(NULL), myGSpaceObj->myOrtho) );
}




/// Trigger the nonlocal computations
void GSpaceDriver::startNonLocalEes(int iteration_loc)
{

    CPcharmParaInfo *sim  = CPcharmParaInfo::get();
    //int natm_nl  = sim->natm_nl;

    if(iteration_loc!=myGSpaceObj->iteration)
        CkAbort("GSpaceDriver::startNonLocalEes - Iteration mismatch between GSpace and someone else who asked to launch NL computations\n");

    /// Set to false, just before I spawn the nonlocal work

    if(natm_nl!=0){
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
/*@}*/
#include "gSpaceDriver.def.h"
