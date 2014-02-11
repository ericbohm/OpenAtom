#include "debug_flags.h"
#include "gSpaceDriver.decl.h"
#include "structureFactor.decl.h"
//#include "RTH.h"
#include "main/CPcharmParaInfoGrp.h"
#include "charm++.h"
#include "uber/Uber.h"

#ifndef GSPACE_DRIVER_H
#define GSPACE_DRIVER_H

// Forward declarations
class CP_State_GSpacePlane;
class CP_State_ParticlePlane;
/** @addtogroup GSpaceState
    @{
*/

/** @brief 2D chare array [\ref nchareG][\ref nstates]
 * Handles flow of control within an instance, always same dimensional cardinality
 * and mapping as \ref CP_Gspace_State_Plane
*/

/** A 2D chare array bound to CP_State_GspacePlane that manages the driver logic and related communication details.
 * 
 * Implements entry methods that are primarily used by other chares for notifying progress on different portions 
 * of the computation. This information is used to orchestrate the work. Holds counters and flags used in the GSpace 
 * driver logic that is currently implemented as RTH code. 
 * 
 * Use starts with the creation of the numInstances bound arrays in init_state_chares(). InstanceController::doneInit() 
 * then invokes startControl() which creates an RTH 'thread' and lets it fall through its logic. The intricacies of the workflow
 * should be evident in GSpacePlane's documentation and in the RTH code comments. RTH invokes different functions and entry 
 * methods and suspends when needed. The RTH 'thread' resumes when someone calls resumeControl().
 * 
 * Moving to other solutions like sdag should mostly result in very localized modifications.  
 */  
class GSpaceDriver: public CBase_GSpaceDriver
{
	public:
		GSpaceDriver_SDAG_CODE
		/// Constructors
		GSpaceDriver() {}
		GSpaceDriver(CkMigrateMessage *msg);
		GSpaceDriver(const UberCollection _thisInstance);
		/// Initializer
		void init();
		/// PUP method
		void pup(PUP::er &p);
		
		/// @entry Creates and invokes the RTH thread that controls GSpace execution. 
 		void startControl();
		/// @entry local. GSpace notifies me that its ready to exit by calling this method
		void readyToExit();
        /// @entry Ortho notifies us that GSpace needs a tolerance update (velocity rotation)
        void needUpdatedPsiV();
		/// @entry Triggers nonlocal energy computations
        void startNonLocalEes(int iteration_loc);
		/// Triggers nonlocal energy computations
		void releaseSFComputeZ(); 					
		
		/// True if this is the first step
		bool isFirstStep;
		///
		int ees_nonlocal;
		///
		int natm_nl;
		///
		int cp_min_opt;
		///
		int gen_wave;
        /// 
        bool isPsiVupdateNeeded;

		/// Pointer to the GSpacePlane object that I am driving (controlling) 
		CP_State_GSpacePlane *myGSpaceObj;
		/// Pointer to the ParticlePlane object that I am driving (controlling) 
		CP_State_ParticlePlane *myParticlePlaneObj;
		/// A handle to the local copy of the config parameters (refresh on migration)
		
	private:
		/// A marker tying this class to a particular instance of interacting chares
		const UberCollection thisInstance;
		/// An RTH runtime thread that executes all the control logic
//		RTH_Runtime* controlThread;
		/// Array section of the structure factor chares that I will be triggering
		CProxySection_StructureFactor sfCompSectionProxy;
};
/*@}*/
#endif // GSPACE_DRIVER_H
