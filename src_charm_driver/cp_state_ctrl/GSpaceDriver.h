#include "debug_flags.h"
#include "gSpaceDriver.decl.h"
#include "structureFactor.decl.h"
#include "RTH.h"
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
		/// @entry local. Called by a compute chare to give control back to this driver logic
		inline void resumeControl() 					{ RTH_Runtime_resume(controlThread); }
		/// @entry local. GSpace notifies me that its ready to exit by calling this method
		void readyToExit();
		/// @entry local. GSpace notifies me that the energy reduction is done by calling this method
		void doneComputingEnergy(const int AtomsGrpIter); 				
		/// @entry local. GSpace notifies me when the atom integration is complete via this method
		void doneMovingAtoms(const int AtomsGrpIter);
		/// @entry local. GSpace notifies me when the nonlocal force computations are done
		void doneNLForces();
        /// @entry Reduction barrier at the end of the Psi loop for all GSpace chares
        void allDonePsi(CkReductionMsg *msg);
        /// @entry Reduction barrier at the end of the Psi write process for all GSpace chares
        void allDoneWritingPsi(CkReductionMsg *msg);
        /// @entry Reduction barrier at the end of the PsiV update loop for all GSpace chares
        void allDonePsiV(CkReductionMsg *msg);
		/// @entry Reduction barrier at the end of the inverse FFT for all GSpace chares
		void allDoneIFFT(CkReductionMsg *msg);
        /// @entry Reduction barrier at the end of the nonlocal computations for all ParticlePlane chares
        void allDoneNLForces(CkReductionMsg *msg);
        /// @entry Ortho notifies us that GSpace needs a tolerance update (velocity rotation)
        void needUpdatedPsiV();
		/// @entry Triggers nonlocal energy computations
        void startNonLocalEes(int iteration_loc);
		/// Triggers nonlocal energy computations
		void releaseSFComputeZ(); 					
		
		// TODO(mikida2): Delete this before committing
		bool tempFlag;

		/// True if this is the first step
		bool isFirstStep;
		///
		int ees_nonlocal;
        /// Indicates if the nonlocal computation loop has completed. Replaces GPP::doneGettingForces
        bool areNLForcesDone;
        /// 
        bool isPsiVupdateNeeded;
		/** Indicates if the atom integration is done and has returned in the current iteration. 
		 * False after I launch, True after return of atoms. Replaces myatom_integrate_flag
		 */
		bool isAtomIntegrationDone;
		/** Indicates if the energy reduction is done and has returned in the current iteration.
		 *  False after I launch eke, True after return of energy. Replaces myenergy_reduc_flagg
		 */
		bool isEnergyReductionDone;
		/// True if we've suspended control while waiting for the energy computations in GSpace. Replaces isuspend_energy
		bool waitingForEnergy;
		/// True if we've suspended control while waiting for the completion of atom integration. Replaces isuspend_atms
		bool waitingForAtoms;

		/// Pointer to the GSpacePlane object that I am driving (controlling) 
		CP_State_GSpacePlane *myGSpaceObj;
		/// Pointer to the ParticlePlane object that I am driving (controlling) 
		CP_State_ParticlePlane *myParticlePlaneObj;
		/// A handle to the local copy of the config parameters (refresh on migration)
		
	private:
		/// A marker tying this class to a particular instance of interacting chares
		const UberCollection thisInstance;
		/// An RTH runtime thread that executes all the control logic
		RTH_Runtime* controlThread;
		/// Array section of the structure factor chares that I will be triggering
		CProxySection_StructureFactor sfCompSectionProxy;
};
/*@}*/
#endif // GSPACE_DRIVER_H
