#include "debug_flags.h"
#include "gSpaceDriver.decl.h"
#include "structureFactor.decl.h"
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
 * driver logic that is currently implemented as SDAG code. 
 * 
 * Use starts with the creation of the numInstances bound arrays in init_state_chares(). InstanceController::doneInit() 
 * then invokes startControl() which calls init and then driveGSpace.  The latter uses Structured DAGger to manage the control flow. The intricacies of the workflow
 * should be evident in gspace.ci's documentation and in the SDAG code comments. SDAG defined the chain of dependencies between entry methods and invokes different local object functions when their input dependencies have been met. Flow returns naturally to the driveGSpace as those functions return.
 * 
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

    /// @entry Creates and invokes the SDAG that controls GSpace execution. 
    void startControl();
    /// @entry local. GSpace notifies me that its ready to exit by calling this method
    void readyToExit();
    /// Prints timing info and lets the instance contoller know that computation is done
    void finishComputation();

    /// @entry Ortho notifies us that GSpace needs a tolerance update (velocity rotation)
    void needUpdatedPsiV();
    /// @entry Triggers nonlocal energy computations
    void startNonLocalEes(int iteration_loc);
    /// Triggers nonlocal energy computations
    void releaseSFComputeZ(); 					

    ///
    int ees_nonlocal;
    ///
    int natm_nl;
    ///
    int cp_min_opt;
    ///
    int cp_bomd_opt;
    ///
    int gen_wave;
    ///
    int ndump_frq;
    ///
    bool isPsiVupdateNeeded;
    ///
    bool isOutputNeeded;

    /// Pointer to the GSpacePlane object that I am driving (controlling) 
    CP_State_GSpacePlane *myGSpaceObj;
    /// Pointer to the ParticlePlane object that I am driving (controlling) 
    CP_State_ParticlePlane *myParticlePlaneObj;
    /// A handle to the local copy of the config parameters (refresh on migration)

  private:
    /// A marker tying this class to a particular instance of interacting chares
    const UberCollection thisInstance;
    /// Array section of the structure factor chares that I will be triggering
    CProxySection_StructureFactor sfCompSectionProxy;
    CkVec <CkArrayIndex3D> sfelems;
};
/*@}*/
#endif // GSPACE_DRIVER_H
