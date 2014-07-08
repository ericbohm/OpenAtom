#ifndef _ORTHOSTARTUP_H
#define _ORTHOSTARTUP_H

#include "paircalcstartup.h"
	      // Blame Ram for ugly crime against readability.  More
	      // redundant config objects and builders doesn't help
	      // global clarity at all.
	      // EJB: Sequestered in this file to sweep it under the rug

	      cp::ortho::orthoConfig orthoCfg;
	      orthoCfg.isDynamics    = (sim->cp_min_opt==1)? false: true;
	      orthoCfg.isGenWave     = (sim->gen_wave==1)? true: false;
	      orthoCfg.numStates     = config.nstates;
	      orthoCfg.grainSize     = config.orthoGrainSize;
	      orthoCfg.instanceIndex = thisInstance.getPO();
	      orthoCfg.maxTolerance  = sim->tol_norb;
	      orthoCfg.uponToleranceFailure = CkCallback(CkIndex_GSpaceDriver::needUpdatedPsiV(), UgSpaceDriverProxy[thisInstance.getPO()]);

	      // Fill in the paircalc configs that are instance dependent
	      cfgSymmPC.gSpaceAID            = UgSpacePlaneProxy[thisInstance.getPO()].ckGetArrayID();
	      cfgAsymmPC.gSpaceAID           = UgSpacePlaneProxy[thisInstance.getPO()].ckGetArrayID();
	      cfgSymmPC.instanceIndex        = thisInstance.getPO();
	      cfgAsymmPC.instanceIndex       = thisInstance.getPO();
	      // Init the post-init callbacks that the paircalcs will trigger (after ortho<-->PC comm setup)
	      cfgSymmPC.uponSetupCompletion  = CkCallback(CkIndex_InstanceController::doneInit(NULL),CkArrayIndex1D(thisInstance.getPO()),instControllerProxy.ckGetArrayID());
	      cfgAsymmPC.uponSetupCompletion = CkCallback(CkIndex_InstanceController::doneInit(NULL),CkArrayIndex1D(thisInstance.getPO()),instControllerProxy.ckGetArrayID());

	      // Identify who is the owner for this bubble
	      CkCallback pcHandleCB(CkIndex_CP_State_GSpacePlane::acceptPairCalcAIDs(0), UgSpacePlaneProxy[thisInstance.getPO()]);

          // Fill out a structure with all configs needed for PC mapping
	      cp::startup::PCMapConfig pcMapCfg(boxSize, 
						*peList4PCmapping, 
						&GSImaptable[thisInstance.getPO()],
						(config.torusMap == 1),
						(config.fakeTorus == 1),
						mapOffsets[numInst]);

	      // Delegate the actual construction/initialization to a creation manager
	      cp::startup::PCCreationManager pcCreator(cfgSymmPC, cfgAsymmPC, orthoCfg);
	      pcCreator.build(pcHandleCB, pcMapCfg);
#endif
