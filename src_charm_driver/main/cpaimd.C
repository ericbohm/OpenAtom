/** \file cpaimd.C
 *                         cpaimd-charm-driver
 * \brief
 *    This file contains cpaimd-charm-driver main. It creates and
 *    initializes all the arrays and libraries.
 */

/** \mainpage

 * OpenAtom is a **Car-Parrinello Ab-Initio Molecular Dynamics**
 * software developed by the **Parallel Programing Laboratory, UIUC**
 * in collaboration with **IBM, NYU, and Yale**. See also the
 * [OpenAtom webpage](http://Charm.cs.illinois.edu/OpenAtom).
 *
 * ###Introduction###
 * OpenAtom is parallel simulation software for studying atomic and molecular
 * systems based on quantum chemical principles. In contrast to classical
 * computational molecular dynamics based on Newtonian mechanims, it uses
 * the Car-Parrinello Ab Initio Molecular Dynamics (CPAIMD) approach.
 * Instead of using an empirical force function, the CPAIMD algorithm
 * computes the forces acting on each atom as the summation of multiple
 * terms derived from plane-wave density functional theory. This allows
 * OpenAtom to study complex atomic and electronic physics in semiconductor,
 * metallic, biological and other molecular systems.
 *
 * OpenAtom is implemented on top of [Charm++](http://Charm.cs.illinois.edu),
 * which is an over-decomposition based parallel programming framework that provides
 *
 * support for message-driven execution of migratable entities empowered by an
 * adaptive runtime system. Charm++ encourages decomposition of parallel computation
 * using units that are natural to the application domain, instead of dividing
 * data into as many pieces as processors. In particular, OpenAtom decomposes the data
 * and the computation across a number of *chare* objects, whose type and/or number
 * only depend on the CPAIMD algorithm and the desired grainsize. This allows OpenAtom
 * to exploit the underlying mathematics via a seamless mix of both data and functional
 * decompositions resulting in greater expressed parallelism, and several overlapping
 * phases of computation combined with a longer critical path of dependent computations.
 * The current implementation of OpenAtom in Charm++ is highly scalable, and has
 * exhibited portable performance across three generations of the IBM Blue Gene family,
 * apart from other supercomputing platforms.
 *
 * ###Downloading OpenAtom
 * OpenAtom is hosted using git and can be downloaded using the following command:
 *
 *        git clone http://Charm.cs.uiuc.edu/gerrit/openatom.git
 *
 * Recent commit history can be view
 * [here](http://\charm.cs.illinois.edu/cgi-bin/gitweb.cgi?p=openatom.git;a=summary).
 *
 * You will also need to download either a stable version of Charm++ from this
 * [weblink](http://\charm.cs.uiuc.edu/software) or the nightly version using the
 * following command:
 *
 *        git clone http://Charm.cs.uiuc.edu/gerrit/charm.git
 *
 * Sample data sets can be obtained using one of the following commands:
 *
 *        git clone http://Charm.cs.uiuc.edu/gerrit/datasets/openatom/water_32M_10Ry.git
 *        git clone http://Charm.cs.uiuc.edu/gerrit/datasets/openatom/water_32M_70Ry.git
 *        git clone http://Charm.cs.uiuc.edu/gerrit/datasets/openatom/water_64M_70Ry.git
 *        git clone http://Charm.cs.uiuc.edu/gerrit/datasets/openatom/water_128_70Ry.git
 *
 * ###Compilation###
 * Before OpenAtom is compiled, one needs to get access to a compiled version of
 * Charm++. Detailed instructions on compiling Charm++ can be obtained
 * [here](http://\charm.cs.illinois.edu/manuals/html/charm++/A.html). On a typical
 * 64-bit linux machine, Charm++ can be compiled using the following command (executed
 * within Charm++ directory):
 *
 *        ./build charm++ net-linux-x86_64 --with-production -j8 (production version)
 *        ./build charm++ net-linux-x86_64 -j8 -g (debug version)
 *
 * You will also need FFTW library, configured for double precision, to compile OpenAtom.
 *
 * The INSTALL file in OpenAtom provides detailed instruction for its compilation. Here
 * is a quick summary:
 *
 * 1. Copy a machine specific configuration file (*config.MACHINE.mk*) from
 *    the *makefiles* directory to the OpenAtom base directory and rename it to
 *    *config.mk*.
 * 2. Update the values of **CHARMBASE** and **FFT_HOME** in the
 *    beginning of the newly created *config.mk*. If FFTW was built
 *    with --enable-type-prefix, you must enable DUAL_FFTW, otherwise
 *    set DUAL_FFTW_OFF.
 * 3. Customize the config.mk for any desired compilation/link flags etc.
 * 4. Now type "make", which should create a binary called *OpenAtom* in *build*
 *    directory on successful compilation.
 *
 *###Running OpenAtom###
 *
 *Before executing OpenAtom, obtain a dataset using the *git* command
 *mentioned in the download section, and either place it in the *data*
 *directory, or modify the **w3210** variable in *config.mk*.
 *
 * If the dataset uses the old format, you will need to execute
 * *setup* in dataset directory. *setup* is located in *BASEDIR/util*
 * directory.
 *
 * OpenAtom is to be executed as a Charm++ application, which is
 * explained in detail at
 * [this link](http://Charm.cs.illinois.edu/manuals/html/charm++/C.html). The
 * general syntax is as follows:
 *
 *       ./charmrun +p<N> ./OpenAtom <path to cpaimd config> <path to water input> (N
 *       is the number of processors to execute the job on)
 *
 * On most machine layers (excluding net), a Charm++ application can
 * be launched in the same manner as an MPI application. For example,
 * if built on top of MPI, one can launch OpenAtom as follows
 *
 *       mpirun -np <N> ./OpenAtom <path to cpaimd config> <path to water input>
 *
 * In the datasets downloaded using the command listed above, these
 * files are located in *regression* directory. You can also use
 * *tidy* located in the *util* directory to clean up the dataset
 * directory before performing a new run.
 *
 * See \ref GSpaceDriver::driveGSpace for SDAG flow of control.
 *
 * ![Overview Of OpenAtom Control Flow](controlFlowAmongstChareArrays_small.gif)
 */

int numPes;
#include "cpaimd.h"

#include "InstanceController.h"
#include "ENL_EKE_Collector.h"
#include "pcCreationManager.h"
#include "eesCache.h"
#include "AtomsCache.h"
#include "AtomsCompute.h"
#include "energyGroup.h"

#include "cp_state_ctrl/CP_State_Plane.h"
#include "cp_state_ctrl/CP_State_ParticlePlane.h"
#include "fft_slab_ctrl/fftCacheSlab.h"
#include "structure_factor/StructFactorCache.h"
#include "structure_factor/StructureFactor.h"
#include "paircalc/pcMapConfig.h"
#include "load_balance/PeList.h"
#include "utility/MapFile.h"
#include "PIBeadAtoms.h"
#include "utility/util.h"
#include "src_piny_physics_v1.0/include/class_defs/Interface_ctrl.h"
#include "src_piny_physics_v1.0/include/class_defs/PINY_INIT/PhysicsParamTrans.h"
#include "src_piny_physics_v1.0/include/class_defs/PINY_INIT/PhysicsAtomPosInit.h"
#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h"

#include "TopoManager.h"
#include "TimeKeeper.h"

#include "PhysScratchCache.h"
#include <cmath>
#include <unistd.h>
#include <iostream>
#include <sstream>
#include <string>
#include <assert.h>

int TimeKeeperID=0;
std::vector <std::string> TimeKeeperNames;
UberCollection thisInstance;

//============================================================================
/**
 * \addtogroup mapping
 * torus_vars Defining the size of the torus, handy when debugging
 * torus map logic on non torus architectures.
 */
/**@{*/

bool fakeTorus;
extern void initFFTLock(void);

/**@}*/

/**@defgroup piny_vars piny_vars
 * \brief Defining all Charm++ readonly variables for PINY physics
 *
 */
//============================================================================
/** \addtogroup piny_vars
 **@{*/
extern MDINTEGRATE  readonly_mdintegrate;
extern MDATOMS      readonly_mdatoms;
extern MDINTER      readonly_mdinter;
extern MDINTRA      readonly_mdintra;
extern GENERAL_DATA readonly_general_data;
extern CP           readonly_cp;
/**@}*/
//============================================================================

//declares the maptables,  pelists, uber proxies, etc.
#include "mapvariables.h"
CProxy_PlatformSpecific platformSpecificProxy;

//============================================================================
// readonly globals

double globalTimer;
int Ortho_UE_step2;
int Ortho_UE_step3;
int Ortho_UE_error;
bool Ortho_use_local_cb;
bool HartreeFockOn;
int done_init=0;
int planes_per_pe;
extern int numOrthosPerDim;
extern int totalOrthos;
extern int diagonalization;

//============================================================================

/// Multicast manager group that handles many mcast/redns in the code. Grep for info
CkGroupID            mCastGrpId;
CkReduction::reducerType complexVectorAdderType;
//============================================================================

/// The build system should define this macro to be the commit identifier
#ifndef OPENATOM_REVISION
#define OPENATOM_REVISION Unknown
#endif
#define _QUOTEIT(x) #x
#define INQUOTES(x) _QUOTEIT(x)

/// A global constant for use in the code
const char OpenAtomRevision[] = INQUOTES(OPENATOM_REVISION);

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

/**
 * \brief The Main of CPAIMD, it calls all the init functions.
 */
main::main(CkArgMsg *msg) {
  globalTimer = CkWallTimer();
  topoMgr = NULL;
  //============================================================================
  /**
     # Sequential startup within Main */
  /* Check arguments : Verbose output about startup procedures */

  done_init = 0;

  if (msg->argc < 3) {
    CkAbort("Usage: cpaimd.x cpaimd_config pinysystem.input");
  }//endif

  CkPrintf("Executing OpenAtom: BINARY - %s\n", msg->argv[0]);
  CkPrintf("Binary produced from source-tree at commit: %s\n",OpenAtomRevision);
  CkPrintf("\n");
  PRINT_LINE_STAR;
  CkPrintf("Starting Cpaimd-Charm-Driver Setup Phase\n");
  PRINT_LINE_DASH;
  CkPrintf("  Cpaimd-Charm-Driver running on %d processors. \n", CkNumPes());

  //============================================================================
  /* check the debug flags for consistency*/

#ifdef _CP_DEBUG_NON_LOCAL_ONLY_
#ifdef _CP_DEBUG_SFNL_OFF_
  CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
  CkPrintf("You can't test non-local by itself and turn it off\n");
  CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
  CkExit();
#endif
#endif

#ifdef _CP_DEBUG_NON_LOCAL_ONLY_
#ifndef _CP_DEBUG_VKS_OFF_
  CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
  CkPrintf("You can't test non-local by itself with vks on\n");
  CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
  CkExit();
#endif
#endif

#ifdef _CP_DEBUG_VKS_ONLY_
#ifdef _CP_DEBUG_VKS_OFF_
  CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
  CkPrintf("You can't test vks by itself and turn it off\n");
  CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
  CkExit();
#endif
#endif

#ifdef _CP_DEBUG_VKS_ONLY_
#ifndef _CP_DEBUG_SFNL_OFF_
  CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
  CkPrintf("You can't test vks by itself with SNFL on\n");
  CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
  CkExit();
#endif
#endif

#ifdef _CP_DEBUG_NON_LOCAL_VKS_ONLY_
#ifdef _CP_DEBUG_SFNL_OFF_
  CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
  CkPrintf("You can't test vks and SNFL only without sfnl \n");
  CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
  CkExit();
#endif
#ifdef _CP_DEBUG_VKS_OFF_
  CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
  CkPrintf("You can't test vks and SNFL only without vks \n");
  CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
  CkExit();
#endif
#endif

  //============================================================================


  /** @defgroup startup startup
   * \brief Startup parses the physics simulation input, parses the parallel driver config file, constructs chares, computes chare placement, coordinates the launch of chares in parallel, the reading of input, and initiates computation when everything is initialized and ready to start timestepping.


   1)  Read PINY config and the rest of the physical system parameter
   files to set the constants which determine problem size.  These are
   stored in the CPcharmParaInfo class object named sim (for
   simulation).*/
  /**@{*/

  /* Invoke PINY input class */
  CkPrintf("At the beginning of the run user mem %lf MB\n",
      (CmiMemoryUsage()/(1024.0*1024.0)));
  CkPrintf("  Reading Physics input from %s\n",msg->argv[2]);

  CkCallback piny_callback (CkCallback::ignore);
  Interface_ctrl piny_interface (msg->argv[2], piny_callback);

  CPcharmParaInfo *sim  = CPcharmParaInfo::get();
  PhysicsParamTransfer::ParaInfoInit(sim);
  int ibinary_opt    = sim->ibinary_opt;
  int natm_nl        = sim->natm_nl;
  int ees_eext_opt   = sim->ees_eext_on;
  int natm_typ       = sim->natm_typ;
  int fftopt         = sim->fftopt;
  /**@}*/
  //============================================================================
  /* Invoke parallel driver input class */
  /** \addtogroup startup
      2)  Read the cpaimd_config file which determines parallel
      decomposition.  This is mostly stored in the config object, but a
      bunch of readonly globals also get instantiated with this
      information.
  */
  /**@{*/

  CkPrintf("  Reading Driver input from %s\n",msg->argv[1]);

  double Timer = CmiWallTimer();
  double phase1start = Timer;
  numPes = CkNumPes();
  int minimization_steps;
  if (sim->cp_bomd_opt) {
    minimization_steps = sim->btime;
  } else {
    minimization_steps = sim->ntime;
  }
  config.readConfig(msg->argv[1],sim->nstates,sim->sizeX,sim->sizeY,sim->sizeZ,
		    minimization_steps,ibinary_opt,natm_nl,fftopt,numPes,natm_typ,
		    ees_eext_opt,sim->gen_wave,sim->ncoef, sim->cp_min_opt,
                    sim->ngrid_eext_c, sim->doublepack,sim->pi_beads,sim->nkpoint,
                    sim->ntemper,sim->nspin);

  fakeTorus        = config.fakeTorus > 0;

  if(fakeTorus) {
    CkAbort("Fake torus based runs are no longer supported\n");
  } else if (CkNumPes() != config.numPes) {
    numPes=config.numPes;
    CkPrintf("numpes set to %d by config file\n",numPes);
  }

  CkPrintf("for numInstances %d numPes %d numPesPerInstance is %d \n",
        config.numInstances, config.numPes, config.numPesPerInstance);

  mapOffsets = new inttriple[config.numInstances];
  int numSfGrps    = config.numSfGrps;  // local copies are nice
  int doublePack   = config.doublePack;
  sim->nchareG     = config.nchareG;

  double newtime= CmiWallTimer();

  PRINT_LINE_DASH;
  CkPrintf("Cpaimd-Charm-Driver input completed in %g\n",newtime-Timer);
  PRINT_LINE_STAR; CkPrintf("\n");
  Timer=newtime;

  // user event trace setup
  setTraceUserEvents();
  // TODO timekeeper registerees will need to distinguish by instance
  // timekeeper itself doesn't care.
  TimeKeeperProxy = CProxy_TimeKeeper::ckNew();
  // Create a multicast manager group that will handle many mcast/redns

  mCastGrpId = CProxy_CkMulticastMgr::ckNew(config.numMulticastMsgs);

  // Create a paircalc config object for the symmetric PC instance
  pc::pcConfig cfgSymmPC;
  // Create a paircalc config object for the asymmetric PC instance
  pc::pcConfig cfgAsymmPC;

  /* choose whether ortho should use local callback */
  Ortho_use_local_cb = true;
  HartreeFockOn = false;
  paircalcstartup(&cfgSymmPC, &cfgAsymmPC, sim, doublePack);
  cp::ortho::orthoConfig orthoCfg;
  //============================================================================
  // Compute structure factor grp parameters and static map for chare arrays
  /**@}*/
  sim->numSfGrps   = numSfGrps;
  int natm_nl_grp_max;
  PhysicsParamTransfer::get_Sfgrp_max(natm_nl,config.numSfGrps,
				      &natm_nl_grp_max);
  sim->natm_nl_grp_max = natm_nl_grp_max;
  /* @} */

  CkPrintf("Before create_line_decomp_descriptor user mem %lf MB\n", 
    (CmiMemoryUsage()/(1024.0*1024.0)));
  create_line_decomp_descriptor(sim);

  CkPrintf("Before control_new_mapping_function user mem %lf MB\n", 
    (CmiMemoryUsage()/(1024.0*1024.0)));
  PhysicsParamTransfer::control_new_mapping_function(sim, doublePack);

  CkPrintf("Before make_rho_runs user mem %lf MB\n", 
    (CmiMemoryUsage()/(1024.0*1024.0)));
  //make_rho_runs(sim);


#include "initializeUber.C"

  // make one controller chare per instance

  GENERAL_DATA *general_data = GENERAL_DATA::get();

  char *output_directory=general_data->gentempering_ctrl.output_directory;

  if(output_directory==NULL)
    {
      output_directory="TEMPER_OUT";
    }
  CkPrintf("tempering output dir %s\n",output_directory);
  
  char *historyfile=general_data->gentempering_ctrl.history_name;


  mainProxy=thishandle;

  // make one controller chare per instance

  config.temperCycle=general_data->gentempering_ctrl.switch_steps;
  CkPrintf("Temperature exchange frequency set to %d\n",config.temperCycle);

  if(historyfile==NULL)
    {
      historyfile="temperature_trace.out";
    }



  int numFFTinstances = 0;
  numFFTinstances += 3; //3 ffts for divRhos
  numFFTinstances += 1; //one for hartExt
  if(sim->ees_eext_on) {
    numFFTinstances += config.nchareHartAtmT + 1; //these for atmSF
  }

  instControllerProxy= CProxy_InstanceController::ckNew(numFFTinstances,
      config.numInstances);
  instControllerProxy.doneInserting();

  diagonalizerBridgeProxy = CProxy_DiagonalizerBridge::ckNew();

  // make one controller temper
  temperControllerProxy= CProxy_TemperController::ckNew(1,sim->temper_t_ext,sim->ntemper, sim->seed, std::string(historyfile), std::string(output_directory),1);
  temperControllerProxy.doneInserting();
  // make one collector per uberKmax
  CkArrayOptions enlopts(config.UberKmax);

  CkPrintf("expected energies %d\n",config.UberImax * 
	   config.UberJmax * config.UberMmax);
  ENLEKECollectorProxy= CProxy_ENL_EKE_Collector::ckNew(config.UberImax * 
	config.UberJmax * config.UberMmax , config.UberKmax, config.UberImax,
        std::string(output_directory),  enlopts); 

  ENLEKECollectorProxy.doneInserting();

  /**@}*/
  //============================================================================
  /* Create the multicast/reduction manager for array sections
     Create the parainfo group from sim
     Initialize chare arrays for real and g-space of states
  */
  /** \addtogroup startup
      3) Topological map setup : we initialize several structures which will be useful for using network topology for more optimal chare placement */
  /**@{*/

  platformSpecificProxy=CProxy_PlatformSpecific::ckNew();

  CkPrintf("Initializing TopoManager\n");
  topoMgr = new TopoManager();

  if(config.torusMap==1) {
    CkAbort("Specialized torusMap are no longer supported. use simpleTopo and simpleTopoCentroid\n");
  }

  CkPrintf("After Init topoManager user mem %lf MB\n", 
    (CmiMemoryUsage()/(1024.0*1024.0)));

  int l = config.Gstates_per_pe;
  int m, pl, pm;
  pl = sim->nstates / l;
  pm = config.numPesPerInstance / pl;
  if(pm == 0) {
    CkPrintf("Choose a larger Gstates_per_pe than %d such that { (no. of "
        "processors [%d] / no. of Instances [%d]) / (no. of states [%d] / "
        "Gstates_per_pe [%d]) } is > 0 \n", l, config.numPes,
        config.numInstances, sim->nstates, l);
    assert(pm > 0);
  }
  m = config.nchareG / pm;

  planes_per_pe=m;
  if(planes_per_pe <= 0) {
    CkPrintf("Choose a smaller Gstates_per_pe than %d such that { (no. of "
        "processors [%d] / no. of Instances [%d]) / (no. of states [%d] / "
        "Gstates_per_pe [%d]) } is > 0 and config.nchareG [%d] / pm [%d] {"
        "where pm = config.numPesPerInstance [%d]/ pl [%d] }  > 0\n",
        l, config.numPes, config.numInstances, sim->nstates, l, config.nchareG,
        pm, config.numPesPerInstance, pl);
    assert(m > 0);
  }

  CkPrintf("Initializing PeList\n");

  PeList *gfoo=NULL;
  PeList *rfoo=NULL;
  int x, y, z;
  PeListFactory *peList4PCmapping;
  
  if(config.simpleTopo) {
    boxSize = config.numPesPerInstance / sim->nchareG;
  } else {
    gfoo = new PeList(1, 0, config.numPesPerInstance);                            // heap it
    peList4PCmapping = new PeListFactory(config.numPesPerInstance);
  }

  if(!config.simpleTopo)
    rfoo = new PeList(1, 0, config.numPesPerInstance);                         // heap it

  CkPrintf("After Init rfoo list user mem %lf MB\n", 
    (CmiMemoryUsage()/(1024.0*1024.0)));
  computeMapOffsets();

  newtime = CmiWallTimer();
  /* these really don't need to be different */
  if(!config.simpleTopo) {
    availGlobG = rfoo;
    availGlobR = gfoo;
    CkPrintf("Pelist initialized in %g with %d elements\n", newtime-Timer, availGlobG->count());
  }
  newtime = CmiWallTimer();
  Timer = newtime;
  CkPrintf("Intermediate creation done in %g\n",newtime-Timer);
  PRINT_LINE_STAR; CkPrintf("\n");
  diagonalization = sim->cp_min_diagonalize;
  /**@}*/
  /*
    ===============================================================================
    Per Instance startup BEGIN
    ===============================================================================
  */
  // used as a signal during startup, to handle things, such as map
  // creation, which will only be done for the first instance.

  // maps will have a transform function to compute the placement
  // for instances after the first.

  /** \addtogroup startup
      4) Parallel Object Proxy Creation
      We loop through a four deep nested launcher by integral, kpoint,
      temper, and spin to construct the appropriate chare arrays for each
      instance.
      + We setup the State chares (GS, RS, PP, RPP) along with their maps.
      + We setup the Paircalc and ortho chares (SymPC, AsymmPC, ortho, CLA)
      + We setup the Density chares (rho, rho*hart)
      + We setup the non-local chares for Ees (if configured)

      In each case we are constructing the proxy and calling for parallel
      object construction of the elements in each of those arrays.  Note that while we are inside Main the chares for PE 0 will be constructed inline in program order.  No assumptions should be made about chares constructed on other PEs until we exit main and pass control to the Charm++ scheduler for parallel launch.
  */
  /**@{*/
  CkPrintf("NumInstances %d: Tempers %d * Beads %d  * Kpoints %d * Spin %d\n",
      config.numInstances, config.UberKmax, config.UberImax, config.UberJmax,
      config.UberMmax);
  for(int spin=0; spin< config.UberMmax; spin++) {
    for(int temper=0; temper< config.UberKmax; temper++) {
      for(int kpoint=0; kpoint< config.UberJmax; kpoint++) {
	for(int integral=0; integral< config.UberImax; integral++) {
                if(config.simpleTopo) {
                  Timer = newtime;
                  int ndims;
                  TopoManager_getDimCount(&ndims);
                  int bdims[10];
                  bdims[0] = bdims[1] = bdims[2] = bdims[3] = bdims[4] = 4;
                  gfoo = new PeList(ndims, bdims, numInst);
                  peList4PCmapping = new PeListFactory(ndims, bdims, numInst);
                  rfoo = new PeList(1, 0, *gfoo);
                  availGlobG = rfoo;
                  availGlobR = gfoo;
                  newtime = CmiWallTimer();
                  CkPrintf("Pelist initialized in %g with %d elements\n", newtime-Timer, availGlobG->count());
                }

		// for each new instance we need a new Uber Index

		CkVec  <int>  peUsedBySF;
		CkVec  <int>  peUsedByNLZ;
		CkVec  <int>  planeUsedByNLZ;

		UberIndex thisInstanceIndex(integral, kpoint, temper, spin); // Internal labels{x,y,z,s}
		thisInstance = UberCollection(thisInstanceIndex);

		UberAlles.push_back(thisInstance);// collection of proxies for all instances
		CkPrintf("Making offset %d: * Temper %d * Bead %d * Kpoint %d * Spin %d, idx(%d,%d,%d,%d) calc %d\n",
			 thisInstance.proxyOffset,
			 temper, integral, kpoint, spin, thisInstance.idxU.x, thisInstance.idxU.y, thisInstance.idxU.z, thisInstance.idxU.s, thisInstance.calcPO());


		//============================================================================
		// We will need a different one of these per instance
		// Transfer parameters from physics to driver
		//    read in atoms : create atoms group
		control_physics_to_driver(thisInstance, sim);

		//============================================================================

		if(config.UberImax>1)	      // handle Path Integrals
		  init_PIBeads(sim, thisInstance);

                //--------------------------------------------------------------------------------
                UsfCacheProxy.push_back( CProxy_StructFactCache::ckNew(config.numSfGrps, 
                    sim->natm_nl, sim->natm_nl_grp_max, thisInstance));
                CkPrintf("created sfcache proxy\n");
                UsfCompProxy.push_back(CProxy_StructureFactor::ckNew());
                CkPrintf("created sfcomp proxy\n");

		// and then we make the usual set of chares to which we pass
		// the Uber Index.
		//Need our maps and groups to exist before anyone tries to use them
                Timer = newtime;
		if(firstInstance || config.simpleTopo)
		  build_all_maps(sim, thisInstance);
		else
		  build_uber_maps(sim, thisInstance);
                newtime = CmiWallTimer();
                CkPrintf("All maps initialized in %g s\n", newtime-Timer);

                Timer = newtime;
		init_state_chares(natm_nl,natm_nl_grp_max,numSfGrps,doublePack,sim, thisInstance);
                newtime = CmiWallTimer();
                CkPrintf("All states created in %g s\n", newtime-Timer);

		//============================================================================
		// Create a paircalc/ortho bubble (symm and asymm pcs, ortho and related frills)

                Timer = newtime;

		// set eigenvector filename based on uber
		char basename[1024];
		GENERAL_DATA *general_data = GENERAL_DATA::get();
		char *ksname=general_data->genfilenames.ksname;
		snprintf(basename,1024, "%s/Spin.%d_Kpt.%d_Bead.%d_Temper.%d/%s",config.dataPathOut, spin, kpoint, integral, temper, ksname);
		orthoCfg.eigenFileName=std::string(basename);
		orthostartup(&orthoCfg, &cfgSymmPC, &cfgAsymmPC, sim, peList4PCmapping);
                newtime = CmiWallTimer();
                CkPrintf("Ortho/PC created in %g s\n", newtime-Timer);

		// Create mapping classes for Paircalcular


		//============================================================================
		// Initialize the density chare arrays
                Timer = newtime;
		init_rho_chares(sim, thisInstance);
                newtime = CmiWallTimer();
                CkPrintf("Rho created in %g s\n", newtime-Timer);

		//============================================================================
		// Initialize commlib strategies for later association and delegation
                Timer = newtime;
		if(sim->ees_nloc_on)
		  init_eesNL_chares( natm_nl, natm_nl_grp_max, doublePack, excludePes, sim, thisInstance);
                newtime = CmiWallTimer();
                CkPrintf("NL created in %g s\n", newtime-Timer);
		firstInstance=false;


                if(config.simpleTopo) {
                  delete rfoo;
                  delete gfoo;
                  delete excludePes;
		}
		// now safe to init atom bead commanders
		UatomsComputeProxy[numInst].init();
		numInst++;
	  }
      }
    }
  } // end of per instance init
  //============================================================================


  if (HartreeFockOn) {
    HFCalculatorProxy = CProxy_HFCalculator::ckNew();
    HFCalculatorProxy.run();
  }

  TimeKeeperProxy.init();

  //============================================================================
  // clean up

  delete msg;
  if(!config.simpleTopo) {
    delete rfoo;
    delete gfoo;
    delete excludePes;
  }

  //============================================================================
  CkPrintf("After all deletes  user mem %lf MB\n", (CmiMemoryUsage()/(1024.0*1024.0)));

  newtime=CmiWallTimer();
  PRINT_LINE_DASH;
  CkPrintf("Cpaimd-Charm-Driver setup phase completed in %g \n",newtime-phase1start);
  CkPrintf("Total time to execute set up in the main chare %g \n",newtime-globalTimer);
  PRINT_LINE_STAR; CkPrintf("\n");
  PRINT_LINE_STAR;
  PRINT_LINE_DASH;CkPrintf("\n");
  CkPrintf("user mem %.2lf MB\n",CmiMemoryUsage()/(1024.0*1024));
  /**@}*/
  //============================================================================
  traceBegin();
}// end Main
//============================================================================

//============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * Cleanup stuff in the hopes of getting clean valgrind
 */
main::~main(){
}//end routine
//============================================================================


//============================================================================
//============================================================================
//Create the structure factor in non-EES case
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void init_SF_non_EES(int natm_nl, int natm_nl_grp_max, int numSfGrps,
		       CPcharmParaInfo *sim,  UberCollection thisInstance)
{
  // Determine placement of structure factors
  CkPrintf("Making SF non-EES\n");
  int nchareG      = config.nchareG;
  int nstates      = config.nstates;
  int *nsend       = new int[nchareG];
  int **listpe     = new int * [nchareG];
  int numproc      = config.numPes;
  int *gspace_proc = new int [numproc];
  for(int i =0;i<numproc;i++){gspace_proc[i]=0;}
  for(int j=0;j<nchareG;j++){
    listpe[j]= new int[nstates];
    nsend[j]=0;
    for(int i=0;i<nstates;i++){
      listpe[j][i] = gsprocNum(sim, i, j, thisInstance.getPO());
      gspace_proc[listpe[j][i]]+=1;
    }//endfor
    lst_sort_clean(nstates, &nsend[j], listpe[j]);
  }//endfor

  FILE *fp = fopen("gspplane_proc_distrib.out","w");
  for(int i=0;i<numproc;i++){
    fprintf(fp,"%d %d\n",i,gspace_proc[i]);
  }//endfor
  fclose(fp);
  delete [] gspace_proc;

  int minsend=nstates;
  int maxsend=0;
  double avgsend=0.0;
  int chareG_use=0;
  PRINT_LINE_STAR;
  CkPrintf("Structure factor chareG dests\n");
  CkPrintf("Number of g-space chares : %d\n",nchareG);
  PRINT_LINE_DASH;
  for(int lsi=0;lsi<nchareG;lsi++){
    chareG_use++;
    CkPrintf("chareG [%d] nsend %d\n",lsi,nsend[lsi]);
    if(nsend[lsi]>maxsend)
      maxsend=nsend[lsi];
    if(nsend[lsi]<minsend)
      minsend=nsend[lsi];
    avgsend+=(double) nsend[lsi];
#ifdef SF_SEND_LIST_OUT
    for(int lsj=0;lsj<nsend[lsi];lsj++){
      CkPrintf("[%d %d] pe %d\n",lsi,lsj,listpe[lsi][lsj]);
    }//endfor
#endif
  }//endfor
  PRINT_LINE_DASH;
  CkPrintf("SFSends min %d max %d avg %g\n",minsend,maxsend,avgsend/(double)chareG_use);
  PRINT_LINE_STAR;

  //--------------------------------------------------------------------------------
  // Insert the objects into the StructureFactor array

  int dupmax=maxsend;  // there is no point in ever having more than that
  if(config.numSfDups<dupmax){dupmax=config.numSfDups;}
  config.numSfDups = dupmax;
  int numSfDups    = dupmax;
  CkPrintf("real numSfdups is %d based on maxsend of %d\n",numSfDups, maxsend);
  CkVec  <int>  peUsedBySF;
  for (int dup=0; dup<dupmax; dup++){
    for (int x = 0; x < nchareG; x += 1){
      int num_dup, istart, iend;
      get_grp_params( nsend[x],  numSfDups,  dup, x ,&num_dup,  &istart, &iend);
      int pe_ind=istart;
      if(x%2==0)
	pe_ind=iend;
      for (int AtmGrp=0; AtmGrp<numSfGrps; AtmGrp++){
	UsfCompProxy[thisInstance.proxyOffset](AtmGrp, x, dup).insert(numSfGrps,numSfDups,
								      natm_nl_grp_max,  num_dup, &(listpe[x][istart]),
								      thisInstance,atmGrpMap(istart, num_dup, nsend[x], listpe[x],AtmGrp,dup,x));
	peUsedBySF.push_back(atmGrpMap(istart, num_dup, nsend[x], listpe[x],   AtmGrp, dup,x));

	pe_ind++;
	if(pe_ind>nsend[x]){ pe_ind=0;}
      }//endfor : AtmGrp
    }//endfor : chareG
  }//endfor : Dups
  UsfCompProxy[thisInstance.proxyOffset].doneInserting();
  UpeUsedBySF.push_back(peUsedBySF);
  for(int j=0;j<nchareG;j++){delete [] listpe[j];}
  delete [] listpe;
  delete [] nsend;
}

void fillInPeUsedBy(CPcharmParaInfo *sim, UberCollection thisInstance) {
  if(!sim->ees_nloc_on) {
    init_SF_non_EES(sim->natm_nl, sim->natm_nl_grp_max, config.numSfGrps, sim, thisInstance);
  }

  CkVec  <int>  peUsedBySF;
  CkVec  <int>  peUsedByNLZ;
  CkVec  <int>  planeUsedByNLZ;

  //============================================================================
  // compute the location for the non-local Z reduction roots for each plane
  // this can then be used in exclusion mapping to avoid overloading them
  int *usedProc= new int[config.numPes];
  memset(usedProc, 0, sizeof(int)*config.numPes);
  int charperpe = sim->nstates/(config.numPesPerInstance);
  if(sim->nstates % config.numPesPerInstance != 0)  charperpe++;
  if(charperpe<1) charperpe=1;
  for(int state=0; state<sim->nstates; state++) {
    int plane = config.nchareG-1;
    while(plane >= 0) {
      bool used = false;
      int thisstateplaneproc = GSImaptable[thisInstance.getPO()].get(state,plane);
      if(usedProc[thisstateplaneproc]>charperpe) {
        used=true;
      }
      if(!used || plane==0) {
        peUsedByNLZ.push_back(thisstateplaneproc);
        planeUsedByNLZ.push_back(plane);
        usedProc[thisstateplaneproc]++;
        plane=-1;
      }
      plane--;
    }
  }
  peUsedByNLZ.quickSort();
  delete [] usedProc;
  UpeUsedByNLZ.push_back(peUsedByNLZ);
  UplaneUsedByNLZ.push_back(planeUsedByNLZ);
  // CkPrintf("UplaneUsedByNLZ length now %d\n",UplaneUsedByNLZ.length());
}

//============================================================================
//============================================================================
//Create the chare array element to PE maps for initial instance of all CkArrays
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void build_all_maps(CPcharmParaInfo *sim, UberCollection thisInstance)
{
  double newtime=CmiWallTimer();
  double Timer = newtime;
  availGlobG->reset();
  GSImaptable[numInst].buildMap(config.nstates, config.nchareG);
  int maploaded = 0;
  MapFile *mf=NULL;
  int size[2]={config.nstates, config.nchareG};
  if(config.loadMapFiles) {

     mf = new MapFile("GSMap", 2, size , config.numPes, "TXYZ", 2, 1, 1, 1);
     maploaded=mf->loadMap("GSMap", &GSImaptable[numInst]);
  }

  if(!maploaded) GSMapTable gsTable = GSMapTable(&GSImaptable[0], &GSImaptable[numInst],
					availGlobG, config.nchareG,
					config.nstates, config.Gstates_per_pe,
					config.useCuboidMap, numInst);
  newtime=CmiWallTimer();
  CkPrintf("GSMap created in %g\n", newtime-Timer);
  if(config.dumpMapFiles)
    {
      if(!mf)
	mf = new MapFile("GSMap", 2, size , config.numPes, "TXYZ", 2, 1, 1, 1);
      mf->dumpMap(&GSImaptable[thisInstance.getPO()], thisInstance.getPO());
    }
  if(mf) delete mf;
  mf=NULL;
  if(config.dumpMapCoordFiles) {
    // if someone wants to dump nontopo maps, they'll need a manager
    if(topoMgr == NULL)  topoMgr = new TopoManager();
    MapFile *mfc = new MapFile("GSMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
    mfc->dumpMapCoords(&GSImaptable[thisInstance.getPO()], thisInstance.getPO());
    delete mfc;
  }

  /* Realspace */

  RSImaptable[numInst].buildMap(config.nstates, sim->sizeZ);
  Timer=CmiWallTimer();
  availGlobR->reset();
  maploaded=0;
  size[0] = config.nstates; size[1] = sim->sizeZ;
  if(config.loadMapFiles) {
      mf = new MapFile("RSMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
      maploaded = mf->loadMap("RSMap", &RSImaptable[numInst]);
  }
  if(!maploaded) {
    RSMapTable RStable= RSMapTable(&RSImaptable[0], &RSImaptable[numInst], availGlobR,
				   config.nstates, sim->sizeZ, config.Rstates_per_pe, config.useCuboidMapRS,
				   &GSImaptable[numInst], config.nchareG, numInst);
    }
  newtime=CmiWallTimer();
  CkPrintf("RSMap created in %g\n", newtime-Timer);
  Timer=newtime;
  if(config.dumpMapFiles) {
    MapFile *mf = new MapFile("RSMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
    mf->dumpMap(&RSImaptable[thisInstance.getPO()], thisInstance.getPO());
    delete mf;
  }

  if(config.dumpMapCoordFiles) {
    MapFile *mf = new MapFile("RSMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
    mf->dumpMapCoords(&RSImaptable[thisInstance.getPO()], thisInstance.getPO());
    delete mf;
  }

  availGlobR->reset();

  Timer=newtime;
  fillInPeUsedBy(sim, thisInstance);
  newtime=CmiWallTimer();
  CkPrintf("FillIn created in %g\n", newtime-Timer);

  /* Rho */

  PeList *RhoAvail=NULL;
  if(!config.loadMapFiles)
    {
      availGlobR->reset();
      // make RhoAvail based on the pes used by RS
      RhoAvail=new PeList(1, 1, *availGlobR);
    }
  //------------------------------------------------------------------------
  // subtract processors used by other nonscaling chares (i.e., non local reduceZ)
  excludePes = new PeList(1, 0, 0);
  if(config.excludePE0 && !config.loadMapFiles)
    excludePes->checkAndAdd(0);

  // avoid co-mapping with the particle plane non-local reduction roots
  if(config.useReductionExclusionMap && !config.loadMapFiles)
  {
    if( config.nchareRhoR_x * config.nchareRhoR_y +
        UpeUsedByNLZ[thisInstance.proxyOffset].size() < RhoAvail->count()){

      CkPrintf("subtracting %d NLZ nodes from %d for RhoR Map\n",
          UpeUsedByNLZ[thisInstance.proxyOffset].size(),RhoAvail->count());
      //       nlz.dump();
      PeList nlz(0, 0, UpeUsedByNLZ[thisInstance.proxyOffset]);
      RhoAvail->deleteList(nlz, 0, 1); //unary minus operator defined in PeList.h
      CkPrintf("Leaving %d for RhoR Map\n",RhoAvail->count());
    }//endif

    //------------------------------------------------------------------------
    // subtract processors used by other nonscaling chares

    if(sim->ees_nloc_on==0){
      if( config.nchareRhoRHart_x * config.nchareRhoRHart_y +
          UpeUsedBySF[thisInstance.proxyOffset].size() < RhoAvail->count()){
        CkPrintf("subtracting %d SF nodes from %d for RhoR Map\n",
            UpeUsedBySF[thisInstance.proxyOffset].size(),RhoAvail->count());
        PeList sf(0, 0, UpeUsedBySF[thisInstance.proxyOffset]);
        RhoAvail->deleteList(sf, 0, 1);
      }//endif
    }//endif
  }
  if(!config.loadMapFiles && RhoAvail->count()>2 ) {
    RhoAvail->reset();
  }

  maploaded=0;
  size[0] = config.nchareRhoR_x; size[1] = config.nchareRhoR_y;
  RhoRSImaptable[thisInstance.getPO()].buildMap(size[0], size[1]);
  Timer=newtime;
  if(config.loadMapFiles) {
    mf = new MapFile("RhoRSMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
    maploaded = mf->loadMap("RhoRSMap", &RhoRSImaptable[thisInstance.getPO()]);
  }
  if(!maploaded) {
    RhoRSMapTable RhoRStable(&RhoRSImaptable[thisInstance.getPO()], RhoAvail,
    size[0], size[1], config.nstates, config.useCentroidMapRho,
    &RSImaptable[thisInstance.getPO()], excludePes);
  }
  newtime=CmiWallTimer();
  CkPrintf("RhoRSMap created in %g\n", newtime-Timer);
  if(config.dumpMapFiles) {
    if(!mf)
      mf = new MapFile("RhoRSMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
    mf->dumpMap(&RhoRSImaptable[thisInstance.getPO()], thisInstance.getPO());
  }
  if(mf) { delete mf; mf=NULL;}
  if(config.dumpMapCoordFiles) {
    MapFile *mfc = new MapFile("RhoRSMap_coord", 2, size, config.numPes, "TXYZ",
        2, 1, 1, 1);
    mfc->dumpMapCoords(&RhoRSImaptable[thisInstance.getPO()], thisInstance.getPO());
    delete mfc;
  }

  if(!config.loadMapFiles && config.nchareRhoG > RhoAvail->count()) {
    CkPrintf("Refreshing avail list count %d less than rhog %d\n",
    RhoAvail->count(), config.nchareRhoG);
    RhoAvail->reset();
  }

  maploaded = 0;
  size[0] = config.nchareRhoG;
  RhoGSImaptable[thisInstance.getPO()].buildMap(size[0]);
  Timer=newtime;
  if(config.loadMapFiles) {
    mf = new MapFile("RhoGSMap", 1, size, config.numPes, "TXYZ", 2, 1, 1, 1);
    maploaded = mf->loadMap("RhoGSMap", &RhoGSImaptable[thisInstance.getPO()]);
  }
  if(!maploaded) {
    RhoGSMapTable RhoGStable(&RhoGSImaptable[thisInstance.getPO()], RhoAvail,
        size[0], config.useCentroidMapRho, &RhoRSImaptable[thisInstance.getPO()],
        excludePes);
  }
  newtime=CmiWallTimer();
  CkPrintf("RhoGSMap created in %g\n", newtime-Timer);
  if(config.dumpMapFiles) {
    if(!mf) mf = new MapFile("RhoGSMap", 1, size, config.numPes, "TXYZ", 2, 1,
        1, 1);
    mf->dumpMap(&RhoGSImaptable[thisInstance.getPO()], thisInstance.getPO());
  }
  if(mf) { delete mf; mf=NULL;}
  if(config.dumpMapCoordFiles) {
    MapFile *mfc = new MapFile("RhoGSMap_coord", 1, size, config.numPes, "TXYZ",
        2, 1, 1, 1);
    mfc->dumpMapCoords(&RhoGSImaptable[thisInstance.getPO()], thisInstance.getPO());
    delete mfc;
  }
  //---------------------------------------------------------------------------
  // rho RHart
  if(sim->ees_eext_on) {
    // if there aren't enough free procs refresh the avail list;
    if(!config.loadMapFiles && 
        config.nchareRhoRHart_x*config.nchareRhoRHart_y*config.nchareHartAtmT > RhoAvail->count())
      RhoAvail->reset();
    int size[3];
    size[0] = config.nchareRhoRHart_x; size[1] = config.nchareRhoRHart_y;
    size[2] = config.nchareHartAtmT;
    maploaded=0;
    Timer=newtime;
    RhoRHartImaptable[thisInstance.getPO()].buildMap(size[0], size[1], size[2]);
    if(config.loadMapFiles) {

      mf = new MapFile("RhoRHartMap", 3, size, config.numPes, "TXYZ", 2, 1, 1, 1);
      maploaded = mf->loadMap("RhoRHartMap", &RhoRHartImaptable[thisInstance.getPO()]);
    }
    if(!maploaded) {
      RhoRHartMapTable RhoRHarttable(&RhoRHartImaptable[thisInstance.getPO()],
          RhoAvail, size[0], size[1], size[2], excludePes);
    }
    newtime=CmiWallTimer();
    CkPrintf("RhoRHart created in %g\n", newtime-Timer);
    if(config.dumpMapFiles) {
      if(!mf) mf = new MapFile("RhoRHartMap", 3, size, config.numPes, "TXYZ",
            2, 1, 1, 1);
      mf->dumpMap(&RhoRHartImaptable[thisInstance.getPO()], thisInstance.getPO());
    }
    if(mf) { delete mf; mf=NULL;}
    if(config.dumpMapCoordFiles) {
      MapFile *mfc = new MapFile("RhoRHartMap_coord", 3, size, config.numPes,
          "TXYZ", 2, 1, 1, 1);
      mfc->dumpMapCoords(&RhoRHartImaptable[thisInstance.getPO()],
          thisInstance.getPO());
      delete mfc;
    }
  } //endif : ees_ext_on

  size[0] = config.nchareRhoG; size[1] = config.nchareHartAtmT;
  if(!config.loadMapFiles && size[0]*size[1] > RhoAvail->count())
    RhoAvail->reset();
  RhoGHartImaptable[thisInstance.getPO()].buildMap(size[0], size[1]);
  maploaded = 0;
  Timer=newtime;
  if(config.loadMapFiles) {
    mf = new MapFile("RhoGHartMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
    maploaded = mf->loadMap("RhoGHartMap", &RhoGHartImaptable[thisInstance.getPO()]);
  }

  if(!maploaded) {
    MapType3 *RhoRHartImaptablep = NULL;
    MapType2 *RhoRSImaptablep = &RhoRSImaptable[thisInstance.getPO()];
    if(sim->ees_eext_on) { //use the right RHart if it exists
      RhoRHartImaptablep = &RhoRHartImaptable[thisInstance.getPO()];
      RhoRSImaptablep = NULL;
    }
    RhoGHartMapTable RhoGHarttable(&RhoGHartImaptable[thisInstance.getPO()],
        RhoAvail, size[0], size[1], config.useCentroidMapRho, RhoRSImaptablep,
        RhoRHartImaptablep, excludePes);
  }
  newtime=CmiWallTimer();
  CkPrintf("RhoGHart created in %g\n", newtime-Timer);

  if(config.dumpMapFiles) {
    if(!mf) mf = new MapFile("RhoGHartMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
    mf->dumpMap(&RhoGHartImaptable[thisInstance.getPO()], thisInstance.getPO());
  }
  if(mf) { delete mf; mf = NULL;}
  if(config.dumpMapCoordFiles) {
    MapFile *mfc = new MapFile("RhoGHartMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
    mfc->dumpMapCoords(&RhoGHartImaptable[thisInstance.getPO()], thisInstance.getPO());
     delete mfc;
  }
  
  if(sim->ees_nloc_on)
    {
      PeList *nlexcludePes=NULL;
      if(config.useRhoExclusionMap)
	nlexcludePes=excludePes;
      else if(config.useReductionExclusionMap)
	nlexcludePes= new PeList(1, 1, UpeUsedByNLZ[thisInstance.proxyOffset]);
      else
	nlexcludePes= new PeList(1, 1, 0);
      if(config.excludePE0 && ! config.loadMapFiles)
	nlexcludePes->checkAndAdd(0);
      if(!nlexcludePes) nlexcludePes= new PeList(1, 1, 0);
      availGlobG->reset();
      RPPImaptable[thisInstance.getPO()].buildMap(sim->nstates, sim->ngrid_nloc_c);
      size[0]=sim->nstates; size[1]=sim->ngrid_nloc_c;
      delete mf;
      mf=NULL;
      mf = new MapFile("RPPMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
      if(config.loadMapFiles) {
	maploaded = mf->loadMap("RPPMap", &RPPImaptable[thisInstance.getPO()]);
      }
      Timer=newtime;
      if(!maploaded)
	{
	  RPPMapTable RPPtable= RPPMapTable(
						       &RPPImaptable[thisInstance.getPO()],
						       availGlobG, nlexcludePes,
						       config.nstates,  sim->ngrid_nloc_c,
						       config.Rstates_per_pe, boxSize,
						       config.useCuboidMap,
						       config.nchareG,
						       &GSImaptable[thisInstance.getPO()]);
	}
      newtime=CmiWallTimer();
      CkPrintf("RPPTable created in %g\n", newtime-Timer);
      if(config.dumpMapFiles) {
	if(!mf)	mf = new MapFile("RPPMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
	mf->dumpMap(&RPPImaptable[thisInstance.getPO()], thisInstance.getPO());
	delete mf;
	mf=NULL;
      }
      if(config.dumpMapCoordFiles) {
	MapFile *mfc = new MapFile("RPPMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
	mfc->dumpMapCoords(&RPPImaptable[thisInstance.getPO()], thisInstance.getPO());
	delete mfc;
      }
    }


  //create map tables for Y pencils
  //TODO: to be changed later so that we use plane decomposition
  Timer=newtime;
  if(1 || config.nchareRhoInter_x != 1) {
    size[0] = config.nchareRhoInter_x; size[1] = config.nchareRhoInter_z;
    for(int loop_off = 0; loop_off < 3; loop_off++) {
      if(!config.loadMapFiles && size[0]*size[1] > RhoAvail->count())
        RhoAvail->reset();
      maploaded = 0;
      RhoYPencilImaptable[loop_off][thisInstance.getPO()].buildMap(size[0], size[1]);
      char fileName[256];
      sprintf(fileName, "DivRho_%d", loop_off);
      if(config.loadMapFiles) {
        mf = new MapFile(fileName, 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
        maploaded = mf->loadMap(fileName,
            &RhoYPencilImaptable[loop_off][thisInstance.getPO()]);
      }
      if(!maploaded) {
        RhoYPencilMapTable RhoYPenciltable(
            &RhoYPencilImaptable[loop_off][thisInstance.getPO()], RhoAvail,
            size[0], size[1], config.useCentroidMapRho,
            &RhoRSImaptable[thisInstance.getPO()], excludePes, loop_off);
      }
      if(config.dumpMapFiles) {
        if(!mf) {
          mf = new MapFile(fileName, 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
        }
        mf->dumpMap(&RhoYPencilImaptable[loop_off][thisInstance.getPO()],
            thisInstance.getPO());
      }
      if(mf) { delete mf; mf=NULL;}
      if(config.dumpMapCoordFiles) {
        MapFile *mfc = new MapFile(strncat(fileName,"_coord", 256), 2, size,
            config.numPes, "TXYZ", 2, 1, 1, 1);
        mfc->dumpMapCoords(&RhoYPencilImaptable[loop_off][thisInstance.getPO()],
            thisInstance.getPO());
        delete mfc;
      }
    }//loop_off - fft of three divRhos
  }//nchareRhoInter_x != 1 -> plane decomposition
  //TODO: to be changed later so that we use plane decomposition
  if(1 || config.nchareHartInter_x != 1) {
    //map table for hart Y pencil
    size[0] = config.nchareHartInter_x; size[1] = config.nchareHartInter_z;
    if(!config.loadMapFiles && size[0]*size[1] > RhoAvail->count())
      RhoAvail->reset();
    maploaded = 0;
    RhoHartYPencilImaptable[thisInstance.getPO()].buildMap(size[0], size[1]);
    char fileName[256];
    sprintf(fileName, "RhoHartYPencil");
    if(config.loadMapFiles) {
      mf = new MapFile(fileName, 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
      maploaded = mf->loadMap(fileName,
          &RhoHartYPencilImaptable[thisInstance.getPO()]);
    }
    if(!maploaded) {
      RhoYPencilMapTable RhoYPenciltable(
          &RhoHartYPencilImaptable[thisInstance.getPO()], RhoAvail,
          size[0], size[1], config.useCentroidMapRho,
          &RhoRSImaptable[thisInstance.getPO()], excludePes, 3);
    }
    if(config.dumpMapFiles) {
      if(!mf) {
        mf = new MapFile(fileName, 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
      }
      mf->dumpMap(&RhoHartYPencilImaptable[thisInstance.getPO()],
          thisInstance.getPO());
    }
    if(mf) { delete mf; mf=NULL;}
    if(config.dumpMapCoordFiles) {
      MapFile *mfc = new MapFile(strncat(fileName,"_coord", 256), 2, size,
          config.numPes, "TXYZ", 2, 1, 1, 1);
      mfc->dumpMapCoords(&RhoHartYPencilImaptable[thisInstance.getPO()],
          thisInstance.getPO());
      delete mfc;
    }
  }// nchareHartInter_x != 1 -> plane decomposition
  newtime=CmiWallTimer();
  CkPrintf("RhoYPencils created in %g\n", newtime-Timer);
  if(sim->ees_eext_on) {
    Timer=newtime;
    if(1 || config.nchareAtmSFInter_x != 1) {
      size[0] = config.nchareAtmSFInter_x; size[1] = config.nchareAtmSFInter_z;
      for(int loop_off = 0; loop_off < config.nchareHartAtmT + 1; loop_off++) {
        if(!config.loadMapFiles && size[0]*size[1] > RhoAvail->count())
          RhoAvail->reset();
        maploaded = 0;
        AtmSFYPencilImaptable[loop_off][thisInstance.getPO()].buildMap(size[0], size[1]);
        char fileName[256];
        sprintf(fileName, "AtmSF_%d", loop_off);
        if(config.loadMapFiles) {
          mf = new MapFile(fileName, 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
          maploaded = mf->loadMap(fileName,
              &AtmSFYPencilImaptable[loop_off][thisInstance.getPO()]);
        }
        if(!maploaded) {
          RhoYPencilMapTable RhoYPenciltable(
              &AtmSFYPencilImaptable[loop_off][thisInstance.getPO()], RhoAvail,
              size[0], size[1], config.useCentroidMapRho,
              &RhoRSImaptable[thisInstance.getPO()], excludePes, 4 + loop_off);
        }
        if(config.dumpMapFiles) {
          if(!mf) {
            mf = new MapFile(fileName, 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
          }
          mf->dumpMap(&AtmSFYPencilImaptable[loop_off][thisInstance.getPO()],
              thisInstance.getPO());
        }
        if(mf) { delete mf; mf=NULL;}
        if(config.dumpMapCoordFiles) {
          MapFile *mfc = new MapFile(strncat(fileName,"_coord", 256), 2, size,
              config.numPes, "TXYZ", 2, 1, 1, 1);
          mfc->dumpMapCoords(&AtmSFYPencilImaptable[loop_off][thisInstance.getPO()],
              thisInstance.getPO());
          delete mfc;
        }
      }//loop_off - fft of nchareHartAtmT + 1
    }//nchareAtmSFInter_x != 1 -> plane decomposition
    newtime=CmiWallTimer();
    CkPrintf("AtmSFYPencils created in %g\n", newtime-Timer);
  }//end of if ees_eext_on
}

//============================================================================
//============================================================================
//Create the maps for the subsequent uber instances
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void build_uber_maps(CPcharmParaInfo *sim, UberCollection thisInstance)
{
  int x=mapOffsets[numInst].getx();
  int y=mapOffsets[numInst].gety();
  int z=mapOffsets[numInst].getz();

  if(CkNumPes() == 1) {
    GSImaptable[numInst]=GSImaptable[0];
    RSImaptable[numInst]=RSImaptable[0];
    if(sim->ees_nloc_on) RPPImaptable[numInst]=RPPImaptable[0];
    RhoRSImaptable[numInst]=RhoRSImaptable[0];
    RhoGSImaptable[numInst]=RhoGSImaptable[0];
    if(sim->ees_eext_on)  RhoRHartImaptable[numInst]=RhoRHartImaptable[0];
    RhoGHartImaptable[numInst]=RhoGHartImaptable[0];
    RhoYPencilImaptable[0][numInst]=RhoYPencilImaptable[0][0];
    RhoYPencilImaptable[1][numInst]=RhoYPencilImaptable[1][0];
    RhoYPencilImaptable[2][numInst]=RhoYPencilImaptable[2][0];
    RhoHartYPencilImaptable[numInst]=RhoHartYPencilImaptable[0];
    if(sim->ees_eext_on)  {
      for(int loop_off = 0; loop_off < config.nchareHartAtmT + 1; loop_off++) {
        AtmSFYPencilImaptable[loop_off][numInst]=AtmSFYPencilImaptable[loop_off][0];
      }
    }
  } else {
    GSImaptable[numInst].translate(&GSImaptable[0], x,y,z,
        config.torusMap==1);
    RSImaptable[numInst].translate(&RSImaptable[0], x,y,z,
        config.torusMap==1);
    if(sim->ees_nloc_on) RPPImaptable[numInst].translate(&RPPImaptable[0],
        x,y,z, config.torusMap==1);
    RhoRSImaptable[numInst].translate(&RhoRSImaptable[0], x,y,z,
        config.torusMap==1);
    RhoGSImaptable[numInst].translate(&RhoGSImaptable[0], x,y,z,
        config.torusMap==1);
    if(sim->ees_eext_on) RhoRHartImaptable[numInst].translate(
        &RhoRHartImaptable[0], x,y,z, config.torusMap==1);
    RhoGHartImaptable[numInst].translate(&RhoGHartImaptable[0], x,y,z,
        config.torusMap==1);
    for(int loop_off = 0; loop_off < 3; loop_off++) {
      RhoYPencilImaptable[loop_off][numInst].translate(
          &RhoYPencilImaptable[loop_off][0], x, y, z, config.torusMap==1);
    }
    RhoHartYPencilImaptable[numInst].translate(
        &RhoHartYPencilImaptable[0], x, y, z, config.torusMap==1);
    if(sim->ees_eext_on)  {
      for(int loop_off = 0; loop_off < config.nchareHartAtmT + 1; loop_off++) {
        AtmSFYPencilImaptable[loop_off][numInst].translate(
            &AtmSFYPencilImaptable[loop_off][0], x, y, z, config.torusMap==1);
      }
    }
  }
  fillInPeUsedBy(sim, thisInstance);
}


//============================================================================
//============================================================================
//Create the PIBeadAtoms array
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void init_PIBeads(CPcharmParaInfo *sim, UberCollection thisInstance)
{

  if(thisInstance.idxU.x>0)
    { // the set of chares being created is for a non-zero PI all PI
      // use the same PIBeadAtoms we simply direct the proxyoffset here to
      // the one for the 0th bead.
      CkPrintf("Constructing PIMD Bead proxies for non zero instances\n");
      UberCollection zeroBeadInstance=thisInstance;
      zeroBeadInstance.idxU.x=0;
      int proxyOffset=zeroBeadInstance.setPO();
      UPIBeadAtomsProxy.push_back(UPIBeadAtomsProxy[proxyOffset]);
    }
  else
    {
      int natm = sim->natm_tot;
      CkPrintf("Constructing PIMD Bead array\n");
      CkArrayOptions opts(config.numBeadAtomChares);
      UPIBeadAtomsProxy.push_back( CProxy_PIBeadAtoms::ckNew(thisInstance,config.UberImax,natm,opts));
    }

  //TODO: we should have a map for this array, the default scheme would
  //result in these chares being within the 0th Bead partition which is
  //not an optimal choice
}
//============================================================================

//Create the array elements for the GSpace, Particle and Real Space planes
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void init_state_chares(int natm_nl,int natm_nl_grp_max,int numSfGrps,
		       int doublePack, CPcharmParaInfo *sim,
		       UberCollection thisInstance)
//============================================================================
{ //begin routine
  //============================================================================
  /*
   * Set up the map variables used to control the location of the
   * 2d-array elements over processors
   */
  //============================================================================
  // Useful Local variables


  int ngrida          = sim->sizeX;
  int ngridb          = sim->sizeY;
  int ngridc          = sim->sizeZ;

  int ngridaNl        = sim->ngrid_nloc_a;
  int ngridbNl        = sim->ngrid_nloc_b;
  int ngridcNl        = sim->ngrid_nloc_c;
  int ees_nonlocal_on = sim->ees_nloc_on;

  int ngridaEext      = sim->ngrid_eext_a;
  int ngridbEext      = sim->ngrid_eext_b;
  int ngridcEext      = sim->ngrid_eext_c;
  int ees_eext_on     = sim->ees_eext_on;

  int nchareG         = sim->nchareG;
  int nchareR         = sim->sizeZ;
  int nchareRPP       = ngridcNl;

  int numIterNL       = sim->nlIters;

  int nchareRhoG      = config.nchareRhoG;
  int nchareGHart     = nchareRhoG;
  int nchareRHart     = config.nchareRhoRHart_x * config.nchareRhoRHart_y;

  int Rstates_per_pe  = config.Rstates_per_pe;
  int Gstates_per_pe  = config.Gstates_per_pe;
  int sGrainSize      = config.sGrainSize;
  int numChunks       = config.numChunks;
  int nkpoint         = sim->nkpoint;

  //============================================================================

  double Timer = CmiWallTimer();;
  if(thisInstance.idxU.y > 0 || thisInstance.idxU.s > 0) {
    // the set of chares being created is for a non-zero kpoint
    // all k-points and spins use the same atoms and energies
    // we simply direct the proxyoffset here to the one for
    // the 0th kpoint

    UberCollection zeroKpointInstance = thisInstance;
    zeroKpointInstance.idxU.y = 0;
    zeroKpointInstance.idxU.s = 0;
    int proxyOffset = zeroKpointInstance.setPO();
    // At some future point we could split the cache proxy so that
    // it could be shared across pretty much any kind of instance.
    // Currently we need different ones for beads and tempers due to
    // some atom bits that are blended in with the rest of the cache
    // stuff due to modularity failure in its design.  If eesCache
    // proxy eats all your memory, complain to Glenn.
    UeesCacheProxy.push_back(UeesCacheProxy[proxyOffset]);
  } else {
    UeesCacheProxy.push_back(CProxy_eesCache::ckNew(nchareRPP, nchareG,
          sim->nstates, nkpoint, thisInstance));
  }

  if(firstInstance) CkPrintf("created eescache proxy in %f\n", 
      CmiWallTimer() - Timer);
  Timer = CmiWallTimer();

  int nchareRRhoTot  = config.nchareRhoR_x * config.nchareRhoR_y;
  int nchareRHartTot = config.nchareRhoRHart_x * config.nchareRhoRHart_y;

  int *numGState     = sim->nlines_per_chareG;
  int *numGNL        = sim->nlines_per_chareG;
  int *numGRho       = NULL; //sim->nlines_per_chareRhoG;
  int *numGEext      = NULL; //sim->nlines_per_chareRhoGEext;

  if(firstInstance) {
    int *numRXState      = new int [nchareR];
    int *numRYState      = new int [nchareR];
    int *numRYStateLower = new int [nchareR];
    int *numRXNL         = new int [nchareRPP];
    int *numRYNL         = new int [nchareRPP];
    int *numRYNLLower    = new int [nchareRPP];
    int nplane_x_use     = sim->nplane_x;
    if(!config.doublePack){nplane_x_use = (nplane_x_use+1)/2;}
    for(int i=0;i<nchareR;i++){
      numRXState[i]      = sim->sizeY;
      numRYState[i]      = nplane_x_use;
      numRYStateLower[i] = nplane_x_use - 1;
    }//endif
    for(int i=0;i<nchareRPP;i++){
      numRXNL[i]         = ngridbNl;
      numRYNL[i]         = nplane_x_use;
      numRYNLLower[i]    = nplane_x_use - 1;
    }//endfor
    int *numRXRho      = NULL; //new int [nchareRRhoTot];
    int *numRYRho      = NULL; //new int [nchareRRhoTot];
    int *numRXEext     = NULL; //new int [nchareRHartTot];
    int *numRYEext     = NULL; //new int [nchareRHartTot];
    int *numSubGx      = sim->numSubGx;
    /*create_Rho_fft_numbers(nchareRRhoTot, nchareRHart, 1,
        sim->nplane_rho_x, sim->sizeY, ngridbEext,
        numRXRho, numRYRho, numRXEext, numRYEext, numSubGx);*/
    UfftCacheProxy.push_back(CProxy_FFTcache::ckNew(
          sim->sizeX,sim->sizeY,sim->sizeZ,
          ngridaEext,ngridbEext,ngridcEext,ees_eext_on,
          ngridaNl,  ngridbNl,  ngridcNl,  ees_nonlocal_on,
          sim->nlines_max, sim->nlines_max_rho,
          config.nchareG,nchareR,
          config.nchareG,nchareRPP,
          0,    nchareR,    0,
          0,   nchareRHart, 0,
          numGState,     numRXState, numRYState,numRYStateLower,
          numGNL,        numRXNL,    numRYNL, numRYNLLower,
          numGRho,       numRXRho,   numRYRho,
          numGEext,      numRXEext,  numRYEext,
          config.fftopt,config.fftprogresssplitReal,config.fftprogresssplit,
          1, thisInstance));
    CkPrintf("created fftcache proxy in %lf\n", CmiWallTimer() - Timer);
    Timer = CmiWallTimer(); 
    delete [] numRXState;
    delete [] numRYState;
    delete [] numRYStateLower;
    delete [] numRXNL;
    delete [] numRYNL;
    delete [] numRYNLLower;
    delete [] numRXRho;
    delete [] numRYRho;
    delete [] numRXEext;
    delete [] numRYEext;
  } else { // direct to 0th proxy
    UfftCacheProxy.push_back(UfftCacheProxy[0]);
  }

  //============================================================================
  // Instantiate the Chares with placement determined by the maps

  //---------------------------------------------------------------------------
  // state g-space

  PRINT_LINE_STAR;
  PRINTF("Building G-space (%d %d) and R-space (%d %d/%d) state Chares\n",
	 sim->nstates, nchareG, sim->nstates, nchareR, nchareRPP);
  PRINT_LINE_DASH;printf("\n");

  // there is only one IntMap per chare type, but each instance has
  // its own map group
  CProxy_GSMap gsMap = CProxy_GSMap::ckNew(thisInstance);
  //  CkArrayOptions gSpaceOpts(nstates,nchareG);
  /**@}*/
  /** addtogroup GSpaceState */
  /**@{*/
  CkArrayOptions gSpaceOpts(sim->nstates,nchareG);
  std::string forwardname("GSpaceForward");
  std::ostringstream fwdstrm;
  fwdstrm << forwardname << "." << thisInstance.idxU.x << "." << thisInstance.idxU.y << "." << thisInstance.idxU.z;
  int gforward=keeperRegister(fwdstrm.str());
  std::string backwardname("GSpaceBackward");
  std::ostringstream bwdstrm;
  bwdstrm << backwardname << "." << thisInstance.idxU.x << "." << thisInstance.idxU.y << "." << thisInstance.idxU.z;
  int gbackward=keeperRegister(bwdstrm.str());
  gSpaceOpts.setMap(gsMap);
  gSpaceOpts.setAnytimeMigration(false);
  gSpaceOpts.setStaticInsertion(true);
  UgSpacePlaneProxy.push_back(CProxy_CP_State_GSpacePlane::ckNew(sim->sizeX, 1, 1,
        sGrainSize, gforward, gbackward, thisInstance, gSpaceOpts));
  UgSpacePlaneProxy[thisInstance.proxyOffset].doneInserting();
  CkPrintf("UgSpacePlaneProxy created in %f\n", CmiWallTimer() - Timer);
  Timer = CmiWallTimer(); 
  /**@}*/
  //--------------------------------------------------------------------------------
  // Bind the GSpaceDriver array to the GSpacePlane array so that they migrate together
  CkArrayOptions gspDriverOpts(sim->nstates,nchareG);
  gspDriverOpts.setAnytimeMigration(false);
  gspDriverOpts.setStaticInsertion(true);

  gspDriverOpts.bindTo(UgSpacePlaneProxy[thisInstance.proxyOffset]);
  UgSpaceDriverProxy.push_back( CProxy_GSpaceDriver::ckNew(thisInstance,gspDriverOpts) );
  UgSpaceDriverProxy[thisInstance.proxyOffset].doneInserting();
  CkPrintf("UgSpaceDriverProxy created in %f\n", CmiWallTimer() - Timer);
  /**@}*/
  //--------------------------------------------------------------------------------
  // We bind the particlePlane array to the gSpacePlane array migrate together
  /** addtogroup Particle */
  /**@{*/
  //  CkArrayOptions particleOpts(nstates,nchareG);
  CkArrayOptions particleOpts(sim->nstates,nchareG);
  particleOpts.setAnytimeMigration(false);
  particleOpts.setStaticInsertion(true);
  particleOpts.setMap(gsMap); // the maps for both the arrays are the same
  particleOpts.bindTo(UgSpacePlaneProxy[thisInstance.proxyOffset]);
  UparticlePlaneProxy.push_back(CProxy_CP_State_ParticlePlane::ckNew(
        nchareG, sim->sizeY, sim->sizeZ,ngridaNl,ngridbNl,ngridcNl,
        1, numSfGrps, natm_nl, natm_nl_grp_max, sim->nstates,
        nchareG, Gstates_per_pe, numIterNL, ees_nonlocal_on,
        thisInstance,  particleOpts));
  UparticlePlaneProxy[thisInstance.proxyOffset].doneInserting();
  CkPrintf("UparticlePlaneProxy created in %f\n", CmiWallTimer() - Timer);
  Timer = CmiWallTimer(); 

  //---------------------------------------------------------------------------
  // state r-space

  CProxy_RSMap rsMap= CProxy_RSMap::ckNew(thisInstance);
  CkArrayOptions realSpaceOpts(sim->nstates,nchareR);
  realSpaceOpts.setMap(rsMap);
  realSpaceOpts.setAnytimeMigration(false);
  realSpaceOpts.setStaticInsertion(true);

  /**@}*/
  int rforward=keeperRegister(std::string("RealSpaceForward"));
  int rbackward=keeperRegister(std::string("RealSpaceBackward"));
  /** \addtogroup RealSpaceState */
  /**@{*/
  UrealSpacePlaneProxy.push_back( CProxy_CP_State_RealSpacePlane::ckNew(1, 1,
  ngrida, ngridb, ngridc, rforward, rbackward, thisInstance, realSpaceOpts));
  UrealSpacePlaneProxy[thisInstance.proxyOffset].doneInserting();
  CkPrintf("UrealSpacePlaneProxy created in %f\n", CmiWallTimer() - Timer);
  Timer = CmiWallTimer(); 

  //--------------------------------------------------------------------------------
  // state r-particleplane

  PRINT_LINE_DASH;
  PRINTF("Completed G-space/R-space state chare array build\n");
  PRINT_LINE_STAR;printf("\n");

  //----------------------------------------------------------------------------
}//end routine
//============================================================================

/** /addtogroup Particle */
/**@{*/
//============================================================================
// Creating arrays CP_StateRealParticlePlane
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void init_eesNL_chares(int natm_nl,int natm_nl_grp_max,
		       int doublePack, PeList *exclusion, CPcharmParaInfo *sim,
		       UberCollection thisInstance)
//============================================================================
{ //begin routine
  //============================================================================
  /*
   * Set up the map variables used to control the location of the
   * 2d-array elements over processors
   */
  //============================================================================
  // Useful Local variables


  int ngridaNl        = sim->ngrid_nloc_a;
  int ngridbNl        = sim->ngrid_nloc_b;
  int ngridcNl        = sim->ngrid_nloc_c;
  int ees_nonlocal_on = sim->ees_nloc_on;

  int nchareG         = sim->nchareG;
  int nchareRPP       = ngridcNl;

  int numIterNL       = sim->nlIters;
  int zmatSizeMax     = sim->nmem_zmat_max;

  PeList *nlexcludePes;

  int Rstates_per_pe  = config.Rstates_per_pe;
  availGlobG->reset();
  double newtime=CmiWallTimer();
  CProxy_RPPMap rspMap= CProxy_RPPMap::ckNew(thisInstance);
  CkArrayOptions pRealSpaceOpts(sim->nstates,ngridcNl);
  pRealSpaceOpts.setMap(rspMap);
  pRealSpaceOpts.setAnytimeMigration(false);
  pRealSpaceOpts.setStaticInsertion(true);
  UrealParticlePlaneProxy.push_back(CProxy_CP_State_RealParticlePlane::ckNew(
									     ngridaNl,ngridbNl,ngridcNl,
									     numIterNL,zmatSizeMax,Rstates_per_pe,
									     nchareG,ees_nonlocal_on, thisInstance,
									     pRealSpaceOpts));
  UrealParticlePlaneProxy[thisInstance.proxyOffset].doneInserting();
  printf("\n");
  PRINT_LINE_DASH;
  PRINTF("Completed RealParticle chare array build\n");
  PRINT_LINE_STAR;printf("\n");

}
/**@}*/

/** \addtogroup Density */
/**@{*/
//============================================================================
// Creating arrays CP_Rho_GSpacePlane, CP_Rho_GSpacePlaneHelper
// and CP_Rho_RealSpacePlane
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
int init_rho_chares(CPcharmParaInfo *sim, UberCollection thisInstance)
//============================================================================
{//begin routine
  //============================================================================
  /*
   * create the array for real-space densities (two-dimensional chare array)
   */
  //============================================================================
  //  Chare array sizes and offsets
  int ees_eext_on     = sim->ees_eext_on;
  int ees_nonlocal_on = sim->ees_nloc_on;
  int natmTyp         = sim->natm_typ;
  int nchareRhoG      = config.nchareRhoG; //use this as the number of RhoG / XPencils
  int nchareRhoR_x    = config.nchareRhoR_x; //dims for ZPencils/RhoR chares
  int nchareRhoR_y    = config.nchareRhoR_y;
  int nchareRhoInter_x   = config.nchareRhoInter_x; //dims for YPencils
  int nchareRhoInter_z   = config.nchareRhoInter_z;
  int nchareHartAtmT  = config.nchareHartAtmT;
  int nchareRhoRHart_x  = config.nchareRhoRHart_x;
  int nchareRhoRHart_y = config.nchareRhoRHart_y;
  int nchareHartInter_x = config.nchareHartInter_x;
  int nchareHartInter_z = config.nchareHartInter_z;
  int nchareAtmSFInter_x = config.nchareAtmSFInter_x;
  int nchareAtmSFInter_z = config.nchareAtmSFInter_z;

  if(thisInstance.idxU.y > 0)
  { // the set of chares being created is for a non-zero kpoint
    // all k-points use the same rho, therefore we do not make new
    // chares, we simply direct the proxyoffset here to the one for
    // the 0th kpoint and return

    UberCollection zeroKpointInstance = thisInstance;
    zeroKpointInstance.idxU.y = 0;
    int proxyOffset = zeroKpointInstance.setPO();
    UrhoRealProxy.push_back(UrhoRealProxy[proxyOffset]);
    UrhoGProxy.push_back(UrhoGProxy[proxyOffset]);
    UrhoGHartExtProxy.push_back(UrhoGHartExtProxy[proxyOffset]);
    if(ees_eext_on){
      UrhoRHartExtProxy.push_back(UrhoRHartExtProxy[proxyOffset]);
    }
    Urho_fft_xProxy.push_back(Urho_fft_xProxy[proxyOffset]);
    Urho_fft_yProxy.push_back(Urho_fft_yProxy[proxyOffset]);
    Urho_fft_zProxy.push_back(Urho_fft_zProxy[proxyOffset]);
    Urho_fft_hartProxy.push_back(Urho_fft_hartProxy[proxyOffset]);
    UlsRhoGProxy.push_back(UlsRhoGProxy[proxyOffset]);
    UlsRhoRealProxy.push_back(UlsRhoRealProxy[proxyOffset]);
    if(ees_eext_on){
      for(int loop_off = 0; loop_off < config.nchareHartAtmT + 1; loop_off++) {
        Urho_fft_atmSFProxy[loop_off].push_back(
            Urho_fft_atmSFProxy[loop_off][proxyOffset]);
        if(config.nchareHartAtmT > 1) {
          Urhart_sectionProxy[loop_off].push_back(
              Urhart_sectionProxy[loop_off][proxyOffset]);
          Ughart_sectionProxy[loop_off].push_back(
              Ughart_sectionProxy[loop_off][proxyOffset]);
        }
      }
    }
    return 1;
  }

  //============================================================================
  // Output to the screen
  CkGroupID rho_fft_maps[3][3]; //map groups for fft
  CkGroupID hart_fft_maps[3];
  CkGroupID atmSF_fft_maps[config.nchareHartAtmT + 1][3];

  CkPrintf("Creating density RhoR(%d,%d), RhoG(%d) RhoGHart(%d,%d), "
      "RhoRHart(%d,%d,%d)\n", nchareRhoR_x, nchareRhoR_y, nchareRhoG,
      nchareRhoG, nchareHartAtmT, nchareRhoRHart_x, nchareRhoRHart_y,
      nchareHartAtmT);
  /* Map for RhoR and ZPencil */
  CProxy_RhoRSMap rhorsMap = CProxy_RhoRSMap::ckNew(thisInstance);
  CkArrayOptions rhorsOpts(nchareRhoR_x, nchareRhoR_y);
  rhorsOpts.setMap(rhorsMap);
  rhorsOpts.setAnytimeMigration(false);
  rhorsOpts.setStaticInsertion(true);
  //RhoR has to be colocated with the Z pencil
  rho_fft_maps[0][0] = rho_fft_maps[1][0] = rho_fft_maps[2][0] =  rhorsMap;
  hart_fft_maps[0] = rhorsMap;

  /* Map for RhoG and XPencil */
  CProxy_RhoGSMap rhogsMap = CProxy_RhoGSMap::ckNew(thisInstance, 0);
  CkArrayOptions rhogsOpts(nchareRhoG);
  rhogsOpts.setMap(rhogsMap);
  rhogsOpts.setAnytimeMigration(false);
  rhogsOpts.setStaticInsertion(true);
  CProxy_RhoGSMap rhogsMap_pencil = CProxy_RhoGSMap::ckNew(thisInstance, 1);
  rho_fft_maps[0][2] = rho_fft_maps[1][2] = rho_fft_maps[2][2] = rhogsMap_pencil;

  //---------------------------------------------------------------------------
  // rho GHart

  CProxy_RhoGHartMap rhogHartMap = CProxy_RhoGHartMap::ckNew(thisInstance,
      0, 1, -1);
  CkArrayOptions rhoghartOpts(nchareRhoG, nchareHartAtmT);
  rhoghartOpts.setMap(rhogHartMap);
  rhoghartOpts.setAnytimeMigration(false);
  rhoghartOpts.setStaticInsertion(true);
  int fixIndex = 0;
  if(ees_eext_on && config.nchareHartAtmT > 1) {
    fixIndex = 1;
  }
  CProxy_RhoGHartMap rhogHartMap_pencil = CProxy_RhoGHartMap::ckNew(thisInstance,
      1, 0, fixIndex);
  hart_fft_maps[2] = rhogHartMap_pencil;

  //---------------------------------------------------------------------------
  // rho RHart
  CkArrayOptions rhorhartOpts(nchareRhoRHart_x, nchareRhoRHart_y, nchareHartAtmT);
  if(ees_eext_on) {
    CProxy_RhoRHartMap rhorHartMap = CProxy_RhoRHartMap::ckNew(thisInstance, 0,
        1, -1);
    rhorhartOpts.setAnytimeMigration(false);
    rhorhartOpts.setStaticInsertion(true);
    rhorhartOpts.setMap(rhorHartMap);

    for(int type = 0; type < config.nchareHartAtmT; type++) {
      CProxy_RhoRHartMap rhorHartMap_pencil = CProxy_RhoRHartMap::ckNew(
          thisInstance, 0, 1, type);
      CProxy_RhoGHartMap rhogHartMap_pencil = CProxy_RhoGHartMap::ckNew(
          thisInstance, 1, 0, type);
      CProxy_RhoYPencilMap rhoYPencilMap = CProxy_RhoYPencilMap::ckNew(
          thisInstance, 2, type);
      atmSF_fft_maps[type][0] = rhorHartMap_pencil;
      atmSF_fft_maps[type][1] = rhoYPencilMap;
      atmSF_fft_maps[type][2] = rhogHartMap_pencil;
    }

    CProxy_RhoYPencilMap rhoYPencilMap = CProxy_RhoYPencilMap::ckNew(
        thisInstance, 2, config.nchareHartAtmT);
    atmSF_fft_maps[config.nchareHartAtmT][0] = atmSF_fft_maps[0][0];
    atmSF_fft_maps[config.nchareHartAtmT][1] = rhoYPencilMap;
    atmSF_fft_maps[config.nchareHartAtmT][2] = atmSF_fft_maps[0][2];
  }

  //---------------------------------------------------------------------------
  // Y pencils
  CProxy_RhoYPencilMap rhoYPencilMap_x = CProxy_RhoYPencilMap::ckNew(thisInstance,
      0, 0);
  CProxy_RhoYPencilMap rhoYPencilMap_y = CProxy_RhoYPencilMap::ckNew(thisInstance,
      0, 1);
  CProxy_RhoYPencilMap rhoYPencilMap_z = CProxy_RhoYPencilMap::ckNew(thisInstance,
      0, 2);
  CProxy_RhoYPencilMap rhoYPencilMap_hart = CProxy_RhoYPencilMap::ckNew(thisInstance,
      1, 0);
  rho_fft_maps[0][1] = rhoYPencilMap_x;
  rho_fft_maps[1][1] = rhoYPencilMap_y;
  rho_fft_maps[2][1] = rhoYPencilMap_z;
  hart_fft_maps[1] = rhoYPencilMap_hart;

  //============================================================================
  // Instantiate the chares

  bool dummy = true;

  //--------------------------------------------------------------------------
  // insert rhoreal
  int rhokeeper= keeperRegister(std::string("Density"));
  UrhoRealProxy.push_back(CProxy_CP_Rho_RealSpacePlane::ckNew( rhokeeper,
        thisInstance, rhorsOpts));
  UrhoRealProxy[thisInstance.proxyOffset].doneInserting();

  UrhoRealProxy[thisInstance.proxyOffset].ckSetReductionClient(new
  CkCallback(CkIndex_InstanceController::printEnergyEexc(NULL),
  CkArrayIndex1D(thisInstance.proxyOffset),instControllerProxy));

  //--------------------------------------------------------------------------
  // insert rhog
  UrhoGProxy.push_back(CProxy_CP_Rho_GSpacePlane::ckNew( thisInstance,
        rhogsOpts));
  UrhoGProxy[thisInstance.proxyOffset].doneInserting();

  //--------------------------------------------------------------------------
  // insert rhoghart
  UrhoGHartExtProxy.push_back(CProxy_CP_Rho_GHartExt::ckNew(thisInstance,
							    rhoghartOpts));
  UrhoGHartExtProxy[thisInstance.proxyOffset].doneInserting();
  UrhoGHartExtProxy[thisInstance.proxyOffset].ckSetReductionClient(new
  CkCallback(CkIndex_InstanceController::printEnergyHart(NULL),
  CkArrayIndex1D(thisInstance.proxyOffset),instControllerProxy));

  //create sections if needed
  if(ees_eext_on && config.nchareHartAtmT > 1) {
    CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
    for(int type = 0; type < config.nchareHartAtmT; type++) {
      CProxySection_CP_Rho_GHartExt proxy =
        CProxySection_CP_Rho_GHartExt::ckNew(
            UrhoGHartExtProxy[thisInstance.proxyOffset].ckGetArrayID(),
            0, nchareRhoG - 1, 1, type, type, 1);
      proxy.ckSectionDelegate(mCastGrp);
      Ughart_sectionProxy[type].push_back(proxy);
    }
  }

  //--------------------------------------------------------------------------
  // insert rhoRhart
  if(ees_eext_on){
    UrhoRHartExtProxy.push_back(CProxy_CP_Rho_RHartExt::ckNew(thisInstance,
        rhorhartOpts));
    UrhoRHartExtProxy[thisInstance.proxyOffset].doneInserting();
    //create sections if needed
    if(config.nchareHartAtmT > 1) {
      CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
      for(int type = 0; type < config.nchareHartAtmT; type++) {
        CProxySection_CP_Rho_RHartExt proxy =
          CProxySection_CP_Rho_RHartExt::ckNew(
              UrhoRHartExtProxy[thisInstance.proxyOffset].ckGetArrayID(),
              0, nchareRhoRHart_x - 1, 1, 0, nchareRhoRHart_y - 1, 1, type, type, 1);
        proxy.ckSectionDelegate(mCastGrp);
        Urhart_sectionProxy[type].push_back(proxy);
      }
    }
  }//endif

  //===========================================================================
  //create fft instances
  double tpi, *hmati_first, cutoff = 8*simReadOnly.ecut;
  CPXCFNCTS::CP_fetch_hmati(&hmati_first, &tpi);
  double hmati[9];
  hmati[0] = hmati_first[3];
  hmati[3] = hmati_first[6];
  hmati[6] = hmati_first[9];
  hmati[1] = hmati_first[2];
  hmati[4] = hmati_first[5];
  hmati[7] = hmati_first[8];
  hmati[2] = hmati_first[1];
  hmati[5] = hmati_first[4];
  hmati[8] = hmati_first[7];
  CkCallback fft_callback(CkIndex_InstanceController::doneFFTCreation(NULL),
      CkArrayIndex1D(thisInstance.getPO()), instControllerProxy.ckGetArrayID());

  Urho_fft_xProxy.push_back(Charm_createFFT(sim->sizeZ, sim->sizeY, sim->sizeX,
        nchareRhoR_x, nchareRhoR_y, nchareRhoInter_x, nchareRhoInter_z,
        nchareRhoG, cutoff, hmati, RC, rho_fft_maps[0], fft_callback));
  Urho_fft_yProxy.push_back(Charm_createFFT(sim->sizeZ, sim->sizeY, sim->sizeX,
        nchareRhoR_x, nchareRhoR_y, nchareRhoInter_x, nchareRhoInter_z,
        nchareRhoG, cutoff, hmati, RC, rho_fft_maps[1], fft_callback));
  Urho_fft_zProxy.push_back(Charm_createFFT(sim->sizeZ, sim->sizeY, sim->sizeX,
        nchareRhoR_x, nchareRhoR_y, nchareRhoInter_x, nchareRhoInter_z,
        nchareRhoG, cutoff, hmati, RC, rho_fft_maps[2], fft_callback));
  Urho_fft_hartProxy.push_back(Charm_createFFT(sim->sizeZ, sim->sizeY,
        sim->sizeX, nchareRhoR_x, nchareRhoR_y, nchareHartInter_x,
        nchareHartInter_z, nchareRhoG, cutoff, hmati, RC, hart_fft_maps,
        fft_callback));

  if(ees_eext_on) {
    int miniGrid[3];
    miniGrid[0] = sim->sizeZ;
    miniGrid[1] = sim->sizeY;
    miniGrid[2] = sim->sizeX;
    for(int type = 0; type < config.nchareHartAtmT + 1; type++) {
      Urho_fft_atmSFProxy[type].push_back(Charm_createFFT(sim->ngrid_eext_c,
            sim->ngrid_eext_b, sim->ngrid_eext_a, nchareRhoRHart_x, nchareRhoRHart_y,
            nchareAtmSFInter_x, nchareAtmSFInter_z, nchareRhoG, cutoff,
            hmati, RC, atmSF_fft_maps[type], fft_callback, miniGrid));
    }
  }

  // Output to the screen
  // need to add maps for these, for now just let em default
  // IF some condition which triggers QMMM
  int nchareRhoGLSP = 1;
  int nchareRhoRealLSP = 1;
  int nchareRhoRealLSPsubplanes = 1;
  CkArrayOptions lspgspOpts(nchareRhoGLSP);
  CkArrayOptions lsprealOpts(nchareRhoRealLSP, nchareRhoRealLSPsubplanes);
  UlsRhoGProxy.push_back(CProxy_CP_LargeSP_RhoGSpacePlane::ckNew(thisInstance,lspgspOpts));
  UlsRhoRealProxy.push_back(CProxy_CP_LargeSP_RhoRealSpacePlane::ckNew(thisInstance,lsprealOpts));
  printf("\n");
  PRINT_LINE_DASH;
  PRINTF("Completed G-space/R-space Rho chare array build\n");
  PRINT_LINE_STAR;printf("\n");

  return 1;
  //===========================================================================
}//end routine
//============================================================================
/**@}*/


//============================================================================
// Get the atoms and the parainfo
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void control_physics_to_driver(UberCollection thisInstance, CPcharmParaInfo *sim){
  //============================================================================
  // make a group : create a proxy for the atom class and also a reduction client


  GENERAL_DATA *general_data = GENERAL_DATA::get();
  char *output_directory=general_data->gentempering_ctrl.output_directory;
  if(output_directory==NULL)
    {
      output_directory="TEMPER_OUT";
    }

  // Make  groups for the atoms and energies 

  if(thisInstance.idxU.y>0|| thisInstance.idxU.s>0)
    { // the set of chares being created is for a non-zero kpoint
      // all k-points use the same atoms and energies
      // we simply direct the proxyoffset here to the one for
      // the 0th kpoint
      // same holds true for spin
      UberCollection zeroKpointInstance=thisInstance;
      zeroKpointInstance.idxU.y=0;
      zeroKpointInstance.idxU.s=0;
      int proxyOffset=zeroKpointInstance.setPO();
      UatomsCacheProxy.push_back(UatomsCacheProxy[proxyOffset]);
      UatomsComputeProxy.push_back(UatomsComputeProxy[proxyOffset]);
      UpScratchProxy.push_back(UpScratchProxy[proxyOffset]);
      UegroupProxy.push_back(UegroupProxy[proxyOffset]);
    }
  else
    {
      int ibead         = thisInstance.idxU.x;
      int itemper       = thisInstance.idxU.z;
      PhysicsAtomPosInit *PhysicsAtom  = new PhysicsAtomPosInit(ibead,itemper);
      int natm          = PhysicsAtom->natm_tot;
      int natm_nl       = PhysicsAtom->natm_nl;
      int len_nhc       = PhysicsAtom->len_nhc;
      int iextended_on  = PhysicsAtom->iextended_on;
      int cp_min_opt    = PhysicsAtom->cp_min_opt;
      int cp_wave_opt   = PhysicsAtom->cp_wave_opt;
      int cp_bomd_opt   = PhysicsAtom->cp_bomd_opt;
      int isokin_opt    = PhysicsAtom->isokin_opt;
      int cp_grimme     = PhysicsAtom->cp_grimme;
      double kT         = PhysicsAtom->kT;
      // kT is the temperature/BOLTZ
      Atom *atoms       = new Atom[natm];
      AtomNHC *atomsNHC = new AtomNHC[natm];
      CkPrintf("[%d] Temperature is %g\n",thisInstance.proxyOffset, kT*BOLTZ);

      // every instance with its own atoms (such as beads and tempers) needs its own
      // physcratchcache, because it is a horrible shared memory thing
      // that needs to be walled off
      UpScratchProxy.push_back(CProxy_PhysScratchCache::ckNew());
      PhysicsAtom->DriverAtomInit(natm,atoms,atomsNHC,ibead,itemper);
      UegroupProxy.push_back(CProxy_EnergyGroup::ckNew(thisInstance));
      // FIXME, this needs a real computation
      // also we need a real map. This was a naive, simple mapping scheme. We
      // can do better.
      int nChareAtoms=(config.numPesPerInstance<natm) ? config.numPesPerInstance : natm;

      if (firstInstance || config.simpleTopo) {
	// build the base atom map
	AtomImaptable[numInst].buildMap(nChareAtoms);
	availGlobG->reset();
	AtomMapTable aTable = AtomMapTable(&AtomImaptable[numInst], availGlobG,
					   numInst,nChareAtoms);


      } else {
        int x=mapOffsets[numInst].getx();
        int y=mapOffsets[numInst].gety();
        int z=mapOffsets[numInst].getz();
        AtomImaptable[numInst].translate(&AtomImaptable[0], x,y,z, config.torusMap==1);
      }

      CProxy_AtomComputeMap aMap = CProxy_AtomComputeMap::ckNew(thisInstance);
      CkArrayOptions atomOpts(nChareAtoms);
      atomOpts.setMap(aMap);
      atomOpts.setAnytimeMigration(false);
      atomOpts.setStaticInsertion(true);
      UatomsCacheProxy.push_back( CProxy_AtomsCache::ckNew(natm,natm_nl,
      							   atoms,thisInstance, std::string(output_directory)) );
      UatomsComputeProxy.push_back( CProxy_AtomsCompute::ckNew(natm,natm_nl,
							       len_nhc,
							       iextended_on,
							       cp_min_opt,cp_wave_opt,cp_bomd_opt,isokin_opt,cp_grimme,
							       kT,atoms,
							       atomsNHC,
							       nChareAtoms,
							       thisInstance,
							       atomOpts
							       ));
      delete [] atoms;
      delete [] atomsNHC;
      delete PhysicsAtom;
    }


  //----------------------------------------------------------------------------
}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
static CkReductionMsg *complexSum(int nMsg, CkReductionMsg **msgs){
  //============================================================================
  // get the size of the messages from the first message
  int matrixSize = msgs[0]->getSize() / sizeof(complex);
  complex *matrix = new complex[matrixSize];
  int i, m;
  for (m = 0; m < nMsg; m++) {
    // sanity check
    CkAssert(msgs[m]->getSize() == matrixSize * sizeof(complex));
    complex *data = (complex *) msgs[m]->getData();
    for (i = 0; i < matrixSize; i++)
      matrix[i] += data[i];
  }
  CkReductionMsg *ret = CkReductionMsg::buildNew(matrixSize *
						 sizeof(complex), matrix);
  delete [] matrix;
  return ret;
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void get_grp_params(int natm_nl, int numSfDups, int indexSfGrp, int planeIndex,
		    int *n_ret, int *istrt_ret, int *iend_ret)
//============================================================================
{// begin routine
  //============================================================================

  int n     = (natm_nl/numSfDups);
  int m     = (natm_nl % numSfDups);

  int istrt = n*indexSfGrp;
  if(indexSfGrp>=m){istrt += m;}
  if(indexSfGrp<m) {istrt += indexSfGrp;}
  if(indexSfGrp<m) {n++;}
  int iend  = n+istrt;
  if(numSfDups>natm_nl)
    if(m>=indexSfGrp)
      {
	n=0;
	istrt=natm_nl+1;
	iend=natm_nl;
	CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
	CkPrintf("Redundant DupSF for chare-G %d\n",planeIndex);
	CkPrintf("At present this hangs, so out you go\n");
	CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
	//         CkExit();
      }

  (*n_ret)     = n;
  (*istrt_ret) = istrt;
  (*iend_ret)  = iend;

  //---------------------------------------------------------------------------
}//end routine
//============================================================================


/** \addtogroup mapping */
/**@{*/
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
int atmGrpMap(int istart, int nsend, int listsize, int *listpe, int AtmGrp,
	      int dup, int planeIndex)
//============================================================================
{// begin routine
  //============================================================================

  int numSfDups=config.numSfDups;

  if(listsize <= numSfDups){ //dup list not unique
    return listpe[(dup % listsize)];
  }//endif

  if(nsend <= config.numSfGrps){  //atom group list not unique
    return listpe[(istart + (AtmGrp % nsend))];
  }//endif

  // nsend must be > numSfGrps and its list is unique
  return listpe[(istart + ((AtmGrp + config.numSfGrps * planeIndex)%nsend))];

  //---------------------------------------------------------------------------
}// end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
int gsprocNum(CPcharmParaInfo *sim, int state, int plane, int numInst) {
  int proc;
  proc = GSImaptable[numInst].get(state, plane);
  return(proc);
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void mapOutput()
//============================================================================
{//begin routine
  //============================================================================

  //----------------------------------------------------------------------------
}//end routine
//============================================================================

//============================================================================
/**@}*/

#ifdef CMK_BALANCED_INJECTION_API
#define DEFAULT_LOW_BI_VALUE 32;
uint16_t lowBIValue=DEFAULT_LOW_BI_VALUE;
uint16_t origBIValue;
#endif

void set_GNI_LOW_BI()
{
/* because gemini dies like a fish on a hot beach if you overuse the
 *    network */
#ifdef CMK_BALANCED_INJECTION_API
  uint16_t currval =ck_get_GNI_BIConfig();
  if(currval!=lowBIValue)
     {
      origBIValue=currval;
      ck_set_GNI_BIConfig(lowBIValue);
     }
  if(CmiMyNode()==0)
    CkPrintf("Balanced Injection initial value was %d, is now %d\n",origBIValue,lowBIValue);
#endif
}

void  PlatformSpecific::reset_BI(){
#ifdef CMK_BALANCED_INJECTION_API
    ck_set_GNI_BIConfig(origBIValue);
#endif
  }

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CkReduction::reducerType sumFastDoubleType;
void registersumFastDouble(void){
  sumFastDoubleType=CkReduction::addReducer(sumFastDouble);
}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// sum together matrices of doubles
// possibly faster than sum_double due to minimizing copies
// and using sameer's fastadd
//============================================================================
inline CkReductionMsg *sumFastDouble(int nMsg, CkReductionMsg **msgs){

  int size0=msgs[0]->getSize();
  int size=size0/sizeof(double);

  double *inmatrix;
  //  int progcount=0;

  double *ret=(double *)msgs[0]->getData();

  for(int i=1; i<nMsg;i++)
    {
      inmatrix=(double *) msgs[i]->getData();
      for(int d=0;d<size;d++)
	ret[d]+=inmatrix[d];
    }
  return CkReductionMsg::buildNew(size0,ret);
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void create_Rho_fft_numbers(int nchareR,int nchareRHart,int rhoRsubplanes,
			    int nplane, int sizeY, int ngridbHart,
			    int *numRXRho,int *numRYRho,int *numRXEext,int *numRYEext,
			    int *numSubGx)
//============================================================================
{//begin routine
  //============================================================================
  // Simple case  : everyone is the same size

  if(rhoRsubplanes==1){

    for(int i=0;i<nchareR;i++){
      numRXRho[i]  = sizeY;
      numRYRho[i]  = nplane;
    }//endfor

    for(int i=0;i<nchareRHart;i++){
      numRXEext[i] = ngridbHart;
      numRYEext[i]   = nplane;
    }//endfor

  }//endif

  //============================================================================
  // Subplane case : different sizes with subplanes label

  if(rhoRsubplanes>1){
    int div,rem;

    //---------------------------------------
    // how many x FFTs for Rho
    div  = (sizeY / rhoRsubplanes);
    rem  = (sizeY % rhoRsubplanes);
    for(int j=0;j<rhoRsubplanes;j++){
      int mySizeY = (j < rem ? div+1 : div);
      for(int i=0;i<nchareR;i++){
        int ind       = j*nchareR + i;
        numRXRho[ind] = mySizeY;
      }//endfor
    }//endfor

    //---------------------------------------
    // how many y FFTs for Rho
    for(int j=0;j<rhoRsubplanes;j++){
      int myNplane = numSubGx[j];
      for(int i=0;i<nchareR;i++){
        int ind       = j*nchareR + i;
        numRYRho[ind] = myNplane;
      }//endfor
    }//endfor

    //---------------------------------------
    // how many x FFTs for HartEExt EES
    div  = (ngridbHart / rhoRsubplanes);
    rem  = (ngridbHart % rhoRsubplanes);
    for(int j=0;j<rhoRsubplanes;j++){
      int myNgridb   = (j < rem ? div+1 : div);
      for(int i=0;i<nchareRHart;i++){
        int ind        = j*nchareRHart + i;
        numRXEext[ind] = myNgridb;
      }//endfor
    }//endfor

    //---------------------------------------
    // how many y FFTs for HartEExt EES
    for(int j=0;j<rhoRsubplanes;j++){
      int myNplane = numSubGx[j];
      for(int i=0;i<nchareRHart;i++){
        int ind        = j*nchareRHart + i;
        numRYEext[ind] = myNplane;
      }//endfor
    }//endfor

  }//endif

  //---------------------------------------------------------------------------
}//end routine
//============================================================================


void setTraceUserEvents()
{
  //============================================================================
  // Set user trace events for projections optimizations


  traceRegisterUserEvent("doRealFwFFT", doRealFwFFT_);

  traceRegisterUserEvent("GspaceFwFFT", GspaceFwFFT_);

  traceRegisterUserEvent("fwFFTGtoR0", fwFFTGtoR0_);
  traceRegisterUserEvent("fwFFTGtoRnot0", fwFFTGtoRnot0_);
  traceRegisterUserEvent("OrthoDGEMM1", OrthoDGEMM1_);
  traceRegisterUserEvent("GradCorrGGA", GradCorrGGA_);
  traceRegisterUserEvent("WhiteByrdFFTX", WhiteByrdFFTX_);
  traceRegisterUserEvent("doRealBwFFT", doRealBwFFT_);
  traceRegisterUserEvent("WhiteByrdFFTY", WhiteByrdFFTY_);
  traceRegisterUserEvent("WhiteByrdFFTZ", WhiteByrdFFTZ_);
  traceRegisterUserEvent("PostByrdfwFFTGtoR", PostByrdfwFFTGtoR_);
  traceRegisterUserEvent("RhoRtoGFFT", RhoRtoGFFT_);
  traceRegisterUserEvent("BwFFTRtoG", BwFFTRtoG_);
  traceRegisterUserEvent("OrthoDGEMM2", OrthoDGEMM2_);
  traceRegisterUserEvent("ByrdanddoFwFFTGtoR",ByrdanddoFwFFTGtoR_);
  traceRegisterUserEvent("eesHartExcG",eesHartExcG_);
  traceRegisterUserEvent("eesEwaldG",eesEwaldG_);
  traceRegisterUserEvent("eesAtmForcR",eesAtmForcR_);
  traceRegisterUserEvent("eesAtmBspline",eesAtmBspline_);
  traceRegisterUserEvent("eesZmatR",eesZmatR_);
  traceRegisterUserEvent("eesEnergyAtmForcR",eesEnergyAtmForcR_);
  traceRegisterUserEvent("eesProjG",eesProjG_);
  traceRegisterUserEvent("doNlFFTGtoR",doNlFFTGtoR_);
  traceRegisterUserEvent("doNlFFTRtoG",doNlFFTRtoG_);
  traceRegisterUserEvent("eesPsiForcGspace",eesPsiForcGspace_);
  traceRegisterUserEvent("enlMatrixCalc",enlMatrixCalc_);
  traceRegisterUserEvent("enlAtmForcCalc",enlAtmForcCalc_);
  traceRegisterUserEvent("enlForcCalc",enlForcCalc_);
  traceRegisterUserEvent("doEextFFTRtoG",doEextFFTRtoG_);
  traceRegisterUserEvent("doEextFFTGtoR",doEextFFTGtoR_);

  traceRegisterUserEvent("HartExcVksG",HartExcVksG_);
  traceRegisterUserEvent("divRhoVksGspace",divRhoVksGspace_);
  traceRegisterUserEvent("GspaceBwFFT", GspaceBwFFT_);
  traceRegisterUserEvent("DoFFTContribute", DoFFTContribute_);
  traceRegisterUserEvent("IntegrateModForces", IntegrateModForces_);
  traceRegisterUserEvent("Scalcmap", Scalcmap_);
  traceRegisterUserEvent("AcceptStructFact", AcceptStructFact_);
  traceRegisterUserEvent("doEextFFTGxtoRx", doEextFFTGxtoRx_);
  traceRegisterUserEvent("doEextFFTRytoGy", doEextFFTRytoGy_);
  traceRegisterUserEvent("doRhoFFTRytoGy", doRhoFFTRytoGy_);
  traceRegisterUserEvent("doRhoFFTGxtoRx", doRhoFFTGxtoRx_);

  traceRegisterUserEvent("GSProcnum", 10000);
  traceRegisterUserEvent("RSProcnum", 20000);
  traceRegisterUserEvent("SCProcnum", 30000);
  traceRegisterUserEvent("GHartAtmForcCopy",GHartAtmForcCopy_);
  traceRegisterUserEvent("GHartAtmForcSend",GHartAtmForcSend_);
  Ortho_UE_step2 = traceRegisterUserEvent("Ortho step 2");
  Ortho_UE_step3 = traceRegisterUserEvent("Ortho step 3");
  Ortho_UE_error = traceRegisterUserEvent("Ortho error");

}

void computeMapOffsets()
{
  int x, y, z, x1, y1, z1;
  // correction to accomodate multiple instances
  for(int thisInst=1; thisInst<config.numInstances; thisInst++) {
    x=thisInst*config.numPesPerInstance;
    y=0;
    z=0;
    mapOffsets[thisInst]=inttriple(x,y,z);
  }
}


void paircalcstartup(pc::pcConfig *cfgSymmPC, pc::pcConfig *cfgAsymmPC, CPcharmParaInfo *sim, int doublePack)
{



  // Stuff it with the actual configurations
  cfgSymmPC->isDynamics         = (sim->cp_min_opt==1)? false: true;
  cfgSymmPC->useComplexMath     = false;

  cfgSymmPC->numPlanes          = config.nchareG;
  cfgSymmPC->numStates          = sim->nstates;
  cfgSymmPC->grainSize          = config.sGrainSize;
  cfgSymmPC->orthoGrainSize     = config.orthoGrainSize;

  cfgSymmPC->conserveMemory     = config.conserveMemory;
  cfgSymmPC->isLBon             = config.lbpaircalc;

  cfgSymmPC->areBWTilesCollected= config.PCCollectTiles;
  cfgSymmPC->isBWstreaming      = config.PCstreamBWout;
  cfgSymmPC->isBWbarriered      = config.useBWBarrier;
  cfgSymmPC->shouldDelayBWsend  = config.PCdelayBWSend;
  cfgSymmPC->isInputMulticast   = !config.usePairDirectSend;
  cfgSymmPC->isOutputReduced    = !config.gSpaceSum;
  cfgSymmPC->inputSpanningTreeFactor = config.PCSpanFactor;

  cfgSymmPC->gemmSplitFWk       = config.gemmSplitFWk;
  cfgSymmPC->gemmSplitFWm       = config.gemmSplitFWm;
  cfgSymmPC->gemmSplitBW        = config.gemmSplitBW;



  // Configurations specific to the symmetric PC instance
  cfgSymmPC->isSymmetric        = true;
  // Turning phantoms off doesn't work with the current sparse-array
  // construction used in pcBuilder.C
  cfgSymmPC->arePhantomsOn      = true;
  cfgSymmPC->numChunks          = config.numChunksSym;
  cfgSymmPC->isDoublePackOn     = doublePack;
  cfgSymmPC->inputMsgPriority   = config.psipriority;
  cfgSymmPC->resultMsgPriority  = config.gsfftpriority;

  // Copy baseline parameters from symm, then override for asymm case
  *cfgAsymmPC=*cfgSymmPC;

  // Configurations specific to the asymmetric PC instance
  cfgAsymmPC->isSymmetric        = false;
  cfgAsymmPC->arePhantomsOn      = false;
  cfgAsymmPC->numChunks          = config.numChunksAsym;
  cfgAsymmPC->isDoublePackOn     = 0;
  cfgAsymmPC->inputMsgPriority   = config.lambdapriority;
  cfgAsymmPC->resultMsgPriority  = config.lambdapriority+2;

  // Configure the GSpace entry methods that the PCs will callback
  if(cfgSymmPC->isOutputReduced)
    {
      cfgSymmPC->gSpaceEP        = CkIndex_CP_State_GSpacePlane::acceptNewPsi ((CkReductionMsg*)NULL);
      cfgSymmPC->PsiVEP          = CkIndex_CP_State_GSpacePlane::acceptNewPsiV((CkReductionMsg*)NULL);
    }
  else
    {
      cfgSymmPC->gSpaceEP        = CkIndex_CP_State_GSpacePlane::acceptNewPsi ((partialResultMsg*)NULL);
      cfgSymmPC->PsiVEP          = CkIndex_CP_State_GSpacePlane::acceptNewPsiV((partialResultMsg*)NULL);
    }

  if(cfgAsymmPC->isOutputReduced)
    {
      cfgAsymmPC->gSpaceEP       = CkIndex_CP_State_GSpacePlane::acceptLambda ((CkReductionMsg*)NULL);
      cfgAsymmPC->PsiVEP         = 0;
    }
  else
    {
      cfgAsymmPC->gSpaceEP       = CkIndex_CP_State_GSpacePlane::acceptLambda ((partialResultMsg*)NULL);
      cfgAsymmPC->PsiVEP         = 0;
    }

#ifdef _CP_SUBSTEP_TIMING_
  //symmetric AKA Psi
  cfgSymmPC->forwardTimerID      = keeperRegister("Sym Forward");
  cfgSymmPC->backwardTimerID     = keeperRegister("Sym Backward");
  cfgSymmPC->beginTimerCB        = CkCallback(CkIndex_TimeKeeper::collectStart(NULL),0,TimeKeeperProxy);
  cfgSymmPC->endTimerCB          = CkCallback(CkIndex_TimeKeeper::collectEnd(NULL),0,TimeKeeperProxy);
  //asymmetric AKA Lambda AKA Gamma
  cfgAsymmPC->forwardTimerID     = keeperRegister("Asym Forward");
  cfgAsymmPC->backwardTimerID    = keeperRegister("Asym Backward");
  cfgAsymmPC->beginTimerCB       = CkCallback(CkIndex_TimeKeeper::collectStart(NULL),0,TimeKeeperProxy);
  cfgAsymmPC->endTimerCB         = CkCallback(CkIndex_TimeKeeper::collectEnd(NULL),0,TimeKeeperProxy);
#endif
}

void orthostartup( cp::ortho::orthoConfig *orthoCfg,  pc::pcConfig *cfgSymmPC, pc::pcConfig *cfgAsymmPC, CPcharmParaInfo *sim, PeListFactory *peList4PCmapping)
{
  // Blame Ram for ugly crime against readability.  More
  // redundant config objects and builders doesn't help
  // global clarity at all.
  // EJB: Sequestered in this function call to sweep it under the rug

  orthoCfg->isDynamics    = (sim->cp_min_opt==1)? false: true;
  orthoCfg->isGenWave     = (sim->gen_wave==1)? true: false;
  orthoCfg->numStates     = config.nstates;
  orthoCfg->grainSize     = config.orthoGrainSize;
  orthoCfg->instanceIndex = thisInstance.getPO();
  orthoCfg->maxTolerance  = sim->tol_norb;
  orthoCfg->uponToleranceFailure = CkCallback(CkIndex_GSpaceDriver::needUpdatedPsiV(), UgSpaceDriverProxy[thisInstance.getPO()]);

  // Fill in the paircalc configs that are instance dependent
  cfgSymmPC->gSpaceAID            = UgSpacePlaneProxy[thisInstance.getPO()].ckGetArrayID();
  cfgAsymmPC->gSpaceAID           = UgSpacePlaneProxy[thisInstance.getPO()].ckGetArrayID();
  cfgSymmPC->instanceIndex        = thisInstance.getPO();
  cfgAsymmPC->instanceIndex       = thisInstance.getPO();
  // Init the post-init callbacks that the paircalcs will trigger (after ortho<-->PC comm setup)
  cfgSymmPC->uponSetupCompletion  = CkCallback(CkIndex_InstanceController::doneInit(),CkArrayIndex1D(thisInstance.getPO()),instControllerProxy.ckGetArrayID());
  cfgAsymmPC->uponSetupCompletion = CkCallback(CkIndex_InstanceController::doneInit(),CkArrayIndex1D(thisInstance.getPO()),instControllerProxy.ckGetArrayID());
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
  cp::startup::PCCreationManager pcCreator(*cfgSymmPC, *cfgAsymmPC, *orthoCfg);
  pcCreator.build(pcHandleCB, pcMapCfg, sim->cp_lsda);

}

//============================================================================
#include "CPcharmParaInfo.def.h"
#include "PlatformSpecific.def.h"
#include "timeKeeper.def.h"
#include "startupMessages.def.h"
#include "cp_state_ctrl/CP_State_GSpacePlane.h"
#include "cpaimd.def.h"

//============================================================================
