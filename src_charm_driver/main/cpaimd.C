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

#include "MeshStreamingStrategy.h"
#include "MultiRingMulticast.h"
#include "OneTimeMulticastStrategy.h"
#include "TopoManager.h"
#include "TimeKeeper.h"

#include "charm++.h"
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


#include "mapvariables.h"


//============================================================================
// readonly globals

double Timer;
int nstates; 
int sizeX;
int nchareG;
int Ortho_UE_step2;
int Ortho_UE_step3;
int Ortho_UE_error;
bool Ortho_use_local_cb;
int done_init=0;
int planes_per_pe;

//============================================================================

#include "commlibhandles.h"

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
  topoMgr = NULL;
  //============================================================================
  /**
     # Sequential startup within Main */
  /* Check arguments : Verbose output about startup procedures */

  done_init=0;
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
  CkPrintf("  Reading Physics input from %s\n",msg->argv[2]);
  CkPrintf("  Reading Driver  input from %s\n",msg->argv[1]);

  PRINT_LINE_DASH; CkPrintf("\n");

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

  CkCallback piny_callback (CkCallback::ignore);
  Interface_ctrl piny_interface (msg->argv[2],piny_callback);

  CPcharmParaInfo *sim  = CPcharmParaInfo::get();
  PhysicsParamTransfer::ParaInfoInit(sim);

  int ibinary_opt    = sim->ibinary_opt;
  int natm_nl        = sim->natm_nl;
  int ees_eext_opt   = sim->ees_eext_on;
  int nchareRhoRHart = sim->ngrid_eext_c;
  int fftopt         = sim->fftopt;
  int natm_typ       = sim->natm_typ;
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

  PRINT_LINE_STAR;
  CkPrintf("Cpaimd-Charm-Driver input started \n");
  PRINT_LINE_DASH; CkPrintf("\n");
  Timer=CmiWallTimer();
  double phase1start=Timer;
  //numPes = 2048; 
  numPes=CkNumPes();
  int minimization_steps;
  if (sim->cp_bomd_opt) {
    minimization_steps = sim->btime;
  } else {
    minimization_steps = sim->ntime;
  }
  config.readConfig(msg->argv[1],sim->nstates,sim->sizeX,sim->sizeY,sim->sizeZ,
		    minimization_steps,ibinary_opt,natm_nl,fftopt,numPes,natm_typ,
		    ees_eext_opt,sim->gen_wave,sim->ncoef, sim->cp_min_opt, sim->ngrid_eext_c,
		    sim->doublepack,sim->pi_beads,sim->nkpoint,sim->ntemper,sim->nspin);

  fakeTorus        = config.fakeTorus>0;

  if(fakeTorus)
    {
      numPes=config.torusDimNX * config.torusDimNY * config.torusDimNZ * config.torusDimNT;
      CkPrintf("numpes set to %d by faketorus\n",numPes);
    }
  else if (CkNumPes() != config.numPes)
    {
      numPes=config.numPes;
      CkPrintf("numpes set to %d by config file\n",numPes);
    }
  CkPrintf("for numInstances %d numPes %d numPesPerInstance is %d \n",config.numInstances, config.numPes, config.numPesPerInstance);
  mapOffsets=new inttriple[config.numInstances];
  int numSfGrps    = config.numSfGrps;  // local copies are nice
  int doublePack   = config.doublePack;
  nchareG          = config.nchareG;
  sim->nchareG     = nchareG; 

  nstates          = config.nstates;    // globals : avail on all procs
  sizeX            = sim->sizeX;

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
  create_line_decomp_descriptor(sim);

  PhysicsParamTransfer::control_new_mapping_function(sim,doublePack);

  make_rho_runs(sim);

#include "initializeUber.C"

  pScratchProxy = CProxy_PhysScratchCache::ckNew();


  mainProxy=thishandle;
  // make one controller chare per instance
  instControllerProxy= CProxy_InstanceController::ckNew(config.numInstances);
  instControllerProxy.doneInserting();

  // make one controller temper
  temperControllerProxy= CProxy_TemperController::ckNew(1);
  temperControllerProxy.doneInserting();

  // make one collector per uberKmax
  CkArrayOptions enlopts(config.UberKmax);
  ENLEKECollectorProxy= CProxy_ENL_EKE_Collector::ckNew(config.UberImax*config.UberJmax*config.UberMmax, config.UberKmax, enlopts);
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

  CkPrintf("Initializing TopoManager\n");
  if(config.fakeTorus) {
    CkPrintf("Initializing TopoManager with fakeTorus\n");
    topoMgr = new TopoManager(config.torusDimNX, config.torusDimNY, 
			      config.torusDimNZ, config.torusDimNT);
  }
  else {
    topoMgr = new TopoManager();
  }
  CkPrintf("         Torus %d x %d x %d nodes %d x %d x %d VN %d DimNT %d .........\n", 
	   topoMgr->getDimX(), topoMgr->getDimY(), topoMgr->getDimZ(),
	   topoMgr->getDimNX(), topoMgr->getDimNY(), topoMgr->getDimNZ(),
	   topoMgr->hasMultipleProcsPerNode(), topoMgr->getDimNT());
  if(config.torusMap==1) {
    PRINT_LINE_STAR; CkPrintf("\n");
    CkPrintf("         Topology Sensitive Mapping being done for RSMap, GSMap, ....\n");
    CkPrintf("            ......., PairCalc, RhoR, RhoG and RhoGHart .........\n\n");
    PRINT_LINE_STAR; CkPrintf("\n");
  }

  int l=config.Gstates_per_pe;
  int m, pl, pm;
  pl = nstates / l;
  pm = config.numPesPerInstance / pl;
  if(pm == 0) {
    CkPrintf("Choose a larger Gstates_per_pe than %d such that { (no. of processors [%d] / no. of Instances [%d]) / (no. of states [%d] / Gstates_per_pe [%d]) } is > 0 \n", 
	     l, config.numPes, config.numInstances, nstates, l);
    CkAssert( pm > 0);
  }
  m = config.nchareG / pm;

  planes_per_pe=m;
  if(planes_per_pe <= 0) {
    CkPrintf("Choose a smaller Gstates_per_pe than %d such that { (no. of processors [%d] / no. of Instances [%d]) / (no. of states [%d] / Gstates_per_pe [%d]) } is > 0 and config.nchareG [%d] / pm [%d] {where pm = config.numPesPerInstance [%d]/ pl [%d] }  > 0\n", 
	     l, config.numPes, config.numInstances, nstates, l, config.nchareG, pm, config.numPesPerInstance, pl);
    CkAssert( m > 0);
  }

  // multiple instance mapping breaks if there isn't a topomanager
  CkPrintf("Initializing PeList\n");

  PeList *gfoo=NULL;
  PeList *rfoo=NULL;
  int x, y, z;
  PeListFactory *peList4PCmapping;
  if(!config.loadMapFiles && config.useCuboidMap)
    {
      if( config.numPesPerInstance % config.nchareG != 0)
	{
	  CkPrintf("To use CuboidMap nchareG %d should be chosen as a factor of numprocs %d / numInstances %d = %d\n",config.nchareG, config.numPes, config.numInstances, config.numPesPerInstance);
	  CkExit();
	}
      int procsPerPlane = config.numPesPerInstance / nchareG;
      int bx, by, bz;
      if(config.torusMap == 1) {
	boxSize = procsPerPlane;
	int order;

	// correction to accomodate multiple instances
	int dimNX, dimNY, dimNZ, dimNT;
	int longDim = -1, maxD;

	dimNX = topoMgr->getDimNX();
	dimNY = topoMgr->getDimNY();
	dimNZ = topoMgr->getDimNZ();
	dimNT = topoMgr->getDimNT();
	int  x1 = dimNX, y1 = dimNY, z1 = dimNZ;

	maxD = dimNX;
	longDim = 1;

	if (dimNY > maxD) { maxD = dimNY; longDim = 2; }
	if (dimNZ > maxD) { maxD = dimNZ; longDim = 3; }

	for(int i=0; i<config.numInstances; i++) {
	  switch(longDim) {
	  case 1:
	    x = i*(maxD/config.numInstances); y = 0; z = 0;
	    x1 = dimNX / config.numInstances;
	    break;
	  case 2:
	    x = 0; y = i*(maxD/config.numInstances); z = 0;
	    y1 = dimNY / config.numInstances;
	    break;
	  case 3:
	    x = 0; y = 0; z = i*(maxD/config.numInstances);
	    z1 = dimNZ / config.numInstances;
	    break;
	  }
	  // printf("Origin %d %d %d Size %d %d %d\n", x, y, z, x1, y1, z1);
	}

	CkPrintf("         Box per Instance: nodes %d x %d x %d .........\n", x1, y1, z1);
	if(findCuboid(bx, by, bz, order, x1, y1, z1, topoMgr->getDimNT(), boxSize, topoMgr->hasMultipleProcsPerNode()))
	  {
	    CkPrintf("Using %d, %d, %d dimensions for box %d mapping order %d\n", bx, by, bz, boxSize, order);
	    gfoo = new PeList(bx, by, bz, order, x1, y1, z1, dimNT);	// heap it
	    peList4PCmapping = new PeListFactory(bx,by,bz,order,x1,y1,z1,dimNT);
	  }
	else
	  {
	    peList4PCmapping = new PeListFactory(config.numPes);
	    CkPrintf("no box for %d\n", boxSize);
	    config.useCuboidMap = 0;
	    gfoo = new PeList(config.numPes);				// heap it
	  }
      }
      else
	{
	  peList4PCmapping = new PeListFactory(config.numPes);
	  // just split by numInstances
	  x=numInst*config.numPesPerInstance;
	  y=0;
	  z=0;
	  gfoo = new PeList(config.numPes);				// heap it
	}
    }
  else
    {
      gfoo = new PeList(config.numPes);				// heap it
      peList4PCmapping = new PeListFactory(config.numPes);
    }
  if(!config.loadMapFiles && config.useCuboidMapRS)
    rfoo = new PeList(*gfoo);
  else
    rfoo = new PeList(config.numPes);				// heap it

  computeMapOffsets();
  /* these really don't need to be different */
  availGlobG = rfoo;
  availGlobR = gfoo;
  newtime = CmiWallTimer();
  CkPrintf("Pelist initialized in %g with %d elements\n", newtime-Timer, availGlobG->count());
  // availGlobG->dump();
  Timer = newtime;
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
  CkPrintf("NumInstances %d: Beads %d  * Kpoints %d * Tempers %d * Spin %d\n",config.numInstances, config.UberImax, config.UberJmax, config.UberKmax,config.UberMmax);
  for(int integral=0; integral< config.UberImax; integral++)
    {
      for(int kpoint=0; kpoint< config.UberJmax; kpoint++)
	{
	  for(int temper=0; temper< config.UberKmax; temper++)
	    {
	      for(int spin=0; spin< config.UberMmax; spin++) {
		// for each new instance we need a new Uber Index
		CkVec  <int>  peUsedBySF;
		CkVec  <int>  peUsedByNLZ;
		CkVec  <int>  planeUsedByNLZ;

		UberIndex thisInstanceIndex(integral, kpoint, temper, spin); // Internal labels{x,y,z,s}
		thisInstance=UberCollection(thisInstanceIndex);
		UberAlles.push_back(thisInstance);// collection of proxies for all instances

		//============================================================================    
		// We will need a different one of these per instance
		// Transfer parameters from physics to driver
		//    read in atoms : create atoms group 
		control_physics_to_driver(thisInstance);

		//============================================================================ 

		if(config.UberImax>1)	      // handle Path Integrals
		  init_PIBeads(sim, thisInstance);

		// and then we make the usual set of chares to which we pass
		// the Uber Index.
		init_state_chares(natm_nl,natm_nl_grp_max,numSfGrps,doublePack,sim, thisInstance);

		//============================================================================
		// Create a paircalc/ortho bubble (symm and asymm pcs, ortho and related frills)

		orthostartup(&orthoCfg, &cfgSymmPC, &cfgAsymmPC, sim, peList4PCmapping);


		//============================================================================
		// compute the location for the non-local Z reduction roots for each plane
		// this can then be used in exclusion mapping to avoid overloading them
		int *usedProc= new int[config.numPesPerInstance];
		memset(usedProc, 0, sizeof(int)*config.numPesPerInstance);
		int charperpe = nstates/(config.numPesPerInstance);
		if(nstates % config.numPesPerInstance != 0)  charperpe++;
		if(charperpe<1) charperpe=1;
		for(int state=0; state<nstates; state++) {
		  int plane = nchareG-1;
		  while(plane >= 0) {
		    bool used = false;
		    int thisstateplaneproc = GSImaptable[thisInstance.getPO()].get(state,plane) % config.numPesPerInstance;
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
		// Create mapping classes for Paircalcular

		//============================================================================ 
		// Initialize the density chare arrays
		init_rho_chares(sim, thisInstance);

		//============================================================================ 
		// Initialize commlib strategies for later association and delegation
		if(sim->ees_nloc_on)
		  init_eesNL_chares( natm_nl, natm_nl_grp_max, doublePack, excludePes, sim, thisInstance);
		firstInstance=false;
		numInst++;
	      }
	    } 
	}
      // now safe to init atom bead commanders
      UatomsComputeProxy[thisInstance.getPO()].init();
    } // end of per instance init
  //============================================================================ 
  // Initialize commlib strategies for later association and delegation
  if(config.numInstances>1)
    CkPrintf("WARNING!!! Commlib does not work for multiple instances\n");

  init_commlib_strategies(sim->nchareRhoG, sim->sizeZ,nchareRhoRHart, thisInstance);

  TimeKeeperProxy.init();

  //============================================================================
  // clean up

  delete msg;
  //    delete sim;
  delete rfoo;
  delete gfoo;
  delete excludePes;
  //    delete availGlobR;
  //    delete availGlobG;

  //============================================================================

  newtime=CmiWallTimer();
  PRINT_LINE_DASH;
  CkPrintf("Cpaimd-Charm-Driver setup phase completed in %g \n",newtime-phase1start);
  PRINT_LINE_STAR; CkPrintf("\n");
  PRINT_LINE_STAR; 
  PRINT_LINE_DASH;CkPrintf("\n");
  CkPrintf("user mem %d\n",CmiMemoryUsage());
  Timer=newtime;
  /**@}*/
  //============================================================================
}// end Main
//============================================================================

//============================================================================    
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================    
/**
 * Cleanup stuff in the hopes of getting clean valgrind
 */
main::~main(){
  if (config.useCommlib) {        
    if(config.usePairEtoM){
    }//endif
  }//endif
}//end routine
//============================================================================    


//============================================================================ 
/**
 * Initialize Commlib communication strategies  
 */ 
//============================================================================ 
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void init_commlib_strategies(int numRhoG, int numReal, int numRhoRhart, UberCollection thisInstance){
  //============================================================================

  PRINT_LINE_STAR;
  PRINTF("Building Commlib strategies\n");
  PRINT_LINE_DASH;printf("\n");
  int i = 0;
  int j = 0;
  int rhoGhelpers = config.rhoGHelpers;
  int numRhoGHart = rhoGhelpers*numRhoG;
  int nchareHartAtmT = config.nchareHartAtmT;
  //============================================================================
#ifdef USE_COMLIB
  if (config.useCommlib) {        

    /*
    //StreamingStrategy *cmstrat = new StreamingStrategy(0.1,10);
    MeshStreamingStrategy *cmstrat = new MeshStreamingStrategy(1,5);
    gAsymInstance= ComlibRegister(cmstrat);    

    //StreamingStrategy *csymstrat = new StreamingStrategy(0.1,10);
    MeshStreamingStrategy *csymstrat = new MeshStreamingStrategy(1,5);
    gSymInstance= ComlibRegister(csymstrat);    
    */
    //================================================================
    CkArrayIndexMax *rhoGElements=NULL;
    CkArrayIndexMax *rhoRealElements = NULL;
    if(config.useRInsRhoGP)
      {
	CkPrintf("Making real_strategy with :");
	CkPrintf("src numReal %d dest numRhoG %d and numHartG %d numHartR %d\n",
		 numReal,numRhoG,numRhoGHart,numRhoRhart);
	//--------------------------------------------------------------
	//  For rho(r) to rho(g)
	rhoGElements = new CkArrayIndexMax[numRhoG];
	for (i = 0; i < numRhoG; i++) {
	  //CkArrayIndex2D idx2d(i,0);
	  CkArrayIndex1D idx1d(i);
	  rhoGElements[i] = idx1d;
	}//endfor

	rhoRealElements = new  CkArrayIndexMax[numReal];
	for(i = 0; i < numReal; i++) {
	  //CkArrayIndex2D idx2d(i,0);
	  CkArrayIndex1D idx1d(i);
	  rhoRealElements[i] = idx1d; 
	}//endfor
#ifdef OLD_COMMLIB
	CharmStrategy *real_strat;
#else
	Strategy *real_strat;
#endif      
	real_strat= new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoRealProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoGProxy[thisInstance.proxyOffset].ckGetArrayID(),
	   numReal, rhoRealElements, numRhoG, rhoGElements);
	commRealInstance= ComlibRegister(real_strat);
	delete [] rhoGElements;
	delete [] rhoRealElements;
      }
    //--------------------------------------------------------------
    //  For drho(r)/dx to igx*rho(g)
    if(config.useRInsIGXRhoGP)
      {
	rhoGElements = new CkArrayIndexMax[numRhoG];
	for (i = 0; i < numRhoG; i++) {
	  rhoGElements[i] = CkArrayIndex1D(i);
	}//endfor

	rhoRealElements = new  CkArrayIndexMax[numReal];
	for(i = 0; i < numReal; i++) {
	  rhoRealElements[i] = CkArrayIndex1D(i);
	}//endfor
#ifdef OLD_COMMLIB
	CharmStrategy *real_strat_igx;
#else
	Strategy *real_strat_igx;
#endif
	real_strat_igx= new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoRealProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoGProxy[thisInstance.proxyOffset].ckGetArrayID(),
	   numReal, rhoRealElements, numRhoG,rhoGElements);
	commRealIGXInstance= ComlibRegister(real_strat_igx);
	delete [] rhoGElements;
	delete [] rhoRealElements;
      }
    if(config.useRInsIGYRhoGP)
      {
	//--------------------------------------------------------------
	//  For drho(r)/dy to igy*rho(g)
	rhoGElements = new CkArrayIndexMax[numRhoG];
	for (i = 0; i < numRhoG; i++) {
	  rhoGElements[i] = CkArrayIndex1D(i);
	}//endfor

	rhoRealElements = new  CkArrayIndexMax[numReal];
	for(i = 0; i < numReal; i++) {
	  rhoRealElements[i] = CkArrayIndex1D(i);
	}//endfor
#ifdef OLD_COMMLIB
	CharmStrategy *real_strat_igy;
#else
	Strategy *real_strat_igy;
#endif
	real_strat_igy= new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoRealProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoGProxy[thisInstance.proxyOffset].ckGetArrayID(),
	   numReal, rhoRealElements, numRhoG,rhoGElements);
	commRealIGYInstance= ComlibRegister(real_strat_igy);

	delete [] rhoGElements;
	delete [] rhoRealElements;
      }
    //--------------------------------------------------------------
    //  For drho(r)/dz to igz*rho(g)
    if(config.useRInsIGZRhoGP)
      {
	rhoGElements = new CkArrayIndexMax[numRhoG];
	for (i = 0; i < numRhoG; i++) {
	  rhoGElements[i] = CkArrayIndex1D(i);
	}//endfor

	rhoRealElements = new  CkArrayIndexMax[numReal];
	for(i = 0; i < numReal; i++) {
	  rhoRealElements[i] = CkArrayIndex1D(i);
	}//endfor
#ifdef OLD_COMMLIB        
	CharmStrategy *real_strat_igz;
#else
	Strategy *real_strat_igz;
#endif
	real_strat_igz= new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoRealProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoGProxy[thisInstance.proxyOffset].ckGetArrayID(),
	   numReal, rhoRealElements, numRhoG,rhoGElements);
	commRealIGZInstance= ComlibRegister(real_strat_igz);
	delete [] rhoGElements;
	delete [] rhoRealElements;
      }
    //--------------------------------------------------------------
    // For hartree-Ext(g) to vks(r)
    if(config.useGHartInsRhoRP)
      {
	rhoGElements = new CkArrayIndexMax[numRhoGHart*nchareHartAtmT];
	for (j= 0; j < nchareHartAtmT; j++) {
	  for (i = 0; i < numRhoGHart; i++) {
	    rhoGElements[i+j*numRhoGHart] = CkArrayIndex1D(i+j*numRhoGHart);
	  }//endfor
	}
	rhoRealElements = new  CkArrayIndexMax[numReal];
	for(i = 0; i < numReal; i++) {
	  rhoRealElements[i] = CkArrayIndex1D(i);
	}//endfor
#ifdef OLD_COMMLIB
	CharmStrategy *gstrathart = new EachToManyMulticastStrategy
#else
	  Strategy *gstrathart = new EachToManyMulticastStrategy
#endif	

	  (USE_DIRECT, UrhoGHartExtProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoRealProxy[thisInstance.proxyOffset].ckGetArrayID(), 
	   numRhoGHart, rhoGElements, numReal, rhoRealElements);
	commGHartInstance = ComlibRegister(gstrathart);
	delete [] rhoGElements;
	delete [] rhoRealElements;
      }
    //--------------------------------------------------------------
    // vks(g), igxrho igyrho igzrho and white byrd to g-space
    if(config.useGIns0RhoRP)
      {

	rhoGElements = new CkArrayIndexMax[numRhoG];
	for (i = 0; i < numRhoG; i++) {
	  rhoGElements[i] = CkArrayIndex1D(i);
	}//endfor

	rhoRealElements = new  CkArrayIndexMax[numReal];
	for(i = 0; i < numReal; i++) {
	  rhoRealElements[i] = CkArrayIndex1D(i);
	}//endfor
#ifdef OLD_COMMLIB
	CharmStrategy *gstrat0;
#else
	Strategy *gstrat0;
#endif
	gstrat0 = new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoGProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoRealProxy[thisInstance.proxyOffset].ckGetArrayID(), 
	   numRhoG, rhoGElements, numReal, rhoRealElements);
	commGInstance0 = ComlibRegister(gstrat0);
	delete [] rhoGElements;
	delete [] rhoRealElements;
      }
    //--------------------------------------------------------------
    // rhog sends to rhor : div_x rho
    if(config.useGIns1RhoRP)
      {
	rhoGElements = new CkArrayIndexMax[numRhoG];
	for (i = 0; i < numRhoG; i++) {
	  rhoGElements[i] = CkArrayIndex1D(i);
	}

	rhoRealElements = new  CkArrayIndexMax[numReal];
	for(i = 0; i < numReal; i++) {
	  rhoRealElements[i] = CkArrayIndex1D(i);
	}
#ifdef OLD_COMMLIB
	CharmStrategy *gstrat1;
#else
	Strategy *gstrat1;
#endif

	gstrat1 = new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoGProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoRealProxy[thisInstance.proxyOffset].ckGetArrayID(), 
	   numRhoG, rhoGElements, numReal, rhoRealElements);
	commGInstance1 = ComlibRegister(gstrat1);
	delete [] rhoGElements;
	delete [] rhoRealElements;
      }
    //--------------------------------------------------------------
    // rhog sends to rhor : div_y rho
    if(config.useGIns2RhoRP)
      {
	rhoGElements = new CkArrayIndexMax[numRhoG];
	for (i = 0; i < numRhoG; i++) {
	  rhoGElements[i] = CkArrayIndex1D(i);
	}

	rhoRealElements = new  CkArrayIndexMax[numReal];
	for(i = 0; i < numReal; i++) {
	  rhoRealElements[i] = CkArrayIndex1D(i);
	}
#ifdef OLD_COMMLIB
	CharmStrategy *gstrat2;
#else
	Strategy *gstrat2;
#endif

	gstrat2 = new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoGProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoRealProxy[thisInstance.proxyOffset].ckGetArrayID(), 
	   numRhoG, rhoGElements, numReal, rhoRealElements);
	commGInstance2 = ComlibRegister(gstrat2);
	delete [] rhoGElements;
	delete [] rhoRealElements;
      }
    //--------------------------------------------------------------
    // rhog sends to rhor : div_z rho
    if(config.useGIns3RhoRP)
      {
	rhoGElements = new CkArrayIndexMax[numRhoG];
	for (i = 0; i < numRhoG; i++) {
	  rhoGElements[i] = CkArrayIndex1D(i);
	}

	rhoRealElements = new  CkArrayIndexMax[numReal];
	for(i = 0; i < numReal; i++) {
	  rhoRealElements[i] = CkArrayIndex1D(i);
	}
#ifdef OLD_COMMLIB
	CharmStrategy *gstrat3;
#else
	Strategy *gstrat3;
#endif

	gstrat3 = new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoGProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoRealProxy[thisInstance.proxyOffset].ckGetArrayID(), 
	   numRhoG, rhoGElements, numReal, rhoRealElements);
	commGInstance3 = ComlibRegister(gstrat3);
	delete [] rhoGElements;
	delete [] rhoRealElements;
      }
    //--------------------------------------------------------------
    // For white-byrd :  rhog sends to rhor : white byrd
    if(config.useGByrdInsRhoRBP)
      {
	rhoGElements = new CkArrayIndexMax[numRhoG];
	for (i = 0; i < numRhoG; i++) {
	  rhoGElements[i] = CkArrayIndex1D(i);
	}

	rhoRealElements = new  CkArrayIndexMax[numReal];
	for(i = 0; i < numReal; i++) {
	  rhoRealElements[i] = CkArrayIndex1D(i);
	}
#ifdef OLD_COMMLIB
	CharmStrategy *gstratByrd;
#else
	Strategy *gstratByrd;
#endif

	gstratByrd = new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoGProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoRealProxy[thisInstance.proxyOffset].ckGetArrayID(), 
	   numRhoG, rhoGElements, numReal, rhoRealElements);
	commGByrdInstance = ComlibRegister(gstratByrd);
	delete [] rhoGElements;
	delete [] rhoRealElements;
      }
    //--------------------------------------------------------------
    // For rhogHart send to rhorHart SF atmtyp
    if(config.useGHartInsRHart && false)
      {
	rhoGElements = new CkArrayIndexMax[numRhoGHart];
	for (i = 0; i < numRhoGHart; i++) {
	  rhoGElements[i] = CkArrayIndex1D(i);
	}

	rhoRealElements = new  CkArrayIndexMax[numRhoRhart];
	for(i = 0; i < numRhoRhart; i++) {
	  rhoRealElements[i] = CkArrayIndex3D(i,0,0);
	}
#ifdef OLD_COMMLIB
	CharmStrategy *gstratEext0;
#else
	Strategy *gstratEext0;
#endif

	gstratEext0 = new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoGHartExtProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoRHartExtProxy[thisInstance.proxyOffset].ckGetArrayID(), 
	   numRhoGHart, rhoGElements, numRhoRhart, rhoRealElements);

	commGHartRHartIns0 = ComlibRegister(gstratEext0);
      }
    //--------------------------------------------------------------
    // For rhogHart send to rhoRHart SF tot
    if(config.useGHartInsRHart && false)
      {
	rhoGElements = new CkArrayIndexMax[numRhoGHart];
	for (i = 0; i < numRhoGHart; i++) {
	  rhoGElements[i] = CkArrayIndex1D(i);
	}

	rhoRealElements = new  CkArrayIndexMax[numRhoRhart];
	for(i = 0; i < numRhoRhart; i++) {
	  rhoRealElements[i] = CkArrayIndex1D(i);
	}
#ifdef OLD_COMMLIB
	CharmStrategy *gstratEext1;
#else
	Strategy *gstratEext1;
#endif

	gstratEext1 = new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoGHartExtProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoRHartExtProxy[thisInstance.proxyOffset].ckGetArrayID(), 
	   numRhoGHart, rhoGElements, numRhoRhart, rhoRealElements);

	commGHartRHartIns1 = ComlibRegister(gstratEext1);
      }
    //--------------------------------------------------------------
    // For rhoRHart send to rhoGHart SF
    if(config.useRHartInsGHart && false){
      rhoGElements = new CkArrayIndexMax[numRhoGHart];
      for (i = 0; i < numRhoGHart; i++) {
        rhoGElements[i] = CkArrayIndex1D(i);
      }//endfor

      rhoRealElements = new  CkArrayIndexMax[numRhoRhart];
      for(i = 0; i < numRhoRhart; i++) {
        rhoRealElements[i] = CkArrayIndex1D(i);
      }//endfor
#ifdef OLD_COMMLIB
      CharmStrategy *real_strat_eext;
#else
      Strategy *real_strat_eext;
#endif

      real_strat_eext = new EachToManyMulticastStrategy
        (USE_DIRECT, UrhoRHartExtProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoGHartExtProxy[thisInstance.proxyOffset].ckGetArrayID(),
         numRhoRhart, rhoRealElements, numRhoGHart,rhoGElements);
      commRHartGHartIns  = ComlibRegister(real_strat_eext);
    }

  }//endif : use commlib
#endif


  //============================================================================
  // Real state space to gspace state and particle plane comm.

#ifdef USE_COMLIB
  if (config.useCommlibMulticast) {
#ifdef OLD_COMMLIB
    DirectMulticastStrategy *dstrat = new DirectMulticastStrategy
      (UrealSpacePlaneProxy[thisInstance.proxyOffset].ckGetArrayID(),1);

    RingMulticastStrategy *rstrat = new RingMulticastStrategy
      (UrealSpacePlaneProxy[thisInstance.proxyOffset].ckGetArrayID(),1);

    RingMulticastStrategy *r1strat = new RingMulticastStrategy
      (UparticlePlaneProxy[thisInstance.proxyOffset].ckGetArrayID(),1);

    MultiRingMulticast *mr1strat = new MultiRingMulticast
      (UparticlePlaneProxy[thisInstance.proxyOffset].ckGetArrayID(),1);

    DirectMulticastStrategy *ppdstrat = new DirectMulticastStrategy
      (UrealParticlePlaneProxy[thisInstance.proxyOffset].ckGetArrayID(),1);

    RingMulticastStrategy *pprstrat = new RingMulticastStrategy
      (UrealParticlePlaneProxy[thisInstance.proxyOffset].ckGetArrayID(),1);
    MultiRingMulticast *ppmr1strat = new MultiRingMulticast
      (UrealParticlePlaneProxy[thisInstance.proxyOffset].ckGetArrayID(),1);
#else
    DirectMulticastStrategy *dstrat = new DirectMulticastStrategy();

    OneTimeMulticastStrategy *rstrat = new OneTimeMulticastStrategy();

    RingMulticastStrategy *r1strat = new RingMulticastStrategy();

    RingMulticastStrategy *mr1strat = new RingMulticastStrategy();

    DirectMulticastStrategy *ppdstrat = new DirectMulticastStrategy();

    RingMulticastStrategy *pprstrat = new RingMulticastStrategy();
    RingMulticastStrategy *ppmr1strat = new RingMulticastStrategy();
#endif

    mcastInstance.reserve(config.UberJmax);
    //multiring should be good on large runs, but not on BG/L
    if(CkNumNodes()>64){
      for (int kp = 0; kp < config.UberJmax; kp++){
        mcastInstance.push_back(ComlibRegister(dstrat));
      }
      mcastInstancePP=ComlibRegister(mr1strat);
      mcastInstanceRPP=ComlibRegister(ppdstrat);
      mcastInstancemRPP=ComlibRegister(ppmr1strat);
    }else{
      for (int kp = 0; kp < config.UberJmax; kp++){
        mcastInstance.push_back(ComlibRegister(rstrat));
      }
      mcastInstancePP=ComlibRegister(r1strat);
      mcastInstanceRPP=ComlibRegister(pprstrat);
      mcastInstancemRPP=ComlibRegister(ppmr1strat);
    }//endif
  }
#endif

  //============================================================================
}//end routine
//============================================================================


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
      CkArrayOptions opts(natm);
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
		       int doublePack, CPcharmParaInfo *sim, UberCollection thisInstance)
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

  int nchareRhoG      = sim->nchareRhoG;
  int rhoGHelpers     = config.rhoGHelpers;
  int nchareGHart     = rhoGHelpers*nchareRhoG;
  int nchareRHart     = ngridcEext;

  int Rstates_per_pe  = config.Rstates_per_pe;
  int Gstates_per_pe  = config.Gstates_per_pe;
  int sGrainSize      = config.sGrainSize; 
  int numChunks       = config.numChunks;
  int nkpoint         = sim->nkpoint;

  //Need our maps and groups to exist before anyone tries to use them

  //============================================================================

  //============================================================================
  //need some groups to exist before we kick off the state which use them
  //--------------------------------------------------------------------------------
  // Groups : no placement required 

  UsfCacheProxy.push_back( CProxy_StructFactCache::ckNew(numSfGrps,natm_nl,natm_nl_grp_max, thisInstance));
  if(firstInstance) CkPrintf("created sfcache proxy\n");
  UsfCompProxy.push_back(CProxy_StructureFactor::ckNew());
  if(firstInstance) CkPrintf("created sfcomp proxy\n");

  if(thisInstance.idxU.y>0 || thisInstance.idxU.s >0)
    { // the set of chares being created is for a non-zero kpoint
      // all k-points and spins use the same atoms and energies
      // we simply direct the proxyoffset here to the one for
      // the 0th kpoint

      UberCollection zeroKpointInstance=thisInstance;
      zeroKpointInstance.idxU.y=0;
      zeroKpointInstance.idxU.s=0;
      int proxyOffset=zeroKpointInstance.setPO();
      // At some future point we could split the cache proxy so that
      // it could be shared across pretty much any kind of instance.
      // Currently we need different ones for beads and tempers due to
      // some atom bits that are blended in with the rest of the cache
      // stuff due to modularity failure in its design.  If eesCache
      // proxy eats all your memory, complain to Glenn.
      UeesCacheProxy.push_back(UeesCacheProxy[proxyOffset]);
    }
  else
    {
      UeesCacheProxy.push_back(CProxy_eesCache::ckNew(nchareRPP,nchareG,nchareRHart,nchareGHart,
						      nstates,nchareRhoG, nkpoint, thisInstance));
    }

  if(firstInstance) CkPrintf("created eescache proxy\n");

  int nchareRRhoTot  = nchareR*(config.rhoRsubplanes);
  int nchareRHartTot = nchareRHart*(config.rhoRsubplanes);

  int *numGState     = sim->nlines_per_chareG;
  int *numGNL        = sim->nlines_per_chareG;
  int *numGRho       = sim->nlines_per_chareRhoG;
  int *numGEext      = sim->nlines_per_chareRhoGEext;

  if(firstInstance)
    {
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
      int *numRXRho      = new int [nchareRRhoTot];
      int *numRYRho      = new int [nchareRRhoTot];
      int *numRXEext     = new int [nchareRHartTot];
      int *numRYEext     = new int [nchareRHartTot];
      int *numSubGx      = sim->numSubGx;
      create_Rho_fft_numbers(nchareR,nchareRHart,config.rhoRsubplanes,
			     sim->nplane_rho_x,sim->sizeY,ngridbEext,
			     numRXRho,numRYRho,numRXEext,numRYEext,numSubGx);
      UfftCacheProxy.push_back(CProxy_FFTcache::ckNew(
						      sim->sizeX,sim->sizeY,sim->sizeZ,
						      ngridaEext,ngridbEext,ngridcEext,ees_eext_on,
						      ngridaNl,  ngridbNl,  ngridcNl,  ees_nonlocal_on, 
						      sim->nlines_max, sim->nlines_max_rho,
						      config.nchareG,nchareR,
						      config.nchareG,nchareRPP, 
						      nchareRhoG,    nchareR,    nchareRRhoTot,
						      nchareGHart,   nchareRHart,nchareRHartTot,
						      numGState,     numRXState, numRYState,numRYStateLower,
						      numGNL,        numRXNL,    numRYNL, numRYNLLower,
						      numGRho,       numRXRho,   numRYRho,
						      numGEext,      numRXEext,  numRYEext,
						      config.fftopt,config.fftprogresssplitReal,config.fftprogresssplit,
						      config.rhoRsubplanes, thisInstance));
      CkPrintf("created fftcache proxy\n");
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
    }      
  else
    { // direct to 0th proxy
      UfftCacheProxy.push_back(UfftCacheProxy[0]);
    }

  //============================================================================
  // Instantiate the Chares with placement determined by the maps

  //---------------------------------------------------------------------------
  // state g-space

  PRINT_LINE_STAR;
  PRINTF("Building G-space (%d %d) and R-space (%d %d/%d) state Chares\n",
	 nstates, nchareG, nstates, nchareR, nchareRPP);
  PRINT_LINE_DASH;printf("\n");
  availGlobG->reset();
  double newtime=CmiWallTimer();
  int x=mapOffsets[numInst].getx();
  int y=mapOffsets[numInst].gety();
  int z=mapOffsets[numInst].getz();

  /** addtogroup mapping */
  /**@{*/
  if (firstInstance) {
    GSImaptable[numInst].buildMap(nstates, nchareG);
    int success = 0;
    if(config.loadMapFiles) {
      int size[2];
      size[0] = nstates; size[1] = nchareG;
      MapFile *mf = new MapFile("GSMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
      success = mf->loadMap("GSMap", &GSImaptable[0]);
#else
      success = mf->loadMap("GSMap", &GSmaptable);
#endif
      delete mf;
    }

    if(success == 0) {
      GSMapTable gsTable = GSMapTable(&GSImaptable[0], &GSImaptable[numInst], availGlobG, 
				      nchareG, nstates, Gstates_per_pe, config.useCuboidMap, numInst, x, y, z);
    }

    CkPrintf("GSMap created in %g\n", newtime-Timer);
  } else {
    if((CkNumPes()==1) && (config.fakeTorus!=1))
      GSImaptable[numInst]=GSImaptable[0];
    else
      GSImaptable[numInst].translate(&GSImaptable[0], x,y,z, config.torusMap==1);

  }

  // there is only one IntMap per chare type, but each instance has
  // its own map group
  CProxy_GSMap gsMap = CProxy_GSMap::ckNew(thisInstance);
  //  CkArrayOptions gSpaceOpts(nstates,nchareG);
  /**@}*/
  /** addtogroup GSpaceState */
  /**@{*/
  CkArrayOptions gSpaceOpts(nstates,nchareG);
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
  UgSpacePlaneProxy.push_back(CProxy_CP_State_GSpacePlane::ckNew(sizeX, 1, 1, sGrainSize, gforward, gbackward, thisInstance, gSpaceOpts));
  UgSpacePlaneProxy[thisInstance.proxyOffset].doneInserting();
  // CkPrintf("{%d} main uGSpacePlaneProxy[%d] is %d\n",thisInstance.proxyOffset,thisInstance.proxyOffset,CkGroupID(UgSpacePlaneProxy[thisInstance.proxyOffset].ckGetArrayID()).idx);
  /**@}*/
  //--------------------------------------------------------------------------------
  // Bind the GSpaceDriver array to the GSpacePlane array so that they migrate together
  CkArrayOptions gspDriverOpts(nstates,nchareG);
  gspDriverOpts.setAnytimeMigration(false);
  gspDriverOpts.setStaticInsertion(true);

  gspDriverOpts.bindTo(UgSpacePlaneProxy[thisInstance.proxyOffset]);
  UgSpaceDriverProxy.push_back( CProxy_GSpaceDriver::ckNew(thisInstance,gspDriverOpts) );
  UgSpaceDriverProxy[thisInstance.proxyOffset].doneInserting();
  /**@}*/
  //--------------------------------------------------------------------------------
  // We bind the particlePlane array to the gSpacePlane array migrate together
  /** addtogroup Particle */
  /**@{*/
  //  CkArrayOptions particleOpts(nstates,nchareG);
  CkArrayOptions particleOpts(nstates,nchareG);
  particleOpts.setAnytimeMigration(false);
  particleOpts.setStaticInsertion(true);
  particleOpts.setMap(gsMap); // the maps for both the arrays are the same
  particleOpts.bindTo(UgSpacePlaneProxy[thisInstance.proxyOffset]);
  UparticlePlaneProxy.push_back(CProxy_CP_State_ParticlePlane::ckNew(
								     nchareG, sim->sizeY, sim->sizeZ,ngridaNl,ngridbNl,ngridcNl,
								     1, numSfGrps, natm_nl, natm_nl_grp_max, nstates, 
								     nchareG, Gstates_per_pe, numIterNL, ees_nonlocal_on, 
								     thisInstance,  particleOpts));
  UparticlePlaneProxy[thisInstance.proxyOffset].doneInserting();

  if(config.dumpMapFiles) {
    int size[2];
    size[0] = nstates; size[1] = nchareG;
    MapFile *mf = new MapFile("GSMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMap(&GSImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMap(&GSmaptable);
#endif
    delete mf;
  }

  if(config.dumpMapCoordFiles) {
    // if someone wants to dump nontopo maps, they'll need a manager
    if(topoMgr == NULL)
      topoMgr = new TopoManager();

    int size[2];
    size[0] = nstates; size[1] = nchareG;
    MapFile *mf = new MapFile("GSMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMapCoords(&GSImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMapCoords(&GSmaptable);
#endif
    delete mf;
  }

  //---------------------------------------------------------------------------
  // state r-space

  // correction to accomodate multiple instances
  /**@}*/
  /** addtogroup mapping */
  /**@{*/
  if(firstInstance) {
    RSImaptable[numInst].buildMap(nstates, nchareR);
    Timer=CmiWallTimer();
    availGlobR->reset();

    int success = 0;
    if(config.loadMapFiles) {
      int size[2];
      size[0] = nstates; size[1] = nchareR;
      MapFile *mf = new MapFile("RSMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
      success = mf->loadMap("RSMap", &RSImaptable[0]);
#else
      success = mf->loadMap("RSMap", &RSmaptable);
#endif
      delete mf;
    }
    if(success == 0) {
#ifdef USE_INT_MAP
      RSMapTable RStable= RSMapTable(&RSImaptable[0], &RSImaptable[numInst], availGlobR, 
				     nstates, nchareR, Rstates_per_pe, config.useCuboidMapRS, 
				     &GSImaptable[0], config.nchareG, numInst, x, y, z);
#else
      RSMapTable RStable= RSMapTable(&RSmaptable, availGlobR, nstates, nchareR, 
				     Rstates_per_pe, config.useCuboidMapRS, &GSmaptable, config.nchareG);
#endif
    }
    newtime=CmiWallTimer();
    CkPrintf("RSMap created in %g\n", newtime-Timer);
    Timer=newtime;
  } else {
    if((CkNumPes()==1) && (config.fakeTorus!=1))
      RSImaptable[numInst]=RSImaptable[0];
    else
      RSImaptable[numInst].translate(&RSImaptable[0], x,y,z, config.torusMap==1);
    CkPrintf("RSMap instance %d created in %g\n", numInst, newtime-Timer);
  }

  CProxy_RSMap rsMap= CProxy_RSMap::ckNew(thisInstance);
  //  CkArrayOptions realSpaceOpts(nstates,nchareR);
  CkArrayOptions realSpaceOpts(nstates,nchareR);
  realSpaceOpts.setMap(rsMap);
  realSpaceOpts.setAnytimeMigration(false);
  realSpaceOpts.setStaticInsertion(true);

  /**@}*/
  int rforward=keeperRegister(std::string("RealSpaceForward"));
  int rbackward=keeperRegister(std::string("RealSpaceBackward"));
  /** \addtogroup RealSpaceState */
  /**@{*/
  UrealSpacePlaneProxy.push_back( CProxy_CP_State_RealSpacePlane::ckNew(1, 1, ngrida, ngridb, ngridc, rforward, rbackward, thisInstance, realSpaceOpts));
  UrealSpacePlaneProxy[thisInstance.proxyOffset].doneInserting();  

  if(config.dumpMapFiles) {
    int size[2];
    size[0] = nstates; size[1] = nchareR;
    MapFile *mf = new MapFile("RSMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMap(&RSImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMap(&RSmaptable);
#endif
    delete mf;
  }

  if(config.dumpMapCoordFiles) {
    int size[2];
    size[0] = nstates; size[1] = nchareR;
    MapFile *mf = new MapFile("RSMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMapCoords(&RSImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMapCoords(&RSmaptable);
#endif
    delete mf;

  }
  /**@}*/
  //--------------------------------------------------------------------------------
  // state r-particleplane

  availGlobR->reset();

  //--------------------------------------------------------------------------------
  // Do a fancy dance to determine placement of structure factors
  CkPrintf("Making SF non-EES\n");
  if(!ees_nonlocal_on)
    {
      CkPrintf("Making SF non-EES\n");
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

  //============================================================================
  // Set some com strategy of Sameer
  if(firstInstance) {

#ifdef USE_COMLIB
    if(config.useCommlib) {
      // TODO: do we need to do this for each Uber instance?
      // technically they are operating on mutually exclusive process
      // sets, so modulo insane commlib global state, it should be
      // safe to reuse these across ubers

      CkPrintf("Making State streaming strats useGssRP = %d, useMssGP = %d,  useGssRPP = %d,  useMssGPP = %d\n",(int)config.useGssInsRealP,(int)config.useMssInsGP,(int)config.useGssInsRealPP,(int)config.useMssInsGPP);
      //mstrat->enableShortArrayMessagePacking();
      //rspaceState to gspaceState : gspaceState to rspaceState 

      if(config.useGssInsRealP)
	{
	  StreamingStrategy *gmstrat = new StreamingStrategy(config.gStreamPeriod,
							     config.gBucketSize);
	  gssInstance= ComlibRegister(gmstrat);    
	}
      //mstrat->enableShortArrayMessagePacking();
      //rPPState to gPPState : gPPState to rPPState 
      if(config.useMssInsGP){
	StreamingStrategy *mstrat = new StreamingStrategy(config.rStreamPeriod,
							  config.rBucketSize);

	mssInstance= ComlibRegister(mstrat);    
      }
      if (config.useMssInsGPP){
	StreamingStrategy *rpmstrat = new StreamingStrategy(config.rStreamPeriod,
							    config.rBucketSize);

	mssPInstance= ComlibRegister(rpmstrat);    
      }
      if (config.useGssInsRealPP){
	StreamingStrategy *gpmstrat = new StreamingStrategy(config.gStreamPeriod,
							    config.gBucketSize);
	gssPInstance= ComlibRegister(gpmstrat);    
      }
    }//endif
#endif

  }
  //============================================================================

  printf("\n");

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
  if(config.useRhoExclusionMap)
    nlexcludePes=exclusion;
  else if(config.useReductionExclusionMap)
    nlexcludePes= new PeList(UpeUsedByNLZ[thisInstance.proxyOffset]);   
  else
    nlexcludePes= new PeList(0);
  if(config.excludePE0 && ! config.loadMapFiles)
    nlexcludePes->appendOne(0);
  int Rstates_per_pe  = config.Rstates_per_pe;
  availGlobG->reset();
  double newtime=CmiWallTimer();
#ifdef USE_INT_MAP
  RPPImaptable[thisInstance.getPO()].buildMap(nstates, nchareRPP);
#endif
  if(firstInstance)
    {
      int success = 0;
      if(config.loadMapFiles) {
	int size[2];
	size[0] = nstates; size[1] = nchareRPP;
	MapFile *mf = new MapFile("RPPMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
	success = mf->loadMap("RPPMap", &RPPImaptable[thisInstance.getPO()]);
#else
	success = mf->loadMap("RPPMap", &RPPmaptable);
#endif
	delete mf;
      }

      if(success == 0) {
#ifdef USE_INT_MAP
	RPPMapTable RPPtable= RPPMapTable(&RPPImaptable[thisInstance.getPO()], availGlobG, nlexcludePes, 
					  nstates,  nchareRPP, Rstates_per_pe,
					  boxSize, config.useCuboidMap, 
					  config.nchareG, &GSImaptable[thisInstance.getPO()]);
#else
	RPPMapTable RPPtable= RPPMapTable(&RPPmaptable, availGlobG, nlexcludePes, 
					  nstates,  nchareRPP, Rstates_per_pe,
					  boxSize, config.useCuboidMap, 
					  config.nchareG, &GSmaptable);
#endif
      }
    }
  else
    {
      int x=mapOffsets[numInst].getx();
      int y=mapOffsets[numInst].gety();
      int z=mapOffsets[numInst].getz();
      if((CkNumPes()==1) && (config.fakeTorus!=1))
	RPPImaptable[numInst]=RPPImaptable[0];
      else
	RPPImaptable[numInst].translate(&RPPImaptable[0], x,y,z, config.torusMap==1);
    }
  CProxy_RPPMap rspMap= CProxy_RPPMap::ckNew(thisInstance);
  newtime=CmiWallTimer();
  CkPrintf("RPPMap created in %g\n",newtime-Timer);

  if(config.dumpMapFiles) {
    int size[2];
    size[0] = nstates; size[1] = nchareRPP;
    MapFile *mf = new MapFile("RPPMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMap(&RPPImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMap(&RPPmaptable);
#endif
    delete mf;
  }
  if(config.dumpMapCoordFiles) {
    int size[2];
    size[0] = nstates; size[1] = nchareRPP;
    MapFile *mf = new MapFile("RPPMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMapCoords(&RPPImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMapCoords(&RPPmaptable);
#endif
    delete mf;
  }

  if(config.useRhoExclusionMap)
    {
      nlexcludePes=NULL;
    }
  else
    {
      if(nlexcludePes!=NULL)
	delete nlexcludePes;
    }
  if(config.excludePE0 &&  !config.loadMapFiles)
    {
      if(nlexcludePes==NULL)
	nlexcludePes=new PeList(0);
      nlexcludePes->appendOne(0);
    }
  Timer=newtime;
  CkArrayOptions pRealSpaceOpts(nstates,ngridcNl);
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
  int ngrid_eext_a    = sim->ngrid_eext_a;
  int ngrid_eext_b    = sim->ngrid_eext_b;
  int ngrid_eext_c    = sim->ngrid_eext_c;
  int ees_eext_on     = sim->ees_eext_on;
  int ees_nonlocal_on = sim->ees_nloc_on;
  int natmTyp         = sim->natm_typ;
  int nchareRhoG      = sim->nchareRhoG;
  int nchareRhoR      = sim->sizeZ;
  int rhoGHelpers     = config.rhoGHelpers;
  int nchareHartAtmT  = config.nchareHartAtmT;
  int nchareRhoGHart  = rhoGHelpers*nchareRhoG;
  int nchareRhoRHart  = ngrid_eext_c;


  if(thisInstance.idxU.y>0)
    { // the set of chares being created is for a non-zero kpoint
      // all k-points use the same rho, therefore we do not make new
      // chares, we simply direct the proxyoffset here to the one for
      // the 0th kpoint and return


      UberCollection zeroKpointInstance=thisInstance;
      zeroKpointInstance.idxU.y=0;
      int proxyOffset=zeroKpointInstance.setPO();
      UrhoRealProxy.push_back(UrhoRealProxy[proxyOffset]);      
      UrhoGProxy.push_back(UrhoGProxy[proxyOffset]);      
      UrhoGHartExtProxy.push_back(UrhoGHartExtProxy[proxyOffset]);
      if(ees_eext_on){
	UrhoRHartExtProxy.push_back(UrhoRHartExtProxy[proxyOffset]);
      }
      UlsRhoGProxy.push_back(UlsRhoGProxy[proxyOffset]);
      UlsRhoRealProxy.push_back(UlsRhoRealProxy[proxyOffset]);

      return 1;
    }




  //============================================================================
  // Output to the screen

  PRINT_LINE_STAR;
  CkPrintf("Building RhoR, RhoG, RhoGHartExt, RhoRHartExt Chares %d %d %d %d natmtyp:%d\n",
	   nchareRhoR,nchareRhoG,nchareRhoGHart,nchareRhoRHart,natmTyp);
  PRINT_LINE_DASH;printf("\n");

  //============================================================================
  // Nuke some procs from the list : reset, nuke, reset if you run out


  PeList *RhoAvail=NULL;
  if(!config.loadMapFiles)
    {
      availGlobR->reset();
      // make RhoAvail based on the pes used by RS
      RhoAvail=new PeList(*availGlobR);
    }
  //------------------------------------------------------------------------
  // subtract processors used by other nonscaling chares (i.e., non local reduceZ)
  excludePes= new PeList(0);   
  if(config.excludePE0 && !config.loadMapFiles)
    excludePes->appendOne(0);

  // avoid co-mapping with the particle plane non-local reduction roots
  if(config.useReductionExclusionMap && !config.loadMapFiles)
    {
      if( nchareRhoR*config.rhoRsubplanes+UpeUsedByNLZ[thisInstance.proxyOffset].size() <
	  RhoAvail->count()){

	CkPrintf("subtracting %d NLZ nodes from %d for RhoR Map\n",
		 UpeUsedByNLZ[thisInstance.proxyOffset].size(),RhoAvail->count());
	//       nlz.dump();
	*RhoAvail-*excludePes; //unary minus operator defined in PeList.h
	RhoAvail->reindex();
	CkPrintf("Leaving %d for RhoR Map\n",RhoAvail->count());
      }//endif

      //------------------------------------------------------------------------
      // subtract processors used by other nonscaling chares

      if(ees_nonlocal_on==0){
	if( nchareRhoR*config.rhoRsubplanes+UpeUsedBySF[thisInstance.proxyOffset].size()<RhoAvail->count()){
	  CkPrintf("subtracting %d SF nodes from %d for RhoR Map\n",
		   UpeUsedBySF[thisInstance.proxyOffset].size(),RhoAvail->count());
	  PeList sf(UpeUsedBySF[thisInstance.proxyOffset]);
	  *RhoAvail-sf;
	  RhoAvail->reindex();
	}//endif
      }//endif
    }
  if(!config.loadMapFiles && RhoAvail->count()>2 ) { RhoAvail->reindex(); }

  //============================================================================
  // Maps and options
  //CkPrintf("RhoR map for %d x %d=%d chares, using %d procs\n",nchareRhoR, config.rhoRsubplanes, nchareRhoR*config.rhoRsubplanes, RhoAvail->count());

  //---------------------------------------------------------------------------
  // rho RS 
  if(firstInstance)
    {
      RhoRSImaptable[thisInstance.getPO()].buildMap(nchareRhoR, config.rhoRsubplanes);
      int success = 0;
      if(config.loadMapFiles) {
	int size[2];
	size[0] = nchareRhoR; size[1] = config.rhoRsubplanes;
	MapFile *mf = new MapFile("RhoRSMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
	success = mf->loadMap("RhoRSMap", &RhoRSImaptable[thisInstance.getPO()]);
	delete mf;
      }
      if(success == 0) {
	RhoRSMapTable RhoRStable(&RhoRSImaptable[thisInstance.getPO()], RhoAvail, nchareRhoR, config.rhoRsubplanes, config.nstates, config.useCentroidMapRho, &RSImaptable[thisInstance.getPO()], excludePes);
      }
    }
  else
    {
      int x=mapOffsets[numInst].getx();
      int y=mapOffsets[numInst].gety();
      int z=mapOffsets[numInst].getz();
      if((CkNumPes()==1) && (config.fakeTorus!=1))
	RhoRSImaptable[numInst]=RhoRSImaptable[0];
      else
	RhoRSImaptable[numInst].translate(&RhoRSImaptable[0], x,y,z, config.torusMap==1);
    }

  CProxy_RhoRSMap rhorsMap = CProxy_RhoRSMap::ckNew(thisInstance);
  CkArrayOptions rhorsOpts(nchareRhoR, config.rhoRsubplanes);
  //CkArrayOptions rhorsOpts;
  rhorsOpts.setMap(rhorsMap);
  rhorsOpts.setAnytimeMigration(false);
  rhorsOpts.setStaticInsertion(true);
  if(config.dumpMapFiles) {
    int size[2];
    size[0] = nchareRhoR; size[1] = config.rhoRsubplanes;
    MapFile *mf = new MapFile("RhoRSMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMap(&RhoRSImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMap(&RhoRSmaptable);
#endif
    delete mf;
  }
  if(config.dumpMapCoordFiles) {
    int size[2];
    size[0] = nchareRhoR; size[1] = config.rhoRsubplanes;
    MapFile *mf = new MapFile("RhoRSMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMapCoords(&RhoRSImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMapCoords(&RhoRSmaptable);
#endif
    delete mf;
  }
  //---------------------------------------------------------------------------
  // rho GS 
  // if there aren't enough free procs refresh the RhoAvail list;
  if(!config.loadMapFiles && nchareRhoG>RhoAvail->count()) 
    {
      CkPrintf("refreshing avail list count %d less than rhog %d\n",RhoAvail->count(), nchareRhoG);
      RhoAvail->reset();
    }
  if(firstInstance)
    {
      RhoGSImaptable[thisInstance.getPO()].buildMap(nchareRhoG, 1);
      int success = 0;
      if(config.loadMapFiles) {
	int size[2];
	size[0] = nchareRhoG; size[1] = 1;
	MapFile *mf = new MapFile("RhoGSMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
	success = mf->loadMap("RhoGSMap", &RhoGSImaptable[thisInstance.getPO()]);
	delete mf;
      }
      if(success == 0) {
	RhoGSMapTable RhoGStable(&RhoGSImaptable[thisInstance.getPO()], RhoAvail,nchareRhoG, config.useCentroidMapRho, &RhoRSImaptable[thisInstance.getPO()], excludePes);

      }
    }
  else
    {
      int x=mapOffsets[numInst].getx();
      int y=mapOffsets[numInst].gety();
      int z=mapOffsets[numInst].getz();
      if((CkNumPes()==1) && (config.fakeTorus!=1))
	RhoGSImaptable[numInst]=RhoGSImaptable[0];
      else
	RhoGSImaptable[numInst].translate(&RhoGSImaptable[0], x,y,z, config.torusMap==1);
    }


  CProxy_RhoGSMap rhogsMap = CProxy_RhoGSMap::ckNew(thisInstance);
  CkArrayOptions rhogsOpts(nchareRhoG,1);
  //CkArrayOptions rhogsOpts;
  rhogsOpts.setMap(rhogsMap);
  rhogsOpts.setAnytimeMigration(false);
  rhogsOpts.setStaticInsertion(true);
  if(config.dumpMapFiles) {
    int size[2];
    size[0] = nchareRhoG; size[1] = 1;
    MapFile *mf = new MapFile("RhoGSMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMap(&RhoGSImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMap(&RhoGSmaptable);
#endif
    delete mf;
  }
  if(config.dumpMapCoordFiles) {
    int size[2];
    size[0] = nchareRhoG; size[1] = 1;
    MapFile *mf = new MapFile("RhoGSMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMapCoords(&RhoGSImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMapCoords(&RhoGSmaptable);
#endif
    delete mf;
  }


  //---------------------------------------------------------------------------
  // rho RHart 
  // if there aren't enough free procs refresh the avail list;
  if(!config.loadMapFiles && nchareRhoRHart*nchareHartAtmT > RhoAvail->count())
    RhoAvail->reset();
  CkArrayOptions rhorhartOpts(nchareRhoRHart, config.rhoRsubplanes, nchareHartAtmT);
  //CkArrayOptions rhorhartOpts;
  rhorhartOpts.setAnytimeMigration(false);
  rhorhartOpts.setStaticInsertion(true);    
  if(ees_eext_on) {
    if(firstInstance)
      {
	RhoRHartImaptable[thisInstance.getPO()].buildMap(nchareRhoRHart, config.rhoRsubplanes, nchareHartAtmT);
	int success = 0;
	if(config.loadMapFiles) {
	  int size[3];
	  size[0] = nchareRhoRHart; size[1] = config.rhoRsubplanes;
	  size[2] = nchareHartAtmT;
	  MapFile *mf = new MapFile("RhoRHartMap", 3, size, config.numPes, "TXYZ", 2, 1, 1, 1);
	  success = mf->loadMap("RhoRHartMap", &RhoRHartImaptable[thisInstance.getPO()]);
	  delete mf;
	}
	if(success == 0) {
	  RhoRHartMapTable RhoRHarttable(&RhoRHartImaptable[thisInstance.getPO()], RhoAvail, 
					 nchareRhoRHart, config.rhoRsubplanes, 
					 config.nchareHartAtmT, excludePes);
	}
      }
    else
      {
	int x=mapOffsets[numInst].getx();
	int y=mapOffsets[numInst].gety();
	int z=mapOffsets[numInst].getz();
	if((CkNumPes()==1) && (config.fakeTorus!=1))
	  RhoRHartImaptable[numInst]=RhoRHartImaptable[0];
	else
	  RhoRHartImaptable[numInst].translate(&RhoRHartImaptable[0], x,y,z, config.torusMap==1);
      }


    CProxy_RhoRHartMap rhorHartMap = CProxy_RhoRHartMap::ckNew(thisInstance);
    rhorhartOpts.setMap(rhorHartMap);

    if(config.dumpMapFiles) {
      int size[3];
      size[0] = nchareRhoRHart; size[1] = config.rhoRsubplanes,
				  size[2] = nchareHartAtmT;				  
      MapFile *mf = new MapFile("RhoRHartMap", 3, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
      mf->dumpMap(&RhoRHartImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
      mf->dumpMap(&RhoRHartmaptable);
#endif
      delete mf;
    }
    if(config.dumpMapCoordFiles) {
      int size[3];
      size[0] = nchareRhoRHart; size[1] = config.rhoRsubplanes,
				  size[2] = nchareHartAtmT;				  
      MapFile *mf = new MapFile("RhoRHartMap_coord", 3, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
      mf->dumpMapCoords(&RhoRHartImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
      mf->dumpMapCoords(&RhoRHartmaptable);
#endif
      delete mf;
    }
  } //endif : ees_ext_on
  CkPrintf("RhoRHartMap built %d x %d x %d\n",nchareRhoRHart, config.rhoRsubplanes, config.nchareHartAtmT);
  //---------------------------------------------------------------------------
  // rho GHart 
  // if there aren't enough free procs refresh the avail list;
  if(!config.loadMapFiles && nchareRhoGHart>RhoAvail->count())
    RhoAvail->reset();
  if(firstInstance)
    {
      RhoGHartImaptable[thisInstance.getPO()].buildMap(nchareRhoGHart, nchareHartAtmT);
      int success = 0;
      if(config.loadMapFiles) {
	int size[2];
	size[0] = nchareRhoGHart; size[1] = nchareHartAtmT;
	MapFile *mf = new MapFile("RhoGHartMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
	success = mf->loadMap("RhoGHartMap", &RhoGHartImaptable[thisInstance.getPO()]);
	delete mf;
      }
      if(success == 0) {

	MapType3 *RhoRHartImaptablep=NULL;
	if(ees_eext_on)
	  RhoRHartImaptablep=&RhoRHartImaptable[thisInstance.getPO()];
	RhoGHartMapTable RhoGHarttable(&RhoGHartImaptable[thisInstance.getPO()], RhoAvail, 
				       nchareRhoGHart, config.nchareHartAtmT,
				       config.useCentroidMapRho, 
				       RhoRHartImaptablep, excludePes);

      }
    }
  else
    {
      int x=mapOffsets[numInst].getx();
      int y=mapOffsets[numInst].gety();
      int z=mapOffsets[numInst].getz();
      if((CkNumPes()==1) && (config.fakeTorus!=1))
	RhoGHartImaptable[numInst]=RhoGHartImaptable[0];
      else
	RhoGHartImaptable[numInst].translate(&RhoGHartImaptable[0], x,y,z, config.torusMap==1);
    }


  CProxy_RhoGHartMap rhogHartMap = CProxy_RhoGHartMap::ckNew(thisInstance);
  CkArrayOptions rhoghartOpts(nchareRhoGHart, nchareHartAtmT);
  //  CkArrayOptions rhoghartOpts;
  rhoghartOpts.setMap(rhogHartMap);
  rhoghartOpts.setAnytimeMigration(false);
  rhoghartOpts.setStaticInsertion(true);
  if(config.dumpMapFiles) {
    int size[2];
    size[0] = nchareRhoGHart; size[1] =  nchareHartAtmT;
    MapFile *mf = new MapFile("RhoGHartMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMap(&RhoGHartImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMap(&RhoGHartmaptable);
#endif
    delete mf;
  }
  if(config.dumpMapCoordFiles) {
    int size[2];
    size[0] = nchareRhoGHart; size[1] =  nchareHartAtmT;
    MapFile *mf = new MapFile("RhoGHartMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMapCoords(&RhoGHartImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMapCoords(&RhoGHartmaptable);
#endif
    delete mf;
  }

  //============================================================================
  // Instantiate the chares

  bool dummy = true;

  //--------------------------------------------------------------------------

  // insert rhoreal
  int rhokeeper= keeperRegister(std::string("Density"));
  // rhorsopts contains the nchareRhoR, rhoRSubplanes, the maps, and will make all at once
  UrhoRealProxy.push_back(CProxy_CP_Rho_RealSpacePlane::ckNew(sizeX,dummy, 
							      ees_eext_on, ngrid_eext_c,
							      rhokeeper,
							      thisInstance,
							      rhorsOpts));
  /*    for (int i = 0; i < nchareRhoR; i++){
        for (int j = 0; j < config.rhoRsubplanes; j++){
        UrhoRealProxy[thisInstance.proxyOffset](i,j).insert(sizeX,dummy, ees_eext_on, ngrid_eext_c, rhokeeper);
        } //endfor
        } //endfor
  */
  UrhoRealProxy[thisInstance.proxyOffset].doneInserting();
  /// @todo: valgrind complains of a tiny memleak here. Check if callbacks get destroyed properly. 
  UrhoRealProxy[thisInstance.proxyOffset].ckSetReductionClient( new CkCallback(CkIndex_InstanceController::printEnergyEexc(NULL),CkArrayIndex1D(thisInstance.proxyOffset),instControllerProxy));
  //--------------------------------------------------------------------------
  // insert rhog
  UrhoGProxy.push_back(CProxy_CP_Rho_GSpacePlane::ckNew(sizeX, 1, 
							1, dummy, thisInstance,
							rhogsOpts));
  /*  for (int i = 0; i < nchareRhoG; i++){
      UrhoGProxy[thisInstance.proxyOffset](i,0).insert(sizeX, 1,1,dummy );
      }//endfor
  */
  UrhoGProxy[thisInstance.proxyOffset].doneInserting();
  //--------------------------------------------------------------------------
  // insert rhoghart
  UrhoGHartExtProxy.push_back(CProxy_CP_Rho_GHartExt::ckNew(ngrid_eext_a,ngrid_eext_b,
							    ngrid_eext_c,ees_eext_on,
							    natmTyp, thisInstance,
							    rhoghartOpts));
  /*
    for (int k = 0; k < nchareHartAtmT; k++){
    for (int i = 0; i < nchareRhoGHart; i++){
    UrhoGHartExtProxy[thisInstance.proxyOffset](i,k).insert(ngrid_eext_a,ngrid_eext_b,
    ngrid_eext_c,ees_eext_on,natmTyp);
    }//endfor
    }//endfor
  */
  /// @todo: valgrind complains of a tiny memleak here. Check if callbacks get destroyed properly. 
  UrhoGHartExtProxy[thisInstance.proxyOffset].ckSetReductionClient(new CkCallback(CkIndex_InstanceController::printEnergyHart(NULL),CkArrayIndex1D(thisInstance.proxyOffset),instControllerProxy));
  UrhoGHartExtProxy[thisInstance.proxyOffset].doneInserting();
  //--------------------------------------------------------------------------
  // insert rhoRhart
  if(ees_eext_on){
    UrhoRHartExtProxy.push_back(CProxy_CP_Rho_RHartExt::ckNew(ngrid_eext_a,ngrid_eext_b,
							      ngrid_eext_c,ees_eext_on,
							      natmTyp,thisInstance, 
							      rhorhartOpts));
    /*
      for (int k = 0; k < nchareHartAtmT; k++){
      for (int i = 0; i < nchareRhoRHart; i++){
      for (int j = 0; j < config.rhoRsubplanes; j++){
      UrhoRHartExtProxy[thisInstance.proxyOffset](i,j,k).insert(ngrid_eext_a,ngrid_eext_b,ngrid_eext_c,
      ees_eext_on,natmTyp);
      }//endfor
      }//endfor
      }//endfor
    */
    UrhoRHartExtProxy[thisInstance.proxyOffset].doneInserting();
  }//endif
  //===========================================================================
  // Output to the screen
  // need to add maps for these, for now just let em default
  // IF some condition which triggers QMMM
  int nchareRhoGLSP=1;
  int nchareRhoRealLSP=1;
  int nchareRhoRealLSPsubplanes=1;
  CkArrayOptions lspgspOpts(nchareRhoGLSP);
  CkArrayOptions lsprealOpts(nchareRhoRealLSP, nchareRhoRealLSPsubplanes);
  UlsRhoGProxy.push_back(CProxy_CP_LargeSP_RhoGSpacePlane::ckNew(thisInstance,lspgspOpts));
  UlsRhoRealProxy.push_back(CProxy_CP_LargeSP_RhoRealSpacePlane::ckNew(thisInstance,lsprealOpts));
  printf("\n");
  PRINT_LINE_DASH;
  PRINTF("Completed G-space/R-space Rho chare array build\n");
  PRINT_LINE_STAR;printf("\n");
  if(RhoAvail!=NULL)
    delete RhoAvail;

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
void control_physics_to_driver(UberCollection thisInstance){
  //============================================================================
  // make a group : create a proxy for the atom class and also a reduction client

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
      double kT         = PhysicsAtom->kT;

      Atom *atoms       = new Atom[natm];
      AtomNHC *atomsNHC = new AtomNHC[natm];

      PhysicsAtom->DriverAtomInit(natm,atoms,atomsNHC,ibead,itemper);
      UegroupProxy.push_back(CProxy_EnergyGroup::ckNew(thisInstance)); 
      // FIXME, this needs a real computation
      // also we need a real map. This was a naive, simple mapping scheme. We
      // can do better.
      int nChareAtoms=(config.numPesPerInstance<natm) ? config.numPesPerInstance : natm;

      if (firstInstance) {
	// build the base atom map
	AtomImaptable[numInst].buildMap(nChareAtoms);
	availGlobG->reset();
	AtomMapTable aTable = AtomMapTable(&AtomImaptable[numInst], availGlobG, 
					   numInst,nChareAtoms);


      } else {
	int x=mapOffsets[numInst].getx();
	int y=mapOffsets[numInst].gety();
	int z=mapOffsets[numInst].getz();
	if((CkNumPes()==1) && (config.fakeTorus!=1))
	  AtomImaptable[numInst]=AtomImaptable[0];
	else
	  AtomImaptable[numInst].translate(&AtomImaptable[0], x,y,z, config.torusMap==1);
      }
      CProxy_AtomComputeMap aMap = CProxy_AtomComputeMap::ckNew(thisInstance);
      CkArrayOptions atomOpts(nChareAtoms);
      atomOpts.setMap(aMap);
      atomOpts.setAnytimeMigration(false);
      atomOpts.setStaticInsertion(true);
      UatomsCacheProxy.push_back( CProxy_AtomsCache::ckNew(natm,natm_nl,
							   atoms,thisInstance));
      UatomsComputeProxy.push_back( CProxy_AtomsCompute::ckNew(natm,natm_nl,
							       len_nhc,
							       iextended_on,
							       cp_min_opt,cp_wave_opt,cp_bomd_opt,isokin_opt,
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
  if(config.fakeTorus)
    proc=proc%CkNumPes();
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
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// return the cuboid x,y,z of a subpartition exactly matching that volume
//============================================================================
bool findCuboid(int &x, int &y, int &z, int &order, int maxX, int maxY, int maxZ, int maxT, int volume, int vn){
  //============================================================================
  int maxD=maxX;
  int minD=maxX;
  int minX=maxX;
  int minY=maxY;
  int minZ=maxZ;
  //  vn=0;
  if(config.forceMappingAxis>-1)
    {
      order=config.forceMappingAxis;
      if (config.forceMappingAxis==0) 
	minD=minX;
      if (config.forceMappingAxis==1) 
	minD=maxY;
      if (config.forceMappingAxis==2) 
	minD=maxZ;
    }
  else
    { 
      order=0;
      maxD = (maxY>maxD) ? maxY : maxD;
      maxD = (maxZ>maxD) ? maxZ : maxD;
      if(vn)
	{  // using Y as the prism axis seems to suck
	  maxD = (maxY>maxD) ? maxY : maxD;
	  maxD = (maxZ>maxD) ? maxZ : maxD;
	  minD = (maxY<minD) ? maxY : minD;
	  minD = (maxZ<maxD) ? maxZ : minD;

	}
      else
	{
	  maxD = (maxY>maxD) ? maxY : maxD;
	  maxD = (maxZ>maxD) ? maxZ : maxD;
	  minD = (maxY<minD) ? maxY : minD;
	  minD = (maxZ<minD) ? maxZ : minD;

	}
    }
  // CkPrintf("minD %d maxD %d\n",minD, maxD);
  if(config.useCuboidMapRS)
    {
      CkPrintf("Using long prisms for useCuboidMapRS\n");
    }

  // We were reducing the volume by maxT and then finding the dimensions of the
  // box in terms of the no. of nodes and not processors

  // That worked ok on BG/L, but BG/P maxT of 4 distorts the long axis logic.
  // in TXYZ you have a 32x8x16 for 1 Rack or 32x8x32 for 2 Rack
  // in XYZT you have 8x8x64 for 1 Rack or 8x8x128 for 2 Rack

  // The latter case makes clear the bizarre and unreal distortion created by
  // pretending that maxT is just a multiplier along one dimension when
  // considering actual box size.  What was merely odd on BG/L is now
  // seriously peculiar.  Also default TXYZ mappings and default Order 0
  // mappings to use the X as long axis for prisms fail to work right because
  // we aren't really spanning the X*maxT dimension with these long prisms.

  int redVol = volume / maxT;
  double cubert= cbrt((double) redVol);
  int cubetrunc= (int) cubert;
  x=y=z=cubetrunc;
  if(cubetrunc>minD)
    cubetrunc=minD;
  if(cubetrunc>maxY)
    cubetrunc=maxY;
  if(cubetrunc>maxZ)
    cubetrunc=maxZ;
  if(redVol==x*y*z && !config.useCuboidMapRS)
    return true;
  bool switchSet=false;
  CkAssert(redVol>0);
  switch (redVol) // for the common values we just pick cuboids we like
    {
    case 1: 
      x=1; y=1; z=1; switchSet=true; break;
    case 2: 
      x=2; y=1; z=1; switchSet=true; break;
    case 3:
      x=3; y=1; z=1; switchSet=true; break;
    case 4:
      x=2; y=2; z=1; switchSet=true; break;
    case 5:
      x=5; y=1; z=1; switchSet=true; break;
    case 6:
      x=3; y=2; z=1; switchSet=true; break;
    case 7:
      x=7; y=1; z=1; switchSet=true; break;
    case 8:
      if(config.useCuboidMapRS)
	{
	  if(minD>=2)

	    { x=2; y=2; z=2; switchSet=true; break;}
	  /* no evidence that these other schemes help
	     if(minD==4)
	     { x=4; y=2; z=1; switchSet=true; break;}
	     if(minD>=8)
	     { x=8; y=1; z=1; switchSet=true; break;}
	  */
	}
    case 9:
      x=3; y=3; z=1; switchSet=true; break;
    case 10:
      x=5; y=2; z=1; switchSet=true; break;
    case 12:
      x=2; y=3; z=2; switchSet=true; break;
    case 14:
      x=7; y=2; z=1; switchSet=true; break;
    case 15:
      x=5; y=3; z=1; switchSet=true; break;
    case 16:
      if(config.useCuboidMapRS)
	{
	  if(minD>=8)
	    { x=8; y=2; z=1; switchSet=true; break;}
	}
      x=4; y=2; z=2; switchSet=true; break;
    case 18:
      x=3; y=3; z=2; switchSet=true; break;
    case 20:
      x=5; y=2; z=2; switchSet=true; break;
    case 21:
      x=7; y=3; z=1; switchSet=true; break;
    case 24:
      x=4; y=3; z=2; switchSet=true; break;
    case 25:
      x=5; y=5; z=1; switchSet=true; break;
    case 27:
      x=3; y=3; z=3; switchSet=true; break;
    case 28:
      x=7; y=2; z=2; switchSet=true; break;
    case 30:
      x=5; y=2; z=2; switchSet=true; break;
    case 32:
      if(config.useCuboidMapRS)
	{
	  if(minD==8)
	    { x=8; y=2; z=2; switchSet=true; break;}
	  if(minD==16)
	    { x=16; y=2; z=1; switchSet=true; break;}
	  if(minD>=32)
	    { x=32; y=1; z=1; switchSet=true; break;}

	}
      x=4; y=2; z=4; switchSet=true; break;
    case 35:
      x=7; y=5; z=1; switchSet=true; break;
    case 36:
      x=4; y=3; z=3; switchSet=true; break;
    case 40:
      x=5; y=4; z=2; switchSet=true; break;
    case 42:
      x=7; y=3; z=2; switchSet=true; break;
    case 43:
      x=7; y=3; z=2; switchSet=true; break;
    case 45:
      x=5; y=3; y=3; switchSet=true; break;
    case 48:
      x=4; y=3; z=4; switchSet=true; break;
    case 50:
      x=5; y=5; z=2; switchSet=true; break;
    case 54:
      x=6; y=3; z=3; switchSet=true; break;
    case 56:
      x=7; y=4; z=2; switchSet=true; break;
    case 60:
      x=5; y= 4; z=3; switchSet=true; break;
    case 64:
      if(config.useCuboidMapRS)
	{
	  if(minD==8)
	    { x=8; y=4; z=2; switchSet=true; break;}
	  if(minD==16)
	    { x=16; y=2; z=2; switchSet=true; break;}
	  if(minD==32)
	    { x=32; y=2; z=1; switchSet=true; break;}
	  if(minD>=64)
	    { x=64; y=1; z=1; switchSet=true; break;}
	}
      x=4; y=4; z=4; switchSet=true; break;
    case 128:
      if(config.useCuboidMapRS)
	{
	  if(minD==8)
	    { x=8; y=4; z=4; switchSet=true; break;}
	  if(minD==16)
	    {x=16; y=4; z=2; switchSet=true; break;}
	  if(minD==32)
	    {  x=32; y=2; z=2; switchSet=true; break;	}
	  if(minD==64)
	    {  x=64; y=2; z=1; switchSet=true; break;	}
	  if(minD>=128)
	    {  x=128; y=1; z=1; switchSet=true; break;	}
	}
      x=8; y=4; z=4; switchSet=true; break;
    case 256:
      if(config.useCuboidMapRS)
	{
	  if(minD==8)
	    { x=8; y=8; z=4; switchSet=true; break;}
	  if(minD==16)
	    { x=16; y=4; z=4; switchSet=true; break;}
	  if(minD==32)
	    { x=32; y=4; z=2; switchSet=true; break;}
	  if(minD==64)
	    { x=64; y=2; z=2; switchSet=true; break;}
	  if(minD>=128)
	    { x=64; y=2; z=1; switchSet=true; break;}
	}
      x=8; y=8; z=4; switchSet=true; break;
    case 512:
      if(config.useCuboidMapRS)
	{
	  if(minD==8)
	    { x=8; y=8; z=8; switchSet=true; break;}
	  if(minD==16)
	    { x=16; y=4; z=8; switchSet=true; break;}
	  if(minD==32)
	    { x=32; y=4; z=4; switchSet=true; break;}
	  if(minD==64)
	    { x=64; y=4; z=2; switchSet=true; break;}
	  if(minD==128)
	    { x=128; y=2; z=2; switchSet=true; break;}
	  if(minD>=256)
	    { x=256; y=2; z=1; switchSet=true; break;}
	}
      x=8; y=8; z=8; switchSet=true; break;

    default:
      break;
    }
  if(switchSet && config.forceMappingAxis>-1)
    {

      switch(config.forceMappingAxis)
	{
	case 0:
	  return true; 	  //no change
	case 1:
	  {
	    order=1; //YXZ
	    int swap=x;
	    x=y;
	    y=swap;
	    return true;
	  }
	case 2:
	  {
	    order=2; //ZXY
	    int swap=x;
	    x=z;
	    z=swap;
	    return true;
	  }
	}
    }
  else if (switchSet)
    {
      // now correct the x,y,z to put long prism axis along the
      // smallest torus dimension which will fit.
      if(x==maxX)
	return true;
      if(x==maxY)
	{ // change to Y
	  order=1; //YXZ
	  int swap=x;
	  x=y;
	  y=swap;
	  return true;
	}
      if(x==maxZ)
	{ // change to Z
	  order=2; //ZXY
	  int swap=x;
	  x=z;
	  z=swap;
	  return true;
	}
      // if we're here then we don't have a spanning prism
      // just pick the smallest which will fit.
      if(x<maxX)
	return true;
      if(x<maxY)
	{ // change to Y
	  order=1; //YXZ
	  int swap=x;
	  x=y;
	  y=swap;
	  return true;
	}
      if(x<maxZ)
	{ // change to Z
	  order=2; //ZXY
	  int swap=x;
	  x=z;
	  z=swap;
	  return true;
	}
    }
  // its something weird so try a best fit
  int start=cubetrunc-1;
  if(config.useCuboidMapRS)
    x = (redVol>=maxX) ? maxX : cubetrunc;
  else
    x=cubetrunc;
  for(; x<=maxX;x++)
    {
      for(y=start; y<=maxY;y++)
	{
	  for(z=start; z<=maxZ;z++)
	    {
	      if(redVol==x*y*z)
		return true;
	    }
	}
    }
  return false;


}
//============================================================================
/**@}*/
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
  for(int thisInst=1; thisInst<config.numInstances; thisInst++)
    {
      if(config.torusMap == 1) {
	int dimNX, dimNY, dimNZ;
	int longDim = -1, maxD;

	x1 = dimNX = topoMgr->getDimNX();
	y1 = dimNY = topoMgr->getDimNY();
	z1 = dimNZ = topoMgr->getDimNZ();

	maxD = dimNX;
	longDim = 1;

	if (dimNY > maxD) { maxD = dimNY; longDim = 2; }
	if (dimNZ > maxD) { maxD = dimNZ; longDim = 3; }

	switch(longDim) {
        case 1:
          x = thisInst*(maxD/config.numInstances); y = 0; z = 0;
          x1 = dimNX / config.numInstances;
          break;
        case 2:
          x = 0; y = thisInst*(maxD/config.numInstances); z = 0;
          y1 = dimNY / config.numInstances;
          break;
        case 3:
          x = 0; y = 0; z = thisInst*(maxD/config.numInstances);
          z1 = dimNZ / config.numInstances;
          break;
	}
      }
      else
	{
	  // just split by numInstances
	  x=thisInst*config.numPesPerInstance;
	  y=0;
	  z=0;
	}
      mapOffsets[thisInst]=inttriple(x,y,z);
    }
}


void paircalcstartup(pc::pcConfig *cfgSymmPC, pc::pcConfig *cfgAsymmPC, CPcharmParaInfo *sim, int doublePack)
{



  // Stuff it with the actual configurations
  cfgSymmPC->isDynamics         = (sim->cp_min_opt==1)? false: true;
  cfgSymmPC->useComplexMath     = false;

  cfgSymmPC->numPlanes          = config.nchareG;
  cfgSymmPC->numStates          = nstates;
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
  cfgSymmPC->arePhantomsOn      = config.phantomSym;
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
  cfgSymmPC->uponSetupCompletion  = CkCallback(CkIndex_InstanceController::doneInit(NULL),CkArrayIndex1D(thisInstance.getPO()),instControllerProxy.ckGetArrayID());
  cfgAsymmPC->uponSetupCompletion = CkCallback(CkIndex_InstanceController::doneInit(NULL),CkArrayIndex1D(thisInstance.getPO()),instControllerProxy.ckGetArrayID());
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
  pcCreator.build(pcHandleCB, pcMapCfg);

}

//============================================================================
#include "CPcharmParaInfo.def.h"
#include "timeKeeper.def.h"
#include "startupMessages.def.h"
#include "cp_state_ctrl/CP_State_GSpacePlane.h"
#include "cpaimd.def.h"

//============================================================================
