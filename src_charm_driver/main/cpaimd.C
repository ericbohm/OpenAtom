//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file cpaimd.C
 *                         cpaimd-charm-driver
 *    Software developed by the Parallel Programing Laboratory, UIUC.
 *    in collaboration with IBM and NYU.
 *   
 *    This file contains cpaimd-charm-driver main. It creates and 
 *    initializes all the arrays and libraries. 
 */      
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \mainpage
<TABLE width="100%" border="0" cellspacing="0" cellpadding="5">
<TR bgcolor="#1C097D"><TD><font color="#FFFFFF" size="4"><table width=100%><tr><td><font size=6 color="#FFFFFF">Ab Initio Molecular Dynamics
</font></td><td align=right><img src="../images/image.jpg"></td></tr></table></FONT></TD></TR>
<TR bgcolor="#FFFFFF"><TD><br>Many important problems in material science, chemistry, solid-state physics, and biophysics require a modeling
 approach based on fundamental quantum mechanical principles. A particular approach that has proved to be relatively
 efficient and useful is Car-Parrinello ab initio molecular dynamics (CPAIMD). Parallelization of this approach beyond a
 few hundred processors is challenging, due to the complex dependencies among various subcomputations, which lead to
 complex communication optimization and load balancing problems. We are parallelizing CPAIMD using Charm++. The computation is modeled using 
 a large number of virtual processors, which are mapped flexibly to
 available processors with assistance from the Charm++ runtime system.
 <br>
 <br>
 This project is a large NSF funded collaboration involving us (PPL : Laxmikant Kale)
 and Drs. Roberto Car, Michael Klein, Glenn Martyna, Mark Tuckerman, Nick Nystrom and Josep Torrellas.
 <br>
 
 <center>
 <img src="../images/Slide1.JPG">
 </center>

 
 
<br>&nbsp;</TD></TR>
<TR bgcolor="#1C097D"><TD><FONT color="#FFFFFF" size="4">Software</FONT></TD></TR>
<TR bgcolor="#FFFFFF"><TD>Currently we have a Charm++ implementation of the core of the CP method. You
 can check out the latest build using CVS. The code is available under the
 module name "leanCP". Charm++ and <a href="http://www.fftw.org">FFTW</a>
 are required to run the code.   
</TD></TR>
<TR bgcolor="#1C097D"><TD><FONT color="#FFFFFF" size="4">People</FONT></TD></TR>
<TR bgcolor="#FFFFFF"><TD><UL>
<LI><A href="mailto:ebohm AT uiuc.edu
">Eric Bohm</A></LI>
<LI><A href="mailto:yanshi AT uiuc.edu
">Yan Shi</A></LI>

<LI><A href="mailto:kale AT cs.uiuc.edu
">L. V. Kale</A></LI>
<LI><A href="mailto:kunzman2 AT uiuc.edu
">David Kunzman</A></LI>
<LI><A href="mailto:skumar2 AT uiuc.edu
">Sameer Kumar</A></LI>

</UL>
</TD></TR>
<TR bgcolor="#1C097D"><TD><FONT color="#FFFFFF" size="4">Papers</FONT></TD></TR>
<TR bgcolor="#FFFFFF"><TD><UL>
<LI><A href="/papers/CpaimdTR03.shtml">03-06
</A>&nbsp;&nbsp;&nbsp;Ramkumar Vadali, L. V. Kale, Glenn Martyna, Mark Tuckerman,&nbsp;&nbsp;<B>Scalable Parallelization of Ab Initio Molecular Dynamics</B>,
&nbsp;&nbsp;<I>Technical Report, communicated to SC 2003</I></LI>

</UL>
</TD></TR>
<TR bgcolor="#1C097D"><TD><FONT color="#FFFFFF" size="4">Related Links</FONT></TD></TR>
<TR bgcolor="#FFFFFF"><TD><UL>
<LI><A href="http://homepages.nyu.edu/~mt33/PINY_MD/PINY.html">PINY MD</A></LI>
<LI><a href="http://charm.cs.uiuc.edu/presentations/cpaimd03/cpaimd.ppt">Presentation on latest state of the code</a></LI>

</UL>
</TD></TR>
<tr bgcolor=#1C097D><td><table width=100%><tr bgcolor=#1C097D>
<td><font color=white>This page maintained by  <A href="mailto:ebohm AT uiuc.edu
">Eric Bohm</A>.</font></td>

<td align=right><font color=white>Back to the <A HREF="/research" onclick="document.research.submit();return false;">PPL Research Page</A></font></td>
</tr></table></td></tr>
</TABLE>
</td></tr></table>
 *
 * end of mainpage
 */

 
#include <math.h>
#include "charm++.h"
#include "ckarray.h"
#include "util.h"

#include "groups.h"
#include "cpaimd.h"
#include "ortho.h"
#include "sim_subroutines.h"
#include "StructFactorCache.h"
#include "StructureFactor.h"
#include "CP_State_Plane.h"
#include "MeshStreamingStrategy.h"

#include "../include/CPcharmParaInfo.h"
#include "../../src_piny_physics_v1.0/include/class_defs/Interface_ctrl.h"
#include "../../src_piny_physics_v1.0/include/charm_defs/Interface_ctrl.decl.h"
#include "../../src_piny_physics_v1.0/include/class_defs/PINY_INIT/PhysicsParamTrans.h"
#include "../../src_piny_physics_v1.0/include/class_defs/PINY_INIT/PhysicsAtomPosInit.h"


//============================================================================
/** 
 * \defgroup  piny_vars Defining all Charm++ readonly variables for PINY physics 
 *
 */
//============================================================================
/* @{ */
extern MDINTEGRATE  readonly_mdintegrate;
extern MDATOMS      readonly_mdatoms;
extern MDINTER      readonly_mdinter;
extern MDINTRA      readonly_mdintra;
extern GENERAL_DATA readonly_general_data;
extern CP           readonly_cp; 
/* @} */

#include "MultiRingMulticast.h"

//============================================================================
/**
 * \defgroup proxy_comlib_vars
 * Defining all the Charm++ readonly variables, which include proxies
 * to access the arrays and groups and the Communication Library
 * handles.
 */
//============================================================================

Config config;
PairCalcID pairCalcID1;
PairCalcID pairCalcID2;
CProxy_CP_State_GSpacePlane gSpacePlaneProxy;
CProxy_CP_State_ParticlePlane particlePlaneProxy;
CProxy_CP_State_RealSpacePlane realSpacePlaneProxy;
CProxy_CP_Rho_RealSpacePlane rhoRealProxy;
CProxy_CP_Rho_GSpacePlane rhoGProxy;
CProxy_CP_Rho_GSpacePlaneHelper rhoGHelperProxy;
CProxy_Ortho orthoProxy;
CProxy_CPcharmParaInfoGrp scProxy;
CProxy_AtomsGrp atomsGrpProxy;
CProxy_EnergyGroup egroupProxy;
CProxy_FFTcache fftCacheProxy;
CProxy_StructFactCache sfCacheProxy;
CProxy_StructureFactor sfCompProxy;

int atom_integrate_done;  // not a real global : more like a group of size 1

int nstates;  // readonly globals
int sizeX;
int nchareG;

//============================================================================
// For using the multicast library :  Set some reduction clients

CkGroupID mCastGrpId; 
CkGroupID mCastGrpId2; 
// For using the communication library
ComlibInstanceHandle commInstance;
ComlibInstanceHandle commRealInstance;
ComlibInstanceHandle mcastInstance;
ComlibInstanceHandle ssInstance;
ComlibInstanceHandle mssInstance;
ComlibInstanceHandle mcastInstancePP;
CkReduction::reducerType complexVectorAdderType;

#include "ReductionClients.h"


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** dummy function for using the PairCalculator library
 *
 */
void myFunc(complex a, complex b) {}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** The Main of CPAIMD. It calls all the init functions.
 *
 */
//============================================================================

main::main(CkArgMsg *m) {

//============================================================================
/** Check arguments : Tell people what we are doing */

    if (m->argc < 3) {
      CkAbort("Usage: pgm config_file_charm config_file_piny");
    }//endif

    CkPrintf("\n================================================\n");
    CkPrintf("Starting Cpaimd-Charm-Driver Setup Phase \n");
    CkPrintf("---------------------------------------------------\n");
    CkPrintf("  Cpaimd-Charm-Driver running on %d processors. \n", CkNumPes());
    CkPrintf("  Reading Physics input from %s\n",m->argv[2]);
    CkPrintf("  Reading Driver  input from %s\n",m->argv[1]);
    CkPrintf("---------------------------------------------------\n\n");

//============================================================================    
/* Invoke PINY input class */

    CkCallback piny_callback (CkCallback::ignore);
    Interface_ctrl piny_interface (m->argv[2],piny_callback);

    CPcharmParaInfo *sim  = new CPcharmParaInfo();
    PhysicsParamTransfer::ParaInfoInit(sim);
    int ibinary_opt = sim->ibinary_opt;
    int natm_nl     = sim->natm_nl;

//============================================================================    
/* Invoke Ramkumar input class */

    CkPrintf("\n================================================\n");
    CkPrintf("Cpaimd-Charm-Driver input started \n");
    CkPrintf("------------------------------------------------\n");

    Config::readConfig(m->argv[1],config,sim->nstates,
                       sim->sizeX,sim->sizeY,sim->sizeZ,
                       sim->ntime,ibinary_opt,natm_nl);

    int numSfGrps    = config.numSfGrps;  // local copies are nice
    int doublePack   = config.doublePack;
    size2d sizeYZ    = size2d(sim->sizeY,sim->sizeZ);
    int gSpacePPC    = config.gSpacePPC;  
    int realSpacePPC = config.realSpacePPC;
    int rhoGPPC      = config.rhoGPPC;
    nchareG          = config.nchareG;
    sim->nchareG     = nchareG; 


    nstates          = config.nstates;    // globals on all procs
    sizeX            = sim->sizeX;

    config.print();

    CkPrintf("------------------------------------------------\n");
    CkPrintf("Cpaimd-Charm-Driver input completed \n");
    CkPrintf("================================================\n\n");

//============================================================================
// Set user trace events

     traceRegisterUserEvent("doRealFwFFT", doRealFwFFT_);
     traceRegisterUserEvent("doRealBwFFT", doRealBwFFT_);
     traceRegisterUserEvent("GspaceFwFFT", GspaceFwFFT_);
     traceRegisterUserEvent("GspaceBwFFT", GspaceBwFFT_);
     traceRegisterUserEvent("RhoRtoGxyFFT", RhoRtoGxzFFT_);
     traceRegisterUserEvent("RhoRtoGyFFT", RhoRtoGyFFT_);
     traceRegisterUserEvent("RhoDivRhoXFFT", RhoDivRhoXFFT_);
     traceRegisterUserEvent("RhoDivRhoYFFT", RhoDivRhoYFFT_);
     traceRegisterUserEvent("RhoDivRhoZFFT", RhoDivRhoZFFT_);
     traceRegisterUserEvent("VksofGFFT", VksofGFFT_);
     traceRegisterUserEvent("VksofRFFT", VksofRFFT_);
     traceRegisterUserEvent("DoFFTContribute", DoFFTContribute_);
     traceRegisterUserEvent("IntegrateModForces", IntegrateModForces_);
     traceRegisterUserEvent("Scalcmap", Scalcmap_);
     traceRegisterUserEvent("AcceptStructFact", AcceptStructFact_);
//============================================================================    
// Compute structure factor grp parameters and static map for chare arrays

    sim->numSfGrps   = numSfGrps;
    int natm_nl_grp_max;
    PhysicsParamTransfer::get_Sfgrp_max(natm_nl,config.numSfGrps, 
                                        &natm_nl_grp_max);
    sim->natm_nl_grp_max = natm_nl_grp_max;

    create_line_decomp_descriptor(sim);
    PhysicsParamTransfer::control_new_mapping_function(sim,doublePack);

    scProxy  = CProxy_CPcharmParaInfoGrp::ckNew(*sim);
    
//============================================================================    
// Create the multicast/reduction manager for array sections
// Create the parainfo group from sim
// Initialize chare arrays for real and g-space of states 

    mCastGrpId = CProxy_CkMulticastMgr::ckNew();


    init_state_chares(sizeYZ,natm_nl,natm_nl_grp_max,numSfGrps,doublePack,
                gSpacePPC,realSpacePPC,sim);

//============================================================================    
// Transfer parameters from physics to driver
//    read in atoms/states : create atoms group 
    
    control_physics_to_driver();


    //* planearray for paircalc and ortho
    int indexSize = nchareG;

    int* indexZ = new int[indexSize];
    for(int i=0, count=0; i<nchareG; i++){
        indexZ[count] = i;
        count++;
    }

//============================================================================ 
// Initialize paircalculators for Psi and Lambda

    init_pair_calculators( nstates,  indexSize, indexZ, gSpacePPC, doublePack, sim);

//============================================================================ 
// initialize Ortho

    init_ortho_chares(nstates, indexSize, indexZ);


//============================================================================ 
// Initialize the density chare arrays
    
    init_rho_chares(sizeYZ,gSpacePPC,realSpacePPC,rhoGPPC);

//============================================================================ 
// Initialize commlib strategies for later association and delegation

    init_commlib_strategies(sizeYZ, rhoGPPC, realSpacePPC);



//============================================================================
// clean up

    delete m;
    delete sim;
    delete [] indexZ;

//============================================================================

    CkPrintf("\n------------------------------------------------\n");
    CkPrintf("Cpaimd-Charm-Driver setup phase complete\n");
    CkPrintf("================================================\n\n");
//--------------------------------------------------------------------------
   }// end Main
//============================================================================

/**
 * Initialize paircalc1 Psi (sym) and paircalc2 Lambda (asym)
 */
void init_pair_calculators(int nstates, int indexSize, int *indexZ, int gSpacePPC, int doublePack, CPcharmParaInfo *sim)
{
  PRINT_LINE_STAR;
  PRINTF("Building Psi and Lambda Pair Calculators\n");
  PRINT_LINE_DASH;printf("\n");
  //-------------------------------------------------------------
  // Create mapping classes for Paircalcular
    CProxy_SCalcMap scMap_sym = CProxy_SCalcMap::ckNew(config.nstates,
                   sizeX / gSpacePPC,config.sGrainSize,CmiTrue,sim->nchareG, 
                   sim->lines_per_chareG, sim->pts_per_chareG) ;
    
    CProxy_SCalcMap scMap_asym = CProxy_SCalcMap::ckNew(config.nstates,
  	           sizeX / gSpacePPC,config.sGrainSize, CmiFalse,sim->nchareG, 
                   sim->lines_per_chareG, sim->pts_per_chareG);
    
    CkGroupID scalc_sym_id = scMap_sym.ckGetGroupID();
    CkGroupID scalc_asym_id = scMap_asym.ckGetGroupID();


  //-------------------------------------------------------------
  // Register the PCs

    int gsp_ep =  CkIndex_CP_State_GSpacePlane::__idx_acceptNewPsi_CkReductionMsg;

   //symmetric AKA Psi
    createPairCalculator(true, nstates, config.sGrainSize, indexSize, indexZ,  CkCallback(CkIndex_Ortho::start_calc(NULL), orthoProxy), &pairCalcID1, gsp_ep, gSpacePlaneProxy.ckGetArrayID(), 1, &scalc_sym_id, doublePack, config.conserveMemory,config.lbpaircalc, config.psipriority, mCastGrpId );

    CkArrayIndex2D myindex(0, 0);

    gsp_ep = CkIndex_CP_State_GSpacePlane::__idx_acceptLambda_CkReductionMsg;
    int myPack=0;

    //asymmetric AKA Lambda
    createPairCalculator(false, nstates,  config.sGrainSize, indexSize, indexZ,CkCallback(CkIndex_CP_State_GSpacePlane::acceptAllLambda(NULL), myindex, gSpacePlaneProxy.ckGetArrayID()), &pairCalcID2, gsp_ep, gSpacePlaneProxy.ckGetArrayID(), 1, &scalc_asym_id, myPack, config.conserveMemory,config.lbpaircalc, config.lambdapriority, mCastGrpId);

}


//============================================================================ 
/**
 * Initialize Commlib communication strategies  
 */ 
//                                          
//============================================================================
void init_commlib_strategies(size2d sizeYZ, int rhoGPPC, int realSpacePPC)
{
  PRINT_LINE_STAR;
  PRINTF("Building Commlib strategies\n");
  PRINT_LINE_DASH;printf("\n");
    int i = 0;

    if (config.useCommlib) {        
        int nsrcelements = sizeYZ[1]/rhoGPPC;
        int ndestelements = sizeYZ[0]/realSpacePPC;

        CkArrayIndexMax *srcelements = new CkArrayIndexMax[sizeYZ[1]/rhoGPPC];
        for (i = 0; i < sizeYZ[1]/rhoGPPC; i++) {
            srcelements[i] = CkArrayIndex1D(i);
        }

        CkArrayIndexMax *destelements = new 
            CkArrayIndexMax[sizeYZ[0]/realSpacePPC];
        for(i = 0; i < sizeYZ[0]/realSpacePPC; i++) {
            destelements[i] = CkArrayIndex1D(i);
        }

        CharmStrategy *strat = new EachToManyMulticastStrategy
            (USE_MESH, rhoGProxy.ckGetArrayID(), rhoRealProxy.ckGetArrayID(), 
             nsrcelements, srcelements, ndestelements, destelements);
        
        srcelements = new CkArrayIndexMax[sizeYZ[1]/rhoGPPC];
        for (i = 0; i < sizeYZ[1]/rhoGPPC; i++) 
            srcelements[i] = CkArrayIndex1D(i);
        
        destelements = new CkArrayIndexMax[sizeYZ[0]/realSpacePPC];
        for(i = 0; i < sizeYZ[0]/realSpacePPC; i++)
            destelements[i] = CkArrayIndex1D(i);
        
        CharmStrategy *real_strat = new EachToManyMulticastStrategy
            (USE_MESH, rhoRealProxy.ckGetArrayID(), rhoGProxy.ckGetArrayID(),
             ndestelements, destelements, nsrcelements, srcelements);
        
        commInstance = ComlibRegister(strat);
        commRealInstance= ComlibRegister(real_strat);
    }

    if (config.useCommlibMulticast) {
        DirectMulticastStrategy *dstrat = new DirectMulticastStrategy
	    (realSpacePlaneProxy.ckGetArrayID(),1);
        
        RingMulticastStrategy *rstrat = new RingMulticastStrategy
            (realSpacePlaneProxy.ckGetArrayID(),1);
        
        MultiRingMulticast *mrstrat = new MultiRingMulticast
            (realSpacePlaneProxy.ckGetArrayID(),1);
        
        DirectMulticastStrategy *d1strat = new DirectMulticastStrategy
            (particlePlaneProxy.ckGetArrayID(),1);

        RingMulticastStrategy *r1strat = new RingMulticastStrategy
            (particlePlaneProxy.ckGetArrayID(),1);

        MultiRingMulticast *mr1strat = new MultiRingMulticast
            (particlePlaneProxy.ckGetArrayID(),1);

	if(CkNumNodes()>64) //multiring should be good on large runs, but not on BG/L
	  {
	      mcastInstance=ComlibRegister(dstrat);
	      mcastInstancePP=ComlibRegister(mr1strat);
	    }
	else
	  {
	      mcastInstance=ComlibRegister(rstrat);
	      mcastInstancePP=ComlibRegister(r1strat);
	  }
	
    }        
    // end Sameer's new communication strategies 
}

//============================================================================
/**
 ** Create stuff for ortho which PC invokes by section reduction
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void init_ortho_chares(int nstates, int indexSize, int *indexZ)
{
  //-------------------------------------------------------------
  PRINT_LINE_STAR;
  PRINTF("Building Ortho Chares\n");
  PRINT_LINE_DASH;printf("\n");

    int chunks = (nstates + config.sGrainSize - 1) / config.sGrainSize;
    CProxy_OrthoMap orthoMap = CProxy_OrthoMap::ckNew(chunks);
    CkArrayOptions orthoOpts;
    orthoOpts.setMap(orthoMap);
    orthoProxy = CProxy_Ortho::ckNew(orthoOpts);

    CkCallback *orthoReduction = new CkCallback(CkIndex_Ortho::collect_error(NULL), orthoProxy(0, 0));
    orthoProxy.ckSetReductionClient(orthoReduction);
    
    // extra triangle ortho elements are really a waste of our time
    // and resources, but we don't have a triangular solver for
    // inv_square, so we'll just make do.

    // They need to exist solely so that the inv_sq method can work.
    // So we need to copy their mirror elements data into them.
    // then when complete they need to know not to call finishpaircalc.
    // Because their redundant data has nowhere to go.

    /* create matrix multiplication objects */
    CLA_Matrix_interface matA1, matB1, matC1;
    CLA_Matrix_interface matA2, matB2, matC2;
    CLA_Matrix_interface matA3, matB3, matC3;

    CkCallback ortho_ready_cb = CkCallback(CkIndex_Ortho::all_ready(),
     orthoProxy(0, 0));
    make_multiplier(&matA1, &matB1, &matC1, orthoProxy, orthoProxy, orthoProxy,
     nstates, nstates, nstates, config.sGrainSize, config.sGrainSize,
     config.sGrainSize, ortho_ready_cb, ortho_ready_cb, ortho_ready_cb,
     MM_ALG_2D);
    make_multiplier(&matA2, &matB2, &matC2, orthoProxy, orthoProxy, orthoProxy,
     nstates, nstates, nstates, config.sGrainSize, config.sGrainSize,
     config.sGrainSize, ortho_ready_cb, ortho_ready_cb, ortho_ready_cb,
     MM_ALG_2D);
    make_multiplier(&matA3, &matB3, &matC3, orthoProxy, orthoProxy, orthoProxy,
     nstates, nstates, nstates, config.sGrainSize, config.sGrainSize,
     config.sGrainSize, ortho_ready_cb, ortho_ready_cb, ortho_ready_cb,
     MM_ALG_2D);


    for (int s1 = 0; s1 < nstates; s1 += config.sGrainSize)
      for (int s2 = 0; s2 < nstates; s2 += config.sGrainSize) {
	int indX = s1 / config.sGrainSize;
	int indY = s2 / config.sGrainSize;
	orthoProxy(indX, indY).insert(config.sGrainSize, config.sGrainSize,
         matA1, matB1, matC1, matA2, matB2, matC2, matA3, matB3, matC3);
      }
    orthoProxy.doneInserting();
    orthoProxy.makeSections(indexSize, indexZ);

}

//============================================================================
/**
 *Create the array elements for the GSpace, Particle and Real Space planes
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void init_state_chares(size2d sizeYZ, int natm_nl,int natm_nl_grp_max,int numSfGrps,
                 int doublePack,int gSpacePPC,int realSpacePPC,
                 CPcharmParaInfo *sim)

//============================================================================
   { //begin routine 
//============================================================================
    /*
     * Set up the map variables used to control the location of the
     * 2d-array elements over processors   
     */

  PRINT_LINE_STAR;
  PRINTF("Building G-space and R-space Chares\n");
  PRINT_LINE_DASH;printf("\n");

    CProxy_GSMap gsMap;
    CProxy_RSMap rsMap;

    rsMap = CProxy_RSMap::ckNew();

    CkArrayOptions realSpaceOpts;
    realSpaceOpts.setMap(rsMap);
    size2d sizeRealPlane(sizeYZ[1], sizeX);

    realSpacePlaneProxy = CProxy_CP_State_RealSpacePlane::ckNew(sizeRealPlane,gSpacePPC,
                                                                realSpacePPC, 
                                                                realSpaceOpts);

								
    fftCacheProxy = CProxy_FFTcache::ckNew(sizeRealPlane, realSpacePPC);
    sfCacheProxy = CProxy_StructFactCache::ckNew(numSfGrps,natm_nl,natm_nl_grp_max);
    sfCompProxy = CProxy_StructureFactor::ckNew();
    
    gsMap = CProxy_GSMap::ckNew(sim->nchareG, sim->lines_per_chareG, sim->pts_per_chareG);

    CkArrayOptions gSpaceOpts;
    gSpaceOpts.setMap(gsMap);

    gSpacePlaneProxy = CProxy_CP_State_GSpacePlane::ckNew(sizeX, sizeYZ, gSpacePPC, 
                                                          realSpacePPC, 
                                                          config.sGrainSize, gSpaceOpts);

    // We bind the particlePlane array to the gSpacePlane array migrate together
    CkArrayOptions particleOpts;
    particleOpts.setMap(gsMap); // the maps for both the arrays are the same
    particleOpts.bindTo(gSpacePlaneProxy);
    particlePlaneProxy = CProxy_CP_State_ParticlePlane::ckNew(nchareG, sizeYZ[0], sizeYZ[1],   
		      gSpacePPC,numSfGrps,natm_nl,natm_nl_grp_max,particleOpts);

    /*
     * Insert the planes in the particle plane array, gSpacePlane array
     */
    gSpacePlaneProxy.setReductionClient(doneInit, (void *) NULL);
    realSpacePlaneProxy.setReductionClient(doneInit, (void *) NULL);
    int s,x;
//    CkPrintf("making nstates %d nchareG %d gspace objects\n",nstates,nchareG);
    for (s = 0; s < nstates; s++){
      for (x = 0; x <nchareG; x++){
             gSpacePlaneProxy(s, x).insert(sizeX, sizeYZ, gSpacePPC, 
                                    realSpacePPC,config.sGrainSize);
             particlePlaneProxy(s, x).insert(sizeX, sizeYZ[0], sizeYZ[1],   
				      gSpacePPC,numSfGrps,natm_nl,natm_nl_grp_max);
      }
    }


    gSpacePlaneProxy.doneInserting();
    particlePlaneProxy.doneInserting();
    particlePlaneProxy.setReductionClient(doneCreatingPP, (void *) NULL);
    
    /*
     * Insert the planes in the real space plane array
     */
    int y;
    for (s = 0;  s < nstates; s++)
        for (y = 0; y < sizeYZ[0]; y += realSpacePPC)
            realSpacePlaneProxy(s, y).insert(sizeRealPlane, gSpacePPC, 
                                             realSpacePPC);
    realSpacePlaneProxy.doneInserting();
    // 

    int *nsend   = new int[nchareG];
    int **listpe = new int * [nchareG];
    int numproc  = CkNumPes();
    int *gspace_proc = new int [numproc];

    for(int i =0;i<numproc;i++){gspace_proc[i]=0;}
    for(int j=0;j<nchareG;j++){   
      listpe[j]= new int[nstates];
      nsend[j]=0;
      for(int i=0;i<nstates;i++){
	listpe[j][i]=cheesyhackgsprocNum(sim, i,j);
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
    CkPrintf("============================\n");
    CkPrintf("Structure factor chareG dests\n");
    CkPrintf("Number of g-space chares : %d\n",nchareG);    
    CkPrintf("---------------------------\n");
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
    CkPrintf("---------------------------\n");
    CkPrintf("SFSends min %d max %d avg %g\n",minsend,maxsend,avgsend/(double)chareG_use);
    CkPrintf("============================\n");
    // Insert the objects into the StructureFactor array
    int dupmax=config.nstates;
    if (config.numSfDups<dupmax)
	dupmax=config.numSfDups;
    config.numSfDups=dupmax;
    for (int dup=0; dup<dupmax; dup++)
      for (x = 0; x < nchareG; x += gSpacePPC)
      {
	  int num_dup, istart, iend;
	  get_grp_params( nsend[x],  config.numSfDups,  dup, x ,&num_dup,  &istart, &iend);
	  int pe_ind=istart;
	  if(x%2==0)
	      pe_ind=iend;
	  for (int AtmGrp=0; AtmGrp<numSfGrps; AtmGrp++)
	  {
//		CkPrintf("inserting SF[%d %d %d] to pe %d\n",AtmGrp,x,dup,listpe[x][pe_ind]);
	      sfCompProxy(AtmGrp, x, dup).insert(numSfGrps, config.numSfDups, natm_nl_grp_max,  num_dup, &(listpe[x][istart]),atmGrpMap(istart, num_dup, nsend[x], listpe[x], AtmGrp, dup,x));	      
	      pe_ind++;
	      if(pe_ind>nsend[x])
		  pe_ind=0;
	  }
      }
    sfCompProxy.doneInserting();

    for(int j=0;j<nchareG;j++)
      delete [] listpe[j];
    delete [] listpe;
    if(config.useCommlib) {
        // Set some com strategy of Sameer
        ssInstance = CkGetComlibInstance();
        StreamingStrategy *strat = new StreamingStrategy(0.2,10);
        //strat->enableShortArrayMessagePacking();
        ssInstance.setStrategy(strat);
        
        // Set some com strategy of Sameer
        StreamingStrategy *mstrat = new StreamingStrategy(0.2,5);
        //mstrat->enableShortArrayMessagePacking();
        mssInstance= ComlibRegister(mstrat);    
    }

  printf("\n");
  PRINT_LINE_DASH;
  PRINTF("Completed g-space chare build\n");
  PRINT_LINE_STAR;printf("\n");

//----------------------------------------------------------------------------
    }//end routine
//============================================================================


//============================================================================
// Creating arrays CP_Rho_GSpacePlane, CP_Rho_GSpacePlaneHelper 
// and CP_Rho_RealSpacePlane
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void init_rho_chares(size2d sizeYZ, int gSpacePPC, int realSpacePPC, int rhoGPPC)

//============================================================================
    {//begin routine
//============================================================================
/*
 * create the array for real-space densities (one-dimensional chare array)
 * and hartree energy computation
 */    
//============================================================================
    

        CkAssert(config.rhoGPPC == 1);
        rhoGHelperProxy = CProxy_CP_Rho_GSpacePlaneHelper::ckNew();
	rhoGHelperProxy.setReductionClient(printEnergyHart, NULL);

        int z, y;
        int helperSize = sizeYZ[0]/config.rhoGHelpers;

        if( (sizeYZ[0] % config.rhoGHelpers) !=0 ){
         CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
         CkPrintf("Helper size must be a mod of %d.\n",sizeYZ[0]);
         CkPrintf("Please fix your cpaimd_config.\n");
         CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
         CkExit();
	}//endif

        int pe = 0;
        for (z = 0; z < sizeYZ[1]; z++){
            for (y = 0; y < sizeYZ[0]; y += helperSize) {
                rhoGHelperProxy(z,y).insert(sizeX, sizeYZ, y, pe, 0);
                pe = (pe + 1)%CkNumPes();
            }
        }
        rhoGHelperProxy.doneInserting();

    
    int pestride =  CkNumPes() / (sizeYZ[0]/realSpacePPC);
    if(pestride < 1)
      pestride = 1;
    
    if (sizeYZ[0]/realSpacePPC == sizeYZ[1]/rhoGPPC ) {
      rhoGProxy = CProxy_CP_Rho_GSpacePlane::ckNew();
      rhoRealProxy = CProxy_CP_Rho_RealSpacePlane::ckNew();
     
      int peg = 0;
      int per = pestride/2;
      int i;

      bool fftuseCommlib = config.fftuseCommlib;
      ComlibInstanceHandle fftcommInstance;
      if (fftuseCommlib) {        
	  int period_in_ms = 1, nmsgs = 1000;
	  StreamingStrategy * strat = new StreamingStrategy(period_in_ms, nmsgs);
	  fftcommInstance = CkGetComlibInstance();
	  fftcommInstance.setStrategy(strat);
      }
     for (i = 0; i < sizeYZ[1]/rhoGPPC; i++) {

	if (peg >= CkNumPes())peg = 0;
	if (per >= CkNumPes())per = pestride/2;

	rhoGProxy[i].insert(sizeX, sizeYZ, realSpacePPC, rhoGPPC, fftuseCommlib, 
                            fftcommInstance, peg);
	rhoRealProxy[i].insert(sizeX, sizeYZ, realSpacePPC,
			       rhoGPPC, fftuseCommlib, fftcommInstance,per);
	
	peg += pestride;
	per += pestride;
                                                                                   
      }
                                                                                   
      rhoGProxy.doneInserting();
      rhoRealProxy.doneInserting();
    }else{ 
      CkAbort("sizeYZ[0]/realSpacePPC != sizeYZ[1]/rhoGPPC");
    }
    
    rhoRealProxy.setReductionClient(printEnergyEexc, 0);

//===========================================================================
}//end routine
//============================================================================



//============================================================================
// Get the atoms and the parainfo
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void control_physics_to_driver(){

//============================================================================
// make a group : create a proxy for the atom class and also a reduction client

    PhysicsAtomPosInit *PhysicsAtom  = new PhysicsAtomPosInit();
    int natm    = PhysicsAtom->natm_tot;
    int natm_nl = PhysicsAtom->natm_nl;
    Atom *atoms = new Atom[natm];
    PhysicsAtom->DriverAtomInit(atoms);

    atomsGrpProxy = CProxy_AtomsGrp::ckNew(natm,natm_nl,atoms);
    delete [] atoms;
    delete PhysicsAtom;

//=====================================================================
// Make a group for the energies

    egroupProxy = CProxy_EnergyGroup::ckNew(); 

//----------------------------------------------------------------------------
    }
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
static CkReductionMsg *complexSum(int nMsg, CkReductionMsg **msgs)
{
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
void get_grp_params(int natm_nl, int numSfGrps, int indexSfGrp, int planeIndex,
		    int *n_ret, int *istrt_ret, int *iend_ret)
{

   int n     = (natm_nl/numSfGrps);
   int m     = (natm_nl % numSfGrps);

   int istrt = n*indexSfGrp;
   if(indexSfGrp>=m){istrt += m;}
   if(indexSfGrp<m) {istrt += indexSfGrp;}
   if(indexSfGrp<m) {n++;}
   int iend  = n+istrt;
   if(numSfGrps>natm_nl)
     if(m>=indexSfGrp)
       {
	 n=0;
	 istrt=natm_nl+1;
	 iend=natm_nl;
         CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
	 CkPrintf("Redundant DupSF for chare-G %d\n",planeIndex);
         CkPrintf("At present this hangs, so out you go\n");
         CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
         CkExit();
       }
  
   (*n_ret)     = n;
   (*istrt_ret) = istrt;
   (*iend_ret)  = iend;
}
//============================================================================



//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
int atmGrpMap(int istart, int nsend, int listsize, int *listpe, int AtmGrp, 
              int dup, int planeIndex)
{

  int numSfDups=config.numSfDups;
  if(listsize <= numSfDups)
    { //dup list not unique
      return listpe[(dup % listsize)];
    }
  if(nsend <= config.numSfGrps)
    {  //atom group list not unique
      return listpe[(istart + (AtmGrp % nsend))];
    }
  // nsend must be > numSfGrps and its list is unique
  return listpe[(istart + ((AtmGrp + config.numSfGrps * planeIndex)%nsend))];

}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
int cheesyhackgsprocNum(CPcharmParaInfo *sim,int state, int plane) {
//============================================================================

  CkArrayIndex2D idx2d(state,plane);

  int pe        = 0;
  int numChareG = 0;
  
  if(config.doublePack) {
     numChareG = nchareG;
  }else{
     CkAbort("not doublepack is broken\n");
  }//endif
  
  double state_load=0.0;
  if(state_load <= 0.0 ){
    for(int x = 0; x  < numChareG; x++) {
      double curload = 0.0;
      double sload = 0.0;
      hackGSpacePlaneLoad(sim, x, &curload, &sload);
      state_load += curload;
    }//endfor
  }//endif

  int pes_per_state = config.GpesPerState;
  int np            = CkNumPes()/pes_per_state;
  
  if(np < 1){np = 1;}
  
  int partition_nstates = config.nstates / np;
  if(config.nstates % np != 0){partition_nstates ++;}
  
  int start_pe    = (idx2d.index[0]/partition_nstates) * pes_per_state;
  int start_state = (idx2d.index[0]/partition_nstates) * partition_nstates;
  
  double cum_load     = 0.0;
  double average_load = state_load * config.nstates / CkNumPes();

  for(int x = 0; x < numChareG; x++) {
    for(int s = 0; s < partition_nstates; s++) {
      double curload = 0.0;
      double sload = 0.0;
      
      hackGSpacePlaneLoad(sim,x, &curload, &sload);

      cum_load += curload;

      if((idx2d.index[0] == s + start_state) && idx2d.index[1] == x) {
      
	double dpe = 0.0;
	dpe = cum_load / average_load;

	pe = (int)dpe;
	pe = pe % pes_per_state;
	pe += start_pe;

	return pe % CkNumPes();
      }//endif
    }//endfor : s
  }//endfor : x

  //  CkPrintf("Warning pe not found for index [%d, %d]\n",idx2d.index[0],idx2d.index[1]);
  return (idx2d.index[0]*1037+idx2d.index[1])%CkNumPes();

//============================================================================
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void hackGSpacePlaneLoad(CPcharmParaInfo *sim,int idx, double *line_load, 
                         double *pt_load){
//============================================================================

  int nchareG             = sim->nchareG;
  double *lines_per_chareG = sim->lines_per_chareG;
  double *pts_per_chareG   = sim->pts_per_chareG;

  if( (idx < 0) || (idx >= nchareG) ){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("GSpace plane index %d out of range: 0 < idx > %d\n",idx,nchareG);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  line_load[0] = lines_per_chareG[idx];
  pt_load[0]   =   pts_per_chareG[idx];

//============================================================================
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void mapOutput()
//============================================================================
    {//begin routine
//============================================================================

	int i;
        i = 1;

//============================================================================
  }//end routine
//============================================================================


//============================================================================
#include "cpaimd.def.h"
//============================================================================
