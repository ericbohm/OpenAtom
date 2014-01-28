//#include "CP_State_GSpacePlane.h"
//
//#ifndef GSPACE_RTH_H
//#define GSPACE_RTH_H
//
//extern Config config;
///** @defgroup ControlFlow ControlFlow
//    @{
//
//*/
///// This routine defines the control flow of GSpacePlane
//RTH_Routine_locals(GSpaceDriver,driveGSpace)
//RTH_Routine_code(GSpaceDriver,driveGSpace)
//    /// Indefinite loop. Exits only when one of the conditions succeed/fail
//while(1)
//  { 
//    ///________________________________________________________________________________
//    /** # (I) Compute the forces and coef evolution code block
//     *
//     * If this is minimization (& gen_wave?) OR if this is not the first step  
//     */
//    if( (CPcharmParaInfo::get()->cp_min_opt==1 && CPcharmParaInfo::get()->gen_wave==0) || c->isFirstStep==false)
//      {
//	c->isEnergyReductionDone = false;
//	c->isAtomIntegrationDone = false;
//	/// If psi has not been turned off, do everything; else, move only the atoms  
//
//#ifndef _CP_DEBUG_PSI_OFF_  
//	//------------------------------------------------------------------------
//	/// ## (A) Start new iteration, reset the counters
//	c->myGSpaceObj->thisProxy(c->thisIndex.x,c->thisIndex.y).startNewIter();
//	//------------------------------------------------------------------------
//	/// ## (B) If nonlocal energy computations are allowed, start SF/computeZ of NL forces.
//#ifndef _CP_DEBUG_SFNL_OFF_
//	/// If EES is turned off, trigger nonlocal computations
//	if (c->ees_nonlocal==0 && CPcharmParaInfo::get()->natm_nl!=0){
//	  /// Launch the structure factor and the Z matrix computations
//	  c->releaseSFComputeZ();
//	  /// If NL computations have been barriered, then wait (GSpaceDriver::doneNLForces or GSpaceDriver::allDoneNLForces resumes)
//#ifdef BARRIER_CP_PARTICLEPLANE_NONLOCAL
//	  RTH_Suspend();
//#endif
//	}
//#endif
//	//------------------------------------------------------------------------
//	/// If VKS forces are allowed, do the forward and back FFTs (The nonlocal ees methods will be launched elsewhere to overlap appropriately with other work)
//#ifndef _CP_DEBUG_VKS_OFF_ 
//	/// ## (D) FFT psi(gx,gy,gz)->psi(gx,gy,z)
//	c->myGSpaceObj->thisProxy(c->thisIndex.x,c->thisIndex.y).doFFT();
//	/// Send psi to RealSpace
//	c->myGSpaceObj->sendFFTData();
//	/// Wait for (psi*vks) = F[gx,gy,z] to arive from RealSpace (CP_State_GSpacePlane::acceptIFFT resumes)
//	RTH_Suspend();
//	/// ## (F) As the Psi forces have come back to us, do inverse FFT
//	c->myGSpaceObj->thisProxy(c->thisIndex.x,c->thisIndex.y).doIFFT();
//	/// If barriered, then wait for all GSpace chares to finish inverse FFT (GSpaceDriver::allDoneIFFT resumes)
//#ifdef BARRIER_CP_GSPACE_IFFT
//	RTH_Suspend();
//#endif
//	/// else, if VKS forces are not turned on, we are responsible for triggering ees nonlocals
//#else
//	/// If nonlocals are turned on, 
//#ifndef _CP_DEBUG_SFNL_OFF_
//	/// If ees methods are turned on, then trigger the ees computations
//	if(c->ees_nonlocal==1 && CPcharmParaInfo::get()->natm_nl!=0)
//	  c->startNonLocalEes(c->myGSpaceObj->iteration);
//#endif
//	/// Even if we're not doing the FFT loop, just set the flag to fool ourselves if/when we check later
//	c->myGSpaceObj->doneDoingIFFT = true;
//#endif
//	//------------------------------------------------------------------------
//	/// ## (G) If NL computations are allowed
//#ifndef _CP_DEBUG_SFNL_OFF_
//	/// If NL computations are not ALL done (all channels/atm types), suspend (GSpaceDriver::doneNLForces or GSpaceDriver::allDoneNLForces resume)
//	if (!c->areNLForcesDone)
//	  RTH_Suspend();
//#endif
//	//------------------------------------------------------------------------
//	/// ## (H) Added up all force contributions
//	c->myGSpaceObj->combineForcesGetEke();
//#endif // only move the atoms
//	//------------------------------------------------------------------------
//	/// ## (I) The atoms can't go until cp forces are completely finished (all chares)
//	// However, the atoms can overlap with all this lambda, psi stuff.
//	c->myGSpaceObj->launchAtoms();
//	//------------------------------------------------------------------------
//	/// ## (G) Add contraint forces (rotate forces to non-orthogonal frame)
//#ifndef _CP_DEBUG_PSI_OFF_  // you are moving everything
//	/// Reset ParticlePlane's flag to false
//	if(CPcharmParaInfo::get()->natm_nl!=0)c->areNLForcesDone = false;
//	c->myGSpaceObj->sendLambda();
//#ifndef _CP_DEBUG_ORTHO_OFF_
//	RTH_Suspend(); // wait for forces to be fixed up  
//	// CP_State_GSpacePlane::acceptLambda resumes
//#endif
//	//------------------------------------------------------------------------
//	/// ## (H) Get sum sq forces (even under dynamics its good to have) : also cg thingy
//	c->myGSpaceObj->computeCgOverlap();
//	RTH_Suspend(); /// wait for cg reduction : CP_State_GSpacePlane::psiCgOvlap resumes
//	//------------------------------------------------------------------------
//	/// ## (I) If user wants state output, and this is dynamics (>2nd iteration, as nothing moves in dynamics iter 1) ... then write the states to file at the correct intervals or at the end of the simulation. Then wait for the data reduction and writes to finish (GSpaceDriver::allDoneWritingPsi resumes)
//	if (config.stateOutput==1 && CPcharmParaInfo::get()->cp_min_opt==0 && c->myGSpaceObj->iteration>1)
//	  if ( (c->myGSpaceObj->iteration-1)%CPcharmParaInfo::get()->ndump_frq==0 || c->myGSpaceObj->iteration==config.maxIter || c->myGSpaceObj->outputFlag==1 ) 
//	    {
//	      c->myGSpaceObj->writeStateDumpFile();
//	      RTH_Suspend(); 
//	    }
//	//------------------------------------------------------------------------
//	/// ## (J) Evolve the electrons to the next step
//	c->myGSpaceObj->integrateModForce();
//	c->myGSpaceObj->sendRedPsi();  // Sync Redundant psi entries
//	if(c->myGSpaceObj->iRecvRedPsi==0)
//	  RTH_Suspend();  /// CP_State_GSpacePlane::acceptRedPsi resumes
//	c->myGSpaceObj->doneRedPsiIntegrate(); // after integrate AND acceptRedPsi
//#endif  // you are moving everyting
//      }// endif determine entry point
//        
//    ///________________________________________________________________________________        
//    /// # (II) Orthogonalization and output code block
//
//#ifndef _CP_DEBUG_PSI_OFF_
//    //------------------------------------------------------------------------
//    /// ## (A) Orthogonalize
//    /// Send psi to the PairCalculators 
//    c->myGSpaceObj->sendPsi();
//    /// If ortho has not been turned off, wait for new Psi (CP_State_GSpacePlane::doNewPsi resumes)
//#ifndef _CP_DEBUG_ORTHO_OFF_
//    RTH_Suspend();  
//#endif
//    //------------------------------------------------------------------------
//    /// ## (B) If user wants state output, and this is minimization ... then write the states to file at the correct intervals or at the end of the simulation. Then wait for the data reduction and writes to finish (GSpaceDriver::allDoneWritingPsi resumes)
//    if (config.stateOutput==1 && CPcharmParaInfo::get()->cp_min_opt==1 && c->myGSpaceObj->iteration>0)
//      if ( c->myGSpaceObj->iteration%CPcharmParaInfo::get()->ndump_frq==0 || c->myGSpaceObj->iteration==config.maxIter || c->myGSpaceObj->exitFlag==1 ) 
//	{
//	  c->myGSpaceObj->writeStateDumpFile();
//	  RTH_Suspend();
//	}
//    //------------------------------------------------------------------------
//    /// ## (C) If this is dynamics, check if we need a tolerance update (Velocity Norb rotation)
//    if(CPcharmParaInfo::get()->cp_min_opt == 0 && c->myGSpaceObj->iteration>1 && c->myGSpaceObj->iteration<= config.maxIter)
//      {
//	/// If Ortho told us that we need a tolerance update
//	if(c->isPsiVupdateNeeded)
//	  {
//	    /// 
//	    c->myGSpaceObj->acceptedVPsi = false;
//	    /// Sync PsiV at time t
//	    c->myGSpaceObj->sendRedPsiV();
//	    /// If I have not received all my PsiV CPcharmParaInfo::get()->, wait for it (CP_State_GSpacePlane::acceptRedPsiV resumes)
//	    if(c->myGSpaceObj->iRecvRedPsiV==0)
//	      RTH_Suspend();
//	    /// Tuck away your g=0 plane velocities
//	    c->myGSpaceObj->doneRedPsiVIntegrate();
//	    /// Rotate yourself 
//	    c->myGSpaceObj->sendPsiV();
//	    /// Wait for new PsiV (CP_State_GSpacePlane::doNewPsiV resumes)
//	    RTH_Suspend();
//	    /// The PsiV update is done. Reset the flag
//	    c->isPsiVupdateNeeded = false;
//	  }
//      }
//    //------------------------------------------------------------------------
//    /// ## (D) Check for Energy reduction completion : should just be a safety
//    if(c->isEnergyReductionDone==false && c->myGSpaceObj->iteration>0)
//      c->waitingForEnergy = true;
//#endif
//    //------------------------------------------------------------------------
//    /// ## (E) Check for atom integration : should just be a safety
//    if(c->isAtomIntegrationDone==false && c->myGSpaceObj->iteration>0)
//      c->waitingForAtoms = true;
//
//    //------------------------------------------------------------------------
//    /// ## (F) If the atom or energy stuff is slow, relax for a bit (doneMovingAtoms or GSpaceDriver::doneComputingEnergy resumes)
//    if(c->waitingForAtoms || c->waitingForEnergy)
//      {
//#ifdef _CP_DEBUG_WARN_SUSPEND_
//	CkPrintf("GSpace[%d,%d] Suspending atm/energy on proc %d : atom:%d energy:%d\n",c->thisIndex.x,c->thisIndex.y,CkMyPe(),c->waitingForAtoms,c->waitingForEnergy);
//#endif
//	RTH_Suspend();
//      }
//    //------------------------------------------------------------------------
//    /// ## (G) If you have triggered an exit condition just chill until ckexit
//    if(c->myGSpaceObj->cleanExitCalled==1)
//      RTH_Suspend();
//    /// Once you've been through this loop once, its not the first step anymore!
//    c->isFirstStep = false;   
//    /*
//    // send orthoT now
//    if(CPcharmParaInfo::get()->cp_min_opt == 0)
//    c->myGSpaceObj->launchOrthoT();
//    */
//  }
//RTH_Routine_end(GSpaceDriver,driveGSpace)
///*@}*/
//#endif // GSPACE_RTH_H
