#include "charm++.h"
#include "util.h"
#include "cpaimd.h"
#include "groups.h"
#include "InstanceController.h"
extern int nstates;

/* ostensibly the InstanceController may need to know about everything */
extern CkVec < CkVec <int> > UplaneUsedByNLZ;
extern CProxy_CPcharmParaInfoGrp         scProxy;
extern CkVec <CProxy_CP_State_GSpacePlane>       UgSpacePlaneProxy;
extern CkVec <CProxy_CP_State_ParticlePlane>     UparticlePlaneProxy;
extern CkVec <CProxy_CP_State_RealParticlePlane> UrealParticlePlaneProxy;
extern CkVec <CProxy_CP_State_RealSpacePlane>    UrealSpacePlaneProxy;
extern CkVec <CProxy_CP_Rho_RealSpacePlane>      UrhoRealProxy;
extern CkVec <CProxy_CP_Rho_GSpacePlane>         UrhoGProxy;
extern CkVec <CProxy_CP_Rho_RHartExt>            UrhoRHartExtProxy;
extern CkVec <CProxy_CP_Rho_GHartExt>            UrhoGHartExtProxy;
extern CkVec <CProxy_Ortho>                      UorthoProxy;
extern CkVec <CProxy_AtomsGrp>                   UatomsGrpProxy;
extern CkVec <CProxy_EnergyGroup>                UegroupProxy;
extern CkVec <CProxy_FFTcache>                   UfftCacheProxy;
extern CkVec <CProxy_StructFactCache>            UsfCacheProxy;
extern CkVec <CProxy_StructureFactor>            UsfCompProxy;
extern CkVec <CProxy_eesCache>                   UeesCacheProxy;
extern CkVec <CProxy_OrthoHelper>                UorthoHelperProxy;
extern CProxy_TimeKeeper                 TimeKeeperProxy;

extern CkVec <UberCollection>         UberAlles;

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void InstanceController::doneInit(CkReductionMsg *msg){
//============================================================================
  CkPrintf("{%d} Done_init for %d userflag %d\n",thisIndex, (int)((int *)msg->getData())[0],msg->getUserFlag());
  CkAssert(msg->getUserFlag()==thisIndex);
  delete msg;
    double newtime=CmiWallTimer();
    CkAssert(done_init<5);

    if(done_init<4){
      CkPrintf("{%d} Completed chare instantiation phase %d in %g\n",thisIndex,done_init+1,newtime-Timer);
      Timer=newtime;
    }else{
      CkPrintf("{%d} Completed chare data acquisition phase %d in %g\n",thisIndex, done_init+1,newtime-Timer);
      //      PRINT_LINE_DASH;
      CkPrintf("{%d} Chare array launch and initialization complete       \n",thisIndex);
      //      PRINT_LINE_STAR; printf("\n");
      Timer=newtime;
    }//endif
    if (done_init==1)
      { // kick off post constructor inits

	UrhoRealProxy[thisIndex].init();
	UrhoGProxy[thisIndex].init();
	UrhoGHartExtProxy[thisIndex].init();
      if(scProxy.ckLocalBranch()->cpcharmParaInfo->ees_eext_on)
	{UrhoRHartExtProxy[thisIndex].init();}


      }
    if (done_init == 3){
      // 2nd to last, we do this after we know gsp, pp, and rp exist
      if(scProxy.ckLocalBranch()->cpcharmParaInfo->ees_nloc_on==1)
	{UrealParticlePlaneProxy[thisIndex].init();}
      // its completion triggers the final phase

      // kick off file reading in gspace
      CkPrintf("{%d} Initiating import of states\n",thisIndex);
      CkPrintf("{%d} IC uGSpacePlaneProxy[%d] is %d\n",thisIndex,thisIndex, CkGroupID(UgSpacePlaneProxy[thisIndex].ckGetArrayID()).idx);
      for(int s=0;s<nstates;s++) {
        UgSpacePlaneProxy[thisIndex](s,UplaneUsedByNLZ[thisIndex][s]).readFile();
      } //endfor

      /* for(int s=0;s<nstates;s++){ ifndef USE_TOPOMAP
        gSpacePlaneProxy(s,0).readFile();
      }//endfor */

    }//endif
    if (done_init >= 4) {
      if (done_init == 4){ 
	//          PRINT_LINE_STAR;
          if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt==1){
            CkPrintf("{%d} Running Open Atom CP Minimization: \n",thisIndex);
	  }else{
            CkPrintf("{%d} Running Open Atom CP Dynamics: \n",thisIndex);
	  }//endif
	  //          PRINT_LINE_STAR; CkPrintf("\n");
	  //          PRINT_LINE_STAR;
	  UgSpacePlaneProxy[thisIndex].run();
      }//endif
    }
    done_init++;
}
//============================================================================

void InstanceController::printEnergyHart(CkReductionMsg *msg){
  //  double ehart = 0, eext = 0.0, ewd = 0.0;
  void *data=msg->getData();
  double ehart = ((double *)data)[0];
  double eext = ((double *)data)[1];
  double ewd  = ((double *)data)[2];
  
  CkPrintf("{%d} EHART       = %5.8lf\n", thisIndex, ehart);
  CkPrintf("{%d} EExt        = %5.8lf\n", thisIndex, eext);
  CkPrintf("{%d} EWALD_recip = %5.8lf\n", thisIndex, ewd);

  UgSpacePlaneProxy[thisIndex](0, 0).computeEnergies(ENERGY_EHART, ehart);
  UgSpacePlaneProxy[thisIndex](0, 0).computeEnergies(ENERGY_EEXT, eext);
  UgSpacePlaneProxy[thisIndex](0, 0).computeEnergies(ENERGY_EWD, ewd);
  delete msg;
}

void InstanceController::printEnergyEexc(CkReductionMsg *msg)
{
  double eexc = 0;
  double egga = 0;
  void *data=msg->getData();
  eexc += ((double *)data)[0];
  egga += ((double *)data)[1];
  
  CkPrintf("{%d} EEXC        = %5.8lf\n", thisIndex, eexc);
  CkPrintf("{%d} EGGA        = %5.8lf\n", thisIndex, egga);
  CkPrintf("{%d} EEXC+EGGA   = %5.8lf\n", thisIndex, eexc+egga);
      
  UgSpacePlaneProxy[thisIndex](0, 0).computeEnergies(ENERGY_EEXC, eexc);
  UgSpacePlaneProxy[thisIndex](0, 0).computeEnergies(ENERGY_EGGA, egga);

#define _GLENN_STUFF_OFF_
#ifdef _GLENN_STUFF_
  CkPrintf("exiting in printEnergyEexc\n");CkExit();
#endif

}
//============================================================================

//============================================================================
// Print out Quantum KE and put all the energies into the message group 
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
  void InstanceController::printEnergyEke(CkReductionMsg *m){
  
  double d = ((double *)m->getData())[0];
  delete m;
#ifdef _CP_DEBUG_SFNL_OFF_
  CkPrintf("ENL         = OFF FOR DEBUGGING\n");
#endif
  CkPrintf("{%d} EKE         = %5.8lf\n", thisIndex, d);
  UgSpacePlaneProxy[thisIndex](0,0).computeEnergies(ENERGY_EKE, d);
}
//============================================================================


//============================================================================
// Print out Fict CP KE and send it to the energy group
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
  void InstanceController::printFictEke(CkReductionMsg *m){
  
  double d0   = ((double *)m->getData())[0];
  double d1   = ((double *)m->getData())[1];
  double d2   = ((double *)m->getData())[2];
  double d3   = ((double *)m->getData())[3];
  double d4   = ((double *)m->getData())[4];
  delete m;

  if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt==0){
    CkPrintf("{%d} Fict Temp   =  %.10g K\n", thisIndex, d0); // per g-chare temp
    CkPrintf("{%d} Fict Eke    =  %.10g K\n", thisIndex, d2); // total kinetic energy
    CkPrintf("{%d} Fict TempNHC=  %.10g K\n", thisIndex, d1); // per g-chare tempNHC
    CkPrintf("{%d} Fict EkeNHC =  %.10g K\n", thisIndex, d3); // total NHC kinetic energy
    CkPrintf("{%d} Fict PotNHC =  %.10g K\n", thisIndex, d4); // total potNHC
    CkPrintf("{%d} Fict EConv  =  %.10g K\n", thisIndex, d2+d3+d4);
  }//endif
  UgSpacePlaneProxy[thisIndex](0,0).computeEnergies(ENERGY_FICTEKE, d0);  

}
//============================================================================


//============================================================================
//  When ALL the cp forces are done, you can integrate the atoms.
//  Could be made slicker by using knowledge of which gpsace chares
//  are on each processor and coordinating that with atom launch.
//  cklocalbranch->count_psi_forc++; zeroed in atms grp constructor
//  and incremented could be incremented in launch atoms.
//  If(cklocalbranch->count_psi_forc==numgspacethisproc){
//    atomsgrpproxy(myid).startrealspaceforces();
//  }
//  That is, you only need to know that all the forces are complete
//  on each PROCECESSOR and you can go ahead to the atoms. The atoms
//  are a group which needs to contribute atom forces to a reduction.
//  If all the forces aren't in there, you can't start the reduction
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void InstanceController::allDoneCPForces(CkReductionMsg *m){
  delete m;
  CkPrintf("All done CP forces\n");
  UatomsGrpProxy[thisIndex].StartRealspaceForces();
}
//============================================================================
