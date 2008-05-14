//============================================================================
/** \file ReductionClients.h
 *
 */
//============================================================================

#ifndef _ReductionClients_
#define _ReductionClients_

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/*
 * Reduction clients for the various reduction operations
 */
void doneCreatingPP(void *param, int dataSize, void *data) {
  //    fileReaderProxy.readFiles();
}
//============================================================================
	

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void printEnergyHart(void *param, int dataSize, void *data){
  //  double ehart = 0, eext = 0.0, ewd = 0.0;
  
  double ehart = ((double *)data)[0];
  double eext = ((double *)data)[1];
  double ewd  = ((double *)data)[2];
  
  CkPrintf("EHART       = %5.8lf\n", ehart);
  CkPrintf("EExt        = %5.8lf\n", eext);
  CkPrintf("EWALD_recip = %5.8lf\n", ewd);

  gSpacePlaneProxy(0, 0).computeEnergies(ENERGY_EHART, ehart);
  gSpacePlaneProxy(0, 0).computeEnergies(ENERGY_EEXT, eext);
  gSpacePlaneProxy(0, 0).computeEnergies(ENERGY_EWD, ewd);
}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void printEnergyEexc(void *param, int dataSize, void *data){
  double eexc = 0;
  double egga = 0;
  
  eexc += ((double *)data)[0];
  egga += ((double *)data)[1];
  
  CkPrintf("EEXC        = %5.8lf\n", eexc);
  CkPrintf("EGGA        = %5.8lf\n", egga);
  CkPrintf("EEXC+EGGA   = %5.8lf\n", eexc+egga);
  
  gSpacePlaneProxy(0, 0).computeEnergies(ENERGY_EEXC, eexc);
  gSpacePlaneProxy(0, 0).computeEnergies(ENERGY_EGGA, egga);

#define _GLENN_STUFF_OFF_
#ifdef _GLENN_STUFF_
  CkPrintf("exiting in reductclients.h\n");CkExit();
#endif

}
//============================================================================




#endif
