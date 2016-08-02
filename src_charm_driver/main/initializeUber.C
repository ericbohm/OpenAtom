#ifndef _INITIALIZEUBER_C
#define _INITIALIZEUBER_C
/** \addtogroup Uber */
/**@{*/
// bump all the INT_MAPs to the right size
AtomImaptable.resize(config.numInstances);
//EnergyCommMgrImaptable.resize(config.numInstances);
PIBImaptable.resize(config.numInstances);
GSImaptable.resize(config.numInstances);
RSImaptable.resize(config.numInstances);
RPPImaptable.resize(config.numInstances);
RhoGSImaptable.resize(config.numInstances);
RhoRSImaptable.resize(config.numInstances);
RhoGHartImaptable.resize(config.numInstances);
RhoRHartImaptable.resize(config.numInstances);
RhoYPencilImaptable.resize(3);
for(int pencil = 0; pencil < 3; pencil++) {
  RhoYPencilImaptable[pencil].resize(config.numInstances);
}
RhoHartYPencilImaptable.resize(config.numInstances);
if(simReadOnly.ees_eext_on) {
  AtmSFYPencilImaptable.resize(config.nchareHartAtmT + 1);
  for(int type = 0; type < config.nchareHartAtmT + 1; type++) {
    AtmSFYPencilImaptable[type].resize(config.numInstances);
  }
}


UberPes.resize(config.numInstances);


for(int inst; inst < config.numInstances; inst++)
  {
    UberPes[inst].reserve(config.numPesPerInstance);
  }

// bump all our proxy vecs to the right size

UPIBeadAtomsProxy.reserve(config.numInstances);
UpScratchProxy.reserve(config.numInstances);
UgSpacePlaneProxy.reserve(config.numInstances);
UgSpaceDriverProxy.reserve(config.numInstances);
UparticlePlaneProxy.reserve(config.numInstances);
UrealParticlePlaneProxy.reserve(config.numInstances);
UrealSpacePlaneProxy.reserve(config.numInstances);
UrhoRealProxy.reserve(config.numInstances);
UrhoGProxy.reserve(config.numInstances);
UrhoRHartExtProxy.reserve(config.numInstances);
UrhoGHartExtProxy.reserve(config.numInstances);
UatomsComputeProxy.reserve(config.numInstances);
UatomsCacheProxy.reserve(config.numInstances);
UegroupProxy.reserve(config.numInstances);
//UeCommProxy.reserve(config.numInstances);
UfftCacheProxy.reserve(config.numInstances);
UsfCacheProxy.reserve(config.numInstances);
UsfCompProxy.reserve(config.numInstances);
UeesCacheProxy.reserve(config.numInstances);
UplaneUsedByNLZ.reserve(config.numInstances);
UlsRhoRealProxy.reserve(config.numInstances);
UlsRhoGProxy.reserve(config.numInstances);
UavailProcs.reserve(config.numInstances);
excludePes=NULL;
Urho_fft_xProxy.reserve(config.numInstances);
Urho_fft_yProxy.reserve(config.numInstances);
Urho_fft_zProxy.reserve(config.numInstances);
Urho_fft_hartProxy.reserve(config.numInstances);
if(simReadOnly.ees_eext_on) {
  Urho_fft_atmSFProxy.resize(config.nchareHartAtmT + 1);
  for(int type = 0; type < config.nchareHartAtmT + 1; type++) {
    Urho_fft_atmSFProxy[type].reserve(config.numInstances);
  }
  if(config.nchareHartAtmT > 1) {
    Urhart_sectionProxy.resize(config.nchareHartAtmT);
    Ughart_sectionProxy.resize(config.nchareHartAtmT);
    for(int type = 0; type < config.nchareHartAtmT; type++) {
      Urhart_sectionProxy[type].reserve(config.numInstances);
      Ughart_sectionProxy[type].reserve(config.numInstances);
    }
  }
}
#endif
