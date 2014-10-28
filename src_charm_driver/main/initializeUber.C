#ifndef _INITIALIZEUBER_C
#define _INITIALIZEUBER_C
/** \addtogroup Uber */
/**@{*/
// bump all the INT_MAPs to the right size
AtomImaptable.resize(config.numInstances);
PIBImaptable.resize(config.numInstances);
GSImaptable.resize(config.numInstances);
RSImaptable.resize(config.numInstances);
RPPImaptable.resize(config.numInstances);
RhoGSImaptable.resize(config.numInstances);
RhoRSImaptable.resize(config.numInstances);

RhoGHartImaptable.resize(config.numInstances);
RhoRHartImaptable.resize(config.numInstances);

// bump all our proxy vecs to the right size
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
UfftCacheProxy.reserve(config.numInstances);
UsfCacheProxy.reserve(config.numInstances);
UsfCompProxy.reserve(config.numInstances);
UeesCacheProxy.reserve(config.numInstances);
UplaneUsedByNLZ.reserve(config.numInstances);
UlsRhoRealProxy.reserve(config.numInstances);
UlsRhoGProxy.reserve(config.numInstances);
UavailProcs.reserve(config.numInstances);
excludePes=NULL;
#endif
