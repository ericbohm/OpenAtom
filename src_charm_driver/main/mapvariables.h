#ifndef _MAPVARIABLES_H
#define _MAPVARIABLES_H

/** \file mapvariables.h 
    \brief collection point for global scope map, proxy, and Uber read-only
    variables
*/

//============================================================================
/**
 * @defgroup proxy_vars proxy_vars
 * Defining all the Charm++ readonly variables, which include proxies
 * 
 * to access the arrays and groups and the Communication Library
 * handles.
 */
//============================================================================
/**@{*/
bool firstInstance = true;  
int numInst=0;
// INT_MAPs are the ones actually used so
// these are being changed to CkVec's for beads
CkVec <int> PIBImaptable;
CkVec <MapType1> AtomImaptable;
CkVec <MapType2> GSImaptable;
CkVec <MapType2> RSImaptable;
CkVec <MapType2> RPPImaptable;
CkVec <MapType2> RhoGSImaptable;
CkVec <MapType2> RhoRSImaptable;
CkVec <MapType2> RhoGHartImaptable;
CkVec <MapType3> RhoRHartImaptable;
//CkVec <MapType1> LSPRhoGSImaptable;
CkVec <MapType2> LSPRhoRSImaptable;

#ifndef USE_INT_MAP
CkHashtableT<intdual, int> GSmaptable(10000,0.25);
CkHashtableT<intdual, int> RSmaptable(10000,0.25);
CkHashtableT<intdual, int> RPPmaptable(10000,0.25);
CkHashtableT<intdual, int> AsymScalcmaptable(10000,0.25);
CkHashtableT<intdual, int> SymScalcmaptable(10000,0.25);
#else
CkHashtableT<intdual, int> GSmaptable;
CkHashtableT<intdual, int> RSmaptable;
CkHashtableT<intdual, int> RPPmaptable;
CkHashtableT<intdual, int> AsymScalcmaptable;
CkHashtableT<intdual, int> SymScalcmaptable;
#endif
CkHashtableT<intdual, int> RhoGSmaptable;
CkHashtableT<intdual, int> RhoRSmaptable;
CkHashtableT<intdual, int> RhoGHartmaptable;
CkHashtableT<inttriple, int> RhoRHartmaptable;
CkHashtableT<intdual, int> Orthomaptable;
CkHashtableT<intdual, int> OrthoHelpermaptable;


CProxy_main                       mainProxy;
CProxy_PhysScratchCache           pScratchProxy;
Config                            config;
CProxy_TimeKeeper                 TimeKeeperProxy;
CProxy_InstanceController         instControllerProxy;
CProxy_TemperController         temperControllerProxy;
CProxy_ENL_EKE_Collector          ENLEKECollectorProxy;
CPcharmParaInfo simReadOnly;
/**@}*/


//============================================================================
/** @defgroup Uber Uber
 * \brief Ubers provide a multidimensional collection of CkArray proxies such that a complete instance of all objects necessary for a simulation are accessible at each unique tuple of indices. 
 * 
 *  Uber proxies for all the things which change per step 
 *  
 *  Indexed by PathIntegral Bead.  Each Bead has its own set of
 *  proxies.  Charm driver startup will construct a different set of
 *  arrays for each bead.  
 *
 *  There is a small flexibility vs performance tradeoff.  If we
 *  assume beads are never co-mapped, then the old readonly can be
 *  overwritten locally on each processor
 *  (e.g. gSpacePlaneProxy=UgSpacePlaneProxy[mybead];)
 *  Otherwise we force a slight indirection penalty to lookup 
 *  U*Proxy[mybead] for every send.  In practice this is probably
 *  noise compared to the real expense of sending a message, but it
 *  does seem a little silly for the default case where there is only
 *  one bead.
 */


/** \addtogroup Uber */
/**@{*/
CkVec <CProxy_PIBeadAtoms>       UPIBeadAtomsProxy;
CkVec <CProxy_CP_State_GSpacePlane>       UgSpacePlaneProxy;
CkVec <CProxy_GSpaceDriver>               UgSpaceDriverProxy;
CkVec <CProxy_CP_State_ParticlePlane>     UparticlePlaneProxy;
CkVec <CProxy_CP_State_RealParticlePlane> UrealParticlePlaneProxy;
CkVec <CProxy_CP_State_RealSpacePlane>    UrealSpacePlaneProxy;
CkVec <CProxy_CP_Rho_RealSpacePlane>      UrhoRealProxy;
CkVec <CProxy_CP_Rho_GSpacePlane>         UrhoGProxy;
CkVec <CProxy_CP_Rho_RHartExt>            UrhoRHartExtProxy;
CkVec <CProxy_CP_Rho_GHartExt>            UrhoGHartExtProxy;
CkVec <CProxy_AtomsCompute>               UatomsComputeProxy;
CkVec <CProxy_AtomsCache>                 UatomsCacheProxy;
CkVec <CProxy_EnergyGroup>                UegroupProxy;
CkVec <CProxy_FFTcache>                   UfftCacheProxy;
CkVec <CProxy_StructFactCache>            UsfCacheProxy;
CkVec <CProxy_StructureFactor>            UsfCompProxy;
CkVec <CProxy_eesCache>                   UeesCacheProxy;
CkVec <CProxy_CP_LargeSP_RhoGSpacePlane>      UlsRhoGProxy;
CkVec <CProxy_CP_LargeSP_RhoRealSpacePlane>      UlsRhoRealProxy;

CkVec <UberCollection>			  UberAlles;
CkVec < PeList * >                        UavailProcs;

/**@}*/


CkVec < CkVec <int> > UpeUsedBySF;
CkVec < CkVec <int> > UpeUsedByNLZ;
CkVec < CkVec <int> > UplaneUsedByNLZ;

PeList *availGlobG=NULL;
PeList *availGlobR=NULL;
PeList *excludePes=NULL;
int boxSize;
TopoManager *topoMgr=NULL;
inttriple *mapOffsets=NULL;


#endif
