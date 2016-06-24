/** \file atomsCache.h
 */
#ifndef ATOMSCACHE_H
#define ATOMSCACHE_H
#include "Atoms.h"
struct EnergyStruct;
#include "Atoms.decl.h"
#include "atomMessages.h"
#include "CPcharmParaInfo.decl.h"


#include "uber/Uber.h"


/** \file AtomsCache.h 
 * @defgroup Atoms Atoms
 *
 * \brief Handles coordinate and force integration for Atoms, keeping a distributed coordinate cache in AtomsCache and computing updates in AtomsCompute 
 *
 *
 * AtomsCache is a mostly passive structure for handling data
 * replication of the atom coordinates and write-back collation of
 * force contributions.  
 *
 * various chares use CkLocal to access the atoms:
 * Rho_RHart RhoGHart StructureFactor RealParticlePlane
 * ParticlePlane 
 * -> fastAtoms read coords write forces
 * RhoRHart
 * -> natm read
 * -> fastAtoms read coords write forces
 * -> iteration read
 * RhoGHart
 * -> natm read
 * -> iteration read coords write forces
 * -> fastAtoms read
 * RealParticlePlane
 * -> temperScreenFile* read
 * -> iteration read
 * -> fastAtoms read coords write forces
 * StructureFactor
 * -> iteration read
 * -> fastAtoms read coords write forces
 * InstanceController
 * -> iteration read
 * -> temperScreenFile* read
 * GSpacePlane
 * -> natm read
 * -> fastAtoms read coords for genwave
 * -> temperScreenFile* read

 * EJB: this could perhaps be better implemented in an MSA, once they
 * shake performance bugs out of MSA, and pay attention to
 * partitioning considerations, we'll revisit it.  


 CkCache could be considered here, but our use case is so simple that
 CkCache is wildly over engineered for our purposes.  The number of
 atoms is always relatively small, so we have no eviction
 requirements.  Similarly the number of clients is static and all
 clients need access to all the atom coordinates every iteration, so
 there is no segmentation and no window of opportunity for clever
 just in time delivery.  

 Nodegroup?  AtomsCache is an excellent candidate for Nodegroup.  All
 force updates are a purely additive operation, each computation will
 use the coordinates and its own chare local data to produce its
 force contribution.  Force updates are handled via by reference
 passage into the PINY compute routines.  Updates are order
 independent, so as long as the memory subsystem doesn't lose its
 mind maintaining associative additive consistency (almost the
 weakest consistency one could ask for in a shared memory context)
 then correctness is guaranteed without locks or explicit fencing.
 False sharing shouldn't be an issue given that each force component
 is in its own vector as are the coordinates.  Having only one of
 these per node gives us a slight memory footprint reduction.  More
 importantly it requires numprocs/node fewer update messages.  If you
 think this gives us grief just replace nodegroup with Group and
 recompile, the nodegroup advantages are all implicit.

 However, AtomsCache directly uses the registration in the eesCache to
 execute releaseGSP.  So, AtomsCache's nodegroupness needs to be 1:1
 with eesCache's nodegroupness, or an alternate launch scheme put
 in place for releaseGSP to break that dependency.
 * @addtogroup Atoms
* @{ 
  */

    class AtomsCache: public Group {
      public:
        const UberCollection thisInstance;
        int natm;
        int natm_nl;
        int iteration;
        FastAtoms fastAtoms;
        FILE *temperScreenFile;
	int pimdchaincount;
        std::string output_directory;

        AtomsCache(int,int,Atom *,UberCollection thisInstance, std::string output_directory);
        AtomsCache(CkMigrateMessage *m) {}
        ~AtomsCache();
        void contributeforces(ContribForcesMsg *);
        void contributeforcesSectBcast();
        void atomsDone(AtomsDoneMsg *);
        void atomsDoneSectBcast(CkReductionMsg *);
        void atomsDoneSectBcast();
        void acceptAtoms(AtomMsg *);  // entry method
        void acceptAtomsSectBCast(AtomMsg *);  // entry method
        void releaseGSP();
	void acceptChainContribution(double PIMDChain);
	void secDone(CkReductionMsg *m){CkPrintf("secdone\n");}
        CProxySection_AtomsCache secProxy;
	CkSectionInfo forcecookie;
	void initForceCookie(ContribForcesMsg *);
	void createSpanningSection();
        void zeroforces() {
          double *fx = fastAtoms.fx;
          double *fy = fastAtoms.fy;
          double *fz = fastAtoms.fz;
          for(int i=0; i<natm; i++){
            fx[i]       = 0;
            fy[i]       = 0;
            fz[i]       = 0;
          }//endfor
        }//end routine
    };

  /*@}*/
#endif // ATOMSCACHE_H
