#include "paircalc/ckPairCalculator.h"

#ifndef PC_SECTION_MANAGER_H
#define PC_SECTION_MANAGER_H

/// Forward declarations
class Ortho; 
class PairCalcID;


namespace cp {
    namespace ortho {

/// Class that manages the paircalc sections that each Ortho chare communicates with
class PCSectionManager
{
    friend class ::Ortho;

    public:
        /// Constructor
        PCSectionManager() {}
        /// An initializer method that fills this with data
        void init(const CkIndex2D orthoIdx, const PairCalcID &pcid, const int orthoGrainSize, CkGroupID oMCastGID, CkGroupID oRedGID);
        /// PUP serializer
        void pup(PUP::er &p);
        /// Creates a paircalc array section given the necessary info. Replaces initOneRedSect()
        void setupArraySection(int numZ, int* z, CkCallback cb, CkCallback synccb, int s1, int s2, bool arePhantomsOn, bool direct, bool commlib);
        /// Sends out the results to the paircalc section. Replaces finishPairCalcSection()
        void sendResults(int n, double *ptr1, double *ptr2, int orthoX, int orthoY, int actionType, int priority);
        /// Used to send OrthoT to the asymm instance. Replaces sendMatrix()
        void sendMatrix(int n, double *ptr1, double *ptr2, int orthoX, int orthoY, int actionType, int priority);

    private:
        /// Paircalc section proxy
        CProxySection_PairCalculator pcSection;

        /// Number of states in this simulation
        int numStates;
        /// The statewise decomposition grain size for the paircalcs
        int pcGrainSize;
        /// The statewise decomposition grain size for the ortho chares
        int orthoGrainSize;
        /// The number of chunks (4th dimension decomposition) of paircalcs
        int numChunks;

        /// Paircalc Array ID
        CkArrayID pcArrayID;
        /// Is this paircalc array a symmetric or asymmetric instance
        bool isSymmetric;
        /// If symmetric, does this paircalc array have phantoms turned on?
        bool arePhantomsOn;

        /// The index of the calling Ortho chare
        CkIndex2D orthoIndex;
        /// The multicast and reduction groups that handle comm
        CkGroupID orthomCastGrpID, orthoRedGrpID;
        /// The priority to use for messages to PC
        int msgPriority;
};

    } // end namespace ortho
} // end namespace cp

#endif // PC_SECTION_MANAGER

