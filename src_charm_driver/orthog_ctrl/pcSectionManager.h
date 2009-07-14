#include "paircalc/ckPairCalculator.h"
#include "configure.h"

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
        void init(const CkIndex2D orthoIdx, const PairCalcID &pcid, const Config &cfg, CkGroupID oMCastGID, CkGroupID oRedGID);
        /// PUP serializer
        void pup(PUP::er &p);
        /// Creates a paircalc array section given the necessary info. Replaces initOneRedSect()
        void setupArraySection(CkCallback cb, CkCallback synccb, bool arePhantomsOn, bool useComlibForOrthoToPC);
        /// Sends out the results to the paircalc section. Replaces finishPairCalcSection()
        void sendResults(int n, double *ptr1, double *ptr2, int orthoX, int orthoY, int actionType, int priority);
        /// Used to send OrthoT to the asymm instance. Replaces sendMatrix()
        void sendMatrix(int n, double *ptr1, double *ptr2, int orthoX, int orthoY, int actionType, int priority);
        /// Identify the state indices of the Paircalc chares this ortho chare needs to talk to
        CkIndex2D computePCStateIndices(const int orthoX, const int orthoY);

    private:
        /// Create a paircalc section containing all chares with the specified two state indices
        void createPCsection(const int s1, const int s2);

        /// Number of planes that GSpace is decomposed into
        int numPlanes;
        /// Number of states in this simulation
        int numStates;
        /// The number of chunks (4th dimension decomposition) of paircalcs
        int numChunks;
        /// The statewise decomposition grain size for the paircalcs
        int pcGrainSize;
        /// The statewise decomposition grain size for the ortho chares
        int orthoGrainSize;

        /// The array ID of the paircalc instance that I will manage comm with
        CkArrayID pcArrayID;
        /// Is this paircalc array a symmetric or asymmetric instance
        bool isSymmetric;
        /// The section of the array that my owner ortho chare will be talking to
        CProxySection_PairCalculator pcSection;

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

