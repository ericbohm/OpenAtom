#include "paircalc/ckPairCalculator.h"

#ifndef PC_SECTION_MANAGER_H
#define PC_SECTION_MANAGER_H

/// Forward declarations
class Ortho; 


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
        void init(const CkIndex2D orthoIdx, const pc::pcConfig &pcCfg, CkArrayID pcAID, CkGroupID oMCastGID, CkGroupID oRedGID);
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
        /// Overloaded version that uses the stored ortho indices to compute the PC state indices
        inline CkIndex2D computePCStateIndices() { return computePCStateIndices(orthoIndex.x,orthoIndex.y); }

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



/**
 * ortho and paircalc grainsizes do not complicate this discussion a whole lot because of the restriction that 
 * ortho grainsize = multiple of paircalc grain size. Because of this equal or exact multiple clause, ortho grains
 * will line up perfectly inside a paircalc grain and, hence, every ortho chare will hold a bunch of states that 
 * will all get delivered to the same paircalc section.
 *
 * paircalcs on the other hand will have to chop up their data along the ortho tile boundaries and contribute to 
 * multiple reductions that end up at the respective ortho chares. Refer PairCalculator::contributeSubTiles.
 *
 */
inline CkIndex2D PCSectionManager::computePCStateIndices(const int orthoX, const int orthoY)
{
    CkIndex2D pc;
    pc.x = orthoX * orthoGrainSize;
    pc.y = orthoY * orthoGrainSize;
    // Do something clever if the grainsizes are not the same
    if(orthoGrainSize != pcGrainSize)
    {
        int maxpcstateindex = (numStates/pcGrainSize - 1) * pcGrainSize;
        pc.x = pc.x / pcGrainSize * pcGrainSize;
        pc.y = pc.y / pcGrainSize * pcGrainSize;
        pc.x = (pc.x>maxpcstateindex) ? maxpcstateindex :pc.x;
        pc.y = (pc.y>maxpcstateindex) ? maxpcstateindex :pc.y;
    }
    return pc;
}

    } // end namespace ortho
} // end namespace cp

#endif // PC_SECTION_MANAGER

