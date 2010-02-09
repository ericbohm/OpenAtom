#ifndef PC_CONFIG_H
#define PC_CONFIG_H

namespace cp {
    namespace paircalc {

/**
 * Dumb structure that holds all the configuration inputs required for paircalc 
 * instantiation, functioning and interaction.
 */
class pcConfig
{
    public:
        /// Is this a minimization or dynamics run
        bool isDynamics; ///< @note: This could turn into an enum as more run modes are introduced
        /// Is this a symmetric or asymmetric paircalc instance
        bool isSymmetric;
        /// If this is a symmetric instance, should it use phantom chares to balance the BW path
        bool arePhantomsOn; ///< @note: Originally, phantomSym

        /// The total number of planes in the system
        int numPlanes; ///< @warning: Not used in the original code. WHY?
        /// The total number of states in the system
        int numStates;
        /// The number of chunks (4th dimension of decomposition)
        int numChunks; ///< @note: Originally, blkSize

        /// The grain size along the states dimensions (plural) (number of states per PC chare)
        int grainSize;
        /** The grain size along the states dimensions for Ortho chares
         *
         * ie, the number of states per ortho chare. Ortho chares along the right and bottom
         * chare array boundaries may have to deal with remainders and may hence 
         * have a larger orthoGrainSize
         */
        int orthoGrainSize;

        /// The array ID of the GSpace chare array this instance talks to
        CkArrayID gSpaceAID; ///< @note: Originally, final_callbackid
        /// The entry point to which this instance should send results to
        int gSpaceEP; ///< @note: Originally,  final_callback_ep
        /// The entry point to which this instance should send PsiV tolerance update results to
        int PsiVEP; ///< @note: Originally,  callback_ep_tol

        /** The mem footprint vs performance setting for the paircalcs (tribool)
         *
         * -1: dont conserve memory at all
         *  0: default (balance mem usage with performance)
         * +1: conserve memory aggresively and free up matrices when not in use
         */
        int conserveMemory;
        /// Should the paircalcs worry about load balancing
        bool isLBon; ///< @note: Originally,  lbpaircalc
        /// What mechanism should be used for reductions @note: defunct
        redtypes reduce;
        /// Is double-packing on?
        bool isDoublePackOn;

        /// Should this instance collect result fragments and send them out together in the BW path?
        bool areBWTilesCollected; ///< @note: Originally, collectAllTiles
        /// Should this instance stream the result fragments in the BW path as they become ready?
        bool isBWstreaming; ///< @note: Originally, PCstreamBWout
        /// Should we impose a hard barrier in the BW path to sync all PC chares?
        bool isBWbarriered; ///< @note: Originally, useBWBarrier
        /// Should we tweak msg priority to delay msgs to GSpace carrying the results?
        bool shouldDelayBWsend; ///< @note: Originally, delaybw
        /// Will the input data be multicast to PC sections or sent directly (p2p)
        bool isInputMulticast; ///< @note: Originally, usePairDirectSend
        /// Should the results from each PC chare be reduced or delivered individually to GSpace?
        bool isOutputReduced; ///< @note: Originally, gSpaceSum
        /// If shouldDelayBWsend, what priority should this instance use for the result msgs
        int resultMsgPriority; ///< @note: Originally,  gpriority
        /// The proxyOffset value of thisInstance of OpenAtom computations
        int instance;

        /**{
         * BGL's painful NIC forces us to split long computations. Configure GEMM splitting here 
         * @note: PURELY for the BGL
         */
        int gemmSplitFWk;
        int gemmSplitFWm;
        int gemmSplitBW;
        ///}
};

    } // end namespace paircalc
} // end namespace cp

/// PUP this class as raw bytes (for now)
PUPbytes(cp::paircalc::pcConfig);

#endif // PC_CONFIG_H
