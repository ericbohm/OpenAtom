#include "debug_flags.h"
#include "charm++.h"

#ifndef ORTHO_CONFIG_H
#define ORTHO_CONFIG_H

namespace cp {
    namespace ortho {

/**
 * Configuration settings for the ortho world. Most of these are invariant
 * post-instantiation.
 */
class orthoConfig
{
    public:
        /// Is this a minimization or dynamics run
        bool isDynamics; ///< @note: This could turn into an enum as more run modes are introduced
        /// If, this is a minimization run, is it for generating the system wave functions?
        bool isGenWave; ///< @todo: Used in ortho only for shifting step number by 1. Eliminate
        //-------------------- Vars that indicate problem size / decomposition --------------------
        /// The number of states in the simulation (the dimension of the input square matrix)
        int numStates;
        /// The block size for parallelization
        int grainSize;

        //------------------- Instance, Array, Group identification etc. --------------------------
        int instanceIndex;
        /// The tolerance threshold for the S->T iterations in Ortho at which to trigger a PsiV update
        double maxTolerance;
        /// Callback to notify bubble owner (GSpace) that a tolerance update is needed
        CkCallback uponToleranceFailure;

        /// Is an OrthoHelper chare array available to perform the step 2 computations
        bool isStep2HelperOn;
        /// Pander to the BGL NIC and split the GEMMs in ortho
        int gemmSplit;



        void pup(PUP::er &p)
        {
            p|isDynamics;
            p|isGenWave;
            p|numStates;
            p|grainSize;
            p|instanceIndex;
            p|maxTolerance;
            p|uponToleranceFailure;
            p|isStep2HelperOn;
            p|gemmSplit;
        }
};

    } // end ortho
} // end namespace cp

#endif // ORTHO_CONFIG_H

