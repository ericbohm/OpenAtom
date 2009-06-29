#ifndef PC_SECTION_MANAGER_H
#define PC_SECTION_MANAGER_H

class PairCalcID;

namespace cp {
    namespace paircalc {

/// Class that manages the paircalc sections that each Ortho chare communicates with
class SectionManager
{
    public:
        /// Constructor
        //SectionManager(const CkIndex2D idx, const int grSize);
        SectionManager() {}
        /// initOneRedSect
        void setupArraySection(int numZ, int* z, int numChunks,  PairCalcID* pcid, CkCallback cb, CkCallback synccb, int s1, int s2, int orthoX, int orthoY, int orthoGrainSize, bool phantom, bool direct, bool commlib);
        /// setGredProxy
        /// finishPairCalcSection
        /// sendMatrix

    public:
        /// A paircalc section proxy
        CProxySection_PairCalculator pcSection;
        /// The index of the calling Ortho chare
        //CkIndex2D orthoIndex;
        /// The ortho grainsize
        //int orthoGrainSize;
};

    } // end namespace paircalc
} // end namespace cp

#endif // PC_SECTION_MANAGER

