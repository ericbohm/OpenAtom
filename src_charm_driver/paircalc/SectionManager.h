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
        /// PUP serializer
        void pup(PUP::er &p);
        /// Creates a paircalc array section given the necessary info. Replaces initOneRedSect()
        void setupArraySection(int numZ, int* z, int numChunks,  PairCalcID* pcid, CkCallback cb, CkCallback synccb, int s1, int s2, int orthoX, int orthoY, int orthoGrainSize, bool phantom, bool direct, bool commlib);
        /// Sends out the results to the paircalc section. Replaces finishPairCalcSection()
        void sendResults(int n, double *ptr1, double *ptr2, PairCalcID *pcid, int orthoX, int orthoY, int actionType, int priority);
        /// Used to send OrthoT to the asymm instance. Replaces sendMatrix()
        void sendMatrix(int n, double *ptr1, double *ptr2, PairCalcID *pcid, int orthoX, int orthoY, int actionType, int priority);

    private:
        /// Paircalc section proxy
        CProxySection_PairCalculator pcSection;
};

    } // end namespace paircalc
} // end namespace cp

#endif // PC_SECTION_MANAGER

