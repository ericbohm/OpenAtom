#include "debug_flags.h"
#include "charm++.h" /// @note: Needed for ckcomplex.h ?
#include "ckcomplex.h"
#include "gParticlePlane.decl.h"
#include "uber/Uber.h"
struct EnergyStruct; /// @warning: Forward declarations of structs seem to choke the ci parser. It doesnt recognize the keyword struct.
#include "paircalc/ckPairCalculator.h"
#include "cpaimd.decl.h"

#ifndef CP_STATE_PARTICLE_PLANE_H
#define CP_STATE_PARTICLE_PLANE_H


class NLDummyMsg: public CMessage_NLDummyMsg 
{
    public:
        int iteration;
};




class EnlCookieMsg : public CkMcastBaseMsg, public CMessage_EnlCookieMsg 
{
    public:
        int foo;
};




class NLFFTMsg: public CMessage_NLFFTMsg 
{
    public:
        int size;
        int senderIndex;
        int step;
        complex *data;
};




class GSPPIFFTMsg: public CMessage_GSPPIFFTMsg 
{
    public:
        int size;
        int offset;
        int iterNL;
        complex *data;
};




class CP_State_ParticlePlane: public CBase_CP_State_ParticlePlane
{
    friend class GSpaceDriver;
    friend class CP_State_GSpacePlane;
    public:
        CP_State_ParticlePlane(CkMigrateMessage *m) {}
        CP_State_ParticlePlane(int ,int ,int ,int ,int ,int ,int ,int ,int ,int ,int ,int ,int ,int ,int, UberCollection );
        ~CP_State_ParticlePlane();
        void pup(PUP::er &);
        /// InitKVectors creates local k vectors for this chare and mallocs some memory.
        void initKVectors();
        void startNLEes(int);
        void lPrioStartNLEes(NLDummyMsg *m);
        /// @entry GSpaceDriver calls this to trigger the launch of the computeZ()s to compute the Z matrix.
        void launchComputeZs();
        /// @entry
        void computeZ(PPDummyMsg *m);
        void setEnlCookie(EnlCookieMsg *m);
        void ResumeFromSync();
        void reduceZ(int, int, complex *,complex *,complex *,complex *);
        void getForces(int, int, complex *);

        void createNLEesFFTdata();
        void FFTNLEesFwd();
        void sendToEesRPP();
        void recvFromEesRPP(GSPPIFFTMsg *msg);
        void FFTNLEesBck();
        void computeNLEesForces();
        void registrationDone(CkReductionMsg *msg);
        void printEnl(CkReductionMsg *msg);
        int myChareG;
        int iteration;
        int iterNL;
        int numNLiter;
        int ees_nonlocal;
        int ngridaNL;
        int ngridbNL;
        int ngridcNL;
        int gSpaceNumPoints;
        int numLines;
        int numFullNL;
        int natm_nl;
        int natm_nl_grp_max;
        int numSfGrps;
        int nstates;
        int nchareG;
        int Gstates_per_pe;
        int countNLIFFT;
        int sendDone;
        int registrationFlag;
    private:
        const UberCollection thisInstance;
        int calcReductionPlaneNum(int);
        complex *myForces, *gspace, *projPsiG;
        complex *zmatrixSum, *zmatrix;
        double *dyp_re,*dyp_im;
        double enl;
        double enl_total;
        double totalEnergy;
        int *haveSFAtmGrp;
        int *count;
        int doneEnl;
        int doneForces;
        int zsize, energy_count;
        int sizeX, sizeY, sizeZ, gSpacePlanesPerChare;
        int reductionPlaneNum;
        complex *zmatrix_fx,*zmatrix_fy,*zmatrix_fz;
        complex *zmatrixSum_fx,*zmatrixSum_fy,*zmatrixSum_fz;
        CkSectionInfo enlCookie;
        CProxySection_CP_State_ParticlePlane particlePlaneENLProxy;
        CProxy_CP_State_RealParticlePlane realPP_proxy;
        #ifdef _CP_GS_DEBUG_COMPARE_VKS_
            complex *savedprojpsiBf;
            complex *savedprojpsiBfsend;
            complex *savedprojpsiGBf;
        #endif
};


#endif // CP_STATE_PARTICLE_PLANE_H
