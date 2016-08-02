#ifndef _INSTANCECONTROLLER_H_
#define _INSTANCECONTROLLER_H_

class ICCookieMsg : public CkMcastBaseMsg, public CMessage_ICCookieMsg {
  public:
    int junk;
};

/****************************************************************************
 * Handles launch synchronization for each instance.
 * 
 * Catch all for per instance synchronization of any kind.  
 *
 ****************************************************************************/
class InstanceController: public CBase_InstanceController {
  public:
    InstanceController_SDAG_CODE
    double Timer;
    InstanceController(int);
    ~InstanceController(){}
    InstanceController(CkMigrateMessage *m){}
    void init();
    void doneInit();
    void doneFFTCreation(idMsg *msg);
    void initDensity();
    void initCookie(ICCookieMsg *msg);
    void printEnergyHart(CkReductionMsg *msg);
    void printEnergyEexc(CkReductionMsg *msg);
    void printFictEke(CkReductionMsg *msg);
    void printEnergyEke(CkReductionMsg *m);
    void allDoneCPForces(int);
    void allDoneCPForcesAllKPoint(int);
    void cleanExit(CkReductionMsg *m);
    void cleanExitAll(CkReductionMsg *m);
    void acceptNewTemperature(double temp);
    void useNewTemperature(double temp);
    void atomsDoneNewTemp(CkReductionMsg *m);
    void gspDoneNewTemp(CkReductionMsg *m);
    void fmagMinTest(CkReductionMsg *m);
    void instancesReady(CkReductionMsg *msg);
    void resumeFromTemper();
    CProxySection_CP_State_GSpacePlane gTemperBeadProxy;
    void doneIteration();
    void allInstancesDoneIteration();
  private:
    int done_init, fft_expected;
    bool done_fft_creation;
    int numKpointforces;
    bool atomsTempDone;
    bool gspTempDone;
    bool printToScreen;
    CkSectionInfo allKPcookie;
};

class DiagonalizerBridge : public CBase_DiagonalizerBridge {
  public:
    DiagonalizerBridge_SDAG_CODE
    DiagonalizerBridge();
    void sendLambdaToDiagonalizer(int x, int y, int n, internalType *lmat);
    void prepareDiagonalizerInput(int x, int y, int n, internalType *lmat);
    void sendLambdaBackToOrtho();
    void integrateLambda(int n, internalType *lmat);
    int x, y;
};

#endif
