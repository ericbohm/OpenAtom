/** \file energyGroup.h
 */
#ifndef ENERGYGROUP_H
#define ENERGYGROUP_H
#include "uber/Uber.h"
#include "energy.h"
#include "EnergyGroupM.decl.h"
#include "CPcharmParaInfo.decl.h"

/** EnergyGroup class.
 */
class EnergyGroup : public CBase_EnergyGroup {

  public:
    const UberCollection thisInstance;
    EnergyStruct estruct;
    EnergyGroup(){};
    EnergyGroup(UberCollection thisInstance);
    EnergyGroup(CkMigrateMessage* msg) {}		
    int iteration_gsp;
    int iteration_atm;
    int kpointEnergyDoneCount;
    void updateEnergiesFromGS(EnergyStructMsg *);
    void updateEnergiesFromGSSectBcast(EnergyStructMsg *);
    void energyDone(CkReductionMsg *m);
    void energyDone();
    void sendToTemper(CkReductionMsg *m);
    void sectionDone(CkReductionMsg *m);
    void resumeFromTemperSectBcast();
    void createSpanningSection();
    inline EnergyStruct getEnergyStruct(){return estruct;}
    void initTemperCookie(ECookieMsg *msg);
    void resumeFromTemper(ECookieMsg *msg);

  private:
    CProxySection_EnergyGroup secProxy;
    int ktemps;
    int countOfEnergies;
    CkSectionInfo temperCookie;
};
/*EnergyStruct GetEnergyStruct();*/
//============================================================================

#endif // ENERGYGROUP
