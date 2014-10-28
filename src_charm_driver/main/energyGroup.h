/** \file energyGroup.h
 */
#ifndef ENERGYGROUP_H
#define ENERGYGROUP_H
#include "uber/Uber.h"
#include "energy.h"
#include "EnergyGroup.decl.h"
#include "CPcharmParaInfo.decl.h"



/** EnergyGroup class.
 */
class EnergyGroup : public Group {

  public:
    const UberCollection thisInstance;
    EnergyStruct estruct;
    EnergyGroup(UberCollection thisInstance);
    int iteration_gsp;
    int iteration_atm;
    int kpointEnergyDoneCount;
    void updateEnergiesFromGS(EnergyStruct &, UberCollection);
    void energyDone(CkReductionMsg *);
    void energyDone();
    void sendToTemper(CkReductionMsg *);
    void resumeFromTemper();

    inline EnergyStruct getEnergyStruct(){return estruct;}
  private:
    int ktemps;
};
/*EnergyStruct GetEnergyStruct();*/
//============================================================================

#endif // ENERGYGROUP
