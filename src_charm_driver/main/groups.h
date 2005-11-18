//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file groups.h
 *
 *
 */
//============================================================================

#ifndef _groups_h_
#define _groups_h_

#include "../../include/Atoms.h"

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \brief AtomsGrp class.
 *
 *
 *
 */
//============================================================================

class AtomsGrp: public Group {
 public:
  AtomsGrp(CkMigrateMessage *m) {}
  AtomsGrp(int n, int n_nl,Atom *);
  ~AtomsGrp();
  double pot_ewd_rs;      // total real space ewald energy
  double pot_ewd_rs_loc;  // local real space ewald energy
  int natm;
  int natm_nl;
  Atom *atoms;

  void positiveDisplace(int localState);
  void negativeDisplace(int localState);
  void contributeforces(double pot_ewald);
  void recvContribute(CkReductionMsg *);
  void StartRealspaceForces();
  void zeroforces() {
    for(int i=0; i<natm; i++){
      atoms[i].fx = 0;
      atoms[i].fy = 0;
      atoms[i].fz = 0;
    }
  }

};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/* Structure that stores all energies */

struct EnergyStruct {
    //non local
    double enl;

    //local external energy
    double eext;
    
    double eke;
    //hartree energy
    double ehart;

    //ion-ion
    double eewald_recip;
    double eewald_real;

    //exchange correlation
    double egga;
    double eexc;

    double totalEnergy;
    double fmagPsi;

    double fictEke;

};
PUPbytes(EnergyStruct);
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \brief EnergyGroup class.
 *
 *
 *
 */
//============================================================================

class EnergyGroup : public Group {
    
 public:
    EnergyStruct estruct;
    EnergyGroup();
    
    void updateEnergies(EnergyStruct &es) {
      estruct.enl          = es.enl;
      estruct.eke          = es.eke;
      estruct.eext         = es.eext;
      estruct.ehart        = es.ehart;
      estruct.eewald_recip = es.eewald_recip;
      estruct.egga         = es.egga;
      estruct.eexc         = es.eexc;
      estruct.fictEke      = es.fictEke;
      estruct.totalEnergy  = es.totalEnergy;
      estruct.fmagPsi      = es.fmagPsi;
#ifdef _DEBUG_ESTRUCT_
       CkPrintf("Energies received %lf, %lf, %lf, %lf, %lf\n", 
                 estruct.enl,estruct.eke,estruct.eext,estruct.ehart, 
                 estruct.eewald_recipt,estruct.egga,estruct.eexc,
                 estruct.fictEke,estruct.totalEnergy);
#endif
    }

    inline EnergyStruct getEnergyStruct(){return estruct;}
};
EnergyStruct GetEnergyStruct();
//============================================================================

#endif // #ifndef _groups_h_
