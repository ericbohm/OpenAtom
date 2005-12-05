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
  int natm;
  int natm_nl;
  int len_nhc;
  int iextended_on;
  int cp_min_opt;
  int cp_wave_opt;
  int iteration;
  double kT;
  double pot_ewd_rs;      // total real space ewald energy
  double pot_ewd_rs_loc;  // local real space ewald energy
  double eKinetic;        // kinetic energy
  double eKineticNhc;     // NHC kinetic energy
  double potNhc;         // NHC potential energy
  Atom *atoms;
  AtomNHC *atomsNHC;

  AtomsGrp(CkMigrateMessage *m) {}
  AtomsGrp(int,int,int, int ,int ,int,double ,Atom *,AtomNHC *);
  ~AtomsGrp();
  void contributeforces(double pot_ewald);
  void recvContribute(CkReductionMsg *);
  void atomsDone(CkReductionMsg *);
  void StartRealspaceForces();
  void outputAtmEnergy();
  void zeroforces() {
    for(int i=0; i<natm; i++){
      atoms[i].fx = 0;
      atoms[i].fy = 0;
      atoms[i].fz = 0;
    }//endfor
  }//end routine
  void zeronhc(){for(int i=0;i<natm;i++){atomsNHC[i].posKT = 0;}}
  void computeEkin(){
    eKinetic = 0.0;
    for(int i=0;i<natm;i++){
      eKinetic += atoms[i].m*(atoms[i].vx*atoms[i].vx+
                              atoms[i].vy*atoms[i].vy+
                              atoms[i].vz*atoms[i].vz);
    }//endfor
    eKinetic *= 0.5;
  }//end routine
  void computeENhc(){
    potNhc = 0.0;
    for(int i=0;i<natm;i++){potNhc  += atomsNHC[i].posKT;}
    eKineticNhc = 0.0;
    for(int i=0;i<natm;i++){
      for(int j=0;j<len_nhc;j++){
        eKineticNhc += atomsNHC[i].m[j]*(atomsNHC[i].vx[j]*atomsNHC[i].vx[j]+
                                         atomsNHC[i].vy[j]*atomsNHC[i].vy[j]+
                                         atomsNHC[i].vz[j]*atomsNHC[i].vz[j]);
      }//endfor
    }//endfor
    eKineticNhc *= 0.5;
  }//end routine
  void computeFNhc(){
    for(int i=0; i<natm; i++){
      atomsNHC[i].fx[0] = (atoms[i].m*atoms[i].vx*atoms[i].vx-kT);
      atomsNHC[i].fy[0] = (atoms[i].m*atoms[i].vy*atoms[i].vy-kT);
      atomsNHC[i].fz[0] = (atoms[i].m*atoms[i].vz*atoms[i].vz-kT);
      for(int j=1,k=0; j<len_nhc; j++,k++){
        atomsNHC[i].fx[j]=(atomsNHC[i].m[k]*atomsNHC[i].vx[k]*atomsNHC[i].vx[k]-kT);
        atomsNHC[i].fy[j]=(atomsNHC[i].m[k]*atomsNHC[i].vy[k]*atomsNHC[i].vy[k]-kT);
        atomsNHC[i].fz[j]=(atomsNHC[i].m[k]*atomsNHC[i].vz[k]*atomsNHC[i].vz[k]-kT);
      }//endfor
      for(int j=0;j<len_nhc;j++){
        atomsNHC[i].fx[j] /= atomsNHC[i].m[j];
        atomsNHC[i].fy[j] /= atomsNHC[i].m[j];
        atomsNHC[i].fz[j] /= atomsNHC[i].m[j];
      }//endfor
    }//endfor
    eKineticNhc *= 0.5;
  }//end routine

};
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/* Structure that stores all energies */

struct EnergyStruct {
    
   int iteration_gsp;       // step upon which energies below computed
    double enl;             // non local
    double eext;            // local external energy
    double eke;             // quantum kinetic energy 
    double ehart;           // hartree energy
    double egga;            // exchange correlation grad corr
    double eexc;            // exchange correlation local
    double fictEke;         // fict KE from cp dynamics
    double fmagPsi;         // coef force magnitude
    double eewald_recip;    // atm(ion)-atm(ion) recip (computed by psi chares)
    double totalElecEnergy; // sum of electronic energies + ewald_recip 
                            // no fict and no ewald_real

   int iteration_atm;       // step upon which energies below computed.
    double eewald_real;     // Real space ewald.
    double eKinetic_atm;    // classical kinetic energy
    double eKineticNhc_atm; // NHC kinetic energy
    double potNhc_atm;      // NHC pot energy
    double fmag_atm;        // magnitude of atm forces
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
    int iteration_gsp;
    int iteration_atm;
    void updateEnergiesFromGS(EnergyStruct &);
    inline EnergyStruct getEnergyStruct(){return estruct;}
};
EnergyStruct GetEnergyStruct();
//============================================================================

#endif // #ifndef _groups_h_
