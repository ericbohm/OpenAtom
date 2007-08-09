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
class AtomMsg: public CMessage_AtomMsg {
 public:
  int nsize;
  int natmStr,natmEnd;
  double *data;  
};
//============================================================================

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
  int isokin_opt;
  int countAtm;
  int nAtmMsgRecv;
  double kT;
  double pot_ewd_rs;      // total real space ewald energy
  double vself, vbgr;
  double eKinetic;        // kinetic energy
  double eKineticNhc;     // NHC kinetic energy
  double potNhc;         // NHC potential energy
  double **px;
  double *ftot; 
  Atom *atoms;
  AtomNHC *atomsNHC;
  FastAtoms fastAtoms;

  AtomsGrp(CkMigrateMessage *m) {}
  AtomsGrp(int,int,int,int, int ,int ,int,double ,Atom *,AtomNHC *);
  ~AtomsGrp();
  void contributeforces();
  void recvContribute(CkReductionMsg *);
  void atomsDone(CkReductionMsg *);
  void atomsDone();
  void sendAtoms(double,double ,double,int,int,int);
  void acceptAtoms(AtomMsg *);
  void StartRealspaceForces();
  void outputAtmEnergy();
  void zeroforces() {
    double *fx = fastAtoms.fx;
    double *fy = fastAtoms.fy;
    double *fz = fastAtoms.fz;
    for(int i=0; i<natm; i++){
      atoms[i].fx = 0;
      atoms[i].fy = 0;
      atoms[i].fz = 0;
      fx[i]       = 0;
      fy[i]       = 0;
      fz[i]       = 0;
    }//endfor
  }//end routine
  void zeronhc(){for(int i=0;i<natm;i++){atomsNHC[i].posKT = 0;}}
  void copySlowToFast(){
    double *x  = fastAtoms.x;
    double *y  = fastAtoms.y;
    double *z  = fastAtoms.z;
    double *fx = fastAtoms.fx;
    double *fy = fastAtoms.fy;
    double *fz = fastAtoms.fz;
    for(int i=0;i<natm;i++){
      x[i]  = atoms[i].x;
      y[i]  = atoms[i].y;
      z[i]  = atoms[i].z;
      fx[i] = atoms[i].fx;
      fy[i] = atoms[i].fy;
      fz[i] = atoms[i].fz;
    }//endfor
  }//end routine
  void copyFastToSlow(){
    double *x  = fastAtoms.x;
    double *y  = fastAtoms.y;
    double *z  = fastAtoms.z;
    double *fx = fastAtoms.fx;
    double *fy = fastAtoms.fy;
    double *fz = fastAtoms.fz;
    for(int i=0;i<natm;i++){
      atoms[i].x  = x[i];
      atoms[i].y  = y[i];
      atoms[i].z  = z[i];
      atoms[i].fx = fx[i];
      atoms[i].fy = fy[i];
      atoms[i].fz = fz[i];
    }//endfor
  }//end routine

};
//============================================================================


#include "energy.h"

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
    void energyDone(CkReductionMsg *);
    void energyDone();
    inline EnergyStruct getEnergyStruct(){return estruct;}
};
EnergyStruct GetEnergyStruct();
//============================================================================


#endif // #ifndef _groups_h_
