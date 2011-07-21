/** \file groups.h
 */
#ifndef GROUPS_H
#define GROUPS_H
#include "Atoms.h"
struct EnergyStruct;
#include "groups.decl.h"
#include "CPcharmParaInfo.decl.h"


#include "uber/Uber.h"

class AtomMsg: public CMessage_AtomMsg 
{
    public:
        int nsize;
        int natmStr,natmEnd;
        double *data;  
};


class AtomXYZMsg: public CMessage_AtomXYZMsg 
{
    public:
        double *x;  
        double *y;  
        double *z;  
	int index;
};



/** AtomsGrp class.
 * various chares use CkLocal to access the atoms:
 * Rho_RHart RhoGHart StructureFactor RealParticlePlane
 * ParticlePlane RealSpacePlane GSpacePlane
 * EJB: this could perhaps be better implemented in an MSA
 */
class AtomsGrp: public Group {
 public:
  const UberCollection thisInstance;
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
  double potPerdCorr;
  double vself, vbgr;
  double eKinetic;        // kinetic energy
  double eKineticNhc;     // NHC kinetic energy
  double potNhc;         // NHC potential energy
  double **px;
  double *ftot; 
  Atom *atoms;
  AtomNHC *atomsNHC;
  FastAtoms fastAtoms;
  PIMD_CM PIMD_CM_Atoms;
  CProxySection_AtomsGrp proxyHeadBeads;
  bool atomsCMrecv,atomsPIMDXrecv;
  bool amBeadRoot, amZerothBead;
  int acceptCountfu;
  int acceptCountX;
  int acceptCountu;
  int PIBeadIndex;
  int TemperIndex;
  FILE *temperScreenFile;
  int numPIMDBeads;
  AtomsGrp(CkMigrateMessage *m) {}
  AtomsGrp(int,int,int,int, int ,int ,int,double ,Atom *,AtomNHC *, UberCollection thisInstance);
  void init();
  ~AtomsGrp();
  void contributeforces();
  void integrateAtoms();
  void accept_PIMD_x(double _x, double _y, double _z, int atomI);
  void accept_PIMD_Fu(double _fxu, double _fyu, double _fzu, int atomI);
  void accept_PIMD_CM(AtomXYZMsg *m);
  void accept_PIMD_u(double _ux, double _uy, double _uz, int atomI);

  void recvContribute(CkReductionMsg *);
  void atomsDone(CkReductionMsg *);
  void send_PIMD_u();
  void send_PIMD_Fx();
  void send_PIMD_x();
  void atomsDone();
  void sendAtoms(double,double ,double,int,int,int);
  void acceptAtoms(AtomMsg *);
  void startRealSpaceForces();
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




#include "energy.h"

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
    inline EnergyStruct getEnergyStruct(){return estruct;}
};
/*EnergyStruct GetEnergyStruct();*/
//============================================================================


#endif // GROUPS_H
