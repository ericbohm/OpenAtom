//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/** \file AtomsCompute.h
 */
//==========================================================================
#ifndef ATOMSCOMPUTE_H
#define ATOMSCOMPUTE_H
#include "Atoms.h"

struct EnergyStruct;

#include "Atoms.decl.h"
#include "atomMessages.h"
#include "CPcharmParaInfo.decl.h"


#include "uber/Uber.h"
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================


/** \file AtomsCompute.h

 * \brief Accepts reduction of forces from AtomsCache,
 * integrates forces to produce new coordinates,
 * distributes new coordinates to AtomsCache.

 * If PIMD, then PIMD coordination is done before producing the
 * updated coordinates.

 * If Tempers, then Temperature exchange is done before producing the
 * updated coordinates.
 * @addtogroup Atoms
 * @{ 
 */
//==========================================================================
class AtomsCompute: public CBase_AtomsCompute {
  //==========================================================================
  public:
    const UberCollection thisInstance;
    /* inherited all the vars from the old AtomsGrp */
    int natm;
    int natm_nl;
    int len_nhc;
    int iextended_on;
    int cp_min_opt;
    int cp_bomd_opt;
    int cp_wave_opt;
    int *iteration;
    int isokin_opt;
    int countAtm;
    int nAtomsPerChare;
    int nChareAtoms;
    int nAtmMsgRecv;
    int PIBeadIndex;
    int TemperIndex;
    int numPIMDBeads;
    int tol_reached;
    unsigned int handleForcesCount;
    double kT;              // temperature
    double pot_ewd_rs;      // total real space ewald energy
    double potPerdCorr;
    double potPIMDChain;    // Path integral chain energy :computed in u space
    double vself, vbgr;
    double eKinetic;        // kinetic energy
    double eKineticNhc;     // NHC kinetic energy
    double potNhc;          // NHC potential energy
    double **px;
    double *ftot; 
    double tau,beta,omega2PIMD;
    double *massPIMDScal;  // veig of PINY : length numPIMDBeads
    CkReductionMsg *contribMsg[2];
    Atom *atoms;
    AtomNHC *atomsNHC;
    FastAtoms fastAtoms;
    PIMD_CM PIMD_CM_Atoms;
    CProxySection_AtomsCompute proxyHeadBeads;
    CProxySection_AtomsCompute proxyAllBeads;
    bool atomsCMrecv,atomsPIMDXrecv;
    bool amBeadRoot, amZerothBead;
    int acceptCountfu;
    int acceptCountX;
    int acceptCountu;

    int ktemps;
    FILE *temperScreenFile;

    AtomsCompute(CkMigrateMessage *m) {}
    AtomsCompute(int, int, int, int, int, int, int ,int, double, Atom *, AtomNHC *, int nChareAtoms, UberCollection thisInstance);
    void init();
    ~AtomsCompute();
    void integrateAtoms();
    void accept_PIMD_x(double _x, double _y, double _z, int atomI);
    void accept_PIMD_Fu(double _fxu, double _fyu, double _fzu, int atomI);
    void accept_PIMD_Fu_and_u(double _fxu, double _fyu, double _fzu, double _xu, double _yu, double _zu, int atomI);
    void accept_PIMD_CM(AtomXYZMsg *m);
    void accept_PIMD_u(double _ux, double _uy, double _uz, int atomI);
    void acceptNewTemperature(double temp);
    void recvContribute(CkReductionMsg *);
    void recvContributeForces(CkReductionMsg *);
    void atomsDone(CkReductionMsg *);
    void send_PIMD_u();
    void send_PIMD_Fx();
    void send_PIMD_Fx_and_x();
    void send_PIMD_x();
    void sendAtoms(double,double ,double,double,int,int,int);
    void acceptAtoms(AtomMsg *);
    void outputAtmEnergy();
    void bcastAtomsToAtomCache();
    void startRealSpaceForces(int);
    void releaseGSP();
    void handleForces();
    //==========================================================================
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
    //==========================================================================
    void zeronhc(){for(int i=0;i<natm;i++){atomsNHC[i].posKT = 0;}}
    //==========================================================================
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
    //==========================================================================
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
    //==========================================================================
    void initPIMD(){
      beta            = 1.0/kT;
      tau             = beta/((double)numPIMDBeads);
      omega2PIMD      = 1.0/(tau*beta);
      massPIMDScal[0] = 1.0;  // vieg[0]=0; careful
      for(int ip=2;ip<=numPIMDBeads;ip++){
        massPIMDScal[(ip-1)] = ((double)(ip))/((double)(ip-1));
      }//endfor
    }//end routine
    //==========================================================================
    void switchPIMDBeadForceMass(int mybead, int natmStr, int natmEnd, double *pot){
#define _CHECK_PIMD_TRANSFORMS_OFF_
#ifdef _CHECK_PIMD_TRANSFORMS_
      CkPrintf("switchPIMDBeadForceMass: I am bead %d strt dumping stuff\n",mybead);
#endif
      for(int i=natmStr;i<natmEnd;i++){
#ifdef _CHECK_PIMD_TRANSFORMS_
        CkPrintf("b=%d atm=%d : %g %g %g : %g %g %g : %g %g %g : %g %g %g\n",
            mybead,i,
            atoms[i].x,atoms[i].y,atoms[i].z,atoms[i].xu,atoms[i].yu,atoms[i].zu,
            atoms[i].fx,atoms[i].fy,atoms[i].fz,atoms[i].fxu,atoms[i].fyu,atoms[i].fzu);
#endif
        atoms[i].fx   = atoms[i].fxu; 
        atoms[i].fy   = atoms[i].fyu; 
        atoms[i].fz   = atoms[i].fzu; 
        atoms[i].xold = atoms[i].x; 
        atoms[i].yold = atoms[i].y; 
        atoms[i].zold = atoms[i].z; 
        atoms[i].x    = atoms[i].xu; 
        atoms[i].y    = atoms[i].yu; 
        atoms[i].z    = atoms[i].zu; 
      }//endfor
#ifdef _CHECK_PIMD_TRANSFORMS_
      CkPrintf("switchPIMDBeadForceMass: I am bead %d end dumping stuff\n",mybead);
#endif
      if(mybead>1){ // mybead goes from 1 to numPIMDBead
        int mybead1 = mybead-1;
        pot[0] = 0.0;
        for(int i=natmStr;i<natmEnd;i++){
          atoms[i].m   *= massPIMDScal[mybead1];
          atoms[i].fx  -= omega2PIMD*atoms[i].m*atoms[i].x;
          atoms[i].fy  -= omega2PIMD*atoms[i].m*atoms[i].y;
          atoms[i].fz  -= omega2PIMD*atoms[i].m*atoms[i].z;
          pot[0]       += omega2PIMD*atoms[i].m*(
              atoms[i].x*atoms[i].x
              +atoms[i].y*atoms[i].y
              +atoms[i].z*atoms[i].z);
        }//endfor
        pot[0] *= 0.5;
      }//endif
    }//end routine
    //==========================================================================
    void unswitchPIMDMass(int mybead, int natmStr, int natmEnd){
      int mybead1 = mybead-1;
      for(int i=natmStr;i<natmEnd;i++){
        atoms[i].m   /= massPIMDScal[mybead1]; // =1 for mybead=0
        atoms[i].xu   = atoms[i].x; 
        atoms[i].yu   = atoms[i].y; 
        atoms[i].zu   = atoms[i].z; 
      }//endfor
    }//end routine
    //==========================================================================

    //==========================================================================
}; //end class
//==========================================================================
/*@}*/
#endif // ATOMSCOMPUTE_H
