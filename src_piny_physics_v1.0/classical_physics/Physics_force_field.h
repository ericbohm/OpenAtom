// Parallel Programming Laboratory
// http://charm.cs.uiuc.edu
//
// University of Illinois at Urbana-Champaign
// Department of Computer Science
//
// LeanMD is a barebones research application to examine performance of
// n-body type algorithms on massively parallel processors.

#ifndef __Physics_h__
#define __Physics_h__

#include <cmath>
#include "leanMD.h"
#include "Physics_defines.h"
#include "group_bond_con.h"
#include "BasePMECellData.h"
#include "BasePMEGSpaceData.h"
#include "BaseInputInterface.h"
#include "PinyInputInterface.h"
#include "RigidBondInfo.h"

#include "InterRoutines.h"

// intramolecular science routines
#include "bond.h"
#include "bend.h"
#include "tors.h"
#include "onfo.h"

void Test(int id1, int id2);

// two (or one) cells must send data here, where interactions are calculated
class Physics {
  public:
    static double getTimeStepLength() { 
      return readonly_general_data.gentimeinfo.dt; 
    }

    // calculate the forces that each atom in the list makes on the other atoms
    // return the potential energy
    static void calcForces (
      int         numAtoms, // number of atoms in the list
      Vector*     position, // size [numAtoms] positions of each atom
      Vector*     force,    // size [numAtoms] forces on each atom (will be 
                            // calc)
      int*        atomID,   // size [numAtoms] id of each atom
      int*        rootID,   // size [numAtoms] tells whether an atom is a root
                            // or a branch
      StepOutput* out
    );
    
    // calculate the forces that the first list of atoms exerts on the second
    // list of atoms
    // return the potential energy
    static void calcForces (
      int         numAtoms1, // 1st list's number of atoms
      Vector*     position1, // size [numAtoms1] 1st list's positions
      Vector*     force1,    // size [numAtoms1] 1st list's forces (will be 
                             // calc)
      int*        atomID1,   // size [numAtoms1] 1st list's atom id's
      int*        rootID1,   // size [numAtoms1] tells whether an atom is a root
                             // or a branch
      int         numAtoms2, // 2nd list's number of atoms
      Vector*     position2, // size [numAtoms2] 2nd list's positions
      Vector*     force2,    // size [numAtoms2] 2nd list's forces (will be 
                             // calc)
      int*        atomID2,   // size [numAtoms2] 2nd list's atom id's
      int*        rootID2,   // size [numAtoms2] tells whether an atom is a root
                             // or a branch
      StepOutput* out                             
    );
    
    // given a set of atom positions, velocities, and forces
    // apply the forces and calculate new positions for the atoms
    // returns the kinetic energy of the new atoms
    // <useHVel> of true means to use half-step Velocity (use on everything
    // except the first timestep
    static void doIntegration (
      int     numAtoms,      // number of atoms to do integration
      int*    atomID,        // size [numAtoms] id of the atoms
      Vector* position,      // size [numAtoms] positions of the atoms
      Vector* position_old,  // size [numAtoms] a buffer to hold old_positions
                             // temporarily
      Vector* velocity,      // size [numAtoms] velocities for each atom
      Vector* hVel,          // size [numAtoms] half-step velocities for each 
                             // atom
      Vector* force,         // size [numAtoms] force to be applied to each atom
      bool    useHVel,       // true if use the special half-velocity
      int     istep_flag,    // tell whether first, intermediate or last 
                             // iteration
      RigidBondInfo*  rBonds,// rigid bond info for shake and rattle
      StepOutput*     out
    );
      
    // force calculations for bonded constraints: (revised from MINDY code)
    // calculate <force1> and <force2> given other info
    static void calcBondForces (
      const Vector& pos1,     // 1st atom's position
      const Vector& pos2,     // 2nd atom's position
      Vector&       force1,   // 1st atom's resultant force
      Vector&       force2,   // 2nd atom's resultant force
      double        x0,
      double        k,
      int           ityp,     // bond type between atoms
      StepOutput*   out
    );

    // force calculations for angle constraints: (revised from MINDY code)
    // calculate <force1>, <force2>, and <force3> given other info
    static void calcAngleForces
    (
      const Vector& pos1,    // 1st atom's position
      const Vector& pos2,    // 2nd atom's position
      const Vector& pos3,    // 3rd atom's position 
      Vector&       force1,  // 1st atom's resultant force
      Vector&       force2,  // 2nd atom's resultant force
      Vector&       force3,  // 3rd atom's resultant force
      double        k, 
      double        theta0, 
      double        k_ub, 
      double        r_ub,
      double        pot_ub_ret,
      int           ityp,     // bond type between atoms
      StepOutput*   out
    );
    // force calculations for dihedral constraints: (revised from MINDY code)
    // calculate <force1>, <force2>, <force3>, <force4> given other info
    // multiplicity (2 or 4), k[4], delta[4], n[4] are parameters of a dihedral
    static void calcDihedralForces
    (
      const Vector& pos1,    // 1st atom's position
      const Vector& pos2,    // 2nd atom's position
      const Vector& pos3,    // 3rd atom's position 
      const Vector& pos4,    // 4th atom's position
      Vector&       force1,  // 1st atom's resultant force
      Vector&       force2,  // 2nd atom's resultant force
      Vector&       force3,  // 3rd atom's resultant force
      Vector&       force4,  // 4th atom's resultant force
      int           multiplicity, 
      double        k[], 
      double        delta[], 
      int           n[],
      StepOutput*   out
    );
    // force calculations for improper constraints: (revised from MINDY code)
    // calculate <force1>, <force2>, <force3>, <force4> given other info
    static void calcImproperForces
    (
      const Vector& pos1,    // 1st atom's position
      const Vector& pos2,    // 2nd atom's position
      const Vector& pos3,    // 3rd atom's position 
      const Vector& pos4,    // 4th atom's position
      Vector&       force1,  // 1st atom's resultant force
      Vector&       force2,  // 2nd atom's resultant force
      Vector&       force3,  // 3rd atom's resultant force
      Vector&       force4,  // 4th atom's resultant force
      int           multiplicity, 
      double        k[], 
      double        delta[], 
      int           n[],
      StepOutput*   out
    );

    static void calcTorsionForces
    (
      const Vector& pos1,    // 1st atom's position
      const Vector& pos2,    // 2nd atom's position
      const Vector& pos3,    // 3rd atom's position 
      const Vector& pos4,    // 4th atom's position
      Vector&       force1,  // 1st atom's resultant force
      Vector&       force2,  // 2nd atom's resultant force
      Vector&       force3,  // 3rd atom's resultant force
      Vector&       force4,  // 4th atom's resultant force
      int           ityp,    // bond type between atoms
      StepOutput*   out
    );
    
    static void calc14Forces
    (
      const Vector& pos1,     // 1st atom's position
      const Vector& pos2,     // 2nd atom's position
      Vector&       force1,   // 1st atom's resultant force
      Vector&       force2,   // 2nd atom's resultant force
      double        q1,
      double        q2,
      int           ityp,
      StepOutput*   out
    );


    static double calcEcorrForces
    (
      const Vector& pos1,     // 1st atom's position
      const Vector& pos2,     // 2nd atom's position
      Vector&       force1,   // 1st atom's resultant force
      Vector&       force2,   // 2nd atom's resultant force
      double        charge1,
      double        charge2,
      double        alp_ewd   //  Ewald parameter : >= 3.5/cutoff
    );

    // find the maximum exclusion distance squared
    static double get_r2max_excl (int     nAtms, // num atoms
                                  double* xPos,  // size [nAtms]
                                  double* yPos,  // size [nAtms]
                                  double* zPos); // size [nAtms]
                                  
    static BasePMECellData* newPMECellData ();
    
    static BasePMEGSpaceData* newPMEGSpaceData ();
    
    static void calcGridSize (double* _hmati, double _alpha_ewald,
                              int& _ngrid_a, int& _ngrid_b, int& _ngrid_c);
                              
    static BaseInputInterface* newInputInterface () {
      return new PinyInputInterface ();
    }
    
  private:
    static const double TIME_DELTA;

//=============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=============================================================================
// calculate the force that atom1 exerts on atom2,
// return the potential energy
//=============================================================================
static void calcForce
(
  const Vector& position1, // 1st atom's position
  Vector&       force1,    // 1st atom's force (will be calculated)
  double        charge1,   // 1st atom's charge
  double        epsilon1,  // 1st atom's well depth
  double        sigma1,    // 1st atom's zero energy pt.
  double        cut1,      // 1st atom's cutoff
  int           ind1_excl, // 1st atom's index in exclusion list
  const Vector& position2, // 2nd atom's position
  Vector&       force2,    // 2nd atom's force (will be calculated)
  double        charge2,   // 2nd atom's charge
  double        epsilon2,  // 2nd atom's well depth
  double        sigma2,    // 2nd atom's zero energy pt.
  double        cut2,      // 2nd atom's cutoff
  int           ind2_excl, // 2nd atom's index in exclusion list
  int&          processed, // flag
  int           iperd,     // periodicity (0 or 3 for now)
  double        alp_ewd,   // ewald sum convergence parameter >= 3.5/cutoff
  double        r2max_excl,// maximum exclusion distance squared
  StepOutput*   out
);

static void doIntegration (
  int         numAtoms,
  int*        atomID,
  Vector*     position,
  Vector*     velocity,
  Vector*     force,
  StepOutput* out
);

static void doIntegration (
  int         numAtoms,
  int*        atomID,
  Vector*     velocity,
  Vector*     force,
  StepOutput* out
);

};
//=============================================================================

inline void Physics::calcForce
(
  const Vector& position1, // 1st atom's position
  Vector&       force1,    // 1st atom's force (will be calculated)
  double        charge1,   // 1st atom's charge
  double        epsilon1,  // 1st atom's well depth
  double        sigma1,    // 1st atom's zero energy pt.
  double        cut1,      // 1st atom's cutoff
  int           ind1_excl, // 1st atom's index in exclusion list
  const Vector& position2, // 2nd atom's position
  Vector&       force2,    // 2nd atom's force (will be calculated)
  double        charge2,   // 2nd atom's charge
  double        epsilon2,  // 2nd atom's well depth
  double        sigma2,    // 2nd atom's zero energy pt.
  double        cut2,      // 2nd atom's cutoff
  int           ind2_excl, // 2nd atom's index in exclusion list
  int&          processed, // flag
  int           iperd,     // periodicity (0 or 3 for now)
  double        alp_ewd,   // ewald sum convergence parameter >= 3.5/cutoff
  double        r2max_excl,// maximum exclusion distance squared
  StepOutput*   out
) {
#ifdef USE_SPLINE_LJ
  InterRoutines::calc_force_spline(position1,force1,charge1,cut1,ind1_excl,
                                   position2,force2,charge2,cut2,ind2_excl,
                                   processed,iperd,r2max_excl,out);
#else
  InterRoutines::calc_force_explicit(position1,force1,charge1,epsilon1,sigma1,
                                     cut1,ind1_excl,position2,force2,charge2,
                                     epsilon2,sigma2,cut2,ind2_excl,
                                     processed,iperd,alp_ewd,r2max_excl,out);  
#endif
}

#endif // __Physics_h__




