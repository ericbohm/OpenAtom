//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                     Lean_MD-PINY_Phyics
//
// Purpose :                                                                
//      Bare bones physics routines. Will be spit into many files eventually.
//
// Authors :
//      Parallel Programming Laboratory, CS, UIUC
//      Tuckerman Group, Chemistry, New York University
//
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================


#include "Physics.h"
#include "Physics_defines.h"
#include "AtomDatabase.h"
#include "PMEGSpaceData.h"
#include "PMECellData.h"

#include "leanMD.decl.h"

// hardwire the time step for now
//const double Physics::TIME_DELTA = 1.0/TIME_CONV;

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// double Physics::calcForces
//     Intermolecular real-space interactions for atoms in the same cell:
//      calculate the forces between atoms in the same cell
//      and return the potential energy
//==============================================================================
  void Physics::calcForces
(
 int         numAtoms,  // number of atoms in the list
 Vector*     position,  // size [numAtoms] positions of each atom
 Vector*     force,     // size [numAtoms] resulting forces on each atom
 int*        atomID,     // size [numAtoms] id of each atom
 int*        rootID,
 StepOutput* out
 )
{
  double q0, q1;
  double m0, m1;
  double epsilon0, epsilon1, sigma0, sigma1;
  int    typ0, typ1;
  int    processed;
  int    nProcessed  = 0;
  int    iperd       = 3;

  AtomInfo* atomInfo = AtomInfoProxy.ckLocalBranch(); 

  double alpha_ewald = atomInfo->_alpha_ewald;
  double cut0, cut1;
  double r2max_excl  = atomInfo->get_r2max_excl ();

  for (int i=0; i<numAtoms; i++) {
    atomInfo->getInfo((int)atomID[i], m0, q0, epsilon0, 
        sigma0, typ0);

    cut0 = atomInfo->getCutoff(typ0);
    for (int j=i+1; j<numAtoms; j++) {
      atomInfo->getInfo((int)atomID[j], m1, q1, epsilon1, 
          sigma1, typ1);
      cut1 = atomInfo->getCutoff(typ1);

      calcForce(position[i], force[i], q0, epsilon0, sigma0, cut0, atomID [i],
          position[j], force[j], q1, epsilon1, sigma1, cut1, atomID [j],
          processed, iperd, alpha_ewald, r2max_excl, out);

      nProcessed += processed;
    }
  }

#ifdef __BLUEGENE__
  BgElapse(27e-8*nProcessed);
#endif
}
//==============================================================================



//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// double Physics::calcForces
//     Intermolecular real-space interactions for atoms in the same cell:
//      calculate the forces between atoms in the different cells
//      return the potential energy
//==============================================================================
  void Physics::calcForces
(
 int         numAtoms1,  // 1st list's number of atoms
 Vector*     position1,  // size [numAtoms1] 1st list's positions
 Vector*     force1,     // size [numAtoms1] 1st list's forces (will be calc)
 int*        atomID1,    // size [numAtoms1] 1st list's atom id's
 int*        rootID1,
 int         numAtoms2,  // 2nd list's number of atoms
 Vector*     position2,  // size [numAtoms2] 2nd list's positions
 Vector*     force2,     // size [numAtoms2] 2nd list's forces (will be calc)
 int*        atomID2,     // size [numAtoms2] 2nd list's atom id's
 int*        rootID2,
 StepOutput* out
 )
{
  double q0, q1;
  double m0, m1;
  double epsilon0, epsilon1, sigma0, sigma1;
  int    processed;
  int    nProcessed  = 0;
  int    iperd       = 3;
  int    typ0, typ1;

  AtomInfo* atomInfo = AtomInfoProxy.ckLocalBranch(); 

  double alpha_ewald = atomInfo->_alpha_ewald;
  double cut0, cut1;
  double r2max_excl  = atomInfo->get_r2max_excl ();

  for (int i=0; i<numAtoms1; i++) {
    atomInfo->getInfo((int)atomID1[i],m0,q0,epsilon0,sigma0,typ0);
    cut0 = atomInfo->getCutoff(typ0);

    for (int j=0; j<numAtoms2; j++) {
      atomInfo->getInfo((int)atomID2[j],m1,q1,epsilon1,sigma1,typ1);
      cut1 = atomInfo->getCutoff(typ1);

      calcForce(position1[i],force1[i],q0,epsilon0,sigma0, cut0, atomID1[i],
          position2[j],force2[j],q1,epsilon1,sigma1, cut1, atomID2[j],
          processed,iperd, 
          alpha_ewald, r2max_excl, out);
      nProcessed += processed;
    }
  }

#ifdef __BLUEGENE__
  BgElapse(27e-8*nProcessed);
#endif
}
//==============================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// double Physics::doIntegration
//       Integration API
//       Given a set of atom positions, velocities, and forces
//       Apply the forces and calculate new positions and velocities for the atoms
//==========================================================================
void Physics::doIntegration (int             numAtoms,
    int*            atomID,
    Vector*         position,
    Vector*         position_old,
    Vector*         velocity,
    Vector*         hVel,
    Vector*         force,
    bool            useHVel,
    int             istep_flag,
    RigidBondInfo*  rBonds,
    StepOutput*     out)// rigid bond info for shake and rattle 
  //==========================================================================
{ // begin 
  //==========================================================================
  int    num_constraint = 0;
  double dt             = getTimeStepLength();

  if (NULL != rBonds) {
    num_constraint += rBonds->num_21;
    num_constraint += rBonds->num_33;
    num_constraint += rBonds->num_46;
  }

  //==========================================================================
  // Perform integration based on istep_flag : 
  //    [ get force, 1st half step, loop over time{2nd (1st,if not last)} ]
  //     0 == 1st time step
  //     1 == nth time step
  //     2 == last time step
  //==========================================================================
  // do 2nd half integration (moves vels and computes Kinetic energy)
  // if its not the very 1st time step

  if (0 != istep_flag){

    doIntegration (numAtoms, atomID, velocity, force, out);
    if(num_constraint>0){ doRattle(numAtoms,
        position,
        position_old,
        velocity,
        atomID,
        dt,
        rBonds); } //rattle constraints amongst numAtoms 

  }// endif

  //==========================================================================
  // do 1st half integration (moves pos and vels)
  // if its not the very last time step

  if (2 != istep_flag){ 

    if(num_constraint>0){
      for(int i=0;i<numAtoms;i++){
        position_old[i].x=position[i].x;
        position_old[i].y=position[i].y;
        position_old[i].z=position[i].z;
      }//end for
    }//endif
    doIntegration (numAtoms, atomID, position, velocity,force, out);
    if(num_constraint>0){ doShake(numAtoms,
        position,
        position_old,
        velocity,
        atomID,
        dt,
        rBonds); } //shake constraints amongst numAtoms 

  }// endif

  //==========================================================================
}// end Physics::doIntegration
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// double Physics::doIntegration  
//       1st half time step integration
//       Given a set of atom positions, velocities, and forces
//       Apply the forces and calculate new positions and velocities for the atoms
//==============================================================================
  void Physics::doIntegration
(
 int         numAtoms,  // number of atoms to do integration
 int*        atomID,    // size [numAtoms] id of the atoms
 Vector*     position,  // size [numAtoms] positions of the atoms
 Vector*     velocity,  // size [numAtoms] velocities for each atom
 Vector*     force,     // size [numAtoms] force to be applied to each atom
 StepOutput* out
 )
  //==============================================================================
{ // begin 
  //==============================================================================
#ifdef __BLUEGENE__
  BgElapse(80e-6);
#endif

  int i;
  double dt             = getTimeStepLength();
  double dt2            = dt/2;

  //--------------------------------------------------------------------------------
  for(i=0;i<numAtoms;i++) {
    double m,c,epsilon, sigma;
    int    typ;
    AtomInfoProxy.ckLocalBranch()->getInfo((int)atomID[i], m, c, epsilon, sigma, typ);
    double dt2m = dt2/m;
    Vector ta   = force[i]*dt2m;
    Vector vel  = velocity[i];
    Vector pos  = position[i];

    // 1st half step integration
    vel   += ta;
    pos   += vel*dt;
    // store half step position and velocity
    velocity[i] = vel;
    position[i] = pos;
  }
  //--------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------
  for(i=0;i<numAtoms;i++){
    force[i].x     = force[i].y = force[i].z = 0.0;
  }
  //--------------------------------------------------------------------------------
}// end Physics::doIntegration
//==============================================================================



//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// double Physics::doIntegration  
//       2nd half time step integration
//       Given a set of atom positions, velocities, and forces
//       Apply the forces and calculate new positions and velocities for the atoms
//==============================================================================
  void Physics::doIntegration
(
 int         numAtoms, // number of atoms to do integration
 int*        atomID,   // size [numAtoms] id of the atoms
 Vector*     velocity, // size [numAtoms] velocities for each atom
 Vector*     force,    // size [numAtoms] force to be applied to each atom
 StepOutput* out
 )
  //==============================================================================
{ // begin
  //==============================================================================

#ifdef __BLUEGENE__
  BgElapse(80e-6);
#endif

  int i;
  double dt             = getTimeStepLength();
  double dt2            = dt/2;
  double Kinetic_Energy = 0;

  //==============================================================================

  for(i=0;i<numAtoms;i++) {
    double m,c,epsilon, sigma;
    int    typ;
    AtomInfoProxy.ckLocalBranch()->getInfo((int)atomID[i], m, c, epsilon, sigma, typ);
    double dt2m = dt2/m;
    Vector ta   = force[i]*dt2m;
    Vector vel  = velocity[i];
    // 2nd half step integration
    vel   += ta;
    // store velocity
    velocity[i]    = vel;
    // kinetic energy
    Kinetic_Energy = Kinetic_Energy + m*(vel.sqrMag())*0.5;
  }

  out->ke += Kinetic_Energy;

  //==============================================================================
}// end Physics::doIntegration
//==============================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// Force calculations for bonded atoms
//    calculate <force1> and <force2> given other info
//    where is the PE?
//==============================================================================
void Physics::calcBondForces (
    const Vector& pos1,     // 1st atom's position
    const Vector& pos2,     // 2nd atom's position
    Vector&       force1,   // 1st atom's resultant force
    Vector&       force2,   // 2nd atom's resultant force
    double        x0,
    double        k,
    int           ityp,     // bond type between atoms
    StepOutput*   out
    )
  //==============================================================================
{
  double* pvten     = out->pt.pvten;    // Pressure tensor [9]
  double* pvten_tot = out->pt.pvten_tot;// Pressure tensor [9]
  // pvten - 1 because piny wants arrays to start from 1  
  bond (pos1, pos2, force1, force2, pvten-1, pvten_tot-1, ityp, out);
}
//==============================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// Force calculations for bends
//    calculate <force1>, <force2>, and <force3> given other info
//    where is PE?
//==============================================================================
  void Physics::calcAngleForces
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
 )
  //==============================================================================
{
  double* pvten     = out->pt.pvten;    // Pressure tensor [9]
  double* pvten_tot = out->pt.pvten_tot;// Pressure tensor [9]

  // pvten - 1 because piny wants arrays to start from 1
  bend (pos1, pos2, pos3, force1, force2, force3, pvten-1, 
      pvten_tot-1, ityp, out);

}
//==============================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// Force calculations for torsions
//   calculate <force1>, <force2>, <force3>, <force4> given other info
//==============================================================================
  void Physics::calcTorsionForces
(
 const Vector& pos1,    // 1st atom's position
 const Vector& pos2,    // 2nd atom's position
 const Vector& pos3,    // 3rd atom's position 
 const Vector& pos4,    // 4th atom's position
 Vector&       force1,  // 1st atom's resultant force
 Vector&       force2,  // 2nd atom's resultant force
 Vector&       force3,  // 3rd atom's resultant force
 Vector&       force4,  // 4th atom's resultant force
 int           ityp,     // bond type between atoms
 StepOutput*   out
 )
  //==============================================================================
{
  double* pvten     = out->pt.pvten;    // Pressure tensor [9]
  double* pvten_tot = out->pt.pvten_tot;// Pressure tensor [9]
  // pvten - 1 because piny wants arrays to start from 1
  tors (pos1, pos2, pos3, pos4,
      force1, force2, force3, force4,
      pvten-1, pvten_tot-1, ityp, out);
}
//==============================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// Force calculations for Dihedral
//   calculate <force1>, <force2>, <force3>, <force4> given other info
//   multiplicity (2 or 4), k[4], delta[4], n[4] are parameters of a dihedral
//   where is the PE
//==============================================================================
  void Physics::calcDihedralForces
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
 )
  //==============================================================================
{
  double& pot = out->pe;
  const Vector r12 = pos1 - pos2;
  const Vector r23 = pos2 - pos3;
  const Vector r34 = pos3 - pos4;

  Vector dcosdA;
  Vector dcosdB;
  Vector dsindC;
  Vector dsindB;
  Vector f1, f2, f3;

  Vector A, B, C;
  A.cross(r12, r23);
  B.cross(r23, r34);
  C.cross(r23, A);

  double rA = A.Mag();
  double rB = B.Mag();
  double rC = C.Mag();

  double cos_phi = (A*B)/(rA*rB);
  double sin_phi = (C*B)/(rC*rB);

  // Normalize B
  rB = 1.0/rB;
  B *= rB;

  double phi = -atan2(sin_phi, cos_phi);

  if (fabs(sin_phi) > 0.1) {
    // Normalize A
    rA = 1.0/rA;
    A *= rA;
    dcosdA = (A*cos_phi-B)*rA;
    dcosdB = (B*cos_phi-A)*rB;
  }
  else {
    // Normalize C
    rC = 1.0/rC;
    C *= rC;
    dsindC = (C*sin_phi-B)*rC;
    dsindB = (B*sin_phi-C)*rB;
  }

  int mult = multiplicity;
  for (int j=0; j<mult; j++) {
    double kk = k[j];
    double nn = n[j];
    double deltaj = delta[j];
    double K, K1;
    if (nn) {
      K = kk * (1.0+cos(nn*phi + deltaj));
      K1 = -nn*kk*sin(nn*phi + deltaj);
    }
    else {
      double diff = phi-deltaj;
      if (diff < -M_PI) diff += 2.0*M_PI;
      else if (diff > M_PI) diff -= 2.0*M_PI;
      K = kk*diff*diff;
      K1 = 2.0*kk*diff;
    }

    // forces
    if (fabs(sin_phi) > 0.1) {
      K1 = K1/sin_phi;
      f1.x += K1*(r23.y*dcosdA.z - r23.z*dcosdA.y);
      f1.y += K1*(r23.z*dcosdA.x - r23.x*dcosdA.z);
      f1.z += K1*(r23.x*dcosdA.y - r23.y*dcosdA.x);

      f3.x += K1*(r23.z*dcosdB.y - r23.y*dcosdB.z);
      f3.y += K1*(r23.x*dcosdB.z - r23.z*dcosdB.x);
      f3.z += K1*(r23.y*dcosdB.x - r23.x*dcosdB.y);

      f2.x += K1*(r12.z*dcosdA.y - r12.y*dcosdA.z
          + r34.y*dcosdB.z - r34.z*dcosdB.y);
      f2.y += K1*(r12.x*dcosdA.z - r12.z*dcosdA.x
          + r34.z*dcosdB.x - r34.x*dcosdB.z);
      f2.z += K1*(r12.y*dcosdA.x - r12.x*dcosdA.y
          + r34.x*dcosdB.y - r34.y*dcosdB.x);
    }
    else {
      //  This angle is closer to 0 or 180 than it is to
      //  90, so use the cos version to avoid 1/sin terms
      K1 = -K1/cos_phi;

      f1.x += K1*((r23.y*r23.y + r23.z*r23.z)*dsindC.x
          - r23.x*r23.y*dsindC.y
          - r23.x*r23.z*dsindC.z);
      f1.y += K1*((r23.z*r23.z + r23.x*r23.x)*dsindC.y
          - r23.y*r23.z*dsindC.z
          - r23.y*r23.x*dsindC.x);
      f1.z += K1*((r23.x*r23.x + r23.y*r23.y)*dsindC.z
          - r23.z*r23.x*dsindC.x
          - r23.z*r23.y*dsindC.y);

      Vector tmpf3;
      tmpf3.cross(dsindB,r23);
      f3 += tmpf3*K1;

      f2.x += K1*(-(r23.y*r12.y + r23.z*r12.z)*dsindC.x
          +(2.0*r23.x*r12.y - r12.x*r23.y)*dsindC.y
          +(2.0*r23.x*r12.z - r12.x*r23.z)*dsindC.z
          +dsindB.z*r34.y - dsindB.y*r34.z);
      f2.y += K1*(-(r23.z*r12.z + r23.x*r12.x)*dsindC.y
          +(2.0*r23.y*r12.z - r12.y*r23.z)*dsindC.z
          +(2.0*r23.y*r12.x - r12.y*r23.x)*dsindC.x
          +dsindB.x*r34.z - dsindB.z*r34.x);
      f2.z += K1*(-(r23.x*r12.x + r23.y*r12.y)*dsindC.z
          +(2.0*r23.z*r12.x - r12.z*r23.x)*dsindC.x
          +(2.0*r23.z*r12.y - r12.z*r23.y)*dsindC.y
          +dsindB.y*r34.x - dsindB.x*r34.y);
    }
  } // end loop over multiplicity
  force1 += f1;
  force2 += f2-f1;
  force3 += f3-f2;
  force4 += -f3;
}
//==============================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// Force calculations for improper dihedrals
//   calculate <force1>, <force2>, <force3>, <force4> given other info
//   WHERE IS PE?
//==============================================================================
  void Physics::calcImproperForces
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
 )
  //==============================================================================
{
  const Vector r12 = pos1 - pos2;
  const Vector r23 = pos2 - pos3;
  const Vector r34 = pos3 - pos4;

  Vector dcosdA;
  Vector dcosdB;
  Vector dsindC;
  Vector dsindB;
  Vector f1, f2, f3;

  Vector A, B, C;
  A.cross(r12, r23);
  B.cross(r23, r34);
  C.cross(r23, A);

  double rA = A.Mag();
  double rB = B.Mag();
  double rC = C.Mag();

  double cos_phi = (A*B)/(rA*rB);
  double sin_phi = (C*B)/(rC*rB);

  // Normalize B
  rB = 1.0/rB;
  B *= rB;

  double phi = -atan2(sin_phi, cos_phi);

  if (fabs(sin_phi) > 0.1) {
    // Normalize A
    rA = 1.0/rA;
    A *= rA;
    dcosdA = (A*cos_phi-B)*rA;
    dcosdB = (B*cos_phi-A)*rB;
  }
  else {
    // Normalize C
    rC = 1.0/rC;
    C *= rC;
    dsindC = (C*sin_phi-B)*rC;
    dsindB = (B*sin_phi-C)*rB;
  }

  int mult = multiplicity;
  for (int j=0; j<mult; j++) {
    double kk = k[j];
    double nn = n[j];
    double deltaj = delta[j];
    double K, K1;
    if (nn) {
      K = kk * (1.0+cos(nn*phi + deltaj));
      K1 = -nn*kk*sin(nn*phi + deltaj);
    }
    else {
      double diff = phi-deltaj;
      if (diff < -M_PI) diff += 2.0*M_PI;
      else if (diff > M_PI) diff -= 2.0*M_PI;
      K = kk*diff*diff;
      K1 = 2.0*kk*diff;
    }

    // forces
    if (fabs(sin_phi) > 0.1) {
      K1 = K1/sin_phi;
      f1.x += K1*(r23.y*dcosdA.z - r23.z*dcosdA.y);
      f1.y += K1*(r23.z*dcosdA.x - r23.x*dcosdA.z);
      f1.z += K1*(r23.x*dcosdA.y - r23.y*dcosdA.x);

      f3.x += K1*(r23.z*dcosdB.y - r23.y*dcosdB.z);
      f3.y += K1*(r23.x*dcosdB.z - r23.z*dcosdB.x);
      f3.z += K1*(r23.y*dcosdB.x - r23.x*dcosdB.y);

      f2.x += K1*(r12.z*dcosdA.y - r12.y*dcosdA.z
          + r34.y*dcosdB.z - r34.z*dcosdB.y);
      f2.y += K1*(r12.x*dcosdA.z - r12.z*dcosdA.x
          + r34.z*dcosdB.x - r34.x*dcosdB.z);
      f2.z += K1*(r12.y*dcosdA.x - r12.x*dcosdA.y
          + r34.x*dcosdB.y - r34.y*dcosdB.x);
    }
    else {
      //  This angle is closer to 0 or 180 than it is to
      //  90, so use the cos version to avoid 1/sin terms
      K1 = -K1/cos_phi;

      f1.x += K1*((r23.y*r23.y + r23.z*r23.z)*dsindC.x
          - r23.x*r23.y*dsindC.y
          - r23.x*r23.z*dsindC.z);
      f1.y += K1*((r23.z*r23.z + r23.x*r23.x)*dsindC.y
          - r23.y*r23.z*dsindC.z
          - r23.y*r23.x*dsindC.x);
      f1.z += K1*((r23.x*r23.x + r23.y*r23.y)*dsindC.z
          - r23.z*r23.x*dsindC.x
          - r23.z*r23.y*dsindC.y);

      Vector tmpf3;
      tmpf3.cross(dsindB,r23);
      f3 += tmpf3*K1;

      f2.x += K1*(-(r23.y*r12.y + r23.z*r12.z)*dsindC.x
          +(2.0*r23.x*r12.y - r12.x*r23.y)*dsindC.y
          +(2.0*r23.x*r12.z - r12.x*r23.z)*dsindC.z
          +dsindB.z*r34.y - dsindB.y*r34.z);
      f2.y += K1*(-(r23.z*r12.z + r23.x*r12.x)*dsindC.y
          +(2.0*r23.y*r12.z - r12.y*r23.z)*dsindC.z
          +(2.0*r23.y*r12.x - r12.y*r23.x)*dsindC.x
          +dsindB.x*r34.z - dsindB.z*r34.x);
      f2.z += K1*(-(r23.x*r12.x + r23.y*r12.y)*dsindC.z
          +(2.0*r23.z*r12.x - r12.z*r23.x)*dsindC.x
          +(2.0*r23.z*r12.y - r12.z*r23.y)*dsindC.y
          +dsindB.y*r34.x - dsindB.x*r34.y);
    }
  }    // end loop over multiplicity
  force1 += f1;
  force2 += f2-f1;
  force3 += f3-f2;
  force4 += -f3;
}
//==============================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// Force calculations for 1-4 atoms
//    calculate <force1> and <force2> given other info
//==============================================================================
  void Physics::calc14Forces
(
 const Vector& pos1,     // 1st atom's position
 const Vector& pos2,     // 2nd atom's position
 Vector&       force1,   // 1st atom's resultant force
 Vector&       force2,   // 2nd atom's resultant force
 double        q1,
 double        q2,
 int           ityp,
 StepOutput*   out
 )
  //==============================================================================
{
  double* pvten     = out->pt.pvten;    // Pressure tensor [9]
  double* pvten_tot = out->pt.pvten_tot;// Pressure tensor [9]

  // pvten - 1 because piny wants arrays to start from 1
  onfo (pos1, pos2, force1, force2, q1, q2,
      pvten-1, pvten_tot-1, ityp, out);

}
//==============================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// Force calculations excluded electrostatic interactions under 3D periodicity
//    calculate <force1> and <force2> given other info
//==============================================================================
  double Physics::calcEcorrForces
(
 const Vector& pos1,     // 1st atom's position
 const Vector& pos2,     // 2nd atom's position
 Vector&       force1,   // 1st atom's resultant force
 Vector&       force2,   // 2nd atom's resultant force
 double        charge1,  
 double        charge2,
 double        alp_ewd   //  Ewald parameter : >= 3.5/cutoff
 )
  //==============================================================================
{
  Vector r12  = pos1 - pos2;
  double r    = r12.Mag();


  Vector f12;
  double pot = 0;

  force1     += f12;
  force2     -= f12;

  return pot;
}
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// find the maximum exclusion distance squared
//==============================================================================
double Physics::get_r2max_excl (int     nAtms, // num atoms
    double* xPos,  // size [nAtms]
    double* yPos,  // size [nAtms]
    double* zPos)  // size [nAtms]
  //==============================================================================                              
{
  // arrays should start from 1
  xPos--; yPos--; zPos--;
  int* num   = readonly_mdintra.mdexcl.num;
  int* j     = readonly_mdintra.mdexcl.j;
  int* j_off = readonly_mdintra.mdexcl.j_off;
  int nAtoms = readonly_mdintra.mdexcl.natm_tot;

  double r2max_excl   = 0;
  double r2curr_excl  = 0;

  const double R2MAX_EXCL_CONST = 1.0;

  for (int i=1; i<=nAtoms; i++) {
    for (int k=1; k<=num[i]; k++) {
      int index = j [j_off [i] + k];

      r2curr_excl = (xPos[i] - xPos[index])*(xPos[i] - xPos[index]) + 
        (yPos[i] - yPos[index])*(yPos[i] - yPos[index]) + 
        (zPos[i] - zPos[index])*(zPos[i] - zPos[index]);

      if (r2curr_excl > r2max_excl) {
        r2max_excl = r2curr_excl;
      }
      //ckout << "id =" << i << ",j_off=" << j_off[i] << ",num=" << num [i] << ",j=" << j[j_off[i]+k] << endl;
    }
  }

  r2max_excl += R2MAX_EXCL_CONST;

  return r2max_excl;
}
//==============================================================================


BasePMECellData* Physics::newPMECellData () {
  return new PMECellData();
}

BasePMEGSpaceData* Physics::newPMEGSpaceData () {
  return new PMEGSpaceData();
}

void Physics::calcGridSize (double* _hmati, double _alpha_ewald,
    int& _ngrid_a, int& _ngrid_b, int& _ngrid_c) {
  PMEGSpaceData::calcGridSize (_hmati, _alpha_ewald, _ngrid_a, _ngrid_b, 
      _ngrid_c);
}
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// Call the appropriate constrain routines
// 
//==============================================================================
//  doShake()  shake constraints amongst numAtoms 
//==============================================================================
//
//  on one process these NEVER change.
//  in parallel, when the patch changes, these puppies follow the atoms
//  in parallel, assume all atoms in the line, triangle and tet arrive together
//
//  num_21_con
//  ind_atm1_21_con
//  ind_atm2_21_con
//  eq_21_con
//
//  loop over lines in this patch
//     shake_21()
//  end
//
//  num_33_con
//  ind_atm1_33_con
//  ind_atm2_33_con
//  eq1_33_con
//  eq2_33_con
//  eq3_33_con
//
//  loop over triagles in this patch
//     shake_33()
//  end
//
//  num_46_con
//  ind_atm1_46_con
//  ind_atm2_46_con
//  eq1_46_con
//  eq2_46_con
//  eq3_46_con
//  eq4_46_con
//  eq5_46_con
//  eq6_46_con
//
//  loop over tets in this patch
//     shake_46()
//  end
//
//==============================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// Call the appropriate constrain routines
// 
//==============================================================================
//  doRattle()  shake constraints amongst numAtoms 
//==============================================================================
//
//  on one process these NEVER change.
//  in parallel, when the patch changes, these puppies follow the atoms
//  in parallel, assume all atoms in the line, triangle and tet arrive together.
//
//  num_21_con
//  ind_atm1_21_con
//  ind_atm2_21_con
//
//  loop over lines in this patch
//     rattle_21()
//  end
//
//  num_33_con
//  ind_atm1_33_con
//  ind_atm2_33_con
//
//  loop over triagles in this patch
//     rattle_33()
//  end
//
//  num_46_con
//  ind_atm1_46_con
//  ind_atm2_46_con
//
//  loop over tets in this patch
//     rattle_46()
//  end
//
//==============================================================================
