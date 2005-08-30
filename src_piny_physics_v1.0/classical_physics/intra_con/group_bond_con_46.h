#ifndef _GROUP_BOND_CON_46_H_
#define _GROUP_BOND_CON_46_H_

#include "leanMD.h"
#include "../mathlib/proto_math.h"

void shake_46(
  Vector& pos1,     // 1st atom's position
  Vector& pos2,     // 2nd atom's position
  Vector& pos3,     // 3rd atom's position
  Vector& pos4,     // 4th atom's position
  Vector& vel1,     // 1st atom's velocity
  Vector& vel2,     // 2nd atom's velocity
  Vector& vel3,     // 3rd atom's velocity
  Vector& vel4,     // 4th atom's velocity
  const Vector& pos1_old, // 1st atom's old position
  const Vector& pos2_old, // 2nd atom's old position
  const Vector& pos3_old, // 3rd atom's old position
  const Vector& pos4_old, // 4th atom's old position
  double mass1,           // 1st atom's mass
  double mass2,           // 2nd atom's mass
  double mass3,           // 3rd atom's mass
  double mass4,           // 4th atom's mass
  double eq1,             // bond length of 1-2 bond
  double eq2,             // bond length of 1-3 bond
  double eq3,             // bond length of 2-3 bond
  double eq4,             // bond length of 2-3 bond
  double eq5,             // bond length of 2-3 bond
  double eq6,             // bond length of 2-3 bond
  double dt,              // time step
  double tol,             // tolerence
  int max_iter        // maximum iterations
 );

void rattle_46(
  const Vector& pos1,     // 1st atom's position
  const Vector& pos2,     // 2nd atom's position
  const Vector& pos3,     // 3rd atom's position
  const Vector& pos4,     // 4th atom's position
  Vector& vel1,     // 1st atom's velocity
  Vector& vel2,     // 2nd atom's velocity
  Vector& vel3,     // 3rd atom's velocity
  Vector& vel4,     // 4th atom's velocity
  double mass1,           // 1st atom's mass
  double mass2,           // 2nd atom's mass
  double mass3,           // 3rd atom's mass
  double mass4            // 4th atom's mass
 );

#endif
