#ifndef _GROUP_BOND_CON_H_
#define _GROUP_BOND_CON_H_

#include "leanMD.h"
#include "RigidBondInfo.h"

#include "group_bond_con_21.h"
#include "group_bond_con_33.h"
#include "group_bond_con_46.h"

void doShake (
  int     numAtoms,
  Vector* position,      // size [numAtoms] positions of the atoms
  Vector* position_old,  // size [numAtoms] positions of the atoms
  Vector* velocity,      // size [numAtoms] velocities for each atom
  int*    atomID,
  double  dt,            // time step
  RigidBondInfo* rBonds
);

void doRattle (
  int     numAtoms,
  Vector* position,      // size [numAtoms] positions of the atoms
  Vector* position_old,  // size [numAtoms] positions of the atoms
  Vector* velocity,      // size [numAtoms] velocities for each atom
  int*    atomID,
  double  dt,            // time step
  RigidBondInfo* rBonds
);

#endif
