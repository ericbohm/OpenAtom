#ifndef _GROUP_BOND_CON_21_H_
#define _GROUP_BOND_CON_21_H_

#include "leanMD.h"

void shake_21(
    Vector& pos1,     // 1st atom's position
    Vector& pos2,     // 2nd atom's position
    Vector& vel1,     // 1st atom's velocity
    Vector& vel2,     // 2nd atom's velocity
    const Vector& pos1_old, // 1st atom's old positions
    const Vector& pos2_old, // 2nd atom's old positions
    double mass1,           // 1st atom's mass
    double mass2,           // 2nd atom's mass
    double eq,              // bond length of r12
    double dt               // time step
    );


void rattle_21(
    const Vector& pos1,     // 1st atom's position
    const Vector& pos2,     // 2nd atom's position
    Vector& vel1,     // 1st atom's velocity
    Vector& vel2,     // 2nd atom's velocity
    double mass1,           // 1st atom's mass
    double mass2            // 2nd atom's mass
    );

#endif
