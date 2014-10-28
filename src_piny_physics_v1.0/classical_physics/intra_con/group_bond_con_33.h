#ifndef _GROUP_BOND_CON_33_H_
#define _GROUP_BOND_CON_33_H_

#include "leanMD.h"

void shake_33(
    Vector& pos1,     // 1st atom's position
    Vector& pos2,     // 2nd atom's position
    Vector& pos3,     // 3rd atom's position
    Vector& vel1,     // 1st atom's velocity
    Vector& vel2,     // 2nd atom's velocity
    Vector& vel3,     // 3rd atom's velocity
    const Vector& pos1_old, // 1st atom's old positions
    const Vector& pos2_old, // 2nd atom's old positions
    const Vector& pos3_old, // 3rd atom's old positions
    double mass1,           // 1st atom's mass
    double mass2,           // 2nd atom's mass
    double mass3,           // 3rd atom's mass
    double eq1,             // bond length of 1-2 bond
    double eq2,             // bond length of 1-3 bond
    double eq3,             // bond length of 2-3 bond
    double dt,              // time step
    double tol,             // tolerence
    int max_iter        // maximum iterations
    );

void rattle_33(
    const Vector& pos1,     // 1st atom's position
    const Vector& pos2,     // 2nd atom's position
    const Vector& pos3,     // 3rd atom's position
    Vector& vel1,     // 1st atom's velocity
    Vector& vel2,     // 2nd atom's velocity
    Vector& vel3,     // 3rd atom's velocity
    double mass1,           // 1st atom's mass
    double mass2,           // 2nd atom's mass
    double mass3            // 3rd atom's mass
    );

#endif
