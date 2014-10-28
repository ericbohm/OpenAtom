#ifndef _BOND_
#define _BOND_

#include "leanMD.h"

#include "../include/class_defs/allclass_gen.h"
#include "../include/class_defs/allclass_mdintra.h"

void bond(const Vector& pos1,                   // 1st atom's position
    const Vector& pos2,                   // 2nd atom's position
    Vector&       force1,                 // 1st atom's force
    Vector&       force2,                 // 2nd atom's force
    double*       pvten,                  // Pressure tensor [9]
    double*       pvten_tot,              // Pressure tensor [9]
    int           ityp,                   // bend typ between atoms   
    StepOutput*   out,
    int           ires_bond = 1);
#endif
