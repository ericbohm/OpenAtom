#ifndef _BEND_
#define _BEND_

#include "leanMD.h"

#include "../include/class_defs/allclass_gen.h"
#include "../include/class_defs/allclass_mdintra.h"

void bend (const Vector& pos1,      // 1st atom's position
           const Vector& pos2,      // 2nd atom's position
           const Vector& pos3,      // 3rd atom's position
           Vector&       force1,    // 1st atom's force
           Vector&       force2,    // 2nd atom's force
           Vector&       force3,    // 3rd atom's force
           double*       pvten,     // Pressure tensor 
           double*       pvten_tot, // Pressure tensor 
           int           ityp,      // bend typ between atoms 
           StepOutput*   out);

#endif
