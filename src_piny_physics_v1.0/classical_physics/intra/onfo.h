#ifndef _ONFO_
#define _ONFO_

#include "leanMD.h"

#include "../include/class_defs/allclass_gen.h"
#include "../include/class_defs/allclass_mdintra.h"

void onfo (const Vector& pos1,      // 1st atom position
           const Vector& pos2,      // 2nd atom position
           Vector&       force1,    // 1st atom force
           Vector&       force2,    // 2nd atom force
           double        q1,        // charge 1st atom
           double        q2,        // charge 2nd atom
           double*       pvten,
           double*       pvten_tot,
           int           ityp,
           StepOutput*   out);

#endif
