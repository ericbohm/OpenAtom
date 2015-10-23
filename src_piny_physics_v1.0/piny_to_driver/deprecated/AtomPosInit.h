#ifndef _ATOMPOSINIT_H_
#define _ATOMPOSINIT_H_

#include "../pentium_par/standard_include.h"

#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_mdintegrate.h"
#include "../class_defs/allclass_mdatoms.h"
#include "../class_defs/allclass_mdinter.h"
#include "../class_defs/allclass_mdintra.h"

#include "../class_defs/ATOMS/class_mdclatoms_pos.h"
#include "../class_defs/INTEGRATE/class_mdtherm_pos.h"

#include "../proto_defs/proto_parse_entry.h"
#include "../proto_defs/proto_coords_local.h"


class AtomPosInit {
  public:
    AtomPosInit ();
    ~AtomPosInit ();

    MDCLATOMS_POS* mdclatoms_pos;
    MDTHERM_POS    therm_class;
    MDTHERM_POS*   therm_bead;
    int            pi_beads;
    int            ip;
    int            natm_tot;
    int            iextended_on;
    int            num_nhc;
    int            len_nhc;

};

#endif // _ATOMPOSINIT_H_
