//===========================================================================

#ifndef _INTERFACE_CTRL_H_
#define _INTERFACE_CTRL_H_

#include "standard_include.h"

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

class Interface_ctrl {
  public:
    Interface_ctrl(char* fname, CkCallback c);

  private:
    // default constructor cannot be used to instantiate Reader object
    Interface_ctrl () {}
};

#endif // _INTERFACE_CTRL_H_


//===========================================================================
