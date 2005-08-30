//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                          
//==========================================================================

#ifndef _PhysicsAtomPosInit_
#define _PhysicsAtomPosInit_

#include "../../../include/Atoms.h"
#include "../class_defs/ATOMS/class_mdclatoms_pos.h"
#include "../class_defs/ATOMS/class_mdclatoms_pos.h"
#include "../class_defs/INTEGRATE/class_mdtherm_pos.h"

class PhysicsAtomPosInit{

 //---------------------------------------------------------------------------
 public:
     PhysicsAtomPosInit ();
    ~PhysicsAtomPosInit ();

    MDCLATOMS_POS* mdclatoms_pos;
    MDTHERM_POS    therm_class;
    MDTHERM_POS*   therm_bead;
    int            pi_beads;
    int            natm_tot;
    int            natm_nl;
    int            iextended_on;
    int            num_nhc;
    int            len_nhc;

 //---------------------------------------------------------------------------
 // functions

    void DriverAtomInit (Atom *);

};
#endif
