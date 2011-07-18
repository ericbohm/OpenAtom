//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#ifndef _PhysicsAtomPosInit_
#define _PhysicsAtomPosInit_

#include "../../../include/Atoms.h"
#include "../class_defs/ATOMS/class_mdclatoms_pos.h"
#include "../class_defs/ATOMS/class_mdclatoms_pos.h"
#include "../class_defs/INTEGRATE/class_mdtherm_pos.h"

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
class PhysicsAtomPosInit{
//---------------------------------------------------------------------------
 public:
     PhysicsAtomPosInit (int , int );
    ~PhysicsAtomPosInit ();
     MDCLATOMS_POS* mdclatoms_pos;
     MDTHERM_POS    therm_class;
     MDTHERM_POS*   therm_bead;
     int            pi_beads_true;
     int            ntemper;
     int            pi_beads;
     int            natm_tot;
     int            natm_nl;
     int            iextended_on;
     int            cp_min_opt;
     int            cp_wave_opt;
     int            num_nhc;
     int            len_nhc;
     int            istart_typ;
     int            isokin_opt;
     int            ibead;
     int            itemper;
     double         kT;
     double       **mass_nhc;
     void DriverAtomInit (int,Atom *,AtomNHC *,int, int);
};
//==========================================================================
#endif
