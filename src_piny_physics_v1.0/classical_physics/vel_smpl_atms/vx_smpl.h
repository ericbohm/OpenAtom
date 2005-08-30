//=================================================================== 
//=================================================================== 
//                                                                    
// Class for velocity resampling                            
//                                                                    
//=================================================================== 
//=================================================================== 

#ifndef _VX_SMPL_
#define _VX_SMPL_

#include "standard_include.h"

#include "../../mathlib/proto_math.h"

#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_mdintegrate.h"
#include "../class_defs/allclass_mdintra.h"

class vx_smpl{
  public: 

     vx_smpl() {printf("constructor no arg \n");}

     vx_smpl(double*        velx,       // x-comp of atom velocities
             double*        vely,       // y-comp of atom velocities
             double*        velz,       // z-comp of atom velocities
             double*        mass,       // Array containting atom masses
             double*        text_atm,   // Temperature of each atom 
             int            natm,       // Number of atoms in system
             GENENSOPTS*    genensopts, // Ensemble options
             MDVEL_SAMP*    vel_samp,   // contains variables for vel sampling
             MDGHOST_ATOMS* ghost_atoms,// contains ghost atom info
             MDCONSTRNT*    mdconstrnt);// contains info about constraints

   void sampl_vx(double* velx,              // x-comp of atom velocities
                 double* vely,              // y-comp of atom velocities
                 double* velz,              // z-comp of atom velocities
                 double* mass,
                 double* text_atm,
                 int     natm_tot,
                 int*    iseed,
                 int*    iseed2,
                 double* qseed);
  private:
  
};

#endif  
