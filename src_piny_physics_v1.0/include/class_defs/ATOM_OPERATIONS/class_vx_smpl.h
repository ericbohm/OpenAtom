//=================================================================== 
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=================================================================== 
//                                                                    
// Class for velocity resampling                            
//                                                                    
//=================================================================== 

#ifndef _VX_SMPL_
#define _VX_SMPL_

#include "standard_include.h"
#include "../../mathlib/proto_math.h"
#include "../../../include/Atoms.h"

//=================================================================== 
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=================================================================== 
class VX_SMPL{
  public: 
    VX_SMPL() {}
    ~VX_SMPL() {}
    static void ctrlSamplAtomVel(
        int     natm,       // Number of atoms in system
        double* velx,       // x-comp of atom velocities
        double* vely,       // y-comp of atom velocities
        double* velz,       // z-comp of atom velocities
        double* mass        // Array containting atom masses
        );
    static void ctrlSamplAtomVel(
        int     natm,       // Number of atoms in system
        Atom*   atoms  
        );
    static void ctrlSamplAtomNhcVel(
        int      natm,        // Number of atoms in system
        AtomNHC *atomsNHC
        );
    static void sampl3DVelMultiT(
        int     natm,
        double* velx,  
        double* vely,  
        double* velz,  
        double* mass,
        double* text_atm,
        long*    iseed,
        long*    iseed2,
        double* qseed
        );
    static void sampl3DVelMultiT(
        int     natm,
        Atom *atoms,
        double* text_atm,
        long*    iseed,
        long*    iseed2,
        double* qseed
        );
    static void sampl3DVelOneT(
        int natm, 
        double* vx, 
        double* vy,
        double* vz,
        double* mass,
        double  kT,
        long* iseed,
        long* iseed2,
        double* qseed
        );
};
//=================================================================== 
#endif  
