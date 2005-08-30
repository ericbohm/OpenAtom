#ifndef _INTERROUTINES_H_
#define _INTERROUTINES_H_

#include "leanMD.h"

#include "Physics_defines.h"

class InterRoutines {
  public:
    static void calc_force_explicit
    (
      const Vector& position1, // 1st atom's position
      Vector&       force1,    // 1st atom's force (will be calculated)
      double        charge1,   // 1st atom's charge
      double        epsilon1,  // 1st atom's well depth
      double        sigma1,    // 1st atom's zero energy pt.
      double        cut1,      // 1st atom's cutoff
      int           ind1_excl, // 1st atom's index in exclusion list
      const Vector& position2, // 2nd atom's position
      Vector&       force2,    // 2nd atom's force (will be calculated)
      double        charge2,   // 2nd atom's charge
      double        epsilon2,  // 2nd atom's well depth
      double        sigma2,    // 2nd atom's zero energy pt.
      double        cut2,      // 2nd atom's cutoff
      int           ind2_excl, // 2nd atom's index in exclusion list
      int&          processed, // flag
      int           iperd,     // periodicity (0 or 3 for now)
      double        alp_ewd,   // ewald sum convergence parameter >= 3.5/cutoff
      double        r2max_excl,// maximum exclusion distance squared  
      StepOutput*   out
    );

    static void calc_force_spline
    (
      const Vector& position1,  // 1st atom's position
      Vector&       force1,     // 1st atom's force (will be calculated)
      double        charge1,    // 1st atom's charge
      double        cut1,       // 1st atom's cutoff
      int           ind1_excl,  // 1st atom's index in exclusion list
      const Vector& position2,  // 2nd atom's position
      Vector&       force2,     // 2nd atom's force (will be calculated)
      double        charge2,    // 2nd atom's charge
      double        cut2,       // 2nd atom's cutoff
      int           ind2_excl,  // 2nd atom's index in exclusion list
      int&          processed,  // flag
      int           iperd,      // periodicity (0 or 3 for now)
      double        r2max_excl, // maximum exclusion distance squared  
      StepOutput*   out
    );

  private:
    
    static double spl_value
    (
      int index0,          // index into look up table
      double del,          // r-rspl[index]
      double swit,         // switching function
      double *spl_data     // look up table
    );

    static void spl_params
    (
      int nsplin,         // number of pts in look up table
      int nsplin_m2,      // number of pts in look up minus 2
      int inter,          // interaction number
      double r,           // distance between particles 1-2
      double *rmin_spl,   // minimum distance stored look up table 
      double *dr_spl,     // window size of look up table rspl[i+1]-rspl[i] all i
      int &ispl,          // index of r in look up table of this interaction
      double &del_r,      // r-rpsl[ispl] 
      double &swit        // switching function to remove if statements
    );
};

#include "calc_forc_spline.h"
#include "calc_forc_explicit.h"

#endif

