#include "group_bond_con.h"

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// double doShake
//       Shake API
//       Given a set of atom positions, velocities, Shake onto the surface
//       of constraint
//==========================================================================

void doShake (
    int     numAtoms,
    Vector* position,      // size [numAtoms] positions of the atoms
    Vector* position_old,  // size [numAtoms] positions of the atoms
    Vector* velocity,      // size [numAtoms] velocities for each atom
    int*    atomID,
    double  dt,            // time step  
    RigidBondInfo* rBonds
    ) 
  //==========================================================================
{ // begin  doShake
  //==========================================================================
  //  Line constraints

  int      ncon_21;
  int*     ind_con21[3];
  int*     jtyp_21;
  double** eq21 = readonly_mdintra.mdgrp_bond_con.eq_21;

  int      ncon_33;
  int*     ind_con33[4];
  int*     jtyp_33;
  double** eq33 = readonly_mdintra.mdgrp_bond_con.eq_33;

  int      ncon_46;
  int*     ind_con46[5];
  int*     jtyp_46;
  double** eq46 = readonly_mdintra.mdgrp_bond_con.eq_46;

  double*  mass;

  double   tolerance;
  int   max_iter;

  // initialize 21 constraints
  ncon_21       = rBonds->num_21;

  // starts from index '0' for 2nd dim and index '1' for first dim
  ind_con21 [1] = rBonds->j1_21;
  ind_con21 [2] = rBonds->j2_21;
  jtyp_21       = rBonds->jtyp_21;

  // initialize 33 constraints
  ncon_33       = rBonds->num_33;

  // starts from index '0' for 2nd dim and index '1' for first dim
  ind_con33 [1] = rBonds->j1_33;
  ind_con33 [2] = rBonds->j2_33;
  ind_con33 [3] = rBonds->j3_33;
  jtyp_33       = rBonds->jtyp_33;

  // initialize 46 constraints
  ncon_46       = rBonds->num_46;

  // starts from index '0' for 2nd dim and index '1' for first dim
  ind_con46 [1] = rBonds->j1_46;
  ind_con46 [2] = rBonds->j2_46;
  ind_con46 [3] = rBonds->j3_46;
  ind_con46 [4] = rBonds->j4_46;
  jtyp_46       = rBonds->jtyp_46;

  // initialize mass, starts from index [1]
  mass          = readonly_mdatoms.mdclatoms_info.mass;

  tolerance     = readonly_mdintra.mdconstrnt.tolshake;
  max_iter      = readonly_mdintra.mdconstrnt.max_iter;

  // main function body
  for(int i=0;i<ncon_21;i++){
    shake_21(position[ind_con21[1][i]],
        position[ind_con21[2][i]],
        velocity[ind_con21[1][i]],
        velocity[ind_con21[2][i]],
        position_old[ind_con21[1][i]],
        position_old[ind_con21[2][i]],
        mass[atomID[ind_con21[1][i]]],
        mass[atomID[ind_con21[2][i]]],
        eq21[1][jtyp_21[i]],dt);
  }/*endfor*/

  //==========================================================================
  //  Triangle constraints

  for(int i=0;i<ncon_33;i++){
    shake_33(position[ind_con33[1][i]],
        position[ind_con33[2][i]],
        position[ind_con33[3][i]],
        velocity[ind_con33[1][i]],
        velocity[ind_con33[2][i]],
        velocity[ind_con33[3][i]],
        position_old[ind_con33[1][i]],
        position_old[ind_con33[2][i]],
        position_old[ind_con33[3][i]],
        mass[atomID[ind_con33[1][i]]],
        mass[atomID[ind_con33[2][i]]],
        mass[atomID[ind_con33[3][i]]],
        eq33[1][jtyp_33[i]],
        eq33[2][jtyp_33[i]],
        eq33[3][jtyp_33[i]],
        dt,
        tolerance,
        max_iter);
  }/*endfor*/

  //==========================================================================
  //  Tet constraints

  for(int i=0;i<ncon_46;i++){
    shake_46(position[ind_con46[1][i]],
        position[ind_con46[2][i]],
        position[ind_con46[3][i]],
        position[ind_con46[4][i]],
        velocity[ind_con46[1][i]],
        velocity[ind_con46[2][i]],
        velocity[ind_con46[3][i]],
        velocity[ind_con46[4][i]],
        position_old[ind_con46[1][i]],
        position_old[ind_con46[2][i]],
        position_old[ind_con46[3][i]],
        position_old[ind_con46[4][i]],
        mass[atomID[ind_con46[1][i]]],
        mass[atomID[ind_con46[2][i]]],
        mass[atomID[ind_con46[3][i]]],
        mass[atomID[ind_con46[4][i]]],
        eq46[1][jtyp_46[i]],
        eq46[2][jtyp_46[i]],
        eq46[3][jtyp_46[i]],
        eq46[4][jtyp_46[i]],
        eq46[5][jtyp_46[i]],
        eq46[6][jtyp_46[i]],
        dt,
        tolerance,
        max_iter);
  }/*endfor*/

  //==========================================================================
} // end  doShake
//==========================================================================




//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
// double doRattle
//       Shake API
//       Given a set of atom positions, velocities, Shake onto the surface
//       of constraint
//==========================================================================

void doRattle (
    int     numAtoms,
    Vector* position,      // size [numAtoms] positions of the atoms
    Vector* position_old,  // size [numAtoms] positions of the atoms
    Vector* velocity,      // size [numAtoms] velocities for each atom
    int*    atomID,
    double  dt,            // time step  
    RigidBondInfo* rBonds
    ) 
  //==========================================================================
{ // begin  doRattle
  //==========================================================================
  //  Line constraints

  int      ncon_21;
  int*     ind_con21[3];
  int*     jtyp_21;

  int      ncon_33;
  int*     ind_con33[4];
  int*     jtyp_33;

  int      ncon_46;
  int*     ind_con46[5];
  int*     jtyp_46;

  double*  mass;

  // initialize 21 constraints
  ncon_21       = rBonds->num_21;

  // starts from index '0' for 2nd dim and index '1' for first dim
  ind_con21 [1] = rBonds->j1_21;
  ind_con21 [2] = rBonds->j2_21;
  jtyp_21       = rBonds->jtyp_21;

  // initialize 33 constraints
  ncon_33       = rBonds->num_33;

  // starts from index '0' for 2nd dim and index '1' for first dim
  ind_con33 [1] = rBonds->j1_33;
  ind_con33 [2] = rBonds->j2_33;
  ind_con33 [3] = rBonds->j3_33;
  jtyp_33       = rBonds->jtyp_33;

  // initialize 46 constraints
  ncon_46       = rBonds->num_46;

  // starts from index '0' for 2nd dim and index '1' for first dim
  ind_con46 [1] = rBonds->j1_46;
  ind_con46 [2] = rBonds->j2_46;
  ind_con46 [3] = rBonds->j3_46;
  ind_con46 [4] = rBonds->j4_46;
  jtyp_46       = rBonds->jtyp_46;   

  // initialize mass, starts from index [1]
  mass          = readonly_mdatoms.mdclatoms_info.mass;

  // main function body
  for(int i=0;i<ncon_21;i++){
    rattle_21(position[ind_con21[1][i]],
        position[ind_con21[2][i]],
        velocity[ind_con21[1][i]],
        velocity[ind_con21[2][i]],
        //                position_old[ind_con21[1][i]],
        //                position_old[ind_con21[2][i]],
        mass[atomID[ind_con21[1][i]]],
        mass[atomID[ind_con21[2][i]]]);
  }/*endfor*/

  //==========================================================================
  //  Triangle constraints

  for(int i=0;i<ncon_33;i++){
    rattle_33(position[ind_con33[1][i]],
        position[ind_con33[2][i]],
        position[ind_con33[3][i]],
        velocity[ind_con33[1][i]],
        velocity[ind_con33[2][i]],
        velocity[ind_con33[3][i]],
        mass[atomID[ind_con33[1][i]]],
        mass[atomID[ind_con33[2][i]]],
        mass[atomID[ind_con33[3][i]]]);
  }/*endfor*/

  //==========================================================================
  //  Tet constraints

  for(int i=0;i<ncon_46;i++){
    rattle_46(position[ind_con46[1][i]],
        position[ind_con46[2][i]],
        position[ind_con46[3][i]],
        position[ind_con46[4][i]],
        velocity[ind_con46[1][i]],
        velocity[ind_con46[2][i]],
        velocity[ind_con46[3][i]],
        velocity[ind_con46[4][i]],
        mass[atomID[ind_con46[1][i]]],
        mass[atomID[ind_con46[2][i]]],
        mass[atomID[ind_con46[3][i]]],
        mass[atomID[ind_con46[4][i]]]);
  }/*endfor*/

  //==========================================================================
} // end  doRattle
//==========================================================================
