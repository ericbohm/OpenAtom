#include "group_bond_con_21.h"

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//  Constraints between 2-atoms not involved in any other constraint : line
//         Example : X-O-H  where O-H bond in constrained but X-O is not
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

void shake_21(
    Vector& pos1,     // 1st atom's position
    Vector& pos2,     // 2nd atom's position
    Vector& vel1,     // 1st atom's velocity
    Vector& vel2,     // 2nd atom's velocity
    const Vector& pos1_old, // 1st atom's old positions
    const Vector& pos2_old, // 2nd atom's old positions
    double mass1,           // 1st atom's mass
    double mass2,           // 2nd atom's mass
    double eq,              // bond length of r12
    double dt               // time step
    )

  //==========================================================================
{//Begin Routine
  //=======================================================================
  // get positions and masses 

  const double dx12  = pos1.x-pos2.x;
  const double dy12  = pos1.y-pos2.y;
  const double dz12  = pos1.z-pos2.z;

  const double dxo12 = pos1_old.x-pos2_old.x;
  const double dyo12 = pos1_old.y-pos2_old.y;
  const double dzo12 = pos1_old.z-pos2_old.z;

  const double rm1   = 1.0/mass1;
  const double rm2   = 1.0/mass2;

  // =============================================================================== 
  // determine value of multiplier using quadratic equation 
  // using the positive root/physical solution

  const double rmu12    = rm1 + rm2;
  const double r12osq   = dxo12*dxo12 + dyo12*dyo12 + dzo12*dzo12;
  const double r12sq    = dx12 *dx12  + dy12 *dy12  + dz12 *dz12;
  const double r12_r12o = dx12 *dxo12 + dy12 *dyo12 + dz12 *dzo12;

  const double a  = r12osq*rmu12*rmu12;
  const double b  = 2.0*r12_r12o*rmu12;
  const double c  = r12sq - eq*eq; 
  const double dd = b*b - 4.0*a*c;

  const double lam1 = (b - sqrt(dd))/(2.0*a);

  // =============================================================================== 
  // position and velocity update 

  const double fx = lam1*dxo12;
  const double fy = lam1*dyo12;
  const double fz = lam1*dzo12;

  pos1.x -= fx*rm1;
  pos1.y -= fy*rm1;
  pos1.z -= fz*rm1;
  vel1.x -= fx*rm1/dt;
  vel1.y -= fy*rm1/dt;
  vel1.z -= fz*rm1/dt;

  pos2.x += fx*rm2;
  pos2.y += fy*rm2;
  pos2.z += fz*rm2;
  vel2.x += fx*rm2/dt;
  vel2.y += fy*rm2/dt;
  vel2.z += fz*rm2/dt;

  //=======================================================================
} // end routine 
//=======================================================================



//=======================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

void rattle_21(
    const Vector& pos1,     // 1st atom's position
    const Vector& pos2,     // 2nd atom's position
    Vector& vel1,     // 1st atom's velocity
    Vector& vel2,     // 2nd atom's velocity
    double mass1,           // 1st atom's mass
    double mass2           // 2nd atom's mass
    )

  //==========================================================================
  //        Begin Routine                                                     
{// Begin routine 
  //=======================================================================
  // masses, positions and velocities

  const double dx12  = pos1.x-pos2.x;
  const double dy12  = pos1.y-pos2.y;
  const double dz12  = pos1.z-pos2.z;

  const double dvx12 = vel1.x-vel2.x;
  const double dvy12 = vel1.y-vel2.y;
  const double dvz12 = vel1.z-vel2.z;

  const double rm1   = 1.0/mass1;
  const double rm2   = 1.0/mass2;

  //=======================================================================
  // compute lagrange multiplier 

  const double r12sq   = dx12*dx12  + dy12*dy12  + dz12*dz12;
  const double r12_v12 = dx12*dvx12 + dy12*dvy12 + dz12*dvz12;
  const double rmu12   = rm1 + rm2;

  const double lam1    = r12_v12/(rmu12*r12sq);

  //=======================================================================
  // velocity update 

  const double fx = lam1*dx12;
  const double fy = lam1*dy12;
  const double fz = lam1*dz12;

  vel1.x -= fx*rm1;
  vel1.y -= fy*rm1;
  vel1.z -= fz*rm1;

  vel2.x += fx*rm1;
  vel2.y += fy*rm1;
  vel2.z += fz*rm1;

  //=======================================================================
} // end routine 
//=======================================================================
