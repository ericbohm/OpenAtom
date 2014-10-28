#include "group_bond_con_33.h"

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//  Constraints between 3-atoms not involved in any other constraint : triangle
//     Example : X-CH2   : CH CH and HH bond constrained but not X-C
//==========================================================================



//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

void shake_33(
    Vector& pos1,     // 1st atom's position
    Vector& pos2,     // 2nd atom's position
    Vector& pos3,     // 3rd atom's position
    Vector& vel1,     // 1st atom's velocity
    Vector& vel2,     // 2nd atom's velocity
    Vector& vel3,     // 3rd atom's velocity
    const Vector& pos1_old, // 1st atom's old positions
    const Vector& pos2_old, // 2nd atom's old positions
    const Vector& pos3_old, // 3rd atom's old positions
    double mass1,           // 1st atom's mass
    double mass2,           // 2nd atom's mass
    double mass3,           // 3rd atom's mass
    double eq1,             // bond length of 1-2 bond
    double eq2,             // bond length of 1-3 bond
    double eq3,             // bond length of 2-3 bond
    double dt,              // time step
    double tol,             // tolerence
    int max_iter        // maximum iterations
    )
  //==========================================================================
{//Begin Routine
  //=======================================================================

  double dxt[4], dyt[4], dzt[4];
  double dx[4],  dy[4],  dz[4];
  double dxn[4], dyn[4], dzn[4];
  double amat[4][4], avec[4], ainv[4][4];
  double lam[4], lam_old[4];
  double det, rdet;

  //=======================================================================
  // get positions, masses and solution vector of linearized equations

  const double rms1 = 1.0/mass1;
  const double rms2 = 1.0/mass2;
  const double rms3 = 1.0/mass3;

  const double rmm11 = -(rms1+rms2); 
  const double rmm12 = -rms1; 
  const double rmm13 =  rms2;
  const double rmm21 = -rms1; 
  const double rmm22 = -(rms1+rms3);
  const double rmm23 = -rms3;
  const double rmm31 =  rms2; 
  const double rmm32 = -rms3;
  const double rmm33 = -(rms2+rms3);

  dxt[1] = pos1.x - pos2.x;
  dyt[1] = pos1.y - pos2.y;
  dzt[1] = pos1.z - pos2.z;

  dxt[2] = pos1.x - pos3.x;
  dyt[2] = pos1.y - pos3.y;
  dzt[2] = pos1.z - pos3.z;

  dxt[3] = pos2.x - pos3.x;
  dyt[3] = pos2.y - pos3.y;
  dzt[3] = pos2.z - pos3.z;

  dx[1] = pos1_old.x - pos2_old.x;
  dy[1] = pos1_old.y - pos2_old.y;
  dz[1] = pos1_old.z - pos2_old.z;

  dx[2] = pos1_old.x - pos3_old.x;
  dy[2] = pos1_old.y - pos3_old.y;
  dz[2] = pos1_old.z - pos3_old.z;

  dx[3] = pos2_old.x - pos3_old.x;
  dy[3] = pos2_old.y - pos3_old.y;
  dz[3] = pos2_old.z - pos3_old.z;

  avec[1] = eq1*eq1 - (dxt[1]*dxt[1] + dyt[1]*dyt[1] + dzt[1]*dzt[1]);
  avec[2] = eq2*eq2 - (dxt[2]*dxt[2] + dyt[2]*dyt[2] + dzt[2]*dzt[2]);
  avec[3] = eq3*eq3 - (dxt[3]*dxt[3] + dyt[3]*dyt[3] + dzt[3]*dzt[3]);

  //=======================================================================
  // Get initial guess for lambda by solving linear eq: amat*lam = vec

  //--------------------------------------------------
  // construct matrix and vector at current positions

  amat[1][1] = 2.0*rmm11*(dxt[1]*dx[1] + dyt[1]*dy[1] + dzt[1]*dz[1]);
  amat[1][2] = 2.0*rmm12*(dxt[1]*dx[2] + dyt[1]*dy[2] + dzt[1]*dz[2]);
  amat[1][3] = 2.0*rmm13*(dxt[1]*dx[3] + dyt[1]*dy[3] + dzt[1]*dz[3]);

  amat[2][1] = 2.0*rmm21*(dxt[2]*dx[1] + dyt[2]*dy[1] + dzt[2]*dz[1]);
  amat[2][2] = 2.0*rmm22*(dxt[2]*dx[2] + dyt[2]*dy[2] + dzt[2]*dz[2]);
  amat[2][3] = 2.0*rmm23*(dxt[2]*dx[3] + dyt[2]*dy[3] + dzt[2]*dz[3]);

  amat[3][1] = 2.0*rmm31*(dxt[3]*dx[1] + dyt[3]*dy[1] + dzt[3]*dz[1]);
  amat[3][2] = 2.0*rmm32*(dxt[3]*dx[2] + dyt[3]*dy[2] + dzt[3]*dz[2]);
  amat[3][3] = 2.0*rmm33*(dxt[3]*dx[3] + dyt[3]*dy[3] + dzt[3]*dz[3]);

  //--------------------------------------------------
  // Get inverse of matrix 

  det = (amat[1][1] * (amat[2][2] * amat[3][3] - amat[3][2] * amat[2][3]) + 
      amat[2][1] * (amat[3][2] * amat[1][3] - amat[1][2] * amat[3][3]) + 
      amat[3][1] * (amat[1][2] * amat[2][3] - amat[2][2] * amat[1][3]));
  rdet = 1.0/det;
  ainv[1][1] = (amat[2][2] * amat[3][3] - amat[3][2] * amat[2][3])*rdet;
  ainv[2][2] = (amat[1][1] * amat[3][3] - amat[3][1] * amat[1][3])*rdet;
  ainv[3][3] = (amat[1][1] * amat[2][2] - amat[2][1] * amat[1][2])*rdet;
  ainv[2][1] = (amat[3][1] * amat[2][3] - amat[2][1] * amat[3][3])*rdet;
  ainv[1][2] = (amat[1][3] * amat[3][2] - amat[1][2] * amat[3][3])*rdet;
  ainv[3][1] = (amat[2][1] * amat[3][2] - amat[3][1] * amat[2][2])*rdet;
  ainv[1][3] = (amat[1][2] * amat[2][3] - amat[1][3] * amat[2][2])*rdet;
  ainv[3][2] = (amat[3][1] * amat[1][2] - amat[3][2] * amat[1][1])*rdet;
  ainv[2][3] = (amat[1][3] * amat[2][1] - amat[2][3] * amat[1][1])*rdet;

  //--------------------------------------------------
  // solve linear equations

  lam[1] = ainv[1][1]*avec[1] + ainv[1][2]*avec[2] + ainv[1][3]*avec[3];
  lam[2] = ainv[2][1]*avec[1] + ainv[2][2]*avec[2] + ainv[2][3]*avec[3];
  lam[3] = ainv[3][1]*avec[1] + ainv[3][2]*avec[2] + ainv[3][3]*avec[3];

  //=======================================================================
  // Iterative loop for multiplier based on initial guess

  double dl1, dl2,dl3;
  double dlmax = tol + 1.0;
  int iter     = 0;

  do {++iter;

    if(iter > max_iter) {
      CkPrintf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
      CkPrintf("Shake 3x3 not converged after %d iterations.\n",max_iter);
      CkPrintf("The present tolerance is %g \n",dlmax);
      CkPrintf("The desired tolerance is %g \n",tol);
      CkPrintf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
      break;
    }//endif

    //--------------------------------------------------
    // Set up new difference vectors 

    dxn[1]=2.0*dxt[1]+(rmm11*lam[1]*dx[1]+rmm12*lam[2]*dx[2]+rmm13*lam[3]*dx[3]);
    dyn[1]=2.0*dyt[1]+(rmm11*lam[1]*dy[1]+rmm12*lam[2]*dy[2]+rmm13*lam[3]*dy[3]);
    dzn[1]=2.0*dzt[1]+(rmm11*lam[1]*dz[1]+rmm12*lam[2]*dz[2]+rmm13*lam[3]*dz[3]);

    dxn[2]=2.0*dxt[2]+(rmm21*lam[1]*dx[1]+rmm22*lam[2]*dx[2]+rmm23*lam[3]*dx[3]);
    dyn[2]=2.0*dyt[2]+(rmm21*lam[1]*dy[1]+rmm22*lam[2]*dy[2]+rmm23*lam[3]*dy[3]);
    dzn[2]=2.0*dzt[2]+(rmm21*lam[1]*dz[1]+rmm22*lam[2]*dz[2]+rmm23*lam[3]*dz[3]);

    dxn[3]=2.0*dxt[3]+(rmm31*lam[1]*dx[1]+rmm32*lam[2]*dx[2]+rmm33*lam[3]*dx[3]);
    dyn[3]=2.0*dyt[3]+(rmm31*lam[1]*dy[1]+rmm32*lam[2]*dy[2]+rmm33*lam[3]*dy[3]);
    dzn[3]=2.0*dzt[3]+(rmm31*lam[1]*dz[1]+rmm32*lam[2]*dz[2]+rmm33*lam[3]*dz[3]);

    //--------------------------------------------------
    // Construct new A-matrix 

    amat[1][1] = rmm11*(dxn[1]*dx[1] + dyn[1]*dy[1] + dzn[1]*dz[1]);
    amat[1][2] = rmm12*(dxn[1]*dx[2] + dyn[1]*dy[2] + dzn[1]*dz[2]);
    amat[1][3] = rmm13*(dxn[1]*dx[3] + dyn[1]*dy[3] + dzn[1]*dz[3]);

    amat[2][1] = rmm21*(dxn[2]*dx[1] + dyn[2]*dy[1] + dzn[2]*dz[1]);
    amat[2][2] = rmm22*(dxn[2]*dx[2] + dyn[2]*dy[2] + dzn[2]*dz[2]);
    amat[2][3] = rmm23*(dxn[2]*dx[3] + dyn[2]*dy[3] + dzn[2]*dz[3]);

    amat[3][1] = rmm31*(dxn[3]*dx[1] + dyn[3]*dy[1] + dzn[3]*dz[1]);
    amat[3][2] = rmm32*(dxn[3]*dx[2] + dyn[3]*dy[2] + dzn[3]*dz[2]);
    amat[3][3] = rmm33*(dxn[3]*dx[3] + dyn[3]*dy[3] + dzn[3]*dz[3]);

    //--------------------------------------------------
    // Get inverse of matrix 

    det = (amat[1][1] * (amat[2][2] * amat[3][3] - amat[3][2] * amat[2][3]) + 
        amat[2][1] * (amat[3][2] * amat[1][3] - amat[1][2] * amat[3][3]) + 
        amat[3][1] * (amat[1][2] * amat[2][3] - amat[2][2] * amat[1][3]));
    rdet = 1.0/det;
    ainv[1][1] = (amat[2][2] * amat[3][3] - amat[3][2] * amat[2][3])*rdet;
    ainv[2][2] = (amat[1][1] * amat[3][3] - amat[3][1] * amat[1][3])*rdet;
    ainv[3][3] = (amat[1][1] * amat[2][2] - amat[2][1] * amat[1][2])*rdet;
    ainv[2][1] = (amat[3][1] * amat[2][3] - amat[2][1] * amat[3][3])*rdet;
    ainv[1][2] = (amat[1][3] * amat[3][2] - amat[1][2] * amat[3][3])*rdet;
    ainv[3][1] = (amat[2][1] * amat[3][2] - amat[3][1] * amat[2][2])*rdet;
    ainv[1][3] = (amat[1][2] * amat[2][3] - amat[1][3] * amat[2][2])*rdet;
    ainv[3][2] = (amat[3][1] * amat[1][2] - amat[3][2] * amat[1][1])*rdet;
    ainv[2][3] = (amat[1][3] * amat[2][1] - amat[2][3] * amat[1][1])*rdet;

    lam_old[1] = lam[1];
    lam_old[2] = lam[2];
    lam_old[3] = lam[3];

    lam[1] = ainv[1][1]*avec[1] + ainv[1][2]*avec[2] + ainv[1][3]*avec[3];
    lam[2] = ainv[2][1]*avec[1] + ainv[2][2]*avec[2] + ainv[2][3]*avec[3];
    lam[3] = ainv[3][1]*avec[1] + ainv[3][2]*avec[2] + ainv[3][3]*avec[3];

    dl1   = fabs(lam[1]-lam_old[1]);
    dl2   = fabs(lam[2]-lam_old[2]);
    dl3   = fabs(lam[3]-lam_old[3]);
    dlmax = MAX3(dl1,dl2,dl3);

  } while(dlmax > tol);

  //=======================================================================
  // position update 

  const double fx1 = -(lam[1]*dx[1] + lam[2]*dx[2])*rms1;
  const double fy1 = -(lam[1]*dy[1] + lam[2]*dy[2])*rms1;
  const double fz1 = -(lam[1]*dz[1] + lam[2]*dz[2])*rms1;

  const double fx2 =  (lam[1]*dx[1] - lam[3]*dx[3])*rms2;
  const double fy2 =  (lam[1]*dy[1] - lam[3]*dy[3])*rms2;
  const double fz2 =  (lam[1]*dz[1] - lam[3]*dz[3])*rms2;

  const double fx3 =  (lam[2]*dx[2] + lam[3]*dx[3])*rms3;
  const double fy3 =  (lam[2]*dy[2] + lam[3]*dy[3])*rms3;
  const double fz3 =  (lam[2]*dz[2] + lam[3]*dz[3])*rms3;

  pos1.x += fx1;
  pos1.y += fy1;
  pos1.z += fz1;
  vel1.x += fx1/dt;
  vel1.y += fy1/dt;
  vel1.z += fz1/dt;

  pos2.x += fx2;
  pos2.y += fy2;
  pos2.z += fz2;
  vel2.x += fx2/dt;
  vel2.y += fy2/dt;
  vel2.z += fz2/dt;

  pos3.x += fx3;
  pos3.y += fy3;
  pos3.z += fz3;
  vel3.x += fx3/dt;
  vel3.y += fy3/dt;
  vel3.z += fz3/dt;

  //=======================================================================
} // end routine 
//=======================================================================



//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

void rattle_33(
    const Vector& pos1,     // 1st atom's position
    const Vector& pos2,     // 2nd atom's position
    const Vector& pos3,     // 3rd atom's position
    Vector& vel1,     // 1st atom's velocity
    Vector& vel2,     // 2nd atom's velocity
    Vector& vel3,     // 3rd atom's velocity
    double mass1,           // 1st atom's mass
    double mass2,           // 2nd atom's mass
    double mass3            // 3rd atom's mass
    )

  //==========================================================================
{//Begin Routine
  //=======================================================================

  double dx[4],  dy[4],  dz[4];
  double dvx[4], dvy[4], dvz[4];
  double rmm[4][4];
  double amat[4][4];
  double ainv[4][4];
  double avec[4];
  double lam[4];
  double det,rdet;

  //=======================================================================
  // set masses, positions and velocities

  const double rms1 = 1.0/mass1;
  const double rms2 = 1.0/mass2;
  const double rms3 = 1.0/mass3;

  rmm[1][1] = -(rms1+rms2); 
  rmm[1][2] = -rms1; 
  rmm[1][3] =  rms2;
  rmm[2][1] = -rms1; 
  rmm[2][2] = -(rms1+rms3); 
  rmm[2][3] = -rms3;
  rmm[3][1] =  rms2; 
  rmm[3][2] = -rms3; 
  rmm[3][3] = -(rms2+rms3);

  dx[1]  = pos1.x-pos2.x;
  dy[1]  = pos1.y-pos2.y;
  dz[1]  = pos1.z-pos2.z;

  dx[2]  = pos1.x-pos3.x;
  dy[2]  = pos1.y-pos3.y;
  dz[2]  = pos1.z-pos3.z;

  dx[3]  = pos2.x-pos3.x;
  dy[3]  = pos2.y-pos3.y;
  dz[3]  = pos2.z-pos3.z;

  dvx[1] = vel1.x-vel2.x;
  dvy[1] = vel1.y-vel2.y;
  dvz[1] = vel1.z-vel2.z;

  dvx[2] = vel1.x-vel3.x;
  dvy[2] = vel1.y-vel3.y;
  dvz[2] = vel1.z-vel3.z;

  dvx[3] = vel2.x-vel3.x;
  dvy[3] = vel2.y-vel3.y;
  dvz[3] = vel2.z-vel3.z;

  //=======================================================================
  // Solve for lambda

  //-----------------------------------------------------------------
  // set up matrix and rhs  A*lam = vec

  amat[1][1] =-rmm[1][1]*(dx[1]*dx[1] + dy[1]*dy[1] + dz[1]*dz[1]);
  amat[1][2] =-rmm[1][2]*(dx[1]*dx[2] + dy[1]*dy[2] + dz[1]*dz[2]);
  amat[1][3] =-rmm[1][3]*(dx[1]*dx[3] + dy[1]*dy[3] + dz[1]*dz[3]);

  amat[2][1] =-rmm[2][1]*(dx[2]*dx[1] + dy[2]*dy[1] + dz[2]*dz[1]);
  amat[2][2] =-rmm[2][2]*(dx[2]*dx[2] + dy[2]*dy[2] + dz[2]*dz[2]);
  amat[2][3] =-rmm[2][3]*(dx[2]*dx[3] + dy[2]*dy[3] + dz[2]*dz[3]);

  amat[3][1] =-rmm[3][1]*(dx[3]*dx[1] + dy[3]*dy[1] + dz[3]*dz[1]);
  amat[3][2] =-rmm[3][2]*(dx[3]*dx[2] + dy[3]*dy[2] + dz[3]*dz[2]);
  amat[3][3] =-rmm[3][3]*(dx[3]*dx[3] + dy[3]*dy[3] + dz[3]*dz[3]);

  avec[1] = dvx[1]*dx[1] + dvy[1]*dy[1] + dvz[1]*dz[1];
  avec[2] = dvx[2]*dx[2] + dvy[2]*dy[2] + dvz[2]*dz[2];
  avec[3] = dvx[3]*dx[3] + dvy[3]*dy[3] + dvz[3]*dz[3];

  //-----------------------------------------------------------------
  // inverse matrix

  det = (amat[1][1] * (amat[2][2] * amat[3][3] - amat[3][2] * amat[2][3]) + 
      amat[2][1] * (amat[3][2] * amat[1][3] - amat[1][2] * amat[3][3]) + 
      amat[3][1] * (amat[1][2] * amat[2][3] - amat[2][2] * amat[1][3]));
  rdet = 1.0/det;

  ainv[1][1] = (amat[2][2] * amat[3][3] - amat[3][2] * amat[2][3])*rdet;
  ainv[2][2] = (amat[1][1] * amat[3][3] - amat[3][1] * amat[1][3])*rdet;
  ainv[3][3] = (amat[1][1] * amat[2][2] - amat[2][1] * amat[1][2])*rdet;
  ainv[2][1] = (amat[3][1] * amat[2][3] - amat[2][1] * amat[3][3])*rdet;
  ainv[1][2] = (amat[1][3] * amat[3][2] - amat[1][2] * amat[3][3])*rdet;
  ainv[3][1] = (amat[2][1] * amat[3][2] - amat[3][1] * amat[2][2])*rdet;
  ainv[1][3] = (amat[1][2] * amat[2][3] - amat[1][3] * amat[2][2])*rdet;
  ainv[3][2] = (amat[3][1] * amat[1][2] - amat[3][2] * amat[1][1])*rdet;
  ainv[2][3] = (amat[1][3] * amat[2][1] - amat[2][3] * amat[1][1])*rdet;

  //-----------------------------------------------------------------
  // solve for multipliers

  lam[1] = ainv[1][1]*avec[1] + ainv[1][2]*avec[2] + ainv[1][3]*avec[3];
  lam[2] = ainv[2][1]*avec[1] + ainv[2][2]*avec[2] + ainv[2][3]*avec[3];
  lam[3] = ainv[3][1]*avec[1] + ainv[3][2]*avec[2] + ainv[3][3]*avec[3];

  //=======================================================================
  // update velocities 

  vel1.x -= (lam[1]*dx[1] + lam[2]*dx[2])*rms1;
  vel1.y -= (lam[1]*dy[1] + lam[2]*dy[2])*rms1;
  vel1.z -= (lam[1]*dz[1] + lam[2]*dz[2])*rms1;

  vel2.x += (lam[1]*dx[1] - lam[3]*dx[3])*rms2;
  vel2.y += (lam[1]*dy[1] - lam[3]*dy[3])*rms2;
  vel2.z += (lam[1]*dz[1] - lam[3]*dz[3])*rms2;

  vel3.x += (lam[2]*dx[2] + lam[3]*dx[3])*rms3;
  vel3.y += (lam[2]*dy[2] + lam[3]*dy[3])*rms3;
  vel3.z += (lam[2]*dz[2] + lam[3]*dz[3])*rms3;

  //=======================================================================
} // end routine 
//=======================================================================
