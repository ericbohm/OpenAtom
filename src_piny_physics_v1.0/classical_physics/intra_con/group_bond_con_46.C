#include "group_bond_con_46.h"

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//  Constraints between 4-atoms not involved in any other constraint : tet
//     Example : X-CH3   : CH CH CH HH HH HH constrained, X-C not
//==========================================================================



//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

void shake_46(
  Vector& pos1,     // 1st atom's position
  Vector& pos2,     // 2nd atom's position
  Vector& pos3,     // 3rd atom's position
  Vector& pos4,     // 4th atom's position
  Vector& vel1,     // 1st atom's velocity
  Vector& vel2,     // 2nd atom's velocity
  Vector& vel3,     // 3rd atom's velocity
  Vector& vel4,     // 4th atom's velocity
  const Vector& pos1_old, // 1st atom's old position
  const Vector& pos2_old, // 2nd atom's old position
  const Vector& pos3_old, // 3rd atom's old position
  const Vector& pos4_old, // 4th atom's old position
  double mass1,           // 1st atom's mass
  double mass2,           // 2nd atom's mass
  double mass3,           // 3rd atom's mass
  double mass4,           // 4th atom's mass
  double eq1,             // bond length of 1-2 bond
  double eq2,             // bond length of 1-3 bond
  double eq3,             // bond length of 2-3 bond
  double eq4,             // bond length of 2-3 bond
  double eq5,             // bond length of 2-3 bond
  double eq6,             // bond length of 2-3 bond
  double dt,              // time step
  double tol,             // tolerence
  int max_iter        // maximum iterations
 )

//==========================================================================
   {//Begin Routine
//=======================================================================
//         Local Variable declarations                                   

   int na=6;
   int job = 1;
   int ipvt[7],info;

   double dxt[7],dyt[7],dzt[7];
   double dx[7],dy[7],dz[7];
   double dxn[7], dyn[7],dzn[7];
   double dij[7];
   double lam[7], lam_old[7];
   double avec[7];
   double rmm[7][7];
   double amat[50];

//======================================================================= 
// I) Collect masses

   dij[1] = eq1;
   dij[2] = eq2;
   dij[3] = eq3;
   dij[4] = eq4;
   dij[5] = eq5;
   dij[6] = eq6;

   const double rms1 = 1.0/mass1;
   const double rms2 = 1.0/mass2;
   const double rms3 = 1.0/mass3;
   const double rms4 = 1.0/mass4;

   rmm[1][1] = -(rms1+rms2); 
   rmm[1][2] = rmm[1][3] = -rms1;
   rmm[1][4] = rmm[1][5] =  rms2; 
   rmm[1][6] = 0.0;

   rmm[2][1] = -rms1; 
   rmm[2][2] = -(rms1+rms3); 
   rmm[2][3] = -rms1;
   rmm[2][4] = -rms3; 
   rmm[2][5] = 0.0;   
   rmm[2][6] = rms3;

   rmm[3][1] = rmm[3][2] = -rms1; 
   rmm[3][3] = -(rms1+rms4);
   rmm[3][4] = 0.0; 
   rmm[3][5] = rmm[3][6] = -rms4;

   rmm[4][1] =  rms2; 
   rmm[4][2] = -rms3; 
   rmm[4][3] =  0.0;
   rmm[4][4] = -(rms2+rms3); 
   rmm[4][5] = -rms2; 
   rmm[4][6] =  rms3;

   rmm[5][1] =  rms2; 
   rmm[5][2] =  0.0; 
   rmm[5][3] = -rms4;
   rmm[5][4] = -rms2; 
   rmm[5][5] = -(rms2+rms4); 
   rmm[5][6] = -rms4;

   rmm[6][1] =  0.0; 
   rmm[6][2] =  rms3; 
   rmm[6][3] = -rms4;
   rmm[6][4] =  rms3; 
   rmm[6][5] = -rms4; 
   rmm[6][6] = -(rms3+rms4);

//======================================================================= 
// Compute difference vectors 

   dxt[1] = pos1.x-pos2.x; 
   dxt[2] = pos1.x-pos3.x; 
   dxt[3] = pos1.x-pos4.x;
   dxt[4] = pos2.x-pos3.x;
   dxt[5] = pos2.x-pos4.x;
   dxt[6] = pos3.x-pos4.x;

   dyt[1] = pos1.y-pos2.y; 
   dyt[2] = pos1.y-pos3.y; 
   dyt[3] = pos1.y-pos4.y;
   dyt[4] = pos2.y-pos3.y;
   dyt[5] = pos2.y-pos4.y;
   dyt[6] = pos3.y-pos4.y;
  
   dzt[1] = pos1.z-pos2.z; 
   dzt[2] = pos1.z-pos3.z; 
   dzt[3] = pos1.z-pos4.z;
   dzt[4] = pos2.z-pos3.z;
   dzt[5] = pos2.z-pos4.z;
   dzt[6] = pos3.z-pos4.z;

   dx[1] = pos1_old.x-pos2_old.x;
   dx[2] = pos1_old.x-pos3_old.x;
   dx[3] = pos1_old.x-pos4_old.x;
   dx[4] = pos2_old.x-pos3_old.x;
   dx[5] = pos2_old.x-pos4_old.x;
   dx[6] = pos3_old.x-pos4_old.x;

   dy[1] = pos1_old.y-pos2_old.y;
   dy[2] = pos1_old.y-pos3_old.y;
   dy[3] = pos1_old.y-pos4_old.y;
   dy[4] = pos2_old.y-pos3_old.y;
   dy[5] = pos2_old.y-pos4_old.y;
   dy[6] = pos3_old.y-pos4_old.y;

   dz[1] = pos1_old.z-pos2_old.z;
   dz[2] = pos1_old.z-pos3_old.z;
   dz[3] = pos1_old.z-pos4_old.z;
   dz[4] = pos2_old.z-pos3_old.z;
   dz[5] = pos2_old.z-pos4_old.z;
   dz[6] = pos3_old.z-pos4_old.z;

//======================================================================= 
// compute vector of rhs

   for(int i=1; i <= 6; i++){
     avec[i] = dij[i]*dij[i] - (dxt[i]*dxt[i] + dyt[i]*dyt[i] + dzt[i]*dzt[i]);
   }//endfor

//======================================================================= 
// Solve for initial lambda  : A lam = vec

   int iii = 0;
   for(int i=1; i <= 6; i++){
     for(int j=1; j <= 6; j++){
       iii++;
       amat[iii] = 2.0*rmm[i][j]*(dxt[i]*dx[j] + dyt[i]*dy[j] + dzt[i]*dz[j]);
     }//endfor
   }//endfor

   for(int i=1; i <= 6; i++){
     lam[i] = dij[i]*dij[i] - (dxt[i]*dxt[i] + dyt[i]*dyt[i] + dzt[i]*dzt[i]);
   }

// Use generic routines here!! Linpack is now stnd

#ifdef IBM_ESSL
   dgef(&(amat[1]),&na,&na,&(ipvt[1]));
#else
   DGEFA(&(amat[1]),&na,&na,&(ipvt[1]),&info);
#endif
#ifdef IBM_ESSL
   dges(&(amat[1]),&na,&na,&(ipvt[1]),&(lam[1]),&job);
#else
   DGESL(&(amat[1]),&na,&na,&(ipvt[1]),&(lam[1]),&job);
#endif

//=======================================================================        
// Iterative loop to convergence 

   double dmax = 1.0+tol;
   int iter = 0;
   do {
     ++iter;
     if(iter > max_iter) {
       CkPrintf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
       CkPrintf("Shake 4x6 not converged after %d iterations.\n",max_iter);
       CkPrintf("The present tolerance is %g \n",dmax);
       CkPrintf("The desired tolerance is %g \n",tol);
       CkPrintf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
       break;
     }//endif

    //---------------------------------------------------------------
    // Set up guess of difference vectors 

     dxn[1] = 2.0*dxt[1]; dyn[1] = 2.0*dyt[1]; dzn[1] = 2.0*dzt[1];
     dxn[2] = 2.0*dxt[2]; dyn[2] = 2.0*dyt[2]; dzn[2] = 2.0*dzt[2];
     dxn[3] = 2.0*dxt[3]; dyn[3] = 2.0*dyt[3]; dzn[3] = 2.0*dzt[3];
     dxn[4] = 2.0*dxt[4]; dyn[4] = 2.0*dyt[4]; dzn[4] = 2.0*dzt[4];
     dxn[5] = 2.0*dxt[5]; dyn[5] = 2.0*dyt[5]; dzn[5] = 2.0*dzt[5];
     dxn[6] = 2.0*dxt[6]; dyn[6] = 2.0*dyt[6]; dzn[6] = 2.0*dzt[6];
     for(int i=1; i <= 6; i++) {
       for(int j=1; j <= 6; j++) {
         dxn[i] += rmm[i][j]*lam[j]*dx[j];
         dyn[i] += rmm[i][j]*lam[j]*dy[j];
         dzn[i] += rmm[i][j]*lam[j]*dz[j];
       }//endfor
     }//endfor

    //---------------------------------------------------------------
    // Construct A-matrix 

     iii = 0;
     for(int i=1; i <= 6; i++) {
       for(int j=1; j <= 6; j++) {
         iii++;
         amat[iii] = rmm[i][j]*(dxn[i]*dx[j]+dyn[i]*dy[j]+dzn[i]*dz[j]);
       }//endfor
     }//endfor
 
     for(int i=1; i <= 6; i++) { 
       lam_old[i] = lam[i];
       lam[i]     = avec[i];
     }// endfor

#ifdef IBM_ESSL
     dgef(&(amat[1]),&na,&na,&(ipvt[1]));
#else
     DGEFA(&(amat[1]),&na,&na,&(ipvt[1]),&info);
#endif

#ifdef IBM_ESSL
     dges(&(amat[1]),&na,&na,&(ipvt[1]),&(lam[1]),&job);
#else
     DGESL(&(amat[1]),&na,&na,&(ipvt[1]),&(lam[1]),&job);
#endif

     dmax = 0.0;
     for(int i=1; i <= 6; i++) { 
       const double delta = (lam[i]-lam_old[i])*(lam[i]-lam_old[i]);
       dmax = MAX(dmax,delta);
     }// endfor
     dmax = sqrt(dmax);

   } while(dmax > tol);

//=======================================================================        
// Position and velocity update 

  double fx,fy,fz;

  fx = -( lam[1]*dx[1] + lam[2]*dx[2] + lam[3]*dx[3])*rms1;
  fy = -( lam[1]*dy[1] + lam[2]*dy[2] + lam[3]*dy[3])*rms1;
  fz = -( lam[1]*dz[1] + lam[2]*dz[2] + lam[3]*dz[3])*rms1;
  pos1.x += fx;  vel1.x += fx/dt;
  pos1.y += fy;  vel1.y += fy/dt;
  pos1.z += fz;  vel1.z += fz/dt;

  fx = -(-lam[1]*dx[1] + lam[4]*dx[4] + lam[5]*dx[5])*rms2;
  fy = -(-lam[1]*dy[1] + lam[4]*dy[4] + lam[5]*dy[5])*rms2;
  fz = -(-lam[1]*dz[1] + lam[4]*dz[4] + lam[5]*dz[5])*rms2;
  pos2.x += fx;  vel2.x += fx/dt;
  pos2.y += fy;  vel2.y += fy/dt;
  pos2.z += fz;  vel2.z += fz/dt;

  fx = -(-lam[2]*dx[2] - lam[4]*dx[4] + lam[6]*dx[6])*rms3;
  fy = -(-lam[2]*dy[2] - lam[4]*dy[4] + lam[6]*dy[6])*rms3;
  fz = -(-lam[2]*dz[2] - lam[4]*dz[4] + lam[6]*dz[6])*rms3;
  pos3.x += fx;  vel3.x += fx/dt;
  pos3.y += fy;  vel3.y += fy/dt;
  pos3.z += fz;  vel3.z += fz/dt;

  fx = -(-lam[3]*dx[3] - lam[5]*dx[5] - lam[6]*dx[6])*rms4;
  fy = -(-lam[3]*dy[3] - lam[5]*dy[5] - lam[6]*dy[6])*rms4;
  fz = -(-lam[3]*dz[3] - lam[5]*dz[5] - lam[6]*dz[6])*rms4;
  pos4.x += fx;  vel4.x += fx/dt;
  pos4.y += fy;  vel4.y += fy/dt;
  pos4.z += fz;  vel4.z += fz/dt;

//==========================================================================
  } // end routine 
//==========================================================================



//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

void rattle_46(
  const Vector& pos1,     // 1st atom's position
  const Vector& pos2,     // 2nd atom's position
  const Vector& pos3,     // 3rd atom's position
  const Vector& pos4,     // 4th atom's position
  Vector& vel1,     // 1st atom's velocity
  Vector& vel2,     // 2nd atom's velocity
  Vector& vel3,     // 3rd atom's velocity
  Vector& vel4,     // 4th atom's velocity
  double mass1,           // 1st atom's mass
  double mass2,           // 2nd atom's mass
  double mass3,           // 3rd atom's mass
  double mass4            // 4th atom's mass
 )

//==========================================================================
//        Begin Routine                                                     
   {//Begin Routine
//=======================================================================
//         Local Variable declarations                                   

   int na=6;
   int job = 1;
   int ipvt[7],info;

   double dx[7],dy[7],dz[7];
   double dvx[7],dvy[7],dvz[7];
   double lam[7];
   double rmm[7][7];
   double amat[50];
  
//=======================================================================
// inverse mass matrix

  const double rms1= 1.0/mass1;
  const double rms2= 1.0/mass2;
  const double rms3= 1.0/mass3;
  const double rms4= 1.0/mass4;

  rmm[1][1] = -(rms1+rms2); rmm[1][2] = rmm[1][3] = -rms1;
  rmm[1][4] = rmm[1][5] = rms2; rmm[1][6] = 0.0;

  rmm[2][1] = -rms1; rmm[2][2] = -(rms1+rms3); rmm[2][3] = -rms1;
  rmm[2][4] = -rms3; rmm[2][5] = 0.0; rmm[2][6] = rms3;

  rmm[3][1] = rmm[3][2] = -rms1; rmm[3][3] = -(rms1+rms4);
  rmm[3][4] = 0.0; rmm[3][5] = rmm[3][6] = -rms4;

  rmm[4][1] = rms2; rmm[4][2] = -rms3; rmm[4][3] = 0.0;
  rmm[4][4] = -(rms2+rms3); rmm[4][5] = -rms2; rmm[4][6] = rms3;

  rmm[5][1] = rms2; rmm[5][2] = 0.0; rmm[5][3] = -rms4;
  rmm[5][4] = -rms2; rmm[5][5] = -(rms2+rms4); rmm[5][6] = -rms4;

  rmm[6][1] = 0.0; rmm[6][2] = rms3; rmm[6][3] = -rms4;
  rmm[6][4] = rms3; rmm[6][5] = -rms4; rmm[6][6] = -(rms3+rms4);

//=======================================================================
// positions and velocity vectors

  dvx[1] = vel1.x-vel2.x; 
  dvx[2] = vel1.x-vel3.x;
  dvx[3] = vel1.x-vel4.x; 
  dvx[4] = vel2.x-vel3.x;
  dvx[5] = vel2.x-vel4.x; 
  dvx[6] = vel3.x-vel4.x;

  dvy[1] = vel1.y-vel2.y; 
  dvy[2] = vel1.y-vel3.y;
  dvy[3] = vel1.y-vel4.y; 
  dvy[4] = vel2.y-vel3.y;
  dvy[5] = vel2.y-vel4.y; 
  dvy[6] = vel3.y-vel4.y;

  dvz[1] = vel1.z-vel2.z; 
  dvz[2] = vel1.z-vel3.z;
  dvz[3] = vel1.z-vel4.z; 
  dvz[4] = vel2.z-vel3.z;
  dvz[5] = vel2.z-vel4.z; 
  dvz[6] = vel3.z-vel4.z;

  dx[1] = pos1.x-pos2.x; 
  dx[2] = pos1.x-pos3.x;
  dx[3] = pos1.x-pos4.x; 
  dx[4] = pos2.x-pos3.x;
  dx[5] = pos2.x-pos4.x; 
  dx[6] = pos3.x-pos4.x;

  dy[1] = pos1.y-pos2.y; 
  dy[2] = pos1.y-pos3.y;
  dy[3] = pos1.y-pos4.y; 
  dy[4] = pos2.y-pos3.y;
  dy[5] = pos2.y-pos4.y; 
  dy[6] = pos3.y-pos4.y;

  dz[1] = pos1.z-pos2.z; 
  dz[2] = pos1.z-pos3.z;
  dz[3] = pos1.z-pos4.z; 
  dz[4] = pos2.z-pos3.z;
  dz[5] = pos2.z-pos4.z; 
  dz[6] = pos3.z-pos4.z;

//=======================================================================
// construct multiplers A * lam = vec

  int iii = 0;
  for(int i=1; i <= na; i++){
    for(int j=1; j <= na; j++){
       iii++;
       amat[iii] =-rmm[i][j]*(dx[i]*dx[j] + dy[i]*dy[j] + dz[i]*dz[j]);
     }//endfor
   }//endfor

   lam[1] = dvx[1]*dx[1] + dvy[1]*dy[1] + dvz[1]*dz[1];
   lam[2] = dvx[2]*dx[2] + dvy[2]*dy[2] + dvz[2]*dz[2];
   lam[3] = dvx[3]*dx[3] + dvy[3]*dy[3] + dvz[3]*dz[3];
   lam[4] = dvx[4]*dx[4] + dvy[4]*dy[4] + dvz[4]*dz[4];
   lam[5] = dvx[5]*dx[5] + dvy[5]*dy[5] + dvz[5]*dz[5];
   lam[6] = dvx[6]*dx[6] + dvy[6]*dy[6] + dvz[6]*dz[6];

#ifdef IBM_ESSL
   dgef(&(amat[1]),&na,&na,&(ipvt[1]));
#else
   DGEFA(&(amat[1]),&na,&na,&(ipvt[1]),&info);
#endif

#ifdef IBM_ESSL
   dges(&(amat[1]),&na,&na,&(ipvt[1]),&(lam[1]),&job);
#else
   DGESL(&(amat[1]),&na,&na,&(ipvt[1]),&(lam[1]),&job);
#endif

//==========================================================================
// Velocity update 

  vel1.x -= ( lam[1]*dx[1] + lam[2]*dx[2] + lam[3]*dx[3])*rms1;
  vel1.y -= ( lam[1]*dy[1] + lam[2]*dy[2] + lam[3]*dy[3])*rms1;
  vel1.z -= ( lam[1]*dz[1] + lam[2]*dz[2] + lam[3]*dz[3])*rms1;

  vel2.x -= (-lam[1]*dx[1] + lam[4]*dx[4] + lam[5]*dx[5])*rms2;
  vel2.y -= (-lam[1]*dy[1] + lam[4]*dy[4] + lam[5]*dy[5])*rms2;
  vel2.z -= (-lam[1]*dz[1] + lam[4]*dz[4] + lam[5]*dz[5])*rms2;

  vel3.x -= (-lam[2]*dx[2] - lam[4]*dx[4] + lam[6]*dx[6])*rms3;
  vel3.y -= (-lam[2]*dy[2] - lam[4]*dy[4] + lam[6]*dy[6])*rms3;
  vel3.z -= (-lam[2]*dz[2] - lam[4]*dz[4] + lam[6]*dz[6])*rms3;

  vel4.x -= (-lam[3]*dx[3] - lam[5]*dx[5] - lam[6]*dx[6])*rms4;
  vel4.y -= (-lam[3]*dy[3] - lam[5]*dy[5] - lam[6]*dy[6])*rms4;
  vel4.z -= (-lam[3]*dz[3] - lam[5]*dz[5] - lam[6]*dz[6])*rms4;

//=======================================================================
   } // end routine 
//=======================================================================
