//=============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=============================================================================
// INLINING BELOW BECAUSE IT IS AT THE BASE OF ALL NESTED LOOPS -JOSHUA
// calculate the force that atom1 exerts on atom2 (and vice versa),
// return the potential energy
//=============================================================================
  inline void InterRoutines::calc_force_explicit
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
 )
  //=============================================================================
{  // begin routine calcForce
  //=============================================================================
  double* hmat = readonly_general_data.gencell.hmat;

  double pot_coul = 0;  // Coulomb energy
  double pot_vdw  = 0;   // vdw energy

  // approximation to error function accurate to 10^{-13} at all r
  const double p  = 0.3614;
  const double e1 = 0.20414220964220030; 
  const double e2 = 0.19975359569614810;
  const double e3 = 0.22131765964055760;
  const double e4 = 0.03360430734640255;
  const double e5 = 0.47325925787217550;
  const double e6 =-0.50907852006973500;
  const double e7 = 0.67726314919476460;
  const double e8 =-0.36991297909221700;
  const double e9 = 0.06965131976970335;
  const double de1 = 1.0*e1; 
  const double de2 = 2.0*e2; 
  const double de3 = 3.0*e3;
  const double de4 = 4.0*e4; 
  const double de5 = 5.0*e5; 
  const double de6 = 6.0*e6;
  const double de7 = 7.0*e7;
  const double de8 = 8.0*e8;
  const double de9 = 9.0*e9;

  //=============================================================================
  // 

  const double x1 = position1.x;
  const double y1 = position1.y;
  const double z1 = position1.z;
  const double x2 = position2.x;
  const double y2 = position2.y;
  const double z2 = position2.z;

  double dx = x2-x1;
  double dy = y2-y1;
  double dz = z2-z1;

  const double cutoff  = cut1+cut2;
  const double cutoff2 = cutoff * cutoff;

  //=============================================================================
  //CmiPrintf ("hmat1=%.6f, hmat5=%.6f, hmat9=%.6f\n", hmat[1], hmat[5], hmat[9]);

  if(iperd==3){
    dx -= NINT(dx/hmat[1])*hmat[1];
    dy -= NINT(dy/hmat[5])*hmat[5];
    dz -= NINT(dz/hmat[9])*hmat[9];
  }// endif
  const double r2  = dx*dx + dy*dy + dz*dz;

  int exclusion=0;

  if (SimParam.DoConstraints || SimParam.DoIntra) {
    int* j;       // exclusion list
    int* j_off;   // exclusion list
    int* nexcl;   // # of exclusions of each atom

    nexcl = readonly_mdintra.mdexcl.num;
    j     = readonly_mdintra.mdexcl.j;
    j_off = readonly_mdintra.mdexcl.j_off;

    if(r2 <= r2max_excl){
      if(ind1_excl > ind2_excl){
        int offset = j_off [ind1_excl];
        for(int i=1;i<=nexcl[ind1_excl];i++){
          if(j[offset+i]==ind2_excl){exclusion=1;}
          if(j[offset+i]>ind2_excl){break;}
        }//endfor
      }else{
        int offset = j_off [ind2_excl];
        for(int i=1;i<=nexcl[ind2_excl];i++){
          if(j[offset+i]==ind1_excl){exclusion=1;}
          if(j[offset+i]>ind1_excl){break;}
        }//endfor
      }//endif
    }//endif
  }//endif

  processed = 0;
#ifdef _PHYSICS_DEBUG_  
  out->numPairsWExcl += 1;
#endif  
  if (r2 > cutoff2 || exclusion == 1){  // Outside cutoff or excluded
    return;
  }/*endif*/
#ifdef _PHYSICS_DEBUG_  
  out->numPairsCalc += 1;
#endif  
  processed = 1;

  //=============================================================================
  // Coulomb energy

  const double q1q2 = charge1*charge2;
  double f          = 0.0;

  if(q1q2*q1q2>0.0){
    const double r = sqrt(r2);
    if(iperd==0){
      pot_coul  += q1q2 / r;
      f         = -q1q2 / (r2*r);
    }else{
      const double palp   = p*alp_ewd;
      const double talp2  = 2.0*alp_ewd*alp_ewd;
      const double ralp   = r * alp_ewd;
      const double eee    = exp(-ralp*ralp);
      const double tt     = 1.0/(1.0+p*ralp);
      const double gerfc  = ((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt
                +e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee;
      const double dgerfc = ((((((((de9*tt+de8)*tt+de7)*tt+de6)*tt+de5)*tt
                +de4)*tt+de3)*tt+de2)*tt+de1)*tt*tt*eee*palp
        +talp2*gerfc*r;
      pot_coul += q1q2 * gerfc/r;
      f        = +(gerfc/r2 + dgerfc/r)*q1q2/r;
    }// endif
  }//endif

  //=============================================================================
  // Vanderwaals forces

  if (SimParam.DoVDW)  {
    const double feps     = 4.0*sqrt(epsilon1*epsilon2);
    const double sig      = 0.5*(sigma1 + sigma2);
    const double sig2     = sig*sig;
    const double sig2r2   = sig2/r2;
    const double sig6r6   = sig2r2*sig2r2*sig2r2;
    const double pre      = feps*sig6r6;
    const double vdwForce = 6.0*pre*(2.0*sig6r6-1.0)/r2;
    pot_vdw  += pre*(sig6r6-1.0);
    f += vdwForce;

#ifdef _NUMERICAL_CHECKER_        
    CmiPrintf("Numerical Force calculation \n");

    double delta = 1.0e-6;
    double dxn = ((x2+delta) - x1);
    double dyn = (y2-y1);
    double dzn = (z2-z1);
    double r2n = (dxn*dxn + dyn*dyn + dzn*dzn);
    double rn  = sqrt(r2n);
    double E1  = feps*( pow(sig/rn,12.) -  pow(sig/rn,6.));

    dxn = (x2 - delta) - x1;
    r2n = (dxn*dxn + dyn*dyn + dzn*dzn);
    rn  = sqrt(r2n);
    double E2  = feps*( pow(sig/rn,12.) -  pow(sig/rn,6.));

    double Fnum_x = (E2 - E1)/(2.0*delta);

    // --------- Numerical Y

    dxn = (x2 - x1);
    dyn = (y2+delta - y1);
    dzn = (z2 - z1);

    r2n = (dxn*dxn + dyn*dyn + dzn*dzn);
    rn  = sqrt(r2n);

    E1 =  feps*( pow(sig/rn,12.) -  pow(sig/rn,6.));

    dyn = (y2-delta - y1);

    r2n = (dxn*dxn + dyn*dyn + dzn*dzn);
    rn  = sqrt(r2n);

    E2 =  feps*( pow(sig/rn,12.) -  pow(sig/rn,6.));
    double Fnum_y = (E2 - E1)/(2.0*delta);

    // --------- Numerical Z

    dxn = (x2 - x1);
    dyn = (y2 - y1);
    dzn = (z2+delta - z1);

    r2n = (dxn*dxn + dyn*dyn + dzn*dzn);
    rn  = sqrt(r2n);

    E1 =  feps*( pow(sig/rn,12.) -  pow(sig/rn,6.));

    dzn = (z2-delta - z1);

    r2n = (dxn*dxn + dyn*dyn + dzn*dzn);
    rn  = sqrt(r2n);

    E2 =  feps*( pow(sig/rn,12.) -  pow(sig/rn,6.));
    double Fnum_z = (E2 - E1)/(2.0*delta);

    CmiPrintf("Numerical Force %e %e %e \n",Fnum_x,Fnum_y,Fnum_z);

#endif //_NUMERICAL_CHECKER_   
  }

  const double fx = f*dx;
  const double fy = f*dy;
  const double fz = f*dz;

  force1.x -= fx;
  force2.x += fx;
  force1.y -= fy;
  force2.y += fy;
  force1.z -= fz;
  force2.z += fz;

  out->pe    += pot_coul+pot_vdw;
  out->vcoul += pot_coul;
  out->vvdw  += pot_vdw;

  //=============================================================================
} // end routine calcForce
//=============================================================================

