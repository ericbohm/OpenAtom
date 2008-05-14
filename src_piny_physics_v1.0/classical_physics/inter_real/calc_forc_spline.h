//=============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=============================================================================
// INLINING BELOW BECAUSE IT IS AT THE BASE OF ALL NESTED LOOPS -JOSHUA
// calculate the force that atom1 exerts on atom2 (and vice versa),
// return the potential energy
//=============================================================================

inline void InterRoutines::calc_force_spline
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
)
//=============================================================================
    {// begin routine calcForce
//=============================================================================
  int**   inter_map_mat; // map to unique interaction index 
  int**   inter_map_mat0;// map to interaction index         
  double* inter_cutoff;  // interaction cutoff <= cut1+cut2   
  int     nsplin;        // number of points in spline lookup 
	
  int     nsplin_m2;     // number of points in spline-2 lookup 
  double* inter_cv0;     // potential spline    
  double* inter_cdv0;    // force spline        
  double* inter_cv0_c;   // ewald couloumb pot spline lookup    
	
  double* inter_cdv0_c;  // ewald couloumb force spline lookup  
  double* inter_rmin_spl;// minimum distance in spline lookup   
  double* inter_dr_spl;  // spacing in spline lookup            
	
  int     ind1_atmtyp;   // 1st atom's type:mdatom_maps->iatm_atm_typ[ind1_excl]
  int     ind2_atmtyp;   // 2nd atom's type:mdatom_maps->iatm_atm_typ[ind2_excl]
  double* hmat;          // cell matrix assumed to be rectangular for now
  double* hmati;         // inverse cell matrix for later

  inter_map_mat = readonly_mdinter.mdinteract.inter_map_mat;
  inter_map_mat0 = readonly_mdinter.mdinteract.inter_map_mat0;
  inter_cutoff = readonly_mdinter.mdinteract.cutoff;
  nsplin = readonly_mdinter.mdinteract.nsplin;
	
  nsplin_m2 = nsplin-2;
  inter_cv0 = readonly_mdinter.mdinteract.cv0;
  inter_cdv0 = readonly_mdinter.mdinteract.cdv0;
  inter_cv0_c = readonly_mdinter.mdinteract.cv0_c;
	
  inter_cdv0_c = readonly_mdinter.mdinteract.cdv0_c;
  inter_rmin_spl = readonly_mdinter.mdinteract.rmin_spl;
  inter_dr_spl = readonly_mdinter.mdinteract.dr_spl;
	
  ind1_atmtyp = readonly_mdatoms.mdatom_maps.iatm_atm_typ[ind1_excl];
  ind2_atmtyp = readonly_mdatoms.mdatom_maps.iatm_atm_typ[ind2_excl];
	
  hmat  = readonly_general_data.gencell.hmat;
  hmati = readonly_general_data.gencell.hmati;

// Get particle displacement

  double dx = position2.x-position1.x;
  double dy = position2.y-position1.y;
  double dz = position2.z-position1.z;

// Fix for general periodicties and general boxes
// Check for the cubic box flag etc.
  if(iperd==3){
    dx -= NINT(dx/hmat[1])*hmat[1];
    dy -= NINT(dy/hmat[5])*hmat[5];
    dz -= NINT(dz/hmat[9])*hmat[9];
  }// endif

  const double r2  = dx*dx + dy*dy + dz*dz;

//=============================================================================
// Check for exclusions 

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

//=============================================================================
// Check to see the interactions is within the cutoff and NOT excluded

  const int inter_index_gen = inter_map_mat0[ind1_atmtyp][ind2_atmtyp];
  const double cutoff_max   = cut1+cut2;
  const double cutoff       = inter_cutoff[inter_index_gen];

#ifdef _CHECK_CUTOFF_
  if(cutoff>cutoff_max){
    CkPrintf("Inconsistent cutoff found for interaction %d %d\n",ind1_atmtyp,ind2_atmtyp);
    CkExit;
  }//endif
#endif
  const double cutoff2     = cutoff * cutoff;

#ifdef _PHYSICS_DEBUG_  
  out->numPairsWExcl += 1;
#endif  

  processed = 0;
  if (r2 > cutoff2 || exclusion == 1){  // Outside cutoff or excluded
    return;
  }//endif
  processed = 1;

#ifdef _PHYSICS_DEBUG_  
  out->numPairsCalc += 1;
#endif  

//=============================================================================
// Prepare the spline lookup, zero potential energy contributions

  const double r        = sqrt(r2);
  const int inter_index = inter_map_mat[ind1_atmtyp][ind2_atmtyp];

  int index_spl;
  double del_r;
  double swit;
  spl_params(nsplin,nsplin_m2,inter_index,r,inter_rmin_spl,inter_dr_spl,
             index_spl,del_r,swit);

//=============================================================================
// Coulomb energy and forces

  double pot_coul   = 0.0;  // Coulomb energy
  double f_coul     = 0.0;
  const double q1q2 = charge1*charge2;

  if(q1q2*q1q2>0.0){

    if(iperd==0){
      pot_coul  =  q1q2 / r;
      f_coul    = -q1q2 / (r2*r);
    }else{
      pot_coul =  q1q2*spl_value(index_spl,del_r,swit,inter_cv0_c);
      f_coul   = -q1q2*spl_value(index_spl,del_r,swit,inter_cdv0_c);;
    }// endif

  }//endif

//=============================================================================
// Vanderwaals energy and forces

  double pot_vdw  = 0.0;   // vdw energy
  double f_vdw    = 0.0;

  if (SimParam.DoVDW)  {
    pot_vdw  = spl_value(index_spl,del_r,swit,inter_cv0);
    f_vdw    = spl_value(index_spl,del_r,swit,inter_cdv0);
  }//endif

//=============================================================================
// Accumulate forces and potentials

  const double f  = f_coul+f_vdw;
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


//==========================================================================
// Spline lookup
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

inline double InterRoutines::spl_value
(
  int index0,          // index into look up table
  double del,          // r-rspl[index]
  double swit,         // switching function
  double *spl_data     // look up table
)

//=======================================================================
    {//Begin Routine
//=======================================================================
//            Local variable declarations                                

   static const double oneth = (1.0/3.0);
   static const double onesi = (1.0/6.0);

//========================================================================
// I) Fit a cubic on the fly and evaluate it 

   const int indexp1 = index0 + 1;
   const int indexm1 = index0 - 1;
   const int indexm2 = index0 - 2;
   const double c0  = spl_data[index0];
   const double fp1 = spl_data[indexp1];
   const double fm1 = spl_data[indexm1];
   const double fm2 = spl_data[indexm2];
   const double c1  = oneth*fp1+0.5*c0-fm1+onesi*fm2;
   const double c2  = 0.5*fp1-c0+0.5*fm1;
   const double c3  = onesi*fp1-0.5*c0+0.5*fm1-onesi*fm2;
   const double spl = swit*(c0 + del*(c1 + del*(c2 + del*c3)));

   return spl;

//--------------------------------------------------------------------------
   }// end routine spl_value
//==========================================================================


//==========================================================================
// Find the spline parameters : del_r and ispl and swit
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

inline void InterRoutines::spl_params
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
)

//=======================================================================
    {//Begin Routine
//=======================================================================


   const double rmin  = rmin_spl[inter];
   const double dr    = dr_spl[inter];
   const double dri   = 1.0/dr;
         int kktemp   = (int) (((r - rmin)*dri)+0.5) + 3;
             kktemp   = MIN(kktemp,nsplin_m2);
             kktemp   = MAX(kktemp,3);
   const double rsp   = ((double)(kktemp-3))*dr+rmin;

   del_r              = (r - rsp)*dri;
   ispl               =  kktemp + (inter-1)*nsplin ;
   swit               = 1.0;

//--------------------------------------------------------------------------
    }// end routine vspl_fetch
//==========================================================================
