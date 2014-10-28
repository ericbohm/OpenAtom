//==========================================================================
//==========================================================================
// This routine computes the energy and forces from                         
// the intramolecular bond potential.                                       
//==========================================================================
//==========================================================================

#include "bond.h"
//==========================================================================
//==========================================================================

#define DEBUGF(x) // CmiPrintf x

void bond(const Vector& pos1,              // 1st atom's position
    const Vector& pos2,              // 2nd atom's position
    Vector&       force1,            // 1st atom's force
    Vector&       force2,            // 2nd atom's force
    double*       pvten,             // Pressure tensor 
    double*       pvten_tot,         // Pressure tensor 
    int           ityp,              // bend typ between atoms   
    StepOutput*   out,
    int           ires_bond)         // = 1
  //=======================================================================
{//begin routine
  //=======================================================================

  //            Local variable declarations                               

  DEBUGF (("x1=%.6f, y1=%.6f, z1=%.6f\n", pos1.x, pos1.y, pos1.z));
  DEBUGF (("x2=%.6f, y2=%.6f, z2=%.6f\n", pos2.x, pos2.y, pos2.z));
  DEBUGF (("Bond type = %d\n", ityp));

  const MDBOND&  bond    = readonly_mdintra.mdbond;
  const GENCELL& gencell = readonly_general_data.gencell;
  int iget_pv_real_inter = readonly_mdinter.mdenergy_ctrl.iget_pv_real_inter;

  int i;
  double r122,r12;                         // (r-r_0)^2 and (r-r_0)  
  double vbond;                            // Bond potential         
  double dvbond;                           // Derivative of bond pot 
  double pre;                              // Force prefactor        
  double rr0;                              // (r-r_0)                
  double wfor;
  double dx12,dy12,dz12;
  double vpot;
  double fx1,fy1,fz1;

  // Define local pointers                                                 

  double bond_eq_pow;

  double bond_c_0        = bond.c_0[ityp];
  double bond_c_1        = bond.c_1[ityp];
  double bond_c_2        = bond.c_2[ityp];
  double bond_c_3        = bond.c_3[ityp];
  double bond_c_4        = bond.c_4[ityp];
  double bond_c_5        = bond.c_5[ityp];
  double bond_c_6        = bond.c_6[ityp];

  double bond_dc_0       = bond.dc_0[ityp];
  double bond_dc_1       = bond.dc_1[ityp];
  double bond_dc_2       = bond.dc_2[ityp];
  double bond_dc_3       = bond.dc_3[ityp];
  double bond_dc_4       = bond.dc_4[ityp];
  double bond_dc_5       = bond.dc_5[ityp];
  double bond_dc_6       = bond.dc_6[ityp];

  double pvten_tmp[10]; 

  //=======================================================================

  if(ires_bond==1){
    bond_eq_pow  = bond.eq_pow[ityp];  // need to create eq_pow_res
  }else{
    bond_eq_pow  = bond.eq_pow[ityp];
  }//endif

  wfor = 1.0;  // Need to assign for respa

  for(i=1;i<=9;i++){(pvten_tmp)[i]=0.0;}

  //--------------------------------------------------------------------
  //  B) Gather positions                                                 

  dx12 = pos1.x - pos2.x;
  dy12 = pos1.y - pos2.y;
  dz12 = pos1.z - pos2.z;

  DEBUGF (("dx12=%.6f, dy12=%.6f, dz=%.6f\n", dx12, dy12, dz12));

  r122 = dx12*dx12 + dy12*dy12 + dz12*dz12;
  r12  = sqrt(r122);

  rr0 = r12 - bond_eq_pow;

  //--------------------------------------------------------------------
  //  F) Get the bonding potential energy using Horner's method         

  vbond =  ((((( bond_c_6
              *rr0 + bond_c_5)
            *rr0 + bond_c_4)
          *rr0 + bond_c_3)
        *rr0 + bond_c_2)
      *rr0 + bond_c_1)
    *rr0 + bond_c_0;

  DEBUGF (("rr0 %.12lg vbond %.12lg \n",rr0,vbond));

  vpot = vbond;

  //---------------------------------------------------------------------
  //  G) Get the force on the atoms using the chain rule 

  dvbond =   (((( bond_dc_6
            *rr0 + bond_dc_5)
          *rr0 + bond_dc_4)
        *rr0 + bond_dc_3)
      *rr0 + bond_dc_2)
    *rr0 + bond_dc_1;

  pre = -dvbond/r12;

  fx1 = dx12*pre;
  fy1 = dy12*pre;
  fz1 = dz12*pre;

  //----------------------------------------------------------------------
  //  H) Sum the potential energy                                        

  //----------------------------------------------------------------------
  //  J) Pressure tensor, etc.                                            

  if(gencell.iperd == 2 || gencell.iperd == 3) {
    pvten_tmp[1] = dx12*fx1;
    pvten_tmp[5] = dy12*fy1;
    pvten_tmp[9] = dz12*fz1;
    pvten_tmp[2] = dx12*fy1;
  }//endif

  if((gencell.iperd) == 3) {
    pvten_tmp[3] = dx12*fz1;
    pvten_tmp[6] = dy12*fz1;
  }//endif

  //----------------------------------------------------------------------

  force1.x += wfor*fx1;
  force1.y += wfor*fy1;
  force1.z += wfor*fz1;

  force2.x -= wfor*fx1;
  force2.y -= wfor*fy1;
  force2.z -= wfor*fz1; 

  //=======================================================================
  // II) Increment the Pressure tensor 

  pvten_tmp[4] = pvten_tmp[2];
  pvten_tmp[7] = pvten_tmp[3];
  pvten_tmp[8] = pvten_tmp[6];

  DEBUGF (("%d pvten_tmp %.12lg \n",1,pvten_tmp[1]));

  for(i=1;i<=9;i++){
    pvten[i] += pvten_tmp[i]*wfor;
  }//endfor

  if(iget_pv_real_inter==1){
    for(i=1;i<=9;i++){
      pvten_tot[i] += pvten_tmp[i];
    }//endfor
  }//endif

  out->pe    += vpot;
  out->vbond += vpot;

  //--------------------------------------------------------------------------
}// end routine 
//==========================================================================


