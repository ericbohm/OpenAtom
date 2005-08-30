//==========================================================================
//==========================================================================

// This routine computes the energy and forces from                          
// the intramolecular onefour potential.                                    

//==========================================================================
//==========================================================================

#include "onfo.h"

//==========================================================================
//==========================================================================

#define DEBUGF(x) // CmiPrintf x

void onfo (const Vector& pos1,       // 1st atom position
           const Vector& pos2,       // 2nd atom position
           Vector&       force1,     // 1st atom force
           Vector&       force2,     // 2nd atom force
           double        q1,         // charge 1st atom
           double        q2,         // charge 2nd atom
           double*       pvten,
           double*       pvten_tot,
           int           ityp,
           StepOutput*   out)
//==========================================================================
 {//begin routine
//==========================================================================
//               Local variable declarations

  double pe         = 0;
  double vonfo_vdw  = 0; // energy return values 
  double vonfo_coul = 0; // energy return values   
  
  // read-only class contains onfo constants
  const MDONFO & onfo    = readonly_mdintra.mdonfo;
  const GENCELL & cell   = readonly_general_data.gencell;
  int iget_pv_real_inter = readonly_mdinter.mdenergy_ctrl.iget_pv_real_inter;
  
  double r122i,r12i;                 // (r-r_0)^2 and (r-r_0)  
  double vonfo,vonfoe;               // Onfo potential         
  double dvonfo,xxx;                 // Derivative of onfo pot 
  double pre;                        // Force prefactor        
  double q12s;
  double wfor;

  double pvten_tmp[10];
  double dx12,dy12,dz12;
  double vpot,vpot_e;
  double fx1,fy1,fz1;

  // Define local pointers             

  double onfo_sc         = onfo.sc[ityp];
  double onfo_s6         = onfo.s6[ityp];
  double onfo_feps       = onfo.feps[ityp];
  
  // local temp variables
  Vector force1_tmp (0,0,0);
  Vector force2_tmp (0,0,0);

//==========================================================================
// I) loop over all the onfos in steps of nlen to save memory               

  wfor   = 1.0; //  Need to assign wght_tra_res 

  int i;
  for(i=1;i<=9;i++){pvten_tmp[i]=0;}

//--------------------------------------------------------------------------

  dx12 = pos1.x - pos2.x;
  dy12 = pos1.y - pos2.y;
  dz12 = pos1.z - pos2.z;

//--------------------------------------------------------------------------
//  E) Periodic boundary conditions                                         

   r122i = 1.0/(dx12*dx12 + dy12*dy12 + dz12*dz12);
   r12i  = sqrt(r122i);

   DEBUGF(("r12i %.12lg \n",r12i));
//--------------------------------------------------------------------------
//  F) Get the onfoing potential energy                                     

    q12s   = q1*q2*onfo_sc;

    xxx    = onfo_s6*r122i*r122i*r122i;
    DEBUGF(("xxx %.12lg \n",xxx));
    vonfo  = onfo_feps*xxx*(xxx-1.0);
    DEBUGF(("onfo_feps %.12lg \n",onfo_feps));
    vonfoe = q12s*r12i;

    vpot   = vonfo+vonfoe;
    vpot_e = vonfoe;

//--------------------------------------------------------------------------
//  G) Get the force on the atoms using the chain rule                      

    dvonfo = 6.0*onfo_feps*xxx*(2.0*xxx-1.0)*r122i
           + q12s*r12i*r122i;

     pre    = dvonfo;

     fx1 = dx12*pre;
     fy1 = dy12*pre;
     fz1 = dz12*pre;

//--------------------------------------------------------------------------
//  H) Sum the potential energy                                             

     pe          += vpot;
     vonfo_coul  += vpot_e;

//--------------------------------------------------------------------------
//  J) Pressure tensor, etc.                                                

    if(cell.iperd == 2 || cell.iperd == 3) {
       pvten_tmp[1] = dx12*fx1;
       pvten_tmp[5] = dy12*fy1;
       pvten_tmp[9] = dz12*fz1;
       pvten_tmp[2] = dx12*fy1;
       pvten_tmp[4] = dy12*fx1;
    }//endif
    if((cell.iperd) == 3) {
         pvten_tmp[3] = dx12*fz1;
         pvten_tmp[6] = dy12*fz1;
    }//endif

//--------------------------------------------------------------------------

    force1_tmp.x += wfor*fx1;
    force1_tmp.y += wfor*fy1;
    force1_tmp.z += wfor*fz1;

    force2_tmp.x -= force1_tmp.x;
    force2_tmp.y -= force1_tmp.y;
    force2_tmp.z -= force1_tmp.z;

    force1 += force1_tmp;
    force2 += force2_tmp;
    
//==========================================================================
// Increment the Pressure tensor 

   pvten_tmp[4] = pvten_tmp[2];
   pvten_tmp[7] = pvten_tmp[3];
   pvten_tmp[8] = pvten_tmp[6];

   for(i=1;i<=9;i++){
    pvten[i]  +=  pvten_tmp[i]*wfor;
   }//endfor

   if(iget_pv_real_inter==1){    
    for(i=1;i<=9;i++){
      pvten_tot[i] += pvten_tmp[i];
    }//endfor
   }//endif

   vonfo_vdw = pe - vonfo_coul;
   
   out->pe         += pe;
   out->vonfo      += vonfo;
   out->vonfo_vdw  += vonfo_vdw;
   out->vonfo_coul += vonfo_coul;

//--------------------------------------------------------------------------
}// end routine
//==========================================================================
