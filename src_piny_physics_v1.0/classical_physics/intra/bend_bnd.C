//==========================================================================
//==========================================================================
// This routine computes the energy and forces from                         
// the intramolecular bend_bnd potential.                                   
//==========================================================================
//==========================================================================

#include "standard_include.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_mdintra.h"

#include "../proto_defs/proto_bend_bnd.h"
#include "../proto_defs/proto_vector.h"

//==========================================================================
//==========================================================================

void bend_bnd( const Vector& pos1,  // 1st atom's position
               const Vector& pos2,  // 2nd atom's position
               const Vector& pos3,  // 3rd atom's position
               Vector& force1,      // 1st atom's force
               Vector& force2,      // 2nd atom's force
               Vector& force3,      // 2nd atom's force
               const MDBEND_BND & bend_bnd, // read-only class with bond constants 
               const GENCELL & gencell, 
               double *pvten,          // Pressure tensor 
               double *pvten_tot,      // Pressure tensor 
               double *vbend_bndt,
               double *vbend_bnd_bend, // contains contribution of bend_bnd_bend to bend energy 
               double *vbend_bnd_bond, // contains contribution of bend_bnd_bond to bond energy 
               int ityp,               // bend_bnd typ between atoms   
               int iget_pv_real_inter)
        

//==========================================================================
{//begin routine
//==========================================================================

  double r122,r12;                 // (r-r_0)^2 and (r-r_0)      
  double r322,r32;
  double cost,sint,seps;
  double vbendc,vbends;            // Bend_bnd potential       
  double dvbendc,dvbends;          // Derivative of bend_bnd pot 
  double vbond, dvbond;
  double sisum;                    // sine sum                 
  double pre;                      // Force prefactor          
  double rpmag;                    // 1/(r12*r32)              
  double cos122,cos322;            // cos/r122, cos/r322       
  double rr0;
  double wfor;

  double fxtmp,fytmp,fztmp;
  double fx1,fy1,fz1;
  double fx2,fy2,fz2;
  double fx3,fy3,fz3;
  double dx12,dy12,dz12;
  double dx23,dy23,dz23;
  double dx43,dy43,dz43;
 
  double vpot;

// Define local pointers                                                

  double bend_bnd_eq_bond = bend_bnd.eq_bond[ityp];

  double bend_bnd_cbend_0 = bend_bnd.cbend_0[ityp];
  double bend_bnd_cbend_1 = bend_bnd.cbend_1[ityp];
  double bend_bnd_cbend_2 = bend_bnd.cbend_2[ityp];
  double bend_bnd_cbend_3 = bend_bnd.cbend_3[ityp];
  double bend_bnd_cbend_4 = bend_bnd.cbend_4[ityp];
  double bend_bnd_cbend_5 = bend_bnd.cbend_5[ityp];
  double bend_bnd_cbend_6 = bend_bnd.cbend_6[ityp];

  double bend_bnd_sbend_0 = bend_bnd.sbend_0[ityp];
  double bend_bnd_sbend_1 = bend_bnd.sbend_1[ityp];
  double bend_bnd_sbend_2 = bend_bnd.sbend_2[ityp];
  double bend_bnd_sbend_3 = bend_bnd.sbend_3[ityp];
  double bend_bnd_sbend_4 = bend_bnd.sbend_4[ityp];
  double bend_bnd_sbend_5 = bend_bnd.sbend_5[ityp];
  double bend_bnd_sbend_6 = bend_bnd.sbend_6[ityp];

  double bend_bnd_dcbend_0= bend_bnd.dcbend_0[ityp];
  double bend_bnd_dcbend_1= bend_bnd.dcbend_1[ityp];
  double bend_bnd_dcbend_2= bend_bnd.dcbend_2[ityp];
  double bend_bnd_dcbend_3= bend_bnd.dcbend_3[ityp];
  double bend_bnd_dcbend_4= bend_bnd.dcbend_4[ityp];
  double bend_bnd_dcbend_5= bend_bnd.dcbend_5[ityp];
  double bend_bnd_dcbend_6= bend_bnd.dcbend_6[ityp];

  double bend_bnd_dsbend_0= bend_bnd.dsbend_0[ityp];
  double bend_bnd_dsbend_1= bend_bnd.dsbend_1[ityp];
  double bend_bnd_dsbend_2= bend_bnd.dsbend_2[ityp];
  double bend_bnd_dsbend_3= bend_bnd.dsbend_3[ityp];
  double bend_bnd_dsbend_4= bend_bnd.dsbend_4[ityp];
  double bend_bnd_dsbend_5= bend_bnd.dsbend_5[ityp];
  double bend_bnd_dsbend_6= bend_bnd.dsbend_6[ityp];

  double bend_bnd_cbond_0 = bend_bnd.cbond_0[ityp];
  double bend_bnd_cbond_1 = bend_bnd.cbond_1[ityp];
  double bend_bnd_cbond_2 = bend_bnd.cbond_2[ityp];
  double bend_bnd_cbond_3 = bend_bnd.cbond_3[ityp];
  double bend_bnd_cbond_4 = bend_bnd.cbond_4[ityp];
  double bend_bnd_cbond_5 = bend_bnd.cbond_5[ityp];
  double bend_bnd_cbond_6 = bend_bnd.cbond_6[ityp];

  double bend_bnd_dcbond_0 = bend_bnd.dcbond_0[ityp];
  double bend_bnd_dcbond_1 = bend_bnd.dcbond_1[ityp];
  double bend_bnd_dcbond_2 = bend_bnd.dcbond_2[ityp];
  double bend_bnd_dcbond_3 = bend_bnd.dcbond_3[ityp];
  double bend_bnd_dcbond_4 = bend_bnd.dcbond_4[ityp];
  double bend_bnd_dcbond_5 = bend_bnd.dcbond_5[ityp];
  double bend_bnd_dcbond_6 = bend_bnd.dcbond_6[ityp];

  double pvten_tmp[10];

//=======================================================================
// 0) Initialize                                                         

  seps = 1.0e-8;
  wfor = 1.0; // Need to set wght_tra_res 

  int i;
  for(i=1;i<=9;i++){(pvten_tmp)[i]=0;}

//=======================================================================
//=======================================================================
//  II) Gather positions                                               
    
      dx12 = pos1.x - pos2.x;
      dy12 = pos1.y - pos2.y;
      dz12 = pos1.z - pos2.z;

      dx23 = pos3.x - pos2.x;
      dy23 = pos3.y - pos2.y;
      dz23 = pos3.z - pos2.z;

      dx43 = pos1.x - pos3.x;
      dy43 = pos1.y - pos3.y;
      dz43 = pos1.z - pos3.z;

//=======================================================================
//=======================================================================
// VI) Get the bend energies                                            

      r122  = (dx12*dx12 + dy12*dy12 + dz12*dz12);
      r12   = sqrt(r122);
      
      r322  = (dx23*dx23 + dy23*dy23 + dz23*dz23);
      r32   = sqrt(r322);
      rpmag = 1.0/(r12*r32);
      
//-------------------------------------------------------------------
//  A) Calculate the cosine and sine of the angle between them       
//     (cos(theta_{123}), sine(theta_{123}))                         
      
      cost = (dx12*dx23 + dy12*dy23 + dz12*dz23)/(r12*r32);
      
      cost   = (cost < 1.0 ? cost:1.0);
      cost   = (cost > -1.0 ? cost:-1.0);
      
      sint   = sqrt(1.0 - cost*cost);
      sint   = (sint > seps ? sint:seps);
      
      cos122 = cost/r122;
      cos322  = cost/r322;

//--------------------------------------------------------------------
//  B) Get the bend_bnding potential energy via cosine and sine power 
//    series using Horner's method                                   
      
      vbendc = ((((((bend_bnd_cbend_6 
                     *cost + bend_bnd_cbend_5)
                    *cost + bend_bnd_cbend_4)
                   *cost + bend_bnd_cbend_3)
                  *cost + bend_bnd_cbend_2)
                 *cost + bend_bnd_cbend_1)
                *cost + bend_bnd_cbend_0);
      
      sisum = (((((bend_bnd_sbend_6 
                  *sint + bend_bnd_sbend_5)
                 *sint + bend_bnd_sbend_4)
                *sint + bend_bnd_sbend_3)
               *sint + bend_bnd_sbend_2)
              *sint);
      
      vbends = sisum*cost + bend_bnd_sbend_1 *sint;
      
      vpot = vbendc + vbends;

//-------------------------------------------------------------------
//  C) Get the force on the atoms using the chain rule and Horner's  
      
      dvbendc =  (((((bend_bnd_dcbend_6 
                      *cost + bend_bnd_dcbend_5)
                     *cost + bend_bnd_dcbend_4)
                    *cost + bend_bnd_dcbend_3)
                   *cost + bend_bnd_dcbend_2)
                  *cost + bend_bnd_dcbend_1);
      
      dvbends =  ((((bend_bnd_dsbend_6 
                     *sint + bend_bnd_dsbend_5)
                    *sint + bend_bnd_dsbend_4)
                   *sint + bend_bnd_dsbend_3)
                  *sint + bend_bnd_dsbend_2);
      
      dvbends = dvbends*cost + bend_bnd_dsbend_1;
      
      pre     = -(dvbendc + sisum - dvbends*cost/sint);

      fx1 = (dx23*rpmag - dx12*cos122)*pre;
      fy1 = (dy23*rpmag - dy12*cos122)*pre;
      fz1 = (dz23*rpmag - dz12*cos122)*pre;
      fx3 = (dx12*rpmag - dx23*cos322)*pre;
      fy3 = (dy12*rpmag - dy23*cos322)*pre;
      fz3 = (dz12*rpmag - dz23*cos322)*pre;
      fx2 =-(fx1 + fx3); 
      fy2 =-(fy1 + fy3);
      fz2 =-(fz1 + fz3);

      pvten_tmp[1] = dx12*fx1 + dx23*fx3;
      pvten_tmp[5] = dy12*fy1 + dy23*fy3;
      pvten_tmp[9] = dz12*fz1 + dz23*fz3;

      pvten_tmp[2] = dx12*fy1 + dx23*fy3;
      pvten_tmp[3] = dx12*fz1 + dx23*fz3;
      pvten_tmp[6] = dy12*fz1 + dy23*fz3;

      *vbend_bnd_bend += vpot;

//====================================================================
// VIII) Get the bond energy                                               

      r122 = dx43*dx43 + dy43*dy43 + dz43*dz43;
      r12  = sqrt(r122);
      rr0 = r12 - bend_bnd_eq_bond ;

//--------------------------------------------------------------------
//  A) Get the bend_bnding potential energy using Horner's method     
      
      vbond =  ((((( bend_bnd_cbond_6 
                    *rr0 + bend_bnd_cbond_5)
                   *rr0 + bend_bnd_cbond_4)
                  *rr0 + bend_bnd_cbond_3)
                 *rr0 + bend_bnd_cbond_2)
                *rr0 + bend_bnd_cbond_1)
        *rr0 + bend_bnd_cbond_0 ;
      
      vpot += vbond;
      
//--------------------------------------------------------------------
//  B) Get the force on the atoms using the chain rule 
      
      dvbond =   (((( bend_bnd_dcbond_6 
                     *rr0 + bend_bnd_dcbond_5)
                    *rr0 + bend_bnd_dcbond_4)
                   *rr0 + bend_bnd_dcbond_3)
                  *rr0 + bend_bnd_dcbond_2)
                 *rr0 + bend_bnd_dcbond_1 ;
      
      pre    = -dvbond/r12;

      fxtmp = dx43*pre;
      fytmp = dy43*pre;
      fztmp = dz43*pre;

      fx1 += fxtmp;
      fy1 += fytmp;
      fz1 += fztmp;
      fx3 -= fxtmp;
      fy3 -= fytmp;
      fz3 -= fztmp;

      pvten_tmp[1] += dx43*fxtmp;
      pvten_tmp[5] += dy43*fytmp;
      pvten_tmp[9] += dz43*fztmp;
      pvten_tmp[2] += dx43*fytmp;
      pvten_tmp[3] += dx43*fztmp;
      pvten_tmp[6] += dy43*fztmp;

//==========================================================================
//  IX) Sum the potential energy                                           
    
      *vbend_bndt += vpot;

//==========================================================================

      force1.x  += wfor*fx1;
      force1.y  += wfor*fy1;
      force1.z  += wfor*fz1;

      force2.x  += wfor*fx2;
      force2.y  += wfor*fy2;
      force2.z  += wfor*fz2;

      force3.x  += wfor*fx3;
      force3.y  += wfor*fy3;
      force3.z  += wfor*fz3;

//======================================================================
// XII) Increment the Pressure tensor 

  pvten_tmp[4] = pvten_tmp[2];
  pvten_tmp[7] = pvten_tmp[3];
  pvten_tmp[8] = pvten_tmp[6];

  for(i=1;i<=9;i++){
    pvten[i] += pvten_tmp[i]*wfor;
  }//endfor

  if(iget_pv_real_inter==1){
   for(i=1;i<=9;i++){
     pvten_tot[i] += pvten_tmp[i];
   }//endfor
  }//endif

  *vbend_bnd_bond = (*vbend_bndt)-(*vbend_bnd_bend);

//----------------------------------------------------------------------
}//end routine
//==========================================================================


