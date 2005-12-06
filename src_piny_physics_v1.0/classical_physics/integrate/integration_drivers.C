//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

#include "standard_include.h"
#include "../../../include/debug_flags.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_mdintegrate.h"
#include "../../../include/Atoms.h"
#include "../class_defs/ATOM_OPERATIONS/class_atomintegrate.h"

//============================================================================
//        NVE Integration : 2nd half
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::integrate_nve_2nd_half(int itime,int natm,Atom *atoms,
                                double *eKinetic)
//============================================================================
    {//begin routine 
//============================================================================
// Local variables

   GENERAL_DATA *general_data = GENERAL_DATA::get();
#include "../class_defs/allclass_strip_gen.h"
   double dt  = gentimeinfo->dt;
   double dt2 = dt*0.5;

//============================================================================
// Evolve the system : if its not the first time step

   if(itime>0){
     for(int i=0;i<natm;i++){
       atoms[i].vx += dt2*atoms[i].fx/atoms[i].m;
       atoms[i].vy += dt2*atoms[i].fy/atoms[i].m;
       atoms[i].vz += dt2*atoms[i].fz/atoms[i].m;
     }//endfor  
   }//endif
 
   computeEkin(natm,atoms,eKinetic);

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//        NVE Integration : 2nd half
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::integrate_nve_1st_half(int natm,Atom *atoms)
//============================================================================
    {//begin routine 
//============================================================================
// Local variables

   GENERAL_DATA *general_data = GENERAL_DATA::get();
#include "../class_defs/allclass_strip_gen.h"
   double dt  = gentimeinfo->dt;
   double dt2 = dt*0.5;

//============================================================================
// Evolve the system 

   for(int i=0;i<natm;i++){
     atoms[i].vx += dt2*atoms[i].fx/atoms[i].m;
     atoms[i].vy += dt2*atoms[i].fy/atoms[i].m;
     atoms[i].vz += dt2*atoms[i].fz/atoms[i].m;
   }//endfor  

   for(int i=0;i<natm;i++){
     atoms[i].x += dt*atoms[i].vx;
     atoms[i].y += dt*atoms[i].vy;
     atoms[i].z += dt*atoms[i].vz;
   }//endfor  

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//        NVT Integration : 2nd half
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::integrate_nvt_2nd_half(int itime,int natm,int len_nhc,
                           Atom *atoms,AtomNHC *atomsNHC,double *eKinetic,
                           double *eKineticNhc,double *potNhc)
//============================================================================
    {//begin routine 
//============================================================================
// Local variables

   MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();
   GENERAL_DATA *general_data = GENERAL_DATA::get();
   MDTHERM_INFO *mdtherm_info = &(mdintegrate->mdtherm_info);
#include "../class_defs/allclass_strip_gen.h"

   double dt  = gentimeinfo->dt;
   double dt2 = dt*0.5;
   int nresp  = mdtherm_info->nres_nhc;
   int nyosh  = mdtherm_info->nyosh_nhc;

//============================================================================
// Evolve the system : if its not the first time step

   if(itime>0){
     for(int i =0;i<natm;i++){
       atoms[i].vx += dt2*atoms[i].fx/atoms[i].m;
       atoms[i].vy += dt2*atoms[i].fy/atoms[i].m;
       atoms[i].vz += dt2*atoms[i].fz/atoms[i].m;
     }//endfor  
     applyNHC(natm,len_nhc,atoms,atomsNHC,dt,nyosh,nresp);
   }//endif

   computeEkin(natm,atoms,eKinetic);
   computeENHC(natm,len_nhc,atomsNHC,eKineticNhc,potNhc);

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//        NVT Integration : 1st half
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::integrate_nvt_1st_half(int natm,int len_nhc,
                                           Atom *atoms,AtomNHC *atomsNHC)
//============================================================================
    {//begin routine 
//============================================================================
// Local variables

   MDINTEGRATE  *mdintegrate  = MDINTEGRATE::get();
   GENERAL_DATA *general_data = GENERAL_DATA::get();
   MDTHERM_INFO *mdtherm_info = &(mdintegrate->mdtherm_info);
#include "../class_defs/allclass_strip_gen.h"

   double dt  = gentimeinfo->dt;
   double dt2 = dt*0.5;
   int nresp  = mdtherm_info->nres_nhc;
   int nyosh  = mdtherm_info->nyosh_nhc;

//============================================================================
// Evolve the system : 

   applyNHC(natm,len_nhc,atoms,atomsNHC,dt,nyosh,nresp);
   for(int i =0;i<natm;i++){
     atoms[i].vx += dt2*atoms[i].fx/atoms[i].m;
     atoms[i].vy += dt2*atoms[i].fy/atoms[i].m;
     atoms[i].vz += dt2*atoms[i].fz/atoms[i].m;
   }//endfor  

   for(int i=0;i<natm;i++){
     atoms[i].x += dt*atoms[i].vx;
     atoms[i].y += dt*atoms[i].vy;
     atoms[i].z += dt*atoms[i].vz;
   }//endfor  

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::computeEkin(int natm, Atom *atoms, double *eKinetic_ret){

   double eKinetic = 0.0;
   for(int i=0;i<natm;i++){
      eKinetic += atoms[i].m*(atoms[i].vx*atoms[i].vx+
                              atoms[i].vy*atoms[i].vy+
                              atoms[i].vz*atoms[i].vz);
   }//endfor
   eKinetic *= 0.5;

   (*eKinetic_ret) = eKinetic;

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::computeENHC(int natm, int len_nhc,AtomNHC *atomsNHC, 
                                double *eKinNHC_ret,double *potNHC_ret)
//============================================================================
  {//begin routine
//============================================================================

   double potNHC = 0.0;
   for(int i=0;i<natm;i++){potNHC += atomsNHC[i].posKT;}

   double eKinNHC = 0.0;
   for(int i=0;i<natm;i++){
     for(int j=0;j<len_nhc;j++){
        eKinNHC += atomsNHC[i].m[j]*(atomsNHC[i].vx[j]*atomsNHC[i].vx[j]+
                                     atomsNHC[i].vy[j]*atomsNHC[i].vy[j]+
                                     atomsNHC[i].vz[j]*atomsNHC[i].vz[j]);
     }//endfor
   }//endfor
   eKinNHC *= 0.5;
   
   (*eKinNHC_ret) = eKinNHC;
   (*potNHC_ret)  = potNHC;

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void ATOMINTEGRATE::applyNHC(int natm,int len_nhc,Atom *atoms,AtomNHC *atomsNHC,
                             double dt,int nyosh,int nresp)
//============================================================================
  {//begin routine
//============================================================================
// Get the yoshida stuff out

  double wdti[20],wdti2[20],wdti4[20],wdti8[20];
  double dti = dt/((double)nresp);
  set_yosh(nyosh,dti,wdti,wdti2,wdti4,wdti8);

//============================================================================
// (I)Compute the forces : Precompute  m*v^2 of atoms

   for(int i=0;i<natm;i++){
     atoms[i].mvx2 = atoms[i].m*atoms[i].vx*atoms[i].vx;
     atoms[i].mvy2 = atoms[i].m*atoms[i].vy*atoms[i].vy;
     atoms[i].mvz2 = atoms[i].m*atoms[i].vz*atoms[i].vz;
   }//endfor
   get_forc_NHC0(natm,atoms,atomsNHC);
   for(int ic=1;ic<len_nhc;ic++){get_forc_NHC(natm,ic,atomsNHC);}

//============================================================================
// (II) Yoshida-Suzuki step yourself to heaven

  for(int iresn=1;iresn<=nresp;iresn++){
  for(int iyosh=1;iyosh<=nyosh;iyosh++){
//--------------------------------------------------------------------------
//  1) Evolve the last therm velocity in each chain                         
    evolve_vNHCM(natm,len_nhc,atomsNHC,wdti4[iyosh]);
//--------------------------------------------------------------------------
//  2) Evolve the last-1 to the first therm velocity in each chain        
    for(int ic=len_nhc-2;ic>=0;ic--){
      evolve_vNHC(natm,ic,atomsNHC,wdti4[iyosh],wdti8[iyosh]);
    }//endfor
//--------------------------------------------------------------------------
//  3) Evolve the particle velocities 
    evolve_vAtmNHC(natm,atoms,atomsNHC,wdti2[iyosh]);
//--------------------------------------------------------------------------
//  4) Evolve the therm positions                                           
    evolve_pNHC(natm,len_nhc,atomsNHC,wdti2[iyosh]);
//--------------------------------------------------------------------------
//  5) Evolve the 1 to last-1 therm velocity in each chain : get forces
    get_forc_NHC0(natm,atoms,atomsNHC);
    for(int ic=0,icp=1;ic<(len_nhc-1);ic++,icp++){
      evolve_vNHC(natm,ic,atomsNHC,wdti4[iyosh],wdti8[iyosh]);
      get_forc_NHC(natm,icp,atomsNHC);
    }//endfor
//--------------------------------------------------------------------------
//  6) Evolve the last therm velocotiy in each chain                        
    evolve_vNHCM(natm,len_nhc,atomsNHC,wdti4[iyosh]);
//--------------------------------------------------------------------------
  }}//endfor: iyosh,iresn

//==========================================================================
  }//end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMINTEGRATE::get_forc_NHC0(int natm,Atom *atoms,AtomNHC *atomsNHC){
  for(int i=0;i<natm;i++){
    atomsNHC[i].fx[0] = (atoms[i].mvx2-atomsNHC[i].kT)/atomsNHC[i].m[0];
    atomsNHC[i].fy[0] = (atoms[i].mvy2-atomsNHC[i].kT)/atomsNHC[i].m[0];
    atomsNHC[i].fz[0] = (atoms[i].mvz2-atomsNHC[i].kT)/atomsNHC[i].m[0];
  }//endfor
}//end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMINTEGRATE::get_forc_NHC(int natm,int ic,AtomNHC *atomsNHC){
  int icm = ic-1;
  for(int i=0;i<natm;i++){
    atomsNHC[i].fx[ic]=(atomsNHC[i].m[icm]*atomsNHC[i].vx[icm]*atomsNHC[i].vx[icm]
                       -atomsNHC[i].kT)/atomsNHC[i].m[ic];
    atomsNHC[i].fy[ic]=(atomsNHC[i].m[icm]*atomsNHC[i].vy[icm]*atomsNHC[i].vy[icm]
                       -atomsNHC[i].kT)/atomsNHC[i].m[ic];
    atomsNHC[i].fz[ic]=(atomsNHC[i].m[icm]*atomsNHC[i].vz[icm]*atomsNHC[i].vz[icm]
                       -atomsNHC[i].kT)/atomsNHC[i].m[ic];
  }//endfor
}//end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMINTEGRATE::evolve_vNHCM(int natm,int len_nhc,AtomNHC *atomsNHC,
                                 double wdti4)
//==========================================================================
{//begin routine
  int ic = len_nhc-1;
  for(int i=0;i<natm;i++){
    atomsNHC[i].vx[ic] += atomsNHC[i].fx[ic]*wdti4;
    atomsNHC[i].vy[ic] += atomsNHC[i].fy[ic]*wdti4;
    atomsNHC[i].vz[ic] += atomsNHC[i].fz[ic]*wdti4;
  }//endfor
}//end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMINTEGRATE::evolve_vNHC(int natm,int ic,AtomNHC *atomsNHC,
                                double wdti4,double wdti8)
//==========================================================================
{//begin routine
  double argx,argy,argz;
  double aax,aay,aaz;
  int icp = ic+1;
  for(int i=0;i<natm;i++){
    argx = -wdti8*atomsNHC[i].vx[icp]; aax = exp(argx);  
    argy = -wdti8*atomsNHC[i].vy[icp]; aay = exp(argy);
    argz = -wdti8*atomsNHC[i].vz[icp]; aaz = exp(argz);
    atomsNHC[i].vx[ic] = atomsNHC[i].vx[ic]*aax*aax + wdti4*atomsNHC[i].fx[ic]*aax;
    atomsNHC[i].vy[ic] = atomsNHC[i].vy[ic]*aay*aay + wdti4*atomsNHC[i].fy[ic]*aay;
    atomsNHC[i].vz[ic] = atomsNHC[i].vz[ic]*aaz*aaz + wdti4*atomsNHC[i].fz[ic]*aaz;
  }//endfor
}//end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMINTEGRATE::evolve_vAtmNHC(int natm,Atom *atoms,AtomNHC *atomsNHC,
                                   double wdti2)
//==========================================================================
{//begin routine
  double argx,argy,argz;
  double aax,aay,aaz;
  for(int i=0;i<natm;i++){
    argx = -wdti2*atomsNHC[i].vx[0]; aax = exp(argx); 
    argy = -wdti2*atomsNHC[i].vy[0]; aay = exp(argy); 
    argz = -wdti2*atomsNHC[i].vz[0]; aaz = exp(argz);
    atoms[i].vx *= aax;              atoms[i].mvx2 *= (aax*aax); 
    atoms[i].vy *= aay;              atoms[i].mvy2 *= (aay*aay); 
    atoms[i].vz *= aaz;              atoms[i].mvz2 *= (aaz*aaz); 
  }//endfor
}//end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMINTEGRATE::evolve_pNHC(int natm,int len_nhc, AtomNHC *atomsNHC,
                                double wdti2)
//==========================================================================
{//begin routine
  for(int i=0;i<natm;i++){
  for(int ic=0;ic<len_nhc;ic++){
    double pre = (atomsNHC[i].kT*wdti2);
    atomsNHC[i].posKT += (pre*(atomsNHC[i].vx[ic]+atomsNHC[i].vy[ic]
                              +atomsNHC[i].vz[ic]));
             
  }}//endfor
}//end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void ATOMINTEGRATE::set_yosh(int nyosh,double dt,double *wdt,double *wdt2,
                             double *wdt4,double *wdt8)
//========================================================================
  {//begin routine
//========================================================================
//             Local variable declarations                                

   double temp,p2,onethird;

//========================================================================
// Find the yoshida steps you want

  switch(nyosh){
    case 1: 
          wdt[1] = 1.0;
         break;
    case 3:
          onethird = 1.0/3.0;
          temp    = pow(2.0,onethird);
          wdt[1] =  1.0/(2.0-temp);
          wdt[2] = -temp/(2.0-temp);
          wdt[3] =  1.0/(2.0-temp);
         break;
    case 5:
          onethird = 1.0/3.0;
          temp = pow(4.0,onethird);
          p2    = 1.0/(4.0-temp);
          wdt[1]  = p2;
          wdt[2]  = p2;
          wdt[3]  = 1.0-4.0*p2;
          wdt[4]  = p2;
          wdt[5]  = p2;
         break;
    case 7:
          wdt[1] =  0.784513610477560;
          wdt[2] =  0.235573213359357;
          wdt[3] = -1.17767998417887;
          wdt[4] =  1.0 - 2.0*(wdt[1]+wdt[2]+wdt[3]);
          wdt[5] = -1.17767998417887;
          wdt[6] =  0.235573213359357;
          wdt[7] =  0.784513610477560;
         break;
    case 9:
          wdt[1] =  0.192;
          wdt[2] =  0.554910818409783619692725006662999;
          wdt[3] =  0.124659619941888644216504240951585;
          wdt[4] = -0.843182063596933505315033808282941;
          wdt[5] =  1.0 - 2.0*(wdt[1]+wdt[2]+wdt[3]+wdt[4]);
          wdt[6] = -0.843182063596933505315033808282941;
          wdt[7] =  0.124659619941888644216504240951585;
          wdt[8] =  0.554910818409783619692725006662999;
          wdt[9] =  0.192;
         break;
  }//switch

//========================================================================
// Scale the yosida steps by the time step

  for(int i=1;i<=nyosh;i++){
    wdt[i]  = wdt[i]*dt;
    wdt2[i] = wdt[i]*0.5;
    wdt4[i] = wdt[i]*0.25;
    wdt8[i] = wdt[i]*0.125;
  }//endfor

//-------------------------------------------------------------------------
   }// end routine
//==========================================================================



