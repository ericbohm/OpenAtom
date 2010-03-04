/** \file groups.C
 *           Processor group class Functions : Atoms and parainfo
 */


#include "groups.h"
#include "eesCache.h"
#include "cp_state_ctrl/CP_State_GSpacePlane.h"
#include "fft_slab_ctrl/fftCacheSlab.h"
#include "utility/util.h"
#include "CPcharmParaInfoGrp.h"
#include "load_balance/IntMap.h"
#include "charm++.h"

#include <cmath>

using namespace std;
//#include "CP_State_Plane.h"

//----------------------------------------------------------------------------
#define CHARM_ON
#include "src_piny_physics_v1.0/include/class_defs/piny_constants.h"
#include "src_piny_physics_v1.0/include/class_defs/ATOM_OPERATIONS/class_atomintegrate.h"
#include "src_piny_physics_v1.0/include/class_defs/ATOM_OPERATIONS/class_atomoutput.h"
#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cprspaceion.h"

//----------------------------------------------------------------------------

extern CkVec <CProxy_PIBeadAtoms>       UPIBeadAtomsProxy;
extern CkVec <IntMap2on2> GSImaptable;
extern CkVec <CProxy_EnergyGroup>          UegroupProxy;
extern Config                      config;
extern CkVec <CProxy_CP_State_GSpacePlane> UgSpacePlaneProxy;
extern CkVec <CProxy_GSpaceDriver>         UgSpaceDriverProxy;
extern CkVec <CProxy_AtomsGrp>             UatomsGrpProxy;
extern CkVec <CProxy_EnergyGroup>          UegroupProxy;
extern CkVec <CProxy_StructFactCache>      UsfCacheProxy;
extern CkVec <CProxy_eesCache>             UeesCacheProxy;
extern CProxy_CPcharmParaInfoGrp   scProxy;

//----------------------------------------------------------------------------

void IntegrationComplete(void *, void *);

//#define _CP_DEBUG_PSI_OFF_
//#define _CP_ENERGY_GRP_VERBOSE_
//#define _CP_DEBUG_ATMS_
//#define _CP_DEBUG_ATMS_EXIT_

//==============================================================================




class GSAtmMsg: public CMessage_GSAtmMsg {
};

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** Constructor
 *
 *
 */
//==============================================================================
AtomsGrp::AtomsGrp(int n, int n_nl, int len_nhc_, int iextended_on_,int cp_min_opt_,
                   int cp_wave_opt_, int isokin_opt_,double kT_, Atom* a, AtomNHC *aNHC, UberCollection _thisInstance) : thisInstance(_thisInstance) {
//==============================================================================
// parameters, options and energies

    natm            = n;
    natm_nl         = n_nl;
    len_nhc         = len_nhc_;
    iextended_on    = iextended_on_;
    cp_min_opt      = cp_min_opt_;
    cp_wave_opt     = cp_wave_opt_;
    isokin_opt      = isokin_opt_;
    iteration       = 0;
    kT              = kT_;
    pot_ewd_rs      = 0.0;
    eKinetic        = 0.0;
    eKineticNhc     = 0.0;
    potNhc          = 0.0;    
    countAtm        = 0;
    acceptCountfu   = 0;
    acceptCountX    = 0;
    acceptCountu    = 0;
    atomsCMrecv=atomsPIMDXrecv=false;
//==============================================================================
// Initial positions, forces, velocities 
    numPIMDBeads    = config.UberImax;
    atoms           = new Atom[natm];
    atomsNHC        = new AtomNHC[natm];
    CmiMemcpy(atoms, a, natm * sizeof(Atom));           // atoms has no vectors
    for(int i=0;i<natm;i++){atomsNHC[i].Init(&aNHC[i]);}

    fastAtoms.natm = natm;
    fastAtoms.q    = (double *)fftw_malloc(natm*sizeof(double));
    fastAtoms.x    = (double *)fftw_malloc(natm*sizeof(double));
    fastAtoms.y    = (double *)fftw_malloc(natm*sizeof(double));
    fastAtoms.z    = (double *)fftw_malloc(natm*sizeof(double));
    fastAtoms.fx   = (double *)fftw_malloc(natm*sizeof(double));
    fastAtoms.fy   = (double *)fftw_malloc(natm*sizeof(double));
    fastAtoms.fz   = (double *)fftw_malloc(natm*sizeof(double));
    PIMD_CM_Atoms.natm = natm;
    PIMD_CM_Atoms.x    = (double *)fftw_malloc(natm*sizeof(double));
    PIMD_CM_Atoms.y    = (double *)fftw_malloc(natm*sizeof(double));
    PIMD_CM_Atoms.z    = (double *)fftw_malloc(natm*sizeof(double));

    ftot           = (double *)fftw_malloc((3*natm+2)*sizeof(double));

    zeroforces();
    if(iextended_on==1 && cp_min_opt==0){
       zeronhc();
    }//endif
    
//==============================================================================
// A copy of the atoms for fast routines


    copySlowToFast();
    double *qq = fastAtoms.q;
    for(int i=0;i<natm;i++){qq[i]=atoms[i].q;}

//==============================================================================
// Number of messages to be received when atoms are moved

    int nproc = CkNumPes();
    int div   = (natm / nproc);
    int rem   = (natm % nproc);

    nAtmMsgRecv = 0;
    for(int myid=0;myid<nproc;myid++){
      int natmNow = div;
      if(myid< rem){natmNow++;}
      if(natmNow>0){nAtmMsgRecv++;}
    }//endfor

//-----------------------------------------------------------------------------
  }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/* Destructor
 *
 *
 */
//==============================================================================
AtomsGrp::~AtomsGrp(){
     delete [] atoms;
     delete [] atomsNHC;
     fftw_free(fastAtoms.x);
     fftw_free(fastAtoms.y);
     fftw_free(fastAtoms.z);
     fftw_free(fastAtoms.fx);
     fftw_free(fastAtoms.fy);
     fftw_free(fastAtoms.fz);
     fftw_free(PIMD_CM_Atoms.x);
     fftw_free(PIMD_CM_Atoms.y);
     fftw_free(PIMD_CM_Atoms.z);



     fftw_free(ftot);
}
//==============================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/** trigger force computation
 * based on real forces available in each processor's chare then contribute to
 * global group reduction -> recvContribute
 */
void AtomsGrp::startRealSpaceForces(){
//==========================================================================
// Get the real space atom forces

   int myid   = CkMyPe();
   int nproc  = CkNumPes();
   pot_ewd_rs = 0.0;
   vself      = 0.0;
   vbgr       = 0.0;
   potPerdCorr= 0.0;

#ifndef  _CP_DEBUG_PSI_OFF_
#ifndef _CP_DEBUG_SCALC_ONLY_ 
   if(myid<natm-1){
     CPRSPACEION::CP_getionforce(natm,&fastAtoms,myid,nproc,&pot_ewd_rs,&vself,&vbgr,&potPerdCorr);
   }//endif
#endif
#endif

#ifdef _CP_DEBUG_ATMS_
   CkPrintf("GJM_DBG: calling contribute atm forces %d\n",myid);
#endif
   contributeforces();
   
//==========================================================================
  }//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void AtomsGrp::contributeforces(){
//==========================================================================

  int i,j;
  int myid     = CkMyPe();

  copyFastToSlow();
  for(i=0,j=0; i<natm; i++,j+=3){
    ftot[j]   = atoms[i].fx;
    ftot[j+1] = atoms[i].fy;
    ftot[j+2] = atoms[i].fz;
  }//endfor
  ftot[3*natm]  =pot_ewd_rs;
  ftot[3*natm+1]=potPerdCorr;

#ifdef _CP_DEBUG_ATMS_
  CkPrintf("GJM_DBG: inside contribute forces %d : %d\n",myid,natm);
#endif
  CkCallback cb(CkIndex_AtomsGrp::recvContribute(NULL), UatomsGrpProxy[thisInstance.proxyOffset]);
  contribute((3*natm+2)*sizeof(double),ftot,CkReduction::sum_double,cb);


//==========================================================================
  }//end routine
//==========================================================================

void AtomsGrp::init()
{
  // you are bead root if you are the processor at the base of the
  // GSP for your Instance, which we don't know until GSP is mapped
  // hence the 2nd phase init since atomsgrp is made before GSP
  const int offset=thisInstance.getPO();
  amBeadRoot      = (CkMyPe()==GSImaptable[offset].get(0,0));
  amZerothBead    = (thisInstance.idxU==0);
  if(numPIMDBeads>1)
    {
            
      // make lists of group IDs and their Commander PE
      
      CkGroupID *atomsgrpids= new CkGroupID[numPIMDBeads];
      // this const stuff is only here because its in the signature
      // for the section constructor.  Good luck initializing a
      // dynamic const ** array without this approach.
      const int **elems_c=  new const int*[numPIMDBeads];
      const int *naelems_c= new int[numPIMDBeads];
      int **elems=const_cast<int **>(elems);
      int *naelems= const_cast<int *>(naelems_c);
      for(int i=0;i<numPIMDBeads;i++)
	{
	  elems[i]= new int[1];
	  UberCollection instance=thisInstance;
	  naelems[i]=1;
	  instance.idxU.x=0;
	  int offset=instance.calcPO();
	  elems[i][0]=GSImaptable[offset].get(0,0);
	  atomsgrpids[i]=UatomsGrpProxy[offset].ckGetGroupID();
	}
      proxyHeadBeads=CProxySection_AtomsGrp(numPIMDBeads,atomsgrpids, elems_c, naelems_c);
    }


}

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/**
 * Apply forces to each processor's copy of the atoms.  This is parallelized so
 * that a subset of the atoms are computed on each processor and their results
 * broadcast to AtmGroup->acceptAtoms().  Move the atoms each processor is
 * responsible for. Set various energyGroup members.  contribute to group
 * reduction ->atomsDone
 */
void AtomsGrp::recvContribute(CkReductionMsg *msg) {
//==========================================================================
// Local pointers

  int i,j;
  double *ftot      = (double *) msg->getData();

  EnergyGroup *eg   = UegroupProxy[thisInstance.proxyOffset].ckLocalBranch();

  int nproc         = CkNumPes();
  int myid          = CkMyPe();
  int output_on     = config.atmOutput;

//============================================================
// Copy out the reduction of energy and forces

#ifdef _CP_DEBUG_ATMS_
  CkPrintf("GJM_DBG: inside recv forces %d : %d\n",myid,natm);
#endif
  double pot_ewd_rs_loc = ftot[3*natm];
  double potPerdCorrLoc = ftot[3*natm+1];
  double fmag = 0.0;
  for(i=0,j=0;i<natm;i++,j+=3){
    atoms[i].fx = ftot[j];    atoms[i].fy = ftot[j+1];
    atoms[i].fz = ftot[j+2];
    fmag += (ftot[j]*ftot[j]+ftot[j+1]*ftot[j+1]+ftot[j+2]*ftot[j+2]);
#ifdef _CP_DEBUG_ATMS_
    if(myid==0){
      CkPrintf("%d : %g %g %g\n",i,atoms[i].fx,atoms[i].fy,atoms[i].fz);
    }//endif
#endif
  }//endfor
  fmag /= (double)(3*natm);
  fmag  = sqrt(fmag);
  delete msg;

#ifdef _CP_DEBUG_ATMS_EXIT_
  if(myid==0){CkExit();}
#endif

//=================================================================
// Model forces 

#ifdef  _CP_DEBUG_PSI_OFF_
  double omega    = (0.0241888/15.0); // 15 fs^{-1}
  double omega2   = omega*omega;
  double pot_harm = 0.0;
  int npts        = 200;
  double sigma    = 1.0/sqrt(315777.0*atoms[0].m*omega2/300.0);
  double dx       = 6.0*sigma/(double)npts;

  if(iteration==0){
    px =  new double *[natm];
    for(i =0;i<natm;i++){
      px[i] = new double[npts];
      for(j=0;j<npts;j++){px[i][j]=0.0;}
    }//endfor
    for(i=0;i<natm;i++){
      atoms[i].x = 0.0;
      atoms[i].y = 0.0;
      atoms[i].z = 0.0;
    }//endfor
  }//endif

  for(i=0;i<natm;i++){
    atoms[i].fx = -omega2*atoms[i].x*atoms[i].m;
    atoms[i].fy = -omega2*atoms[i].y*atoms[i].m;;
    atoms[i].fz = -omega2*atoms[i].z*atoms[i].m;;
    pot_harm   += (atoms[i].m*omega2*(atoms[i].x*atoms[i].x+
                                      atoms[i].y*atoms[i].y+
                                      atoms[i].z*atoms[i].z));
  }//endfor
  pot_harm *= 0.5;
#endif

//=================================================================
// Check the forces

#ifdef _GLENN_CHECK_FORCES
   FILE *fp = fopen("forces.save","w");
   for(i=0;i<natm;i++){
     fprintf(fp,"%.12g %.12g %.12g\n",atoms[i].fx,atoms[i].fy,atoms[i].fz);
   }
   fclose(fp);
   CkExit();
#endif

//============================================================
// Tuck things that can be tucked.


  eg->estruct.eewald_real     = pot_ewd_rs_loc;  
  eg->estruct.fmag_atm        = fmag;

   if(numPIMDBeads>1)
     {
       send_PIMD_fx();
     }
   else
    {
      integrateAtoms();
    }
}


void AtomsGrp::integrateAtoms()
{

//============================================================
// Integrate the atoms

#ifdef _CP_DEBUG_ATMS_
  CkPrintf("GJM_DBG: Before atom integrate %d : %d\n",myid,natm);
#endif
  EnergyGroup *eg       = UegroupProxy[thisInstance.proxyOffset].ckLocalBranch();
  CPcharmParaInfo *sim  = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 

   double eKinetic_loc   = 0.0;
   double eKineticNhc_loc= 0.0;
   double potNhc_loc     = 0.0;
   int iwrite_atm        = 0;
   int myoutput_on       = 0;
   int nproc = CkNumPes();
   int myid          = CkMyPe();
   int div     = (natm / nproc);
   int rem     = (natm % nproc);
   int natmStr = div*myid;
   int natmNow = div;
   if(myid>=rem){natmStr += rem;}
   if(myid< rem){natmStr += myid;}
   if(myid< rem){natmNow++;}
   int natmEnd = natmNow+natmStr;

#ifdef  _CP_DEBUG_SCALC_ONLY_ 
  for(i=0;i<natm;i++){
    fastAtoms.x[i] = atoms[i].x;
    fastAtoms.y[i] = atoms[i].y;
    fastAtoms.z[i] = atoms[i].z;
    atoms[i].fx    = 0.0;
    atoms[i].fy    = 0.0;
    atoms[i].fz    = 0.0;
  }/*endfor*/
#endif

   ATOMINTEGRATE::ctrl_atom_integrate(iteration,natm,len_nhc,cp_min_opt,
                    cp_wave_opt,iextended_on,atoms,atomsNHC,myid,
                    &eKinetic_loc,&eKineticNhc_loc,&potNhc_loc,&iwrite_atm,
                    myoutput_on,natmNow,natmStr,natmEnd);

#ifdef _CP_DEBUG_SCALC_ONLY_ 
  for(i=0;i<natm;i++){
    atoms[i].x = fastAtoms.x[i];
    atoms[i].y = fastAtoms.y[i];
    atoms[i].z = fastAtoms.z[i];
  }/*endfor*/
#endif

//============================================================
  // Debug output

#ifdef  _CP_DEBUG_PSI_OFF_
   double etot_atm;
   if(isokin_opt==0){
     etot_atm = eKinetic_loc+eKineticNhc_loc+potNhc_loc+pot_harm;
     CkPrintf("iteration %d : tot energy %.12g on %d\n",iteration,etot_atm,myid);
   }else{
     etot_atm = eKineticNhc_loc+potNhc_loc;
     CkPrintf("iteration %d : tot energy %.12g %.12g on %d\n",iteration,etot_atm,
                                                    (eKinetic_loc+pot_harm),myid);
   }//endif
   for(i=natmStr;i<natmEnd;i++){
     int i1 = (atoms[i].x+3.0*sigma)/dx;
     int i2 = (atoms[i].y+3.0*sigma)/dx;
     int i3 = (atoms[i].z+3.0*sigma)/dx;
     i1     = (i1 >= 0    ? i1 : 0);
     i2     = (i2 >= 0    ? i2 : 0);
     i3     = (i3 >= 0    ? i3 : 0);
     i1     = (i1 < npts ? i1 : npts-1);
     i2     = (i2 < npts ? i2 : npts-1);
     i3     = (i3 < npts ? i3 : npts-1);
     px[i][i1] += 1.0;
     px[i][i2] += 1.0;
     px[i][i3] += 1.0;
   }//endif
   if(((iteration+1) % 1000)==0){
     for(i=natmStr;i<natmEnd;i++){
       char fname[100];
       sprintf(fname,"atmtest.dat.%d",i);
       double anorm = sqrt(2.0/(sigma*sigma*acos(-1.0)));
       FILE *fp = fopen(fname,"w");
       for(j =1;j<npts-1;j++){
         double x   = ((double)(2*j+1))*0.5*dx - 3.0*sigma;
         double ans = anorm*exp(-0.5*x*x/(sigma*sigma));
         double cnt = 3.0*((double)(iteration+1));
         fprintf(fp,"%g %g %g\n",x,(px[i][j]*anorm/cnt),ans);
       }//endfor       
       fclose(fp);
     }//endfor
   }//endif
#endif
   


//============================================================
// Get ready for the next iteration : 
//      zero forces, outputAtmEnergy, atomsDone

  zeroforces();

  if(cp_wave_opt==0 && cp_min_opt==0){
    if(natmNow>0){
      sendAtoms(eKinetic_loc,eKineticNhc_loc,potNhc_loc,natmNow,natmStr,natmEnd);
    }//endif
  }else{
    eKinetic                    = 0.0;
    eKineticNhc                 = 0.0;
    potNhc                      = 0.0;
    eg->estruct.eKinetic_atm    = 0.0;
    eg->estruct.eKineticNhc_atm = 0.0;
    eg->estruct.potNhc_atm      = 0.0;

    copySlowToFast();
    outputAtmEnergy();

    int i=0;
    CkCallback cb(CkIndex_AtomsGrp::atomsDone(NULL),UatomsGrpProxy[thisInstance.proxyOffset]);
    contribute(sizeof(int),&i,CkReduction::sum_int,cb);
  }//endif

//-------------------------------------------------------------------------
   }//end routine
//==========================================================================


//==========================================================================
// Atom energy output
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void AtomsGrp::outputAtmEnergy() {
//==========================================================================

  EnergyGroup *eg       = UegroupProxy[thisInstance.proxyOffset].ckLocalBranch();
  CPcharmParaInfo *sim  = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
  int myid              = CkMyPe();  
  double eKinetic       = eg->estruct.eKinetic_atm;
  double eKineticNhc    = eg->estruct.eKineticNhc_atm;
  double fmag           = eg->estruct.fmag_atm;
  double pot_ewd_rs_now = eg->estruct.eewald_real;
  double potNhc         = eg->estruct.potNhc_atm;
  double free_atm       = 3*((double)natm);
  int iperd             = sim->iperd;

  if(myid==0){
     if(iperd!=0){
       CkPrintf("{%d} EWALD_REAL  = %5.8lf\n",thisInstance.proxyOffset, pot_ewd_rs_now);
       CkPrintf("{%d} EWALD_SELF  = %5.8lf\n",thisInstance.proxyOffset,vself);
       CkPrintf("{%d} EWALD_BGR   = %5.8lf\n",thisInstance.proxyOffset,vbgr);
       if(iperd!=3){
         CkPrintf("{%d} EWALD_Perd  = %5.8lf\n",thisInstance.proxyOffset,potPerdCorr);
       }//endif
     }else{
       CkPrintf("{%d} ATM_COUL    = %5.8lf\n",thisInstance.proxyOffset,pot_ewd_rs_now);
     }//endif
     if(cp_min_opt==0){
        CkPrintf("{%d} atm eKin    = %5.8lf\n",thisInstance.proxyOffset,eKinetic);
        CkPrintf("{%d} atm Temp    = %5.8lf\n",thisInstance.proxyOffset,(2.0*eKinetic*BOLTZ/free_atm));
        CkPrintf("{%d} atm fmag    = %5.8lf\n",thisInstance.proxyOffset,fmag);
        if(iextended_on==1){
          double free_Nhc;
          if(isokin_opt==0){
            free_Nhc    = free_atm*((double)len_nhc);
          }else{
            free_Nhc    = free_atm*((double)(len_nhc-1));
	  }//endif
          CkPrintf("{%d} atm eKinNhc = %5.8lf\n",thisInstance.proxyOffset,eKineticNhc);
          CkPrintf("{%d} atm TempNHC = %5.8lf\n",thisInstance.proxyOffset,(2.0*eKineticNhc*BOLTZ/free_Nhc));
          CkPrintf("{%d} atm potNHC  = %5.8lf\n",thisInstance.proxyOffset,potNhc);
	}//endif
     }else{
        CkPrintf("{%d} atm fmag    = %5.8lf\n",thisInstance.proxyOffset,fmag);
     }//endif
  }//endif

//-------------------------------------------------------------------------
   }//end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/**
 * Increment the iteration counters in atoms and eesData call
 * gSpaceDriver->doneMovingAtoms() That is done via localbranch to all
 * co-located gSpacePlanes in the localBarrier scheme.  Used to be via
 * messages.  This permits the new step to advance with the new Psi.
 */
  void AtomsGrp::atomsDone(CkReductionMsg *msg) {
//==========================================================================
  delete msg;
  atomsDone();
}
//==========================================================================


//==========================================================================
// Needs to have each proc invoke directly doneMovingAtoms method of the
// GSpaceDrivers which are mapped to it. Without migration, we have that map
// at startup. With migration, one must write an enroll/dropout routine.
// All 
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void AtomsGrp::atomsDone() {
//==========================================================================
// Increment iteration counters

 int myid = CkMyPe();

 EnergyGroup *eg             = UegroupProxy[thisInstance.proxyOffset].ckLocalBranch();
 iteration++;
 eg->estruct.iteration_atm   = iteration;
 eg->iteration_atm           = iteration;

#ifdef _DEBUG_CHECK_ATOM_COMM_
 char fname[100];
 sprintf(fname,"atoms.out.%d.%d",iteration,myid);
 FILE *fp = fopen(fname,"w");
 for(int i=0;i<natm;i++){
   fprintf(fp,"%.10g %.10g %.10g\n",atoms[i].x,atoms[i].y,atoms[i].z);
 }//endfor
 fclose(fp);
#endif

//==========================================================================
// Use the cool new data caching system to say we're done.

 if(1) { // localAtomBarrier
   
   for(int kpoint=0; kpoint< config.UberJmax; kpoint++){ //each
							 //k-point
							 //needs to be
							 //handled
     
     UberCollection thisPoint=thisInstance;
     thisPoint.idxU.y=kpoint; // not at the gamma point
     thisPoint.setPO();
     eesCache *eesData = UeesCacheProxy[thisPoint.proxyOffset].ckLocalBranch ();
     int *indState     = eesData->gspStateInd;
     int *indPlane     = eesData->gspPlaneInd;
     int ngo           = eesData->nchareGSPProcT;
     
     GSAtmMsg *msg = new  GSAtmMsg;
     for(int i=0; i<ngo; i++){
       int iadd = UgSpacePlaneProxy[thisPoint.proxyOffset](indState[i],indPlane[i]).ckLocal()->registrationFlag;
       if(iadd!=1){
	 CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	 CkPrintf("atom : Bad registration cache flag on proc %d %d %d %d\n",
		  myid,iadd,indState[i],indPlane[i]);
	 CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	 CkExit();
       }//endif
       UgSpaceDriverProxy[thisPoint.proxyOffset](indState[i],indPlane[i]).doneMovingAtoms(iteration); 
     }//endfor
   }//endfor

 }
 /*
  else{

   if(myid==0){
      GSAtmMsg *msg = new (8*sizeof(int)) GSAtmMsg;
      CkSetQueueing(msg, CK_QUEUEING_IFIFO);
      *(int*)CkPriorityPtr(msg) = config.sfpriority-10;
      UgSpaceDriverProxy[thisInstance.proxyOffset].doneMovingAtoms(iteration);
   }//endif
 }//endif
 */
}//end routine




//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void AtomsGrp::sendAtoms(double eKinetic_loc,double eKineticNhc_loc,double potNhc_loc,
                         int natmNow,int natmStr,int natmEnd){
//==========================================================================
// Malloc the message

  int nsize    = 9*natmNow+3;
  AtomMsg *msg = new (nsize,8*sizeof(int)) AtomMsg;
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
  *(int*)CkPriorityPtr(msg) = config.sfpriority-10;

  double *atmData = msg->data;
  msg->nsize      = nsize;
  msg->natmStr    = natmStr;
  msg->natmEnd    = natmEnd;

//==========================================================================
// pack atom positions : new for use : old for output
  if(numPIMDBeads>1)
  for(int i=natmStr,j=0;i<natmEnd;i++,j+=9){
    atmData[(j)  ]=atoms[i].x;
    atmData[(j+1)]=atoms[i].y;
    atmData[(j+2)]=atoms[i].z;
    atmData[(j+3)]=atoms[i].xold;
    atmData[(j+4)]=atoms[i].yold;
    atmData[(j+5)]=atoms[i].zold;
    atmData[(j+6)]=atoms[i].vxold;
    atmData[(j+7)]=atoms[i].vyold;
    atmData[(j+8)]=atoms[i].vzold;
  }//endfor
  else
    {
      for(int i=natmStr,j=0;i<natmEnd;i++,j+=9){
	atmData[(j)  ]=atoms[i].xu;
	atmData[(j+1)]=atoms[i].yu;
	atmData[(j+2)]=atoms[i].zu;
	atmData[(j+3)]=atoms[i].xuold;
	atmData[(j+4)]=atoms[i].yuold;
	atmData[(j+5)]=atoms[i].zuold;
	atmData[(j+6)]=atoms[i].vxuold;
	atmData[(j+7)]=atoms[i].vyuold;
	atmData[(j+8)]=atoms[i].vzuold;
      }//endfor

    }
//==========================================================================
// pack the 3 energies

   atmData[(nsize-3)] = eKinetic_loc;
   atmData[(nsize-2)] = eKineticNhc_loc;
   atmData[(nsize-1)] = potNhc_loc;

//==========================================================================
// Send the message to everyone in the group

  UatomsGrpProxy[thisInstance.proxyOffset].acceptAtoms(msg);

//-------------------------------------------------------------------------
  }//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/**
 * Take packed message of a chunk of the atoms with updated positions. Update
 * local copy of atoms.  Update local energyGroup members.  Print atom energies
 * when we have all of them.  Do file output of atoms if desired.
 */
  void AtomsGrp::acceptAtoms(AtomMsg *msg) {
//==========================================================================

  AtomsGrp *ag      = UatomsGrpProxy[thisInstance.proxyOffset].ckLocalBranch();
  double *atmData   = msg->data;
  int    nsize      = msg->nsize;
  int    natmStr    = msg->natmStr;
  int    natmEnd    = msg->natmEnd;

//==========================================================================
// unpack atom position and velocity
  if(numPIMDBeads>1)
  for(int i=natmStr,j=0;i<natmEnd;i++,j+=9){
    atoms[i].x     = atmData[(j)  ];
    atoms[i].y     = atmData[(j+1)];
    atoms[i].z     = atmData[(j+2)];
    atoms[i].xold  = atmData[(j+3)];
    atoms[i].yold  = atmData[(j+4)];
    atoms[i].zold  = atmData[(j+5)];
    atoms[i].vxold = atmData[(j+6)];
    atoms[i].vyold = atmData[(j+7)];
    atoms[i].vzold = atmData[(j+8)];
  }//endfor
  else
    for(int i=natmStr,j=0;i<natmEnd;i++,j+=9){
    atoms[i].xu     = atmData[(j)  ];
    atoms[i].yu     = atmData[(j+1)];
    atoms[i].zu     = atmData[(j+2)];
    atoms[i].xuold  = atmData[(j+3)];
    atoms[i].yuold  = atmData[(j+4)];
    atoms[i].zuold  = atmData[(j+5)];
    atoms[i].vxuold = atmData[(j+6)];
    atoms[i].vyuold = atmData[(j+7)];
    atoms[i].vzuold = atmData[(j+8)];
    }//endfor
  //endif

//==========================================================================
// unpack energy

  countAtm++;

  if(countAtm==1){
    eKinetic    = 0;
    eKineticNhc = 0;
    potNhc      = 0;
  }//endif

  eKinetic    += atmData[(nsize-3)];
  eKineticNhc += atmData[(nsize-2)];  
  potNhc      += atmData[(nsize-1)];

//==========================================================================
// Delete the message

  delete msg;

//==========================================================================
// Copy to the fast vectors and phone home

  if(countAtm==nAtmMsgRecv){
     countAtm = 0;

     EnergyGroup *eg             = UegroupProxy[thisInstance.proxyOffset].ckLocalBranch();
     eg->estruct.eKinetic_atm    = eKinetic;
     eg->estruct.eKineticNhc_atm = eKineticNhc;
     eg->estruct.potNhc_atm      = potNhc;

     copySlowToFast();
     outputAtmEnergy();

     // iteration is time of atoms[i].xold
     // maxIter is 1 more than you need : slightly annoying but livable
     int output_on = config.atmOutput;
     if(output_on==1 && iteration<=config.maxIter-1){ 
       int pi_beads   = 1;
       int iwrite_atm = 0;
       int myid       = CkMyPe();
       ATOMOUTPUT::ctrl_piny_output(iteration,natm,len_nhc,pi_beads,myid,atoms,atomsNHC,
                                    &iwrite_atm,output_on);
       if(myid==0 && iwrite_atm>0){
         CkPrintf("-----------------------------------\n");
         CkPrintf("Writing atoms to disk at time %d\n",iteration);
         CkPrintf("-----------------------------------\n");
       }//endif
     }//endif

     if(numPIMDBeads>1)
       {
	 // transform PIMD U to X
	 if(amBeadRoot)
	   send_PIMD_u();
	 if(amBeadRoot && amZerothBead)
	   {
	     BeadCMMsg *msg = new (natm,natm,natm) BeadCMMsg;
	     for(int atomI=0;atomI<natm;atomI++)
	       {
		 msg->x[atomI]=atoms[atomI].xu;
		 msg->y[atomI]=atoms[atomI].yu;
		 msg->z[atomI]=atoms[atomI].zu;
	       }
	     CkAbort("make this into a message");
	     proxyHeadBeads.accept_PIMD_CM(msg);
	   }
       }
     else
       {
	 //everybody has to have received all the atoms before continuing : not just me
	 int i=0;
	 CkCallback cb(CkIndex_AtomsGrp::atomsDone(NULL),UatomsGrpProxy[thisInstance.proxyOffset]);
	 contribute(sizeof(int),&i,CkReduction::sum_int,cb);
       }
  }//endif

//-------------------------------------------------------------------------
  }//end routine
//==========================================================================

void AtomsGrp::accept_PIMD_CM(BeadCMMsg *msg)
{
  for(int atomnum=0; atomnum<natm; atomnum++)
    {
      PIMD_CM_Atoms.x[atomnum]=msg->x[atomnum];
      PIMD_CM_Atoms.y[atomnum]=msg->y[atomnum];
      PIMD_CM_Atoms.z[atomnum]=msg->z[atomnum];
    }
  delete msg;
  atomsCMrecv=true;
  if(atomsPIMDXrecv)
    {
      int i=0;
      CkCallback cb(CkIndex_AtomsGrp::atomsDone(NULL),UatomsGrpProxy[thisInstance.proxyOffset]);
      contribute(sizeof(int),&i,CkReduction::sum_int,cb);
      atomsCMrecv=atomsPIMDXrecv=false;
    }

}


void AtomsGrp::send_PIMD_u()
{
  for(int atomnum=0;atomnum<natm;atomnum++)
    {
      UPIBeadAtomsProxy[thisInstance.proxyOffset][atomnum].accept_PIMD_u(atoms[atomnum].xu,atoms[atomnum].yu,atoms[atomnum].zu, PIBeadIndex);
    }
}

void AtomsGrp::send_PIMD_fx()
{ 
  for(int atomnum=0;atomnum<natm;atomnum++)
    {
      UPIBeadAtomsProxy[thisInstance.proxyOffset][atomnum].accept_PIMD_Fx(atoms[atomnum].fx,atoms[atomnum].fy,atoms[atomnum].fz, PIBeadIndex);
    }
}

// is broadcast to us
void AtomsGrp::accept_PIMD_fu(double _fxu, double _fyu, double _fzu, int atomI)
{

  atoms[atomI].fxu=_fxu;
  atoms[atomI].fyu=_fyu;
  atoms[atomI].fzu=_fzu;
  fastAtoms.fxu[atomI]=_fxu;
  fastAtoms.fyu[atomI]=_fyu;
  fastAtoms.fzu[atomI]=_fzu;
  acceptCountfu++;
  if(acceptCountfu==natm)
    {
      integrateAtoms();
      acceptCountfu=0;
    }
}

void AtomsGrp::accept_PIMD_x(double _x, double _y, double _z, int atomI)
{
  
  atoms[atomI].x=_x;
  atoms[atomI].y=_y;
  atoms[atomI].z=_z;
  fastAtoms.x[atomI]=_x;
  fastAtoms.y[atomI]=_y;
  fastAtoms.z[atomI]=_z;
  acceptCountX++;
  if(acceptCountX==natm)
    {
      acceptCountX=0;
      atomsPIMDXrecv=true;
      if(atomsCMrecv)
	{
	  int i=0;
	  CkCallback cb(CkIndex_AtomsGrp::atomsDone(NULL),UatomsGrpProxy[thisInstance.proxyOffset]);
	  contribute(sizeof(int),&i,CkReduction::sum_int,cb);
	  atomsCMrecv=atomsPIMDXrecv=false;
	}

    }
}


// done during initialization in 1st iteration
void AtomsGrp::send_PIMD_x()
{
  for(int atomnum=0;atomnum<natm;atomnum++)
    {
      UPIBeadAtomsProxy[thisInstance.proxyOffset][atomnum].accept_PIMD_x(atoms[atomnum].x,atoms[atomnum].y,atoms[atomnum].z, PIBeadIndex);
    }
}

// done during initialization in 1st iteration
void AtomsGrp::accept_PIMD_u(double _xu, double _yu, double _zu, int atomI)
{

  atoms[atomI].xu=_xu;
  atoms[atomI].yu=_yu;
  atoms[atomI].zu=_zu;
  fastAtoms.xu[atomI]=_xu;
  fastAtoms.yu[atomI]=_yu;
  fastAtoms.zu[atomI]=_zu;
  acceptCountu++;
  if(acceptCountu==natm)
    {
      integrateAtoms();
      acceptCountu=0;
    }
}



//==========================================================================
//Energy group that can retrieve the energies from
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
EnergyGroup::EnergyGroup (UberCollection _thisInstance) : thisInstance(_thisInstance) {
    iteration_gsp = 0;
    iteration_atm = 0;
    kpointEnergyDoneCount=0;
    //non local
    estruct.enl = 0;

    //local external energy
    estruct.eext = 0;
    estruct.eke = 0;

    //hartree energy
    estruct.ehart = 0;
    
    //ion-ion
    estruct.eewald_recip = 0;
    estruct.eewald_real = 0;
    
    //exchange correlation
    estruct.egga = 0;
    estruct.eexc = 0;
    estruct.fmagPsi = 0;

    //CP Fict KE
    estruct.fictEke = 0;

    // total electronic part
    estruct.totalElecEnergy = 0; // needs ewald_real to be physical
    estruct.iteration_gsp = 0;

 
    // atm stuff
    estruct.eKinetic_atm    = 0;    // classical kinetic energy
    estruct.eKineticNhc_atm = 0; // NHC kinetic energy
    estruct.potNhc_atm      = 0;      // NHC pot energy
    estruct.fmag_atm        = 0;        // magnitude of atm forces
    estruct.iteration_atm   = 0;

//-------------------------------------------------------------------------
  } //end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/**
 * CP_StateGspacePlane(0,0) calls this to replicate the energies everywhere for
 * consistency and tolerance checking.
 */
void EnergyGroup::updateEnergiesFromGS(EnergyStruct &es, UberCollection sender) {
//==========================================================================

  if(config.UberJmax>1)
    {// need to sum enl and eke across kpoints
      estruct.enl          += es.enl;
      estruct.eke          += es.eke;
      estruct.fictEke      += es.fictEke;
      estruct.fmagPsi      += es.fmagPsi;
    }
  else
    {
      estruct.enl          = es.enl;
      estruct.eke          = es.eke;
      estruct.fictEke      = es.fictEke;
      estruct.fmagPsi      = es.fmagPsi;
    }
  // these other stuff comes from rho and we only have one of them for
  // all kpoints
      estruct.eext         = es.eext;
      estruct.ehart        = es.ehart;
      estruct.eewald_recip = es.eewald_recip;
      estruct.egga         = es.egga;
      estruct.eexc         = es.eexc;
      // these are gspace things that I don't know what to do with in
      // the kpoint>1 case
   
      estruct.totalElecEnergy  = es.totalElecEnergy;
      estruct.fmagPsi      = es.fmagPsi;
      estruct.iteration_gsp= es.iteration_gsp;
      iteration_gsp        = es.iteration_gsp;
#ifdef _DEBUG_ESTRUCT_
       CkPrintf("Energies received %lf, %lf, %lf, %lf, %lf\n", 
                 estruct.enl,estruct.eke,estruct.eext,estruct.ehart, 
                 estruct.eewald_recipt,estruct.egga,estruct.eexc,
                 estruct.fictEke,estruct.totalEnergy);
#endif
    int i=0;
    CkCallback cb(CkIndex_EnergyGroup::energyDone(NULL),UegroupProxy[thisInstance.proxyOffset]);
    contribute(sizeof(int),&i,CkReduction::sum_int,cb);

//-------------------------------------------------------------------------
  }//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/**
 * call GSpaceDriver->doneComputingEnergy() on all co-located gspace chares
 * which allows the new step to advance with new psi
 */
  void EnergyGroup::energyDone(CkReductionMsg *msg) {
//==========================================================================
   delete msg;
   // need to receive one of these for each k-point
   if(++kpointEnergyDoneCount==config.UberJmax)
     {
       kpointEnergyDoneCount=0;
       energyDone();
     }
}
//==========================================================================


//==========================================================================
// Needs to have each proc invoke directly the doneComputingEnergy method of the
// GSpaceDrivers which are mapped to it. Without migration, we have that map
// at startup. With migration, one must write an enroll/dropout routine.
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void EnergyGroup::energyDone(){
//==========================================================================
// Use the cool new data caching system

 int myid          = CkMyPe();
 if(1) { // localAtomBarrier
   for(int kpoint=0; kpoint< config.UberJmax; kpoint++){ //each
							 //k-point
							 //needs to be
							 //handled
     
     UberCollection thisPoint=thisInstance;
     thisPoint.idxU.y=kpoint; // not at the gamma point
     thisPoint.setPO();
     
     eesCache *eesData = UeesCacheProxy[thisPoint.proxyOffset].ckLocalBranch ();
     int *indState     = eesData->gspStateInd;
     int *indPlane     = eesData->gspPlaneInd;
     int ngo           = eesData->nchareGSPProcT;
     GSAtmMsg *msg = new  GSAtmMsg;
     for(int i=0; i<ngo; i++){
       int iadd = UgSpacePlaneProxy[thisPoint.proxyOffset](indState[i],indPlane[i]).ckLocal()->registrationFlag;
       if(iadd!=1){
	 CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	 CkPrintf("Energy : Bad registration cache flag on proc %d %d %d %d\n",
		  myid,iadd,indState[i],indPlane[i]);
	 CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	 CkExit();
       }//endif
       UgSpaceDriverProxy[thisPoint.proxyOffset](indState[i],indPlane[i]).doneComputingEnergy(iteration_atm); 
     }//endfor
   }//endfor
   }
 /*
  else{


   if(myid==0){
      GSAtmMsg *msg = new (8*sizeof(int)) GSAtmMsg;
      CkSetQueueing(msg, CK_QUEUEING_IFIFO);
      *(int*)CkPriorityPtr(msg) = config.sfpriority-10;
      UgSpaceDriverProxy[thisInstance.proxyOffset].doneComputingEnergy(iteration_atm);
   }//endif
 
 }//endif
 */
}//end routine




/*//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
EnergyStruct GetEnergyStruct() {
    return UegroupProxy[thisInstance.proxyOffset].ckLocalBranch()->getEnergyStruct();
}
//==========================================================================
*/

#include "groups.def.h"

