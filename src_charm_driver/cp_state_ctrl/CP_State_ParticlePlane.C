/*
 * Life-cycle of a CP_State_ParticlePlane:
 *
 * The particle-plane array is a shadow of the g-space planes. This means
 * that particlePlaneProxy(s, p) exists on the same processor as
 * gSpacePlaneProxy(s, p). Also,  both chare arrays have the same number of
 * objects.
 *
 * The life-cycle of this object involves computation of a portion of
 * the Z-matrix, reduction of the Z matrix over g-space and computing
 * forces using the particle data.  ComputeZ is triggered by the
 * arrival of the structure factor at the local sfcache.
 *
 * The computation of the Z-matrix is done in the method computeZ().
 * Immediately after this, the matrix is reduced using the reduceZ() method.
 * Through the reduceZ() method, a portion of the reduced matrix is 
 * available at the zeroth plane of each state.
 * 
 * The reduced matrix is sent from the particle-plane 0 to all the particle
 * planes in the state using the getForces() method, which computes
 * forces and adds these forces to the forces compputed using electron
 * density in the g-space planes. 

 * Important information:
 * 1. The dimensionality of the Z-matrix depends upon how may orbital levels
 * are considered (l = 0,1...). The current code assumes that only l = 0
 * is allowed. This has to be changed.
 * 
 * 2. The special point gx=gy=gz=0 is not dealt with.
 */ 
#include "charm++.h"
#include "util.h"
#include "groups.h"
#include "cpaimd.h"
#include "sim_subroutines.h"
#include "CP_State_Plane.h"
#include "StructFactorCache.h"
#include "ckmulticast.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cplocal.h"

extern CProxy_CP_State_GSpacePlane gSpacePlaneProxy;
extern CProxy_AtomsGrp atomsGrpProxy;
extern CProxy_CPcharmParaInfoGrp scProxy;
extern CProxy_CP_State_ParticlePlane particlePlaneProxy;
extern CProxy_StructFactCache sfCacheProxy;
extern int nstates;
extern int nplane_x;
extern Config config;
extern CkGroupID mCastGrpId;
//============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void printEnl(void *param, void *msg){
  CkReductionMsg *m=(CkReductionMsg *)msg;
  double d = ((double *)m->getData())[0];
  CkPrintf("ENL         = %5.8lf\n", d);
  delete m;

  gSpacePlaneProxy(0,0).computeEnergies(ENERGY_ENL, d);  
}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_State_ParticlePlane::CP_State_ParticlePlane(int x, int y, int z, 
                                               int numPlanes, int numSfGrps_in,
                                               int natm_nl_in, int natm_nl_grp_max_in){


	sizeX                = x;
	sizeY                = y;
	sizeZ                = z;
	gSpacePlanesPerChare = numPlanes;
        numSfGrps            = numSfGrps_in;
        natm_nl              = natm_nl_in;
        natm_nl_grp_max      = natm_nl_grp_max_in;
	count = new int[numSfGrps];
	bzero(count,numSfGrps*sizeof(int));
	doneEnl = 0;
	doneForces = 0;
	enl   = 0.0;
	energy_count = 0;
	zmatrixSum = NULL;
	zmatrix = NULL;
	totalEnergy = 0.0;
	reductionPlaneNum = 0;
	CkAssert(reductionPlaneNum % config.gSpacePPC == 0);
	contribute(sizeof(int), &sizeX, CkReduction::sum_int);
        setMigratable(false);
	usesAtSync = CmiFalse;

	//DY NONLOCAL
        zmatrix_fx    = NULL; zmatrix_fy    = NULL; zmatrix_fz    = NULL;
        zmatrixSum_fx = NULL; zmatrixSum_fy = NULL; zmatrixSum_fz = NULL;
	//DY NONLOCAL

        //Specifies the end of the particle plane iteration
        //It is set to 1 in reduceZ.
        //After this only GSpacePlane will begin the next iteration.
  //register with the SFCache
  StructFactCache *sfcache = sfCacheProxy.ckLocalBranch();
  for(int i=0;i<numSfGrps;i++)
    sfcache->registerPP(thisIndex.x, thisIndex.y,i);

  // Create section proxy for ENL reduction ParticlePlane (any state#, reductionPlaneNum)
  if(thisIndex.x==0 && thisIndex.y==reductionPlaneNum)
    {

      CkArrayIndexMax *elems = new CkArrayIndexMax[nstates];

      CkArrayIndex2D idx(0, reductionPlaneNum);  // plane# = this plane#
      for (int j = 0; j < nstates; j++) {
	idx.index[0] = j;
	elems[j] = idx;
      }

      particlePlaneENLProxy = 
	CProxySection_CP_State_ParticlePlane::ckNew(particlePlaneProxy.ckGetArrayID(),
						    elems, nstates);
      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
      particlePlaneENLProxy.ckDelegate(mcastGrp);
      EnlCookieMsg *emsg= new EnlCookieMsg;
      mcastGrp->setSection(particlePlaneENLProxy);
      particlePlaneENLProxy.setEnlCookie(emsg);
    }

//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
// this function is called from CP_State_GSpacePlane::initGSpace()
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void
CP_State_ParticlePlane::initKVectors(GStateSlab *gss)
{
    gss->setKVectors(&gSpaceNumPoints, &k_x, &k_y, &k_z);
    myForces = new complex[gSpaceNumPoints]; // forces on coefs
    //memset(myForces, 0, gSpaceNumPoints * sizeof(complex));
}
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_State_ParticlePlane::~CP_State_ParticlePlane()
{
	delete [] myForces;
	delete [] k_x;
	delete [] k_y;
	delete [] k_z;
	//delete [] zmatrixSum;

        delete [] zmatrix; //FOO-BAR should i delete this?            
        delete [] zmatrix_fx; 
        delete [] zmatrix_fy; 
        delete [] zmatrix_fz; 
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::pup(PUP::er &p)
{
	ArrayElement2D::pup(p);
	p|gSpaceNumPoints;
	p|doneGettingForces;
	p|numSfGrps;
        p|natm_nl;
	p|natm_nl_grp_max;
	int zsize=numSfGrps*natm_nl_grp_max;
	if (p.isUnpacking()) {
		k_x = new int[gSpaceNumPoints];
		k_y = new int[gSpaceNumPoints];
		k_z = new int[gSpaceNumPoints];
		myForces = new complex[gSpaceNumPoints]; // forces on coefs
                memset(myForces, 0, gSpaceNumPoints * sizeof(complex));
		// these 
		zmatrixSum = new complex[zsize];
                zmatrixSum_fx = new complex[zsize];
                zmatrixSum_fy = new complex[zsize];
                zmatrixSum_fz = new complex[zsize];
		zmatrix = new complex[zsize];
                zmatrix_fx = new complex[zsize];
                zmatrix_fy = new complex[zsize];
                zmatrix_fz = new complex[zsize];
		count= new int[numSfGrps];
	}
	p(k_x, gSpaceNumPoints);
	p(k_y, gSpaceNumPoints);
	p(k_z, gSpaceNumPoints);
	p(count,numSfGrps);
	p(zmatrixSum,zsize);
	p(zmatrixSum_fx,zsize);
	p(zmatrixSum_fy,zsize);
	p(zmatrixSum_fz,zsize);
	p(zmatrix,zsize);
	p(zmatrix_fx,zsize);
	p(zmatrix_fy,zsize);
	p(zmatrix_fz,zsize);
	p|reductionPlaneNum;
	// Dont pack zmatrixSum
	p|sizeX;
	p|sizeY;
	p|sizeZ;
	p|gSpacePlanesPerChare;
	p|zsize;
	p|energy_count;
	p|totalEnergy;
	p|enl;
	p|doneEnl;
	p|doneForces;
        p|enl_total;
        p|reductionPlaneNum;
	p|enlCookie;
	p|particlePlaneENLProxy;
	// "gspace" is not pup'ed since it is always assigned to
	

}
//============================================================================


extern void CmiReference(void *blk);



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void
CP_State_ParticlePlane::computeZ(PPDummyMsg *m)
{
    doneGettingForces = false;    
//    doneForces=0;
    CP_State_GSpacePlane *gsp = 
	gSpacePlaneProxy(thisIndex.x, thisIndex.y).ckLocal();
    GStateSlab *gss = &(gsp->gs);
    int atmIndex=m->atmGrp;
    int sfindex=m->sfindex;
    delete m;
  
    StructFactCache *sfcache = sfCacheProxy.ckLocalBranch();

    complex *structureFactor;
    complex *structureFactor_fx;
    complex *structureFactor_fy;
    complex *structureFactor_fz;
    sfcache->getStructFact(thisIndex.y, atmIndex, &structureFactor, &structureFactor_fx, &structureFactor_fy, &structureFactor_fz);
    zsize = 0;
    AtomsGrp *ag = atomsGrpProxy.ckLocalBranch(); // find me the local copy
    zsize = natm_nl_grp_max*numSfGrps;

    if(!zmatrix){
	zmatrix    = new complex[zsize];
	zmatrix_fx = new complex[zsize];
	zmatrix_fy = new complex[zsize];
	zmatrix_fz = new complex[zsize];
	memset(zmatrix, 0, sizeof(complex)*zsize);
	memset(zmatrix_fx, 0, sizeof(complex)*zsize);
	memset(zmatrix_fy, 0, sizeof(complex)*zsize);
	memset(zmatrix_fz, 0, sizeof(complex)*zsize);
    }

    int mydoublePack = config.doublePack;
    int zoffset=natm_nl_grp_max * atmIndex;
    CPNONLOCAL::CP_enl_matrix_calc(gSpaceNumPoints,gss->packedPlaneData, k_x,k_y,k_z, 
				   structureFactor,structureFactor_fx,structureFactor_fy,
				   structureFactor_fz,&zmatrix[zoffset],&zmatrix_fx[zoffset],&zmatrix_fy[zoffset],
				   &zmatrix_fz[zoffset],thisIndex.x,mydoublePack,numSfGrps,atmIndex);

    // reduce zmatrices over the planes of each state
    thisProxy(thisIndex.x, reductionPlaneNum).reduceZ(natm_nl_grp_max, atmIndex, &zmatrix[zoffset],
						      &zmatrix_fx[zoffset],&zmatrix_fy[zoffset],&zmatrix_fz[zoffset]);


//----------------------------------------------------------------------------
}
//============================================================================


//============================================================================
// this function achieves the reduction over gspace for a particular state
//============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::reduceZ(int size, int atmIndex, complex *zmatrix_,
        complex *zmatrix_fx_,complex *zmatrix_fy_,complex *zmatrix_fz_)
{

  // This method is invoked only CP_State_ParticlePlane(*,reductionPlaneNum)
  count[atmIndex]++;
  int i,j;
  int zsize = natm_nl_grp_max * numSfGrps;
  int zoffset= natm_nl_grp_max * atmIndex;
  if (!zmatrixSum) {
      zmatrixSum = new complex[zsize];
      memset(zmatrixSum, 0, zsize * sizeof(complex));
  }

  if (!zmatrixSum_fx) {
      zmatrixSum_fx = new complex[zsize];
      memset(zmatrixSum_fx, 0, zsize * sizeof(complex));
  }

  if (!zmatrixSum_fy) {
      zmatrixSum_fy = new complex[zsize];
      memset(zmatrixSum_fy, 0, zsize * sizeof(complex));
  }

  if (!zmatrixSum_fz) {
      zmatrixSum_fz = new complex[zsize];
      memset(zmatrixSum_fz, 0, zsize * sizeof(complex));
  }
  
  for (i = 0,j=zoffset; i < size; i++,j++){
    zmatrixSum[j]    += zmatrix_[i];
    zmatrixSum_fx[j] += zmatrix_fx_[i];
    zmatrixSum_fy[j] += zmatrix_fy_[i];
    zmatrixSum_fz[j] += zmatrix_fz_[i];
  }

  if ( ((count[atmIndex] == (sizeX/gSpacePlanesPerChare)/2 - 1)&&(!config.doublePack)) || 
       ((count[atmIndex] == config.low_x_size)&&(config.doublePack)) ) { 
    count[atmIndex] = 0;
    if(doneEnl==0)
      enl=0.0;
    doneEnl++;

    // Now we have contributions from all the other gspace planes
    // in this state. Calculate the energies and atoms forces.
    // Normalize ZmatrixSum

    AtomsGrp *ag = atomsGrpProxy.ckLocalBranch();
    Atom *atoms  = ag->atoms;
    int mydoublePack = config.doublePack;
    double myenl;
    CPNONLOCAL::CP_enl_atm_forc_calc(numSfGrps,atmIndex,atoms,&zmatrixSum[zoffset],
	     &zmatrixSum_fx[zoffset],&zmatrixSum_fy[zoffset],&zmatrixSum_fz[zoffset],&myenl,mydoublePack);
    enl+=myenl;

#ifdef _CP_DEBUG_NLMAT_
    for (i = 0; i < zsize; i++){
      CkPrintf("Non-local matrix[%d][%d]={%g, %g}\n",i,thisIndex.x,
	       zmatrix_[i].re,zmatrix_[i].im );}
#endif

    // send the reduced zmatrix to the particleplanes of this state
    // that are non-zero and compute the non-local forces

    for (i = 0; i < nplane_x; i ++){
         thisProxy(thisIndex.x, i).getForces(size,atmIndex, &zmatrixSum[zoffset]);
    }

    // reduce enl over all states
    //    CkPrintf(" PP [%d %d] doneEnl %d\n",thisIndex.x, thisIndex.y,doneEnl);
    if(doneEnl==numSfGrps)
      {
	CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
	CkCallback cb=CkCallback(printEnl, NULL);
	mcastGrp->contribute(sizeof(double),(void*) &enl, 
		   CkReduction::sum_double, enlCookie, cb);
	doneEnl=0;
	enl=0.0;
	memset(zmatrixSum, 0, zsize * sizeof(complex));
	memset(zmatrixSum_fx, 0, zsize * sizeof(complex));
	memset(zmatrixSum_fy, 0, zsize * sizeof(complex));
	memset(zmatrixSum_fz, 0, zsize * sizeof(complex));
      }


  }/*endif*/

//---------------------------------------------------------------------------
   }
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void 
CP_State_ParticlePlane::sumEnergies(double energy_) 
{
        
        enl_total   += energy_;
	totalEnergy += energy_;
	energy_count++;

	if (energy_count == nstates) {
            energy_count = 0;
            totalEnergy = 0;
            enl_total = 0;
	}
//---------------------------------------------------------------------------
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::getForces(int zmatSize, int atmIndex, complex *zmatrixSum_loc)
{
  doneForces++;
  CP_State_GSpacePlane *gsp = 
    gSpacePlaneProxy(thisIndex.x, thisIndex.y).ckLocal();
  GStateSlab *gss = &(gsp->gs);
  AtomsGrp *ag = atomsGrpProxy.ckLocalBranch();
  int state_ind = gss->istate_ind;                

  // compute the non-local forces using the qc_subroutine
  int mydoublePack = config.doublePack;
  StructFactCache *sfcache = sfCacheProxy.ckLocalBranch();
  complex *structureFactor,*structureFactor_fx,*structureFactor_fy,*structureFactor_fz;
  sfcache->getStructFact(thisIndex.y, atmIndex, &structureFactor, &structureFactor_fx, &structureFactor_fy, &structureFactor_fz);

  if (gSpaceNumPoints != 0)
    CPNONLOCAL::CP_enl_force_calc(zmatrixSum_loc,
				  gSpaceNumPoints, k_x, k_y, k_z, structureFactor,
				  myForces, state_ind, mydoublePack, numSfGrps, atmIndex);

  // forces have been computed, send them to the corresponding
  // g-space plane
  if(doneForces==numSfGrps)
    {
      doneGettingForces = true;
      if (gsp->doneDoingIFFT) {
	gsp->getForcesAndIntegrate();
      }
      doneForces=0;
    }
  //----------------------------------------------------------------------------
}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::setEnlCookie(EnlCookieMsg *m){
  CkGetSectionInfo(enlCookie,m);
  delete m;
//----------------------------------------------------------------------------
}
//============================================================================
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_ParticlePlane::ResumeFromSync() 
{
  if(thisIndex.x==0 && thisIndex.y==reductionPlaneNum)
    {
      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();       
      mcastGrp->resetSection(particlePlaneENLProxy);
    }
//----------------------------------------------------------------------------
}
//============================================================================
