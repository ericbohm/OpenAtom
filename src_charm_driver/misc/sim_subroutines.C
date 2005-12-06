//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** \file sim_subroutines.C
 * Add functions to allow application programmers to initialize these and
 * the corresponding functions in Charm++ to call these functions with
 * appropriate parameters 
 */
//==============================================================================

#include "util.h"
#include "../../include/debug_flags.h"
#include <math.h>
#include "cpaimd.h"
#include "groups.h"
#include "sim_subroutines.h"
#include "CP_State_Plane.h"

//==============================================================================

extern CProxy_FFTcache fftCacheProxy;
extern Config config;
extern int nstates;
extern int sizeX;
extern CProxy_AtomsGrp atomsGrpProxy;
extern CProxy_CPcharmParaInfoGrp scProxy;

//==============================================================================

//============================================================================
/**
 * split up an fft call into multiple invocations with cmiprogress
 * calls between them.
 */
//============================================================================
void fft_split(fftw_plan plan, int howmany, fftw_complex *in, int istride,
		 int idist, fftw_complex *out, int ostride, int odist)
{

  
  /*  int sizefft=howmany*idist*istride;
  CkPrintf("plan n %d, howmany %d istride %d idist %d ostride %d odist %d sizefft %d\n",plan->n,howmany, istride, idist, ostride, odist, sizefft);
  // for debugging 
  fftw_complex *scratch1=(fftw_complex *) fftw_malloc(sizefft*sizeof(complex));
  memcpy(scratch1,in,sizeof(complex)*sizefft);
  fftw(plan,         // da plan
       howmany,     // how many 
       scratch1,           //input data
       istride,      // stride betwen elements 
       idist,        // array separation
       out,          // output data (null if inplace)
       ostride,      // stride betwen elements  (0 inplace)
       odist);       //array separation (0 inplace)


  fftw_complex *scratch2=(fftw_complex *) fftw_malloc(sizefft*sizeof(complex));
  memcpy(scratch2,in,sizeof(complex)*sizefft);
  */
  //    fftw_complex *myin=scratch2;
  fftw_complex *myin=in;
  int split=config.fftprogresssplit; 
  int thismany=split;
  int inleft=howmany;
  int numsplits=howmany/split;
  if(numsplits<=0)
    {numsplits=1;}
  if(howmany>split && howmany%split!=0)
    {numsplits++;}
  
  int inoff=0;
  int outoff=0;
  //  CkPrintf("numsplits %d\n",numsplits);
  for(int i=0;i<numsplits;i++)
    {
      thismany=split;
      if(inleft<split)
	{thismany=inleft;}
      //      CkPrintf("split %d thismany %d inleft %d\n",i, thismany,inleft);
      fftw(plan,         // da plan
	   thismany,     // how many 
	   &(myin[inoff]),           //input data
	   istride,      // stride betwen elements 
	   idist,        // array separation
	   (out==NULL) ? NULL : &(out[outoff]),          // output data (null if inplace)
	   ostride,      // stride betwen elements  (0 inplace)
	   odist);       //array separation (0 inplace)

      inoff+=idist*(thismany);  
      outoff+=odist*(thismany);
      inleft-=thismany;      
#ifdef CMK_VERSION_BLUEGENE
      CmiNetworkProgress();
#endif
    }
  /*  for(int j=0;j<sizefft;j++)
    {
      CkPrintf("%g %g\n",scratch1[j].re,scratch1[j].im);
    }
  CkPrintf("----\n");
  for(int j=0;j<sizefft;j++)
    {
      CkPrintf("%g %g\n",scratch2[j].re,scratch2[j].im);
    }
  CkAbort("fftd\n");
  */

}

//============================================================================
/**
 * split up an fft call into multiple invocations with cmiprogress
 * calls between them.
 */
//============================================================================
void rfftwnd_complex_to_real_split(rfftwnd_plan plan, int howmany, fftw_complex *in, int istride,
		 int idist, fftw_real *out, int ostride, int odist)
{
    int split=config.fftprogresssplit; 
  int thismany=split;
  int inleft=howmany;
  int numsplits=howmany/split;
  if(numsplits<=0)
    {numsplits=1;}
  if(howmany>split && howmany%split!=0)
    {numsplits++;}
  
  int inoff=0;
  int outoff=0;
  //  CkPrintf("numsplits %d\n",numsplits);
  for(int i=0;i<numsplits;i++)
    {
      thismany=split;
      if(inleft<split)
	{thismany=inleft;}
      //      CkPrintf("split %d thismany %d inleft %d\n",i, thismany,inleft);
      rfftwnd_complex_to_real(plan,         // da plan
	   thismany,     // how many 
	   &(in[inoff]),           //input data
	   istride,      // stride betwen elements 
	   idist,        // array separation
	   (out==NULL) ? NULL : &(out[outoff]),          // output data (null if inplace)
	   ostride,      // stride betwen elements  (0 inplace)
	   odist);       //array separation (0 inplace)

      inoff+=idist*(thismany);  
      outoff+=odist*(thismany);
      inleft-=thismany;      
#ifdef CMK_VERSION_BLUEGENE
      CmiNetworkProgress();
#endif
    }

}

//============================================================================
/**
 * split up an fft call into multiple invocations with cmiprogress
 * calls between them.
 */
//============================================================================
void  rfftwnd_real_to_complex_split(rfftwnd_plan plan, int howmany, fftw_real *in, int istride,
		 int idist, fftw_complex *out, int ostride, int odist)
{
  int split=config.fftprogresssplit; 
  int thismany=split;
  int inleft=howmany;
  int numsplits=howmany/split;
  if(numsplits<=0)
    {numsplits=1;}
  if(howmany>split && howmany%split!=0)
    {numsplits++;}
  
  int inoff=0;
  int outoff=0;
  //  CkPrintf("numsplits %d\n",numsplits);
  for(int i=0;i<numsplits;i++)
    {
      thismany=split;
      if(inleft<split)
	{thismany=inleft;}
      //      CkPrintf("split %d thismany %d inleft %d\n",i, thismany,inleft);
      rfftwnd_real_to_complex(plan,         // da plan
	   thismany,     // how many 
	   &(in[inoff]),           //input data
	   istride,      // stride betwen elements 
	   idist,        // array separation
	   (out==NULL) ? NULL : &(out[outoff]),          // output data (null if inplace)
	   ostride,      // stride betwen elements  (0 inplace)
	   odist);       //array separation (0 inplace)

      inoff+=idist*(thismany);  
      outoff+=odist*(thismany);
      inleft-=thismany;      
#ifdef CMK_VERSION_BLUEGENE
      CmiNetworkProgress();
#endif
    }

}

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
GStateSlab::~GStateSlab() {

    if(packedPlaneData    !=NULL) delete [] packedPlaneData;
    if(packedPlaneDataTemp!=NULL) delete [] packedPlaneDataTemp;
    if(packedForceData    !=NULL) delete [] packedForceData;
    if(packedPlaneDataScr !=NULL) delete [] packedPlaneDataScr;
    if(packedVelData      !=NULL) delete [] packedVelData;
    if(runs               !=NULL) delete [] runs;
    destroyNHC();

    packedPlaneData     = NULL;
    packedPlaneDataTemp = NULL;
    packedForceData     = NULL;
    packedPlaneDataScr  = NULL;
    packedVelData       = NULL;
    runs                = NULL;

}
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void GStateSlab::pup(PUP::er &p) {
//==============================================================================
// Dont have to pup fftw plans - they live in the fft cache group

        p|numNonZeroPlanes;
	p|numRuns;
	p|numLines;
	p|numPoints;
	p|numFull;
	p|mysizeX; 
	p|planeSize;
        p|fftReqd;
	p|S_grainSize;
	p|xdim;
	p|ydim;
	p|zdim;
	p|iplane_ind;
        p|istate_ind;
        p|ihave_kx0;
        p|kx0_strt;
        p|kx0_end;
        p|nkx0; p|nkx0_uni; p|nkx0_red; p|nkx0_zero;

        p|eke_ret;
        p|fictEke_ret;
	if (p.isUnpacking()) {runs = new RunDescriptor[numRuns];}
	for (int i = 0; i < numRuns; i++){runs[i].pup(p);}

	if (p.isUnpacking()) {
           packedPlaneData     = new complex[numPoints];
	   packedPlaneDataScr  = new complex[numPoints];
	   packedPlaneDataTemp = new complex[numPoints];
	   packedForceData     = new complex[numPoints];
	   packedVelData       = new complex[numPoints];
           packedRedPsi        = new complex[nkx0];
	}//endif
	p((char *) packedPlaneData, numPoints*sizeof(complex));
	p((char *) packedPlaneDataScr, numPoints*sizeof(complex));
	p((char *) packedPlaneDataTemp, numPoints*sizeof(complex));
	p((char *) packedForceData, numPoints*sizeof(complex));
	p((char *) packedVelData, numPoints*sizeof(complex));
	p((char *) packedRedPsi, nkx0*sizeof(complex));

        p|len_nhc_cp;
        p|num_nhc_cp;
        p|kTCP;
        p|tauNHCCP;
        p|xNHC;
        p|mNHC;
        int nsize  = num_nhc_cp*len_nhc_cp;
        double *vt = new double[nsize];
        double *ft = new double[nsize];
	if (!p.isUnpacking()) {
          int iii=0;
          for(int i =0;i<num_nhc_cp;i++){
          for(int j =0;j<len_nhc_cp;j++){
            vt[iii] = vNHC[i][j];
            ft[iii] = fNHC[i][j];
            iii++;
          }}
	}//endif sending
        p(vt,nsize);
        p(ft,nsize);
	if (p.isUnpacking()) {
          initNHC();
          int iii=0;
          for(int i =0;i<num_nhc_cp;i++){
          for(int j =0;j<len_nhc_cp;j++){
            vNHC[i][j] = vt[iii];
            fNHC[i][j] = ft[iii];
            iii++;
          }}
	}//endif receiving 
        delete []vt;
        delete []ft;

//==============================================================================
  }//end routine
//==============================================================================



//=============================================================================
// packed g-space of size numPoints is expanded to numFull =numRuns/2*nfftz
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================

void FFTcache::expandGSpace(complex* data, complex *packedPlaneData, 
                            RunDescriptor *runs, int numRuns, int numFull,
                            int numPoints, int nfftz)

//==============================================================================
  {//begin routine
//==============================================================================
// Expand out lines of z so that an FFT can be performed on them
//   The ``run'' stores each line in two parts
//      0,1,2,3 have offset, joff = r*nfftz 
//     -3,-2,-1 have offset, joff = r*nffz + nfftz-3
//      runs[r].z stores kz if kz>0 and nfftz-kz if kz<0
//   Total size is nlines*nfftz where nlines=numRuns/2

  complex *origdata = packedPlaneData;
  memset(data, 0, sizeof(complex)*numFull);

  int koff = 0;
  for (int r = 0,l=0; r < numRuns; r+=2,l++) {

    int joff = l*nfftz + runs[r].z;
    for (int i=0,j=joff,k=koff; i<runs[r].length; i++,j++,k++) {
      data[j] = origdata[k];
    }//endfor
    koff += runs[r].length;

    int r1=r+1;
    joff = l*nfftz + runs[r1].z;
    for (int i=0,j=joff,k=koff; i<runs[r1].length; i++,j++,k++) {
      data[j] = origdata[k];
    }//endfor
    koff += runs[r1].length;

  }//endfor

  CkAssert(numPoints == koff);

//------------------------------------------------------------------------------
  }//end routine
//==============================================================================


//==============================================================================
// Expanded gspace of size numFull =numRuns/2*nfftz to packed size numPoints
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void GStateSlab::compressGSpace(const complex *expnddata, int type) {
//==============================================================================
// Contract lines of z after the backFFT has been performed.
//   The ``run'' stores each line in two parts
//      0,1,2,3 have offset, joff = r*nfftz 
//     -3,-2,-1 have offset, joff = r*nffz + nfftz-3
//      runs[r].z stores kz if kz>0 and nfftz-kz if kz<0

  complex *origdata;
  if(type == 1){origdata = packedPlaneData;}
  if(type == 2){origdata = packedForceData;}
  if(type!=2){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Unpacking invalid data, aborting %d\n",type);
    CkPrintf("Only the forces should be contracted\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  int nfftz = planeSize[1];

  int koff = 0;
  for (int r = 0,l=0; r < numRuns; r+=2,l++) {

    int joff = l*nfftz + runs[r].z;
    for (int i=0,j=joff,k=koff; i<runs[r].length; i++,j++,k++) {
      origdata[k]= expnddata[j];
    }//endfor
    koff += runs[r].length;

    int r1=r+1;
    joff = l*nfftz + runs[r1].z;
    for (int i=0,j=joff,k=koff; i<runs[r1].length; i++,j++,k++) {
      origdata[k] = expnddata[j];
    }//endfor
    koff += runs[r1].length;

  }//endfor

  CkAssert(numPoints == koff);

//==============================================================================
  }//end routine
//==============================================================================



//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void GStateSlab::addForces(const complex *points, const int *k_x){
//==============================================================================
    int i;
    double wght;

    if(!config.doublePack){

      for(i = 0; i < numPoints; i++){
	 packedForceData[i] += points[i];
      }//endfor

    }else{

      for(i = 0; i < numPoints; i++){
        wght = (k_x[i] == 0 ? 1.0 : 2.0);
#ifdef _CP_DEBUG_NON_LOCAL_ONLY_
	packedForceData[i] = 0.0; 
#endif
#ifdef _CP_DEBUG_VKS_ONLY_
	packedForceData[i] = 0.0; 
#endif
	packedForceData[i] *= wght; 
	packedForceData[i] += points[i];
      }//endfor

    }//endif : doublePack

//==============================================================================
  }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/*
 * In G-space our representation uses wrapping of coordinates.
 * When reading from the files, if x < 0, x += sizeX;
 */
//==============================================================================

void GStateSlab::setKVectors(int *n, int **kk_x, int **kk_y, int **kk_z){

//======================================================================
// Construct the k-vectors

  int *k_x = new int[numPoints];
  int *k_y = new int[numPoints];
  int *k_z = new int[numPoints];
  
  int r, i, dataCovered = 0;
  int x, y, z;
  for (r = 0; r < numRuns; r++) { // 2*number of lines z
    x = runs[r].x;
    if (x > sizeX/2) x -= sizeX;
    y = runs[r].y;
    if (y > planeSize[0]/2) y -= planeSize[0];
    z = runs[r].z;
    if (z > planeSize[1]/2) z -= planeSize[1];

    for (i = 0; i < runs[r].length; i++) { //pts in lines of z
      k_x[dataCovered] = x;
      k_y[dataCovered] = y;
      k_z[dataCovered] = (z+i);
      dataCovered++;
    }//endfor
  }//endfor

  CkAssert(dataCovered == numPoints);

//======================================================================
// Find pts with k_x==0 then check the layout : kx=0 first

  ihave_kx0 = 0;
  nkx0  = 0;
  nkx0_uni=0;
  nkx0_red=0;
  nkx0_zero=0;
  kx0_strt =0;
  for(i=0;i<numPoints;i++){
    if(k_x[i]==0 && k_y[i]>0){nkx0_uni++;}
    if(k_x[i]==0 && k_y[i]<0){nkx0_red++;}
    if(k_x[i]==0 && k_y[i]==0 && k_z[i]>=0){nkx0_uni++;}
    if(k_x[i]==0 && k_y[i]==0 && k_z[i]<0){nkx0_red++;}
    if(k_x[i]==0 && k_y[i]==0 && k_z[i]==0){nkx0_zero++;}
    if(k_x[i]==0){
      if(ihave_kx0==0){kx0_strt=i;}
      ihave_kx0=1;
      nkx0++;
    }//endif
  }//endif
  kx0_end = kx0_strt + nkx0;

  if(kx0_strt!=0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("kx=0 should be stored first | kx_srt !=0\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif
  
  if(nkx0!=nkx0_uni+nkx0_red){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Incorrect count of redundant guys\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  for(i=0;i<nkx0;i++){  
    if(k_x[i]!=0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("kx should be stored consecutively and first\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
  }//endif

  for(i=0;i<nkx0_red;i++){  
    if(k_y[i]>0 || (k_y[i]==0 && k_z[i]>=0)){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("ky <0 should be stored first\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
  }//endfor

  for(i=nkx0_red;i<nkx0_uni;i++){  
    if(k_y[i]<0 || (k_y[i]==0 && k_z[i]<0)){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("ky <0 should be stored first\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
  }//endfor

  packedRedPsi  = new complex[nkx0];
  memset(packedRedPsi, 0, sizeof(complex)*nkx0);

//==============================================================================
// Set the return values

  *n    = numPoints;
  *kk_x = k_x;
  *kk_y = k_y;
  *kk_z = k_z;
  
//==============================================================================
  }//end routine
//==============================================================================



//==============================================================================
// Expand the data set to line form : numLines*nfftz
// Perform the fft along z direction : numLines ffts 
// to get started towards realspace
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
complex* FFTcache::doGSRealFwFFT(complex *packedPlaneData, RunDescriptor *runs, 
                          int numRuns, int numLines,int numFull, int numPoints,
                          int nfftz, bool fftReqd)
//==============================================================================
   {//begin routine
//==============================================================================

    FFTcache *sc   = fftCacheProxy.ckLocalBranch();
    fftw_plan plan = sc->fwdYPlan; // for historic reasons its the `y' plan.
                                   // but sorry, its the z plan.
    expandGSpace(fftData,packedPlaneData,runs,numRuns,numFull,numPoints,nfftz);

    if(fftReqd){
        fft_split(plan,        // direction Z now
	     numLines,    // these many ffts : one for every line of z in the chare
	     (fftw_complex *)fftData, //input data
	     1,           //stride
	     nfftz,       //distance between z-data sets
	     NULL, 0, 0); // junk because input array stores the output (in-place)
    }//endif
    return fftData;

//==============================================================================
   }//end routine
//==============================================================================



//==============================================================================
// Perform the back fft along z direction : numLines ffts
// then compress the data to non-zero points as you are back in g-space
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void GStateSlab::doBwFFT(complex *fftData) {

   FFTcache *sc   = fftCacheProxy.ckLocalBranch();
   fftw_plan plan = sc->bwdYPlan;  // the Y is a misnomer : its Z
   int expandtype = 2;
   int nfftz      = planeSize[1];

   if(fftReqd){
        fft_split(plan,        // direction Z
	     numLines,    // these many ffts : lines of z
	     (fftw_complex *)fftData, //input data
	     1,           //stride
	     nfftz,       //distance between arrays
	     NULL, 0, 0); // junk because input array stores output (in-place)
   }//endfor
   compressGSpace(fftData, expandtype);

//==============================================================================
  }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/* This gets called at the end of the GStatePlane constructor */
//==============================================================================
void initGStateSlab(GStateSlab *gs, int sizeX, size2d size, int gSpaceUnits, 
          int realSpaceUnits, int s_grain, int iplane_ind,int istate_ind) 
//==============================================================================
   {//begin routine
//==============================================================================
// Explanation of organization of data: A point is psi(kx,ky,kz)  
//
//   Each gStateSlab is a collection of pts in lines of constant kx,ky.
//   The maximum number of pts in any line is size[1] = nfftz
//   Some lines are longer/shorter to spherical trunction (|k|<k_cut).
//   Also, kx>=0 when at the Gamma point or doublePack=1, the only
//   thing that has been recently tested. Collections of lines are parallelized. 
//   Each collection does not correspond to a unique plane of kx.
//   The parameter gSpaceUnits is obsolete and must be unity.
//   In order to create the state in real space, psi(x,y,z), an fft 
//   of all kz lines is performed followed by a transpose. The result of the transpose
//   is parallelized by planes of z. This is different than it used to be.
//   The number of points AFTER the fft is nlines*nfftz will be different
//   for each chare array

   if(gSpaceUnits!=1 || realSpaceUnits!=1){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("gspacePPC==1 obsolete. We no longer parallelize g-space by plane\n");
     CkPrintf("realSpacePPC==1 although real space is parallelized by plane\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
   }//endif
   gs->numNonZeroPlanes=1;
   gs->mysizeX   = sizeX;
   gs->planeSize = size;  // fftsizes (sizeY,sizeZ)  sizeX is a global
   gs->numPoints = 0;     // number of packed points (data size before fft)
   gs->numRuns   = 0;     // 2*(number of lines) in the collection
   gs->numLines  = 0;     // (number of lines) in the collection
   gs->fftReqd   = false; // false if this chare has no packed pts
   gs->numFull   = 0;     // number of pts : numLines*nfftz

   gs->S_grainSize = s_grain;   // PC grainsize
   gs->xdim = 1;                // may need some love
   gs->ydim = gs->planeSize[0];
   gs->zdim = gs->planeSize[1];

   gs->iplane_ind   = iplane_ind;
   gs->istate_ind   = istate_ind; 
   gs->initNHC();

//==============================================================================
   }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
RealStateSlab::~RealStateSlab() {
//==============================================================================

  if(planeArr != NULL) {fftw_free(planeArr);planeArr = NULL; }

}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
//FFTcache - sets up a fftcache on each processor
//==============================================================================
FFTcache::FFTcache(size2d planeSIZE, int ArraySize){
//==============================================================================
	
        planeSize = planeSIZE;  //contains nfftz,nfftx

        int nlines_max = scProxy.ckLocalBranch()->cpcharmParaInfo->nlines_max;
        int psize      = nlines_max*planeSize[0];

	fftData = (complex*) fftw_malloc(psize*sizeof(complex)); 

	//Non Double Pack plans : 2 Y planes and 2 X plans.
	fwdZ1DPlan = fftw_create_plan(planeSize[0], FFTW_FORWARD,  
                     FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
	bwdZ1DPlan = fftw_create_plan(planeSize[0], FFTW_BACKWARD, 
                     FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
	fwdX1DPlan = fftw_create_plan(sizeX, FFTW_FORWARD, 
                     FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
	bwdX1DPlan = fftw_create_plan(sizeX, FFTW_BACKWARD, 
                     FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
	int size[3];
	size[0] = planeSize[1]; size[1] = planeSize[0]; size[2] = 1;
	//Double Pack Plans
        fwdZ1DdpPlan = fftw_create_plan(planeSize[0], FFTW_FORWARD, 
                          FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
	
	bwdZ1DdpPlan = fftw_create_plan(planeSize[0], FFTW_BACKWARD, 
                       FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
	fwdX1DdpPlan = rfftwnd_create_plan(1, (const int*)size, FFTW_COMPLEX_TO_REAL, 
                       FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
        bwdX1DdpPlan = rfftwnd_create_plan(1, (const int*)size,FFTW_REAL_TO_COMPLEX,
                        FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
        int sizeZ = planeSize[0];
	fwdYPlan = fftw_create_plan(sizeZ,FFTW_FORWARD, 
                                 FFTW_MEASURE|FFTW_IN_PLACE|FFTW_USE_WISDOM);
	bwdYPlan = fftw_create_plan(sizeZ,FFTW_BACKWARD, 
                                 FFTW_MEASURE|FFTW_IN_PLACE|FFTW_USE_WISDOM);

//==============================================================================
   }//end routine
//==============================================================================

//==============================================================================
// After the transpose : nfftz*nffty array nicely expanded out.
//                       Only nplane_x FFTs along Y need to be performed
//                       due to spherical truncation
//                       data is |psi|^2 : z is chare array index
//                                         x is inner index 
//                                         y is outer index
//                       chare[z].data[x+y*sizeX]
//                       chare[z].planeArr[y+x*sizeY] is input
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================

double* FFTcache::doRealFwFFT(complex *planeArr)

//==============================================================================
   {//begin routine    
//==============================================================================

   double *data;
   int pSize    = planeSize[0] * sizeX;
   int nplane_x = scProxy.ckLocalBranch()->cpcharmParaInfo->nplane_x;

//==============================================================================
// Case doublePack && config.inPlaceFFT
// 

  if(config.doublePack && config.inPlaceFFT) {
      data = new double[pSize];
      int stride = sizeX/2+1;
      // FFT along Y direction : Y moves with stride sizex/2+1 through memory
      fft_split(fwdZ1DdpPlan,    // y-plan (label lies)
  	   nplane_x,        // how many < sizeX/2 + 1
	   (fftw_complex *)(planeArr),//input data
 	   stride,          // stride betwen elements (y is inner)
  	   1,               // array separation (nffty elements)
	   NULL,0,0);       // output data is input data
      // fftw only gives you one sign for real to complex : so do it yourself
      for(int i=0;i<stride*planeSize[0];i++){planeArr[i].im = -planeArr[i].im;}
      rfftwnd_complex_to_real_split(fwdX1DdpPlan,
			      planeSize[0],    // how many
			      (fftw_complex *)planeArr, 
                              1,               // stride (x is inner)
                              stride,          // array separation
    		              NULL,0,0);       // output array is real 
      // x is now the inner index as fftw has transposed for us
      double *realArr = reinterpret_cast<double*> (planeArr);
      for(int i=0,i2=0;i<planeSize[0];i++,i2+=2){
       for(int j=i*sizeX;j<(i+1)*sizeX;j++){
         data[j] = realArr[(j+i2)]*realArr[(j+i2)];
       }//endfor
      }//endfor
    return data;

  }//endif : 

//==============================================================================
// Case !doublePack && config.inPlaceFFT

 if(!config.doublePack && config.inPlaceFFT){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Fix the non-double pack FFT\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
 }//endif

 return data;                

//==============================================================================
  }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RealStateSlab::zeroOutPlanes() {    
  memset(planeArr,0,size*sizeof(complex));
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RealStateSlab::allocate() {
  if(planeArr == NULL) {planeArr = (complex *) fftw_malloc(size*sizeof(complex));}
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RealStateSlab::destroy() {
  if(planeArr != NULL) {fftw_free(planeArr); planeArr=NULL;}
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doRhoRealtoRhoG(double *realArr)
//==============================================================================
 {//begin routine 
   //==============================================================================

   int nplane_rho_x = scProxy.ckLocalBranch()->cpcharmParaInfo->nplane_rho_x;

   //==============================================================================
   // Case : doublePack  and inplace no options thats how it is
   
   rfftwnd_real_to_complex_split(bwdX1DdpPlan,
			   planeSize[0],          // these many 1D ffts
			   realArr,1,(sizeX+2),   // x is inner here
			   NULL,0,0);            
   complex *planeArr = reinterpret_cast<complex*> (realArr);
   int stride = sizeX/2+1;

#ifdef CMK_VERSION_BLUEGENE
   CmiNetworkProgress();
#endif

   for (int i=0; i<stride*planeSize[0];i++){planeArr[i].im = -planeArr[i].im;}

   fft_split(bwdZ1DdpPlan,
	nplane_rho_x, // these many 1D ffts
	(fftw_complex *)planeArr, 
	stride,1,NULL,0,0);

   //==============================================================================
 }//end routine
//==============================================================================


//==============================================================================
//                             vks : z is chare array index
//                                   x is inner index 
//                                   y is outer index
//                                   chare[z].data[x+y*sizeX]
//
//          vks is from a message and should not be modified
//
//                           planeArr : chare[z].data[y+x*sizeY]
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doRealBwFFT(const double *vks, complex *planeArr,
                           int ind_state, int ind_plane)
//==============================================================================
 {//begin routine 
//==============================================================================

   int pSize    = planeSize[0] * sizeX;
   int nplane_x = scProxy.ckLocalBranch()->cpcharmParaInfo->nplane_x;

//==============================================================================
// Output

#ifdef _CP_DEBUG_FFTR_VKSR_
   if(config.doublePack){
     if(ind_state==0 && ind_plane==0){
        double *realArr = reinterpret_cast<double*> (planeArr);
        FILE *fp = fopen("vks_state0_plane0_real.out","w");
          for(int i=0,i2=0;i<planeSize[0];i++,i2+=2){
            for(int j=i*sizeX;j<(i+1)*sizeX;j++){
              fprintf(fp,"%g %g\n",vks[j],realArr[(j+i2)]);
            }//endfor
          }//endfor
        fclose(fp);
     }//endif
   }//endif
#endif

//==============================================================================
// Case : doublePack  and inplace 

  if(config.doublePack && config.inPlaceFFT){
     double *realArr = reinterpret_cast<double*> (planeArr);
     int stride = sizeX/2+1;
     for(int i=0,i2=0;i<planeSize[0];i++,i2+=2){
       for(int j=i*sizeX;j<(i+1)*sizeX;j++){
         realArr[(j+i2)] = realArr[(j+i2)]*vks[j];
       }//endfor
     }//endfor
     rfftwnd_real_to_complex_split(bwdX1DdpPlan,
			planeSize[0],          // these many 1D ffts
			realArr,1,(sizeX+2),   // x is inner here
                        NULL,0,0);            
#ifdef CMK_VERSION_BLUEGENE
   CmiNetworkProgress();
#endif

     for (int i=0; i<stride*planeSize[0];i++){planeArr[i].im = -planeArr[i].im;}
     fft_split(bwdZ1DdpPlan,
	     nplane_x, // these many 1D ffts
	     (fftw_complex *)planeArr, 
  	     stride,1,NULL,0,0);
  }//endif : double pack and in-place

//==============================================================================

  if(!config.doublePack && config.inPlaceFFT){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Fix the non-double pack FFT\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

//==============================================================================
   }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/* This gets called at the end of the CP_State_RealSpacePlane constructor */
//==============================================================================
void initRealStateSlab(RealStateSlab *rs, size2d planeSize, int gSpaceUnits, 
                       int realSpaceUnits, int stateIndex, int planeIndex)
//==============================================================================
   {//begin routine
//==============================================================================
// Explanation of the organization:
//   Each CP_State_RealSpacePlane instance is actually a planes of psi_I(x,y,z).
//   We index into planeArr using the z-coordinate
//   and each chare is an x-y plane.
//   Data arrives from the gspaceplane of total size : nlines_tot*nfftz

   if(gSpaceUnits!=1 || realSpaceUnits!=1){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("gspacePPC==1 obsolete. We no longer parallelize g-space by plane\n");
     CkPrintf("realSpacePPC==1 although real space is parallelized by plane\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
   }//endif

   rs->planeSize  = planeSize;         // FFTSize Y and Z, sizeX is global
   rs->thisState  = stateIndex;        // my state (I)
   rs->thisPlane  = planeIndex;        // my plane (z)
   rs->nsize      = planeSize[0]*sizeX; // when fft is completed
   rs->numPlanesToExpect = scProxy.ckLocalBranch()->cpcharmParaInfo->nchareG;
   int rsize_a    = planeSize[0]*(planeSize[1]/2+1);
   int rsize_b    = planeSize[1]*(planeSize[0]/2+1);
   rs->rsize      = (rsize_a > rsize_b ? rsize_a : rsize_b);

   if(config.doublePack)  {
      rs->size = rs->rsize;
   }else{
      rs->size = rs->nsize;
   }//endif

   rs->planeArr = (complex *) fftw_malloc(rs->size * sizeof(complex));

//==============================================================================    
   }//end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
RhoRealSlab::~RhoRealSlab(){

	//delete [] Vks;  I wonder why this isn't deleted
	delete [] density;
        complex *dummy;
	dummy  = reinterpret_cast<complex*> (doFFTonThis);
	fftw_free((fftw_complex *)dummy);
	dummy  = reinterpret_cast<complex*> (rhoIRX);
	fftw_free((fftw_complex *)dummy);
	dummy  = reinterpret_cast<complex*> (rhoIRY);
	fftw_free((fftw_complex *)dummy);
	dummy  = reinterpret_cast<complex*> (rhoIRZ);
	fftw_free((fftw_complex *)dummy);
	dummy  = reinterpret_cast<complex*> (gradientCorrection);
	fftw_free((fftw_complex *)dummy);
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/* This gets called at the end of the RealSpaceDensity Constructor */
//==============================================================================
void initRhoRealSlab(RhoRealSlab *rho_rs, int xdim, int ydim, int zdim,
                     int numRealSpace, int numRhoG, int myIndexX, int myIndexY)
//==============================================================================
   {//begin routine
//==============================================================================

	rho_rs->sizeX = xdim;
	rho_rs->sizeY = ydim;
	rho_rs->sizeZ = zdim;

	rho_rs->xdim = rho_rs->sizeX;
	rho_rs->ydim = 1;
	rho_rs->zdim = rho_rs->sizeZ;

	rho_rs->startx = 0;
	rho_rs->starty = 0;
	rho_rs->startz = 0; 

	rho_rs->size     = (rho_rs->sizeX+2) * (rho_rs->sizeZ);
	rho_rs->trueSize = (rho_rs->sizeX) * (rho_rs->sizeZ);
        int sizenow      = rho_rs->size;

	rho_rs->Vks     =  new double[sizenow];
	bzero(rho_rs->Vks, sizenow*sizeof(double));
	rho_rs->density =  new double[sizenow];
	bzero(rho_rs->density, sizenow*sizeof(double));
        int csizenow    = sizenow/2;

        complex *dummy;
	dummy               = (complex*) fftw_malloc(csizenow*sizeof(complex));
	rho_rs->doFFTonThis = reinterpret_cast<double*> (dummy);
	dummy               = (complex*) fftw_malloc(csizenow*sizeof(complex));
	rho_rs->rhoIRX      = reinterpret_cast<double*> (dummy);
	dummy               = (complex*) fftw_malloc(csizenow*sizeof(complex));
	rho_rs->rhoIRY      = reinterpret_cast<double*> (dummy);
	dummy               = (complex*) fftw_malloc(csizenow*sizeof(complex));
	rho_rs->rhoIRZ      = reinterpret_cast<double*> (dummy);
	dummy               = (complex*) fftw_malloc(csizenow*sizeof(complex));
	rho_rs->gradientCorrection = reinterpret_cast<double*> (dummy);

//==============================================================================
   }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
RhoGSlab::~RhoGSlab()
{
  fftw_free(Rho);
  fftw_free( divRhoX);
  fftw_free( divRhoY);
  fftw_free( divRhoX);
  fftw_free( packedRho);
  fftw_free( packedVks);
  fftw_free( Vks);
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/*
 * In G-space our representation uses wrapping of coordinates.
 * When reading from the files, if x < 0, x += sizeX;
 */
//==============================================================================

void RhoGSlab::setKVectors(int *n){

//======================================================================
// Construct the k-vectors

  k_x = new int[numPoints];
  k_y = new int[numPoints];
  k_z = new int[numPoints];
  
  int r, i, dataCovered = 0;
  int x, y, z;
  for (r = 0; r < numRuns; r++) { // 2*number of lines z
    x = runs[r].x;
    if (x > sizeX/2) x -= sizeX;
    y = runs[r].y;
    if (y > sizeY/2) y -= sizeY;
    z = runs[r].z;
    if (z > sizeZ/2) z -= sizeZ;

    for (i = 0; i < runs[r].length; i++) { //pts in lines of z
      k_x[dataCovered] = x;
      k_y[dataCovered] = y;
      k_z[dataCovered] = (z+i);
      dataCovered++;
    }//endfor
  }//endfor

  CkAssert(dataCovered == numPoints);

//==============================================================================
// Set the return values

  *n    = numPoints;
  
//==============================================================================
  }//end routine
//==============================================================================

//=============================================================================
// packed g-space of size numPoints is expanded to numFull =numRuns/2*nfftz
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================

void RhoGSlab::expandRhoGSpace(complex* data, complex *packedData)

//==============================================================================
  {//begin routine
//==============================================================================
// Expand out lines of z so that an FFT can be performed on them
//   The ``run'' stores each line in two parts
//      0,1,2,3 have offset, joff = r*nfftz 
//     -3,-2,-1 have offset, joff = r*nffz + nfftz-3
//      runs[r].z stores kz if kz>0 and nfftz-kz if kz<0
//   Total size is nlines*nfftz where nlines=numRuns/2

  memset(data,0,sizeof(complex)*numFull);
  int nfftz = sizeZ;

  int koff = 0;
  for (int r = 0,l=0; r < numRuns; r+=2,l++) {

    int joff = l*nfftz + runs[r].z;
    for (int i=0,j=joff,k=koff; i<runs[r].length; i++,j++,k++) {
      data[j] = packedData[k];
    }//endfor
    koff += runs[r].length;

    int r1=r+1;
    joff = l*nfftz + runs[r1].z;
    for (int i=0,j=joff,k=koff; i<runs[r1].length; i++,j++,k++) {
      data[j] = packedData[k];
    }//endfor
    koff += runs[r1].length;

  }//endfor

  CkAssert(numPoints == koff);

//------------------------------------------------------------------------------
  }//end routine
//==============================================================================

//==============================================================================
// packed g-space of size numPoints is expanded to numFull =numRuns/2*nfftz
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================

void RhoGSlab::divRhoGdot(double *hmati, double tpi){

//==============================================================================

  int nfftz = sizeZ;
  double gx,gy,gz;

//==============================================================================

  memset(divRhoX, 0, sizeof(complex)*numFull);
  memset(divRhoY, 0, sizeof(complex)*numFull);
  memset(divRhoZ, 0, sizeof(complex)*numFull);

  int koff = 0;
  for (int r = 0,l=0; r < numRuns; r+=2,l++) {

    int joff = l*nfftz + runs[r].z;
    for (int i=0,j=joff,k=koff; i<runs[r].length; i++,j++,k++) {
      gx = tpi*(k_x[k]*hmati[1] + k_y[k]*hmati[2] + k_z[k]*hmati[3]);
      gy = tpi*(k_x[k]*hmati[4] + k_y[k]*hmati[5] + k_z[k]*hmati[6]);
      gz = tpi*(k_x[k]*hmati[7] + k_y[k]*hmati[8] + k_z[k]*hmati[9]);
      complex tmp = (packedRho[k].multiplyByi())*(-1.0);
      divRhoX[j] = tmp*gx;
      divRhoY[j] = tmp*gy;
      divRhoZ[j] = tmp*gz;
    }//endfor
    koff += runs[r].length;

    int r1=r+1;
    joff = l*nfftz + runs[r1].z;
    for (int i=0,j=joff,k=koff; i<runs[r1].length; i++,j++,k++) {
      gx = tpi*(k_x[k]*hmati[1] + k_y[k]*hmati[2] + k_z[k]*hmati[3]);
      gy = tpi*(k_x[k]*hmati[4] + k_y[k]*hmati[5] + k_z[k]*hmati[6]);
      gz = tpi*(k_x[k]*hmati[7] + k_y[k]*hmati[8] + k_z[k]*hmati[9]);
      complex tmp = (packedRho[k].multiplyByi())*(-1.0);
      divRhoX[j] = tmp*gx;
      divRhoY[j] = tmp*gy;
      divRhoZ[j] = tmp*gz;
    }//endfor
    koff += runs[r1].length;

  }//endfor

  CkAssert(numPoints == koff);

//------------------------------------------------------------------------------
  }//end routine
//==============================================================================


//==============================================================================
// packed g-space of size numPoints is expanded to numFull =numRuns/2*nfftz
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================

void RhoGSlab::createWhiteByrd(double *hmati, double tpi){

//==============================================================================

  int nfftz = sizeZ;
  double gx,gy,gz;
  complex tmp;

//==============================================================================

  complex *whitebyrd = Rho;
  memset(whitebyrd,0,sizeof(complex)*numFull);

  int koff = 0;
  for (int r = 0,l=0; r < numRuns; r+=2,l++) {

    int joff = l*nfftz + runs[r].z;
    for (int i=0,j=joff,k=koff; i<runs[r].length; i++,j++,k++) {
      gx  = tpi*(k_x[k]*hmati[1] + k_y[k]*hmati[2] + k_z[k]*hmati[3]);
      gy  = tpi*(k_x[k]*hmati[4] + k_y[k]*hmati[5] + k_z[k]*hmati[6]);
      gz  = tpi*(k_x[k]*hmati[7] + k_y[k]*hmati[8] + k_z[k]*hmati[9]);
      tmp = divRhoX[j]*gx + divRhoY[j]*gy + divRhoZ[j]*gz;
      whitebyrd[j] = tmp.multiplyByi()*(-1.0); 
    }//endfor
    koff += runs[r].length;

    int r1=r+1;
    joff = l*nfftz + runs[r1].z;
    for (int i=0,j=joff,k=koff; i<runs[r1].length; i++,j++,k++) {
      gx  = tpi*(k_x[k]*hmati[1] + k_y[k]*hmati[2] + k_z[k]*hmati[3]);
      gy  = tpi*(k_x[k]*hmati[4] + k_y[k]*hmati[5] + k_z[k]*hmati[6]);
      gz  = tpi*(k_x[k]*hmati[7] + k_y[k]*hmati[8] + k_z[k]*hmati[9]);
      tmp = divRhoX[j]*gx + divRhoY[j]*gy + divRhoZ[j]*gz;
      whitebyrd[j] = tmp.multiplyByi()*(-1.0); 
    }//endfor
    koff += runs[r1].length;

  }//endfor

  CkAssert(numPoints == koff);

//------------------------------------------------------------------------------
  }//end routine
//==============================================================================


//==============================================================================
//==============================================================================
// Expanded rhogspace of size numFull =numRuns/2*nfftz to packed size numPoints
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RhoGSlab::compressGSpace(const complex *expnddata, int type) {
//==============================================================================
// Contract lines of z after the backFFT has been performed.
//   The ``run'' stores each line in two parts
//      0,1,2,3 have offset, joff = r*nfftz 
//     -3,-2,-1 have offset, joff = r*nffz + nfftz-3
//      runs[r].z stores kz if kz>0 and nfftz-kz if kz<0

  complex *data;
  switch(type){
    case 0 : data = packedRho; break;
  }//endif

  int nfftz = sizeZ;

  int koff = 0;
  for (int r = 0,l=0; r < numRuns; r+=2,l++) {

    int joff = l*nfftz + runs[r].z;
    for (int i=0,j=joff,k=koff; i<runs[r].length; i++,j++,k++) {
      data[k]= expnddata[j];
    }//endfor
    koff += runs[r].length;

    int r1=r+1;
    joff = l*nfftz + runs[r1].z;
    for (int i=0,j=joff,k=koff; i<runs[r1].length; i++,j++,k++) {
      data[k] = expnddata[j];
    }//endfor
    koff += runs[r1].length;

  }//endfor

  CkAssert(numPoints == koff);

//==============================================================================
  }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// complete the transform of rho(r) to rho(g) : last FFT after tranpose
// 
//==============================================================================
void RhoGSlab::doBwFFTRtoG(int expandtype){
//==============================================================================

   FFTcache *sc   = fftCacheProxy.ckLocalBranch();
   fftw_plan plan = sc->bwdYPlan;  // Y label is a misnomer : its Z
   int nfftz      = sizeZ;

   complex *partlyfftData;
   switch(expandtype){
     case 0: partlyfftData = Rho; break;
     case 1: partlyfftData = divRhoX; break;
     case 2: partlyfftData = divRhoY; break;
     case 3: partlyfftData = divRhoZ; break;
   }//endif

   fft_split(plan,        // direction Z
        numLines,    // these many ffts : lines of z
	(fftw_complex *)partlyfftData, //input data
	1,           //stride
	nfftz,       //distance between arrays
	NULL, 0, 0); // junk because input array stores output (in-place)

   if(expandtype==0){compressGSpace(partlyfftData, expandtype);}

//==============================================================================
  }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RhoGSlab::doFwFFTGtoR(int iopt, int index){
//==============================================================================

   int nfftz      = sizeZ;
   FFTcache *sc   = fftCacheProxy.ckLocalBranch();
   fftw_plan plan = sc->fwdYPlan; // for historic reasons its called `y' plan
                                  // but it contains the `z' plan.

   complex *fftData;
   switch(iopt){
     case 0 : fftData = Rho;     break;
     case 1 : fftData = divRhoX; break;
     case 2 : fftData = divRhoY; break;
     case 3 : fftData = divRhoZ; break;
     case 4 : fftData = Vks;  break;
   default: CkAbort("impossible iopt"); break;
   }//end switch

   fft_split(plan,        // direction Z now
        numLines,    // these many ffts : one for every line of z in the chare
	reinterpret_cast<fftw_complex*>(fftData), //input data
         1,           //stride
         nfftz,       //distance between z-data sets
         NULL, 0, 0); // junk because input array stores the output (in-place)

//==============================================================================
   }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RhoRealSlab::doFwFFTGtoR(int iopt,double probScale){
//==============================================================================

   FFTcache *fsc             = fftCacheProxy.ckLocalBranch();
   fftw_plan    fwdZ1DdpPlan = fsc->fwdZ1DdpPlan; 
   rfftwnd_plan fwdX1DdpPlan = fsc->fwdX1DdpPlan; 

   CPcharmParaInfo *sim   = (scProxy.ckLocalBranch ())->cpcharmParaInfo;      
   int nplane_x           = sim->nplane_rho_x;
   int pSize              = (sizeX+2) * sizeY;
   int stride             = sizeX/2+1;

   double *data;
   switch(iopt){
     case 0: data = Vks; break;
     case 1: data = rhoIRX; break;
     case 2: data = rhoIRY; break;
     case 3: data = rhoIRZ; break;
     case 4: data = doFFTonThis; break;
   }//endif
  complex *planeArr = reinterpret_cast<complex*> (data);

//==============================================================================
// FFT along Y direction : Y moves with stride sizex/2+1 through memory

   fft_split(fwdZ1DdpPlan,    // y-plan (label lies)
	nplane_x,        // how many < sizeX/2 + 1
        (fftw_complex *)(planeArr),//input data
 	stride,          // stride betwen elements (x is inner)
  	1,               // array separation (nffty elements)
	NULL,0,0);       // output data is input data

#ifdef CMK_VERSION_BLUEGENE
   CmiNetworkProgress();
#endif

  // fftw only gives you one sign for real to complex : so do it yourself
  for(int i=0;i<stride*sizeY;i++){planeArr[i].im = -planeArr[i].im;}

//==============================================================================
// FFT along X direction : X moves with stride 1 through memory

  rfftwnd_complex_to_real_split(fwdX1DdpPlan,
		          sizeY,    // how many
			  (fftw_complex *)planeArr, 
                          1,               // stride (x is inner)
                          stride,          // array separation
  		          NULL,0,0);       // output array is real 

//==============================================================================
// copy out to a nice stride in memory

  double *temp = gradientCorrection;
  memcpy(temp,data,pSize*sizeof(double));

  for(int i=0,i2=0;i<sizeY;i++,i2+=2){
    for(int j=i*sizeX;j<(i+1)*sizeX;j++){
      data[j] = temp[(j+i2)]*probScale;
    }//endfor
  }//endfor

//==============================================================================
  }//end if
//==============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void RhoRealSlab::uPackAndScale(double *uPackData, double *PackData, 
                                double scale)
//============================================================================
   {//begin routine
//============================================================================

   for(int i=0,i2=0;i<sizeZ;i++,i2+=2){
     for(int j=i*sizeX;j<(i+1)*sizeX;j++){
       uPackData[(j+i2)] = PackData[j]*scale;
     }//endfor
     uPackData[sizeX]   = 0.0;
     uPackData[sizeX+1] = 0.0;
   }//endfor

//============================================================================
   }
//============================================================================


