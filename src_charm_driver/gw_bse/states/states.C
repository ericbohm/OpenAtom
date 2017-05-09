#include "states.decl.h"
#include "standard_include.h"
#include "states.h"
#include "messages.h"
#include "controller.h"
#include "pmatrix.h"
#include "main.decl.h"
#include "allclass_gwbse.h"
#include <assert.h>
#if CMK_PROJECTIONS_USE_ZLIB
#include "zlib.h"
#endif
#include "fft_routines.h"
#include "fft_controller.h"

States::States() {

  // Set my relevant indices and print out a debug message
  ispin = thisIndex.x, ikpt = thisIndex.y, istate = thisIndex.z;
  //CkPrintf("State chare with spin, kpt, state (%d %d %d) constructed on PE %d.\n", ispin, ikpt, istate, CkMyPe());

  // Get access to the global struct of readonlies
  GWBSE *gwbse = GWBSE::get();

  // set the option for the type of file read
  ibinary_opt = gwbse->gwbseopts.ibinary_opt;
  // set file name to read my state from
  sprintf(fileName, "./STATES_IN/Spin.%d_Kpt.%d_Bead.0_Temper.0/state%d.out", ispin, ikpt, istate+1);
  // if doublePack is true -> gamma point only calculation (only 1 k point)
  // if doublePack is false -> multiple k points
  // double pack is hard coded. FIXME : minjung needs to change this
  doublePack = false;
  // read states from file
  readState(fileName);

  // read states from shifted k grid (for q=0 calculation) 
  // set file name to read my state from  (there is additional "0" right after "Kpt." )
  // we do this only for the occupied states
  nocc = gwbse->gw_epsilon.nocc;
  int qindex = Q_IDX;
  if( istate < nocc && qindex == 0){
    sprintf(fileName,"./STATES_IN/Spin.%d_Kpt.0%d_Bead.0_Temper.0/state%d.out", ispin, ikpt, istate+1);
    readStateShifted(fileName);
  }
  
  // set the size of fft grid
  for(int i=0;i<3;i++){nfft[i] = gwbse->gw_parallel.fft_nelems[i];}

}

// constructor needed for chare object migration (ignore for now)
// note: this constructor does not need to appear in the ".ci" file
States::States(CkMigrateMessage *msg) { }

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// Entry methods for sending the state to either P or the PsiCache as needed.
// These methods will be called by the top level controller that is in charge
// of coordinating the number of stages and the pipeline depth for forming P.
//==============================================================================

//==============================================================================
// Pack up our realspace coefficients and broadcast them to the cache for longer
// term storage on each node.
void States::sendToCache() {
  //CkPrintf("[%i,%i,%i]: Sending psi to node cache...\n", ispin, ikpt, istate);
  int ndata = nfft[0]*nfft[1]*nfft[2];
  PsiMessage* msg = new (ndata, ndata) PsiMessage(ndata, stateCoeffR);
  msg->spin_index = ispin;
  msg->k_index = ikpt;
  msg->state_index = istate;
  msg->shifted = false;
  psi_cache_proxy.receivePsi(msg);

  msg = new (ndata, ndata) PsiMessage(ndata, stateCoeffR_shifted);//TODO:cleanup
  msg->spin_index = ispin;
  msg->k_index = ikpt;
  msg->state_index = istate;
  msg->shifted = true;
  psi_cache_proxy.receivePsi(msg);
}

//==============================================================================
// Pack up our realspace coefficients and broadcast them to the cache to be
// multiplied with each occupied psi to create a set of f vectors.
void States::sendToComputeF() {
  //CkPrintf("[%i,%i,%i]: Sending psi for f-comp...\n", ispin, ikpt, istate);
  int ndata = nfft[0]*nfft[1]*nfft[2];
  //int accept_arr[accept.size()];
  //for(int i=0;i<accept.size();i++)
  //  accept_arr[i] = accept[i];
  PsiMessage* msg = new (ndata, ndata) PsiMessage(ndata, stateCoeffR);
  msg->spin_index = ispin;
  msg->k_index = ikpt;
  msg->state_index = istate;
  msg->shifted = false;
  //msg->accept_size = accept.size();
  psi_cache_proxy.computeFs(msg);
 
}

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// FFTW Routines
//==============================================================================

void States::fftGtoR() {

  // Set up the FFT data structures in the FFTController
  FFTController* fft_controller = fft_controller_proxy.ckLocalBranch();
  int backward = 1;

  fft_controller->setup_fftw_3d(nfft,backward);

  fftw_complex* in_pointer = fft_controller->get_in_pointer();
  fftw_complex* out_pointer = fft_controller->get_out_pointer();
  
   // we need to setup fftidx
  int *g[3]; // put_into_fftbox routine takes 2D g array, so we need to do this
  g[0] = ga;
  g[1] = gb;
  g[2] = gc;

  int **fftidx;
  fftidx = new int *[numCoeff];
  for(int i=0; i<numCoeff;i++){ fftidx[i] = new int [3]; }

  // this routine changes negative g index to be a positive numbers
  // since it is origianlly written with Fortran, fftidx has fortran counting,
  // i.e., if gidx is (0,0,0), then (1,1,1) in fftidx
  gidx_to_fftidx(numCoeff, g, nfft, fftidx);

  // state coefficients are copied to in_pointer
  // put_into_fftbox was originally written for doublePack = 0 (false)
  // for gamma point calculation, put_into_fftbox has been modified from the original version
  put_into_fftbox(numCoeff, stateCoeff, fftidx, nfft, in_pointer, doublePack);
  
  // tell the FFTController to do the fft
  fft_controller->do_fftw();

  // transfer data from out_pointer to stateCoeffR
  // malloc stateRspace first
  int ndata = nfft[0]*nfft[1]*nfft[2];
  stateCoeffR = new complex [ndata];
  double scale = sqrt(1.0/double(ndata)); // IFFT requires normalization
  fftbox_to_array(ndata, out_pointer, stateCoeffR, scale);

  // delete stateCoeff
  delete [] stateCoeff;

  // fft for shifted states (only occupied states)
  int qindex = Q_IDX;
  if( istate < nocc && qindex == 0){
    stateCoeffR_shifted = new complex [ndata];
    put_into_fftbox(numCoeff, stateCoeff_shifted, fftidx, nfft, in_pointer, doublePack);
    fft_controller->do_fftw();
    fftbox_to_array(ndata, out_pointer, stateCoeffR_shifted, scale);
    delete [] stateCoeff_shifted;
  }

  // delete space used for fftidx
  for (int i = 0; i < numCoeff; i++) { delete [] fftidx[i]; }
  delete [] fftidx;

  // tell the controller that the states are ready
  contribute(CkCallback(CkReductionTarget(Controller, fftComplete), controller_proxy));
  
}


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// Scalar routines
//==============================================================================

void States::readState(char *fromFile) 

  //===================================================================================
{//begin routine
  //===================================================================================
  // A little screen output for the fans
    //CkPrintf("Reading state file: %s for chare (%d %d %d), with binary option %d.\n",fromFile,ispin,ikpt,istate,ibinary_opt);


  if(ibinary_opt < 0 || ibinary_opt > 3){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Bad binary option : %d\n",ibinary_opt);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  //===================================================================================
  // First read in the state and k-vectors : allows parsing of doublePack option

#if !CMK_PROJECTIONS_USE_ZLIB
  if(ibinary_opt>1)
  {
    CkPrintf("Attempt to use ZLIB Failed! Please review compilation\n");
    //CkPrintf("Macro cmk-projections-use-zlib  is %d \n", CMK_PROJECTIONS_USE_ZLIB);
    CkExit();
  }
#endif
  int nx,ny,nz;
  int nktot = 0;
  if(ibinary_opt==0){
    FILE *fp=fopen(fromFile,"r");
    if (fp==NULL){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Can't open state file %s\n",fromFile);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
    int numCoeffLoc;
    if(4!=fscanf(fp,"%d%d%d%d",&numCoeffLoc,&nx,&ny,&nz)){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Can't parse size line of file %s\n",fromFile);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
    numCoeff = numCoeffLoc;
    stateCoeff = new complex [numCoeffLoc];
    ga = new int [numCoeffLoc];
    gb = new int [numCoeffLoc];
    gc = new int [numCoeffLoc];
    for(int pNo=0;pNo<numCoeff;pNo++) {
      int x,y,z;
      double re,im;
      if(5!=fscanf(fp,"%lg%lg%d%d%d",&re,&im,&x,&y,&z)){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("Can't parse packed state location %s\n",fromFile);
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit();
      }//endif
      stateCoeff[pNo] = complex(re, im);
      ga[pNo]    = x;
      gb[pNo]    = y;
      gc[pNo]    = z;
      nktot++;
      if(doublePack && x==0 && y==0 && z==0){break;}
    }//endfor
    fclose(fp);
  }
#if CMK_PROJECTIONS_USE_ZLIB
  else if(ibinary_opt==2){
    //      CkPrintf("Using ZLIB to load ascii states\n");
    char bigenough[1000];  //we know our lines are shorter than this
    char localFile[1000]; // fromFile is const
    strcpy(localFile,fromFile);
    strcat(localFile,".gz");
    gzFile zfp=gzopen(localFile,"rb");
    if (zfp==NULL){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Can't open state file %s\n",localFile);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
    int numCoeffLoc;
    if(gzgets(zfp,bigenough,1000)!=Z_NULL)
    {
      if(4!=sscanf(bigenough,"%d%d%d%d",&numCoeffLoc,&nx,&ny,&nz)){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("Can't parse size line of file %s\n",localFile);
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit();
      }//endif
      numCoeff = numCoeffLoc;
      stateCoeff = new complex [numCoeffLoc];
      ga = new int [numCoeffLoc];
      gb = new int [numCoeffLoc];
      gc = new int [numCoeffLoc];
    }
    else
    {
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Can't parse size line of file %s\n",localFile);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }
    for(int pNo=0;pNo<numCoeff;pNo++) {
      int x,y,z;
      double re,im;
      if(gzgets(zfp,bigenough,1000)!=Z_NULL)
      {
        if(5!=sscanf(bigenough,"%lg%lg%d%d%d",&re,&im,&x,&y,&z)){
          CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          CkPrintf("Can't parse packed state location %s\n",localFile);
          CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          CkExit();
        }//endif
      }
      else
      {
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("Can't parse packed state location %s\n",localFile);
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit();
      }
      stateCoeff[pNo] = complex(re, im);
      ga[pNo]    = x;
      gb[pNo]    = y;
      gc[pNo]    = z;
      nktot++;
      if(doublePack && x==0 && y==0 && z==0){break;}
    }//endfor
    gzclose(zfp);

  }
  else if (ibinary_opt==3)
  {
    //	CkPrintf("Using ZLIB to load binary states\n");
    char localFile[1000]; // fromFile is const
    strcpy(localFile,fromFile);
    strcat(localFile,".gz");
    gzFile zfp=gzopen(localFile,"rb");
    if (zfp==NULL){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Can't open state file %s\n",localFile);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }
    int numCoeffLoc;
    int n=1;
    gzread(zfp,&(numCoeffLoc),sizeof(int));
    gzread(zfp,&(nx),sizeof(int));
    gzread(zfp,&(ny),sizeof(int));
    gzread(zfp,&(nz),sizeof(int));
    numCoeff = numCoeffLoc;
    stateCoeff = new complex [numCoeffLoc];
    ga = new int [numCoeffLoc];
    gb = new int [numCoeffLoc];
    gc = new int [numCoeffLoc];

    for(int pNo=0;pNo<numCoeff;pNo++) {
      int x,y,z;
      double re,im;
      gzread(zfp,&(re),sizeof(double));
      gzread(zfp,&(im),sizeof(double));
      gzread(zfp,&(x),sizeof(int));
      gzread(zfp,&(y),sizeof(int));
      gzread(zfp,&(z),sizeof(int));
      stateCoeff[pNo] = complex(re, im);
      ga[pNo]    = x;
      gb[pNo]    = y;
      gc[pNo]    = z;
      nktot++;
      if(doublePack && x==0 && y==0 && z==0){break;}
    }
    gzclose(zfp);

  }
#endif
  else{
    FILE *fp=fopen(fromFile,"rb");
    if (fp==NULL){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Can't open state file %s\n",fromFile);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }
    int numCoeffLoc;
    int n=1;
    assert(fread(&(numCoeffLoc),sizeof(int),n,fp)>0);
    assert(fread(&(nx),sizeof(int),n,fp));
    assert(fread(&(ny),sizeof(int),n,fp));
    assert(fread(&(nz),sizeof(int),n,fp));
    numCoeff = numCoeffLoc;
    stateCoeff = new complex [numCoeffLoc];
    ga = new int [numCoeffLoc];
    gb = new int [numCoeffLoc];
    gc = new int [numCoeffLoc];

    for(int pNo=0;pNo<numCoeff;pNo++) {
      int x,y,z;
      double re,im;
      assert(fread(&(re),sizeof(double),n,fp));
      assert(fread(&(im),sizeof(double),n,fp));
      assert(fread(&(x),sizeof(int),n,fp));
      assert(fread(&(y),sizeof(int),n,fp));
      assert(fread(&(z),sizeof(int),n,fp));
      stateCoeff[pNo] = complex(re, im);
      ga[pNo]    = x;
      gb[pNo]    = y;
      gc[pNo]    = z;
      nktot++;
      if(doublePack && x==0 && y==0 && z==0){break;}
    }//endfor
    fclose(fp);
  }//endif

#ifdef _CP_DEBUG_UTIL_VERBOSE_
  CkPrintf("Done reading state from file: %s\n",fromFile);
#endif

  //===================================================================================
  if (nktot!=numCoeff){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Inconsistent number of coefficients in %s\n",fromFile);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
  }

#ifdef VERIFICATION
  char filename[100];
  sprintf(filename,"state_0_%d_%d", thisIndex.y, thisIndex.z);
  FILE *fp = fopen(filename, "w");
  for(int pNo=0;pNo<numCoeff;pNo++) {
    fprintf(fp, "%lg %lg ", stateCoeff[pNo].re, stateCoeff[pNo].im);
  }
  fclose(fp);
#endif
  //----------------------------------------------------------------------------------
}//end routine
//===================================================================================



////////////////////////////////////////////////////////////
// read shifted states
// shifted states are needed to calculate epsilon when q=0
////////////////////////////////////////////////////////////
void States::readStateShifted(char *fromFile) 

  //===================================================================================
{//begin routine
  //===================================================================================
  // A little screen output for the fans
    //CkPrintf("Reading state file: %s for chare (%d %d %d), with binary option %d.\n",fromFile,ispin,ikpt,istate,ibinary_opt);


  if(ibinary_opt < 0 || ibinary_opt > 3){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Bad binary option : %d\n",ibinary_opt);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  //===================================================================================
  // First read in the state and k-vectors : allows parsing of doublePack option

#if !CMK_PROJECTIONS_USE_ZLIB
  if(ibinary_opt>1)
  {
    CkPrintf("Attempt to use ZLIB Failed! Please review compilation\n");
    //CkPrintf("Macro cmk-projections-use-zlib  is %d \n", CMK_PROJECTIONS_USE_ZLIB);
    CkExit();
  }
#endif
  int nx,ny,nz;
  int nktot = 0;
  if(ibinary_opt==0){
    FILE *fp=fopen(fromFile,"r");
    if (fp==NULL){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Can't open state file %s\n",fromFile);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
    int numCoeffLoc;
    if(4!=fscanf(fp,"%d%d%d%d",&numCoeffLoc,&nx,&ny,&nz)){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Can't parse size line of file %s\n",fromFile);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
    numCoeff = numCoeffLoc;
    stateCoeff_shifted = new complex [numCoeffLoc];
    //ga = new int [numCoeffLoc];
    //gb = new int [numCoeffLoc];
    //gc = new int [numCoeffLoc];
    for(int pNo=0;pNo<numCoeff;pNo++) {
      int x,y,z;
      double re,im;
      if(5!=fscanf(fp,"%lg%lg%d%d%d",&re,&im,&x,&y,&z)){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("Can't parse packed state location %s\n",fromFile);
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit();
      }//endif
      stateCoeff_shifted[pNo] = complex(re, im);
      /* we do not need to save this for shifted states
      ga[pNo]    = x;
      gb[pNo]    = y;
      gc[pNo]    = z;
      */
      nktot++;
      if(doublePack && x==0 && y==0 && z==0){break;}
    }//endfor
    fclose(fp);
  }
#if CMK_PROJECTIONS_USE_ZLIB
  else if(ibinary_opt==2){
    //      CkPrintf("Using ZLIB to load ascii states\n");
    char bigenough[1000];  //we know our lines are shorter than this
    char localFile[1000]; // fromFile is const
    strcpy(localFile,fromFile);
    strcat(localFile,".gz");
    gzFile zfp=gzopen(localFile,"rb");
    if (zfp==NULL){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Can't open state file %s\n",localFile);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
    int numCoeffLoc;
    if(gzgets(zfp,bigenough,1000)!=Z_NULL)
    {
      if(4!=sscanf(bigenough,"%d%d%d%d",&numCoeffLoc,&nx,&ny,&nz)){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("Can't parse size line of file %s\n",localFile);
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit();
      }//endif
      numCoeff = numCoeffLoc;
      stateCoeff_shifted = new complex [numCoeffLoc];
      //ga = new int [numCoeffLoc];
      //gb = new int [numCoeffLoc];
      //gc = new int [numCoeffLoc];
    }
    else
    {
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Can't parse size line of file %s\n",localFile);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }
    for(int pNo=0;pNo<numCoeff;pNo++) {
      int x,y,z;
      double re,im;
      if(gzgets(zfp,bigenough,1000)!=Z_NULL)
      {
        if(5!=sscanf(bigenough,"%lg%lg%d%d%d",&re,&im,&x,&y,&z)){
          CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          CkPrintf("Can't parse packed state location %s\n",localFile);
          CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          CkExit();
        }//endif
      }
      else
      {
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("Can't parse packed state location %s\n",localFile);
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit();
      }
      stateCoeff_shifted[pNo] = complex(re, im);
      /* we do not need to save this for shifted states
      ga[pNo]    = x;
      gb[pNo]    = y;
      gc[pNo]    = z;
      */
      nktot++;
      if(doublePack && x==0 && y==0 && z==0){break;}
    }//endfor
    gzclose(zfp);

  }
  else if (ibinary_opt==3)
  {
    //	CkPrintf("Using ZLIB to load binary states\n");
    char localFile[1000]; // fromFile is const
    strcpy(localFile,fromFile);
    strcat(localFile,".gz");
    gzFile zfp=gzopen(localFile,"rb");
    if (zfp==NULL){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Can't open state file %s\n",localFile);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }
    int numCoeffLoc;
    int n=1;
    gzread(zfp,&(numCoeffLoc),sizeof(int));
    gzread(zfp,&(nx),sizeof(int));
    gzread(zfp,&(ny),sizeof(int));
    gzread(zfp,&(nz),sizeof(int));
    numCoeff = numCoeffLoc;
    stateCoeff_shifted = new complex [numCoeffLoc];
    //ga = new int [numCoeffLoc];
    //gb = new int [numCoeffLoc];
    //gc = new int [numCoeffLoc];

    for(int pNo=0;pNo<numCoeff;pNo++) {
      int x,y,z;
      double re,im;
      gzread(zfp,&(re),sizeof(double));
      gzread(zfp,&(im),sizeof(double));
      gzread(zfp,&(x),sizeof(int));
      gzread(zfp,&(y),sizeof(int));
      gzread(zfp,&(z),sizeof(int));
      stateCoeff_shifted[pNo] = complex(re, im);
      /* we do not need to save this for shifted states
      ga[pNo]    = x;
      gb[pNo]    = y;
      gc[pNo]    = z;
      */
      nktot++;
      if(doublePack && x==0 && y==0 && z==0){break;}
    }
    gzclose(zfp);

  }
#endif
  else{
    FILE *fp=fopen(fromFile,"rb");
    if (fp==NULL){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Can't open state file %s\n",fromFile);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }
    int numCoeffLoc;
    int n=1;
    assert(fread(&(numCoeffLoc),sizeof(int),n,fp)>0);
    assert(fread(&(nx),sizeof(int),n,fp));
    assert(fread(&(ny),sizeof(int),n,fp));
    assert(fread(&(nz),sizeof(int),n,fp));
    numCoeff = numCoeffLoc;
    stateCoeff_shifted = new complex [numCoeffLoc];
    //ga = new int [numCoeffLoc];
    //gb = new int [numCoeffLoc];
    //gc = new int [numCoeffLoc];

    for(int pNo=0;pNo<numCoeff;pNo++) {
      int x,y,z;
      double re,im;
      assert(fread(&(re),sizeof(double),n,fp));
      assert(fread(&(im),sizeof(double),n,fp));
      assert(fread(&(x),sizeof(int),n,fp));
      assert(fread(&(y),sizeof(int),n,fp));
      assert(fread(&(z),sizeof(int),n,fp));
      stateCoeff_shifted[pNo] = complex(re, im);
      /* we do not need to save this for shifted states
      ga[pNo]    = x;
      gb[pNo]    = y;
      gc[pNo]    = z;
      */
      nktot++;
      if(doublePack && x==0 && y==0 && z==0){break;}
    }//endfor
    fclose(fp);
  }//endif

#ifdef _CP_DEBUG_UTIL_VERBOSE_
  CkPrintf("Done reading state from file: %s\n",fromFile);
#endif

  //===================================================================================
  if (nktot!=numCoeff){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Inconsistent number of coefficients in %s\n",fromFile);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
  }

  //----------------------------------------------------------------------------------
}//end routine
//===================================================================================










#include "states.def.h"
