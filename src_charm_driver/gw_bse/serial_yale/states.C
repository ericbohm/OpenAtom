#include "states.decl.h"
#include "standard_include_gwbse.h"
#include "states.h"
#include "main.decl.h"
#include "allclass_gwbse.h"
#include <assert.h>
#if CMK_PROJECTIONS_USE_ZLIB
#include "zlib.h"
#endif
#include "my_fftw.h"
#include "fft_routines.h"

extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_states_occ states_occ_proxy;



// ------------ OCCUPIED STATES MEMBERFUNCTIONS 
states_occ::states_occ() {

  CkPrintf("Hi I'm occupied state %d and I am constructed on processor %d.\n",thisIndex,CkMyPe());
  countdebug = 0;
  GWBSE *gwbse = GWBSE::get();

  // read state file
  ibinary_opt=gwbse->gwbseopts.ibinary_opt;
  // set file name
  ikpt = 0, ispin = 0;
  sprintf(fileName, "./Spin.%d_Kpt.%d_Bead.0_Temper.0/state%d.out", ispin, ikpt, thisIndex+1);
  // now we consider to have only gamma point
  doublePack = true;
  // read states from file
  readState(fileName);

  fft_G_to_R();
  // nothing to do when the states_occ chare object is created.
  //   this is where member variables would be initialized
  //   just like in a c++ class constructor.
}

// constructor needed for chare object migration (ignore for now)
// note: this constructor does not need to appear in the ".ci" file
states_occ::states_occ(CkMigrateMessage *msg) { }

void states_occ::fft_G_to_R(){

  // set fftsize
  set_fftsize();
  
  // set 3D fft plan
  int backward = 1;
  // this routine creates 3D fftw_plan
  setup_fftw_3d(nfft,backward);
  
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
  // in_pointer and out_pointer are set in my_fftw.C
  // put_into_fftbox was originally written for doublePack = 0 (false)
  // for gamma point calculation, put_into_fftbox has been modified from the original version
  put_into_fftbox(numCoeff, stateCoeff, fftidx, nfft, in_pointer, doublePack);
  
  // call fftw
  do_fftw();

  // transfer data from out_pointer to stateCoeffR
  // malloc stateRspace first
  int ndata = nfft[0]*nfft[1]*nfft[2];
  stateCoeffR = new complex [ndata];
  double scale = 1.0; // scale is 1 for now
  fftbox_to_array(ndata, out_pointer, stateCoeffR, scale);

  // delete stateCoeff
  delete [] stateCoeff;
  
}


// set size of fft box. The values are saved in nfft[3] array
void states_occ::set_fftsize(){

  // find the largest and the smallest value of ga, gb, and gc
  // ga >= 0
  int ndata = numCoeff;

  // if ndata=0, there is no data stored. Print error message and exit the program
  if (ndata==0){
    CkPrintf("There seems no data to FFT. I'm exiting the program. Check your states.");
    CkExit();
  }

  int maxga=0, maxgb=0, maxgc=0;
  int minga=0, mingb=0, mingc=0;
  for(int i=0; i<ndata; i++){
    // check if maxga is smaller than ga[i]
    if (maxga == ga[i] || maxga < ga[i]){ maxga = ga[i];}
    if (maxgb == gb[i] || maxgb < gb[i]){ maxgb = gb[i];}
    if (maxgc == gc[i] || maxgc < gc[i]){ maxgc = gc[i];}
    // check if minga is larger than ga[i]
    if (minga == ga[i] || minga > ga[i]){ minga = ga[i];}
    if (mingb == gb[i] || mingb > gb[i]){ mingb = gb[i];}
    if (mingc == gc[i] || mingc > gc[i]){ mingc = gc[i];}
  }

  if (doublePack){
    if (minga!=0){
      CkPrintf("doublePack flag is on, but the minimum index for ga is not zero.\n");
      CkPrintf("Are you sure about this calculation? I'm exiting the program. Check your state.\n");
      CkExit();
    }
  }
  if (doublePack){
    nfft[0] = 2*maxga + 1;
  }
  else{
    nfft[0] = (maxga - minga) + 1;
  }
  nfft[1] = (maxgb - mingb) + 1;
  nfft[2] = (maxgc - mingc) + 1;

  // if nfft is not an even number, make it even
  for (int i=0; i<3; i++){
    if ( nfft[i] % 2 != 0 ){ nfft[i] += 1; }
  }
  
#ifdef DEBUG_FFT
  CkPrintf("----------------\n");
  CkPrintf("size of fft grid: %d %d %d \n",nfft[0], nfft[1], nfft[2]);
#endif
    
}



void states_occ::beamoutMyState(int iteration, int qindex) {

  // Have this chare object say states_occ to the user.
  CkPrintf("\"States_Occ\" from States_Occ chare # %d on "
           "processor %d (for iteration %d and qindex %d).\n",
           thisIndex, CkMyPe(), iteration, qindex);

  // Tell the next chare object in this array of chare objects
  //   to also say states_occ.  If this is the last chare object in
  //   the array of chare objects, then tell the main chare
  //   object to exit the program.
#ifdef _ROUND_ROBIN_
  GWBSE *gwbse = GWBSE::get();
  int nocc = gwbse->gwbseopts.nocc;

  if (thisIndex < (nocc - 1))
    thisProxy[thisIndex + 1].beamoutMyState(0,0);  // message to next chare object
  else
    mainProxy.done();                           // message to main chare object
#endif

  states_occ_proxy[0].exitfordebugging();

}



void states_occ::exitfordebugging(){

  countdebug++;
  GWBSE *gwbse = GWBSE::get();
  int nocc = gwbse->gwbseopts.nocc;
  int nunocc = gwbse->gwbseopts.nunocc;

  CkPrintf("Hi! I am exitfordebugging!! %d %d %d\n",thisIndex, CkMyPe(),countdebug);
  if (countdebug == nocc+nunocc){
    CkPrintf("Finished!\n");
    countdebug = 0;
    mainProxy.done();
  }

}

//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================

void states_occ::readState(char *fromFile) 

  //===================================================================================
{//begin routine
  //===================================================================================
  // A little screen output for the fans
    CkPrintf("Reading state file: %s for chare %d, with binary option %d.\n ",fromFile,thisIndex,ibinary_opt);


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

  //----------------------------------------------------------------------------------
}//end routine
//===================================================================================

















//------------------- UNOCCUPIED STATES MEMBER FUNCTIONS

states_unocc::states_unocc() {

  CkPrintf("Hi I'm unoccupied state %d and I am constructed on processor %d.\n",thisIndex,CkMyPe());

  //fft_G_to_R();
  // nothing to do when the states_occ chare object is created.
  //   this is where member variables would be initialized
  //   just like in a c++ class constructor.
}

// constructor needed for chare object migration (ignore for now)
// note: this constructor does not need to appear in the ".ci" file
states_unocc::states_unocc(CkMigrateMessage *msg) { }

void states_unocc::fft_G_to_R(){
}



void states_unocc::beamoutMyState(int iteration, int qindex) {

  // Have this chare object say states_occ to the user.
  CkPrintf("\"States_unocc\" from States_unocc chare # %d on "
           "processor %d (for iteration %d and qindex %d).\n",
           thisIndex, CkMyPe(), iteration, qindex);

  // Tell the next chare object in this array of chare objects
  //   to also say states_occ.  If this is the last chare object in
  //   the array of chare objects, then tell the main chare
  //   object to exit the program.
#ifdef _ROUND_ROBIN_
  GWBSE *gwbse = GWBSE::get();
  int nunocc = gwbse->gwbseopts.nunocc;
  if (thisIndex < (nunocc - 1))
    thisProxy[thisIndex + 1].beamoutMyState(0,0);  // message to next chare object

  else
    mainProxy.done();                           // message to main chare object
#endif

  states_occ_proxy[0].exitfordebugging();

}






//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================

void states_unocc::readState(char *fromFile) 

  //===================================================================================
{//begin routine
  //===================================================================================
  // A little screen output for the fans

    CkPrintf("Reading state file: %s for chare %d, with binary option %d.\n ",fromFile,thisIndex,ibinary_opt);

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

  //----------------------------------------------------------------------------------
}//end routine
//===================================================================================



#include "states.def.h"
