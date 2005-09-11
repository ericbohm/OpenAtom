#include "util.h"
#include "para_grp_parse.h"
#include "CPcharmParaInfo.h"
#include "../../src_piny_physics_v1.0/friend_lib/proto_friend_lib_entry.h"
extern Config config;
extern int sizeX;


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================

void readStateIntoRuns(int nPacked, complex *arrCP, CkVec<RunDescriptor> &runs, 
                       const char *fromFile,int ibinary_opt,
                       int *nline_tot_ret,int *nplane_ret) {

//===================================================================================
// A little screen output for the fans

#ifdef _CP_UTIL_VERBOSE_
    CkPrintf("Reading state from file: %s\n",fromFile);
#endif
    if(ibinary_opt < 0 || ibinary_opt > 1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Bad binary option : %d %s\n",ibinary_opt,fromFile);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

//===================================================================================
// First read in the state and k-vectors : allows parsing of doublePack option

    int nx,ny,nz;
    int *kx= new int [nPacked];
    int *ky= new int [nPacked];
    int *kz= new int [nPacked];
    int nktot = 0;
    readState(nPacked, arrCP, fromFile, ibinary_opt, nline_tot_ret, 
              nplane_ret, kx, ky, kz, &nx, &ny, &nz);
    int nplane=*nplane_ret;
    int nline_tot=*nline_tot_ret;

//===================================================================================
// Read the state into the rundescriptor puppy dog
	
    if(!config.doublePack){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("The rundescriptor needs some love for the non-double pack\n"); 
      CkPrintf("It is not consistent with new FFT logic due to input data order\n");
      CkPrintf("If the data is just reordered all should be well, %s\n",fromFile);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    int nrun_tot       = 1;
    int run_length_sum = 0;
    int run_length     = 1;
    int curr_x         = kx[0];
    int curr_y         = ky[0];
    int curr_z         = kz[0];
    if(curr_x<0){curr_x+=nx;}
    if(curr_y<0){curr_y+=ny;}
    if(curr_z<0){curr_z+=nz;}
    int tmpz           = curr_z;
    int nline_tot_now  = 1;

    for(int pNo=1;pNo<nPacked;pNo++) {
      int x = kx[pNo];
      int y = ky[pNo];
      int z = kz[pNo];
      if (x<0){x+=nx;}
      if (y<0){y+=ny;}
      if (z<0){z+=nz;}
      // Count the lines of z by twos
      if (z == tmpz + 1 && x == curr_x && y == curr_y) {
        // same half line : keep counting
        run_length++;
        tmpz += 1;
      }else{
        // We have changed half lines so increment run index
        // Each line of z, constant x,y is stored in 2 run descriptors 
        // Example : -3 -2 -1 is a separate ``run of z''
        //            0 1 2 3 is a separte  ``run of z''
        //            for a line with only a 0 add a zero length descriptor
        //            to represent the missing negative part of the line.
        if(kz[pNo]==0 && kz[(pNo-1)]>=0){
          runs.push_back(RunDescriptor(curr_x,curr_y,curr_z,run_length_sum,0,1));
          nrun_tot      +=1;
	}//endif
        runs.push_back(RunDescriptor(curr_x,curr_y,curr_z,run_length_sum,run_length,1));
        nrun_tot      +=1;
        run_length_sum += run_length;
        curr_x          = x;
        curr_y          = y;
        curr_z          = z;
        tmpz            = z;
        run_length      = 1;
      }//endif
      if(kx[pNo]!=kx[(pNo-1)] || ky[pNo]!=ky[(pNo-1)] ){
        nline_tot_now++;
        if( (nrun_tot-1)/2 != nline_tot_now-1){
          CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	  CkPrintf("Broken Run Desscriptor : %d %d %d : %d %d %d: %d %d\n",
		   kx[pNo],ky[pNo],kz[pNo],kx[(pNo-1)],ky[(pNo-1)],kz[(pNo-1)],
                   nrun_tot-1,nline_tot_now-1);
          CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          CkExit();
        }//endif
      }//endif
    }//endfor
    // Record the last run of z.
    runs.push_back(RunDescriptor(curr_x,curr_y,curr_z,run_length_sum,run_length,1));
    run_length_sum += run_length;

    if(run_length_sum!=nPacked){
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("The rundescriptor didn't assign all pts to runs %s\n",fromFile); 
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkExit();
    }//endif

    if((nrun_tot %2)!=0 || nrun_tot != 2*nline_tot){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("The rundescriptor MUST find an even number of half-lines\n");
      CkPrintf("The rundescriptor MUST find the correct number of lines\n");
      CkPrintf("%d %d %d %d %s\n",nrun_tot,nrun_tot/2,nline_tot,
                                  nline_tot_now,fromFile);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    for(int i=0;i<nline_tot;i+=2){
      RunDescriptor *Desi  = &runs[i];
      RunDescriptor *Desi1 = &runs[(i+1)];
      if( (Desi->x != Desi1->x) || (Desi->y != Desi1->y) ){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("The rundescriptor MUST pair up the half-lines\n");
        CkPrintf("or you will not be a happy camper : %s\n",fromFile);
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit();
      }//endfor
    }//endfor

    delete[] kx;
    delete[] ky;
    delete[] kz;

//===================================================================================
// A little output to the screen!

#ifdef _CP_UTIL_VERBOSE_
     CkPrintf("Done reading state from file: %s\n",fromFile);
#endif

//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================

void readState(int nPacked, complex *arrCP, const char *fromFile,int ibinary_opt,
	       int *nline_tot_ret,int *nplane_ret, int *kx, int *ky, int *kz, 
               int *nx_ret, int *ny_ret, int *nz_ret) {

//===================================================================================
// A little screen output for the fans

#ifdef _CP_UTIL_VERBOSE_
    CkPrintf("Reading state from file: %s\n",fromFile);
#endif
    if(ibinary_opt < 0 || ibinary_opt > 1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Bad binary option : %d\n",ibinary_opt);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

//===================================================================================
// First read in the state and k-vectors : allows parsing of doublePack option

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
         int nPackedLoc;
         if(4!=fscanf(fp,"%d%d%d%d",&nPackedLoc,&nx,&ny,&nz)){
             CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
             CkPrintf("Can't parse size line of file %s\n",fromFile);
             CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
             CkExit();
         }//endif
         for(int pNo=0;pNo<nPacked;pNo++) {
           int x,y,z;
           double re,im;
  	   if(5!=fscanf(fp,"%lg%lg%d%d%d",&re,&im,&x,&y,&z)){
              CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
              CkPrintf("Can't parse packed state location %s\n",fromFile);
              CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
              CkExit();
           }//endif
           arrCP[pNo] = complex(re, im);
           kx[pNo]    = x;
           ky[pNo]    = y;
           kz[pNo]    = z;
           nktot++;
           if(config.doublePack && x==0 && y==0 && z==0){break;}
	}//endfor
       fclose(fp);
    }else{
       FILE *fp=fopen(fromFile,"rb");
         if (fp==NULL){
            CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
            CkPrintf("Can't open state file %s\n",fromFile);
            CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
            CkExit();
         }
         int nPackedLoc;
         int n=1;
         fread(&(nPackedLoc),sizeof(int),n,fp);
         fread(&(nx),sizeof(int),n,fp);
         fread(&(ny),sizeof(int),n,fp);
         fread(&(nz),sizeof(int),n,fp);
         for(int pNo=0;pNo<nPacked;pNo++) {
           int x,y,z;
           double re,im;
           fread(&(re),sizeof(double),n,fp);
           fread(&(im),sizeof(double),n,fp);
           fread(&(x),sizeof(int),n,fp);
           fread(&(y),sizeof(int),n,fp);
           fread(&(z),sizeof(int),n,fp);
           arrCP[pNo] = complex(re, im);
           kx[pNo]    = x;
           ky[pNo]    = y;
           kz[pNo]    = z;
           nktot++;
           if(config.doublePack && x==0 && y==0 && z==0){break;}
	 }//endfor
       fclose(fp);
    }//endif

//===================================================================================
// Add the bottom half of plane zero because the code likes to have it.
// Eventually add reordering for non-doublePack case.

    if(config.doublePack){
       int n_ret;
       ParaGrpParse::flip_data_set(nktot,&n_ret,kx,ky,kz,arrCP);
       if(n_ret!=nPacked){
          CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          CkPrintf("Bad num pts in readState() %s: %d %d \n",nktot,n_ret,fromFile); 
          CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          CkExit();
       }//endif
    }//endif

    int nline_tot = 1;
    int nplane    = 1;
    for(int i = 1;i<nPacked;i++){
      if(kx[i]!=kx[(i-1)]){nplane++;}
      if(kx[i]!=kx[(i-1)] || ky[i]!=ky[(i-1)]){nline_tot++;}
      if(kx[i]<kx[(i-1)]){
         CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
         CkPrintf("Bad x-flip in readState() %s\n",fromFile); 
         CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
         CkExit();
      }
      if(kx[i]==kx[(i-1)]){
       if(ky[i]<ky[(i-1)]){
         CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
         CkPrintf("Bad y-flip in readState() %s\n",fromFile); 
         CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
         CkExit();
       }
      }
      if(kx[i]==kx[(i-1)] && ky[i]==ky[(i-1)]){
        if(kz[i]!=(kz[(i-1)]+1)){
         CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
         CkPrintf("Bad y-flip in readState() %s\n",fromFile); 
         CkPrintf("  %d %d %d\n",kx[i],ky[i],kz[i]);
         CkPrintf("  %d %d %d\n",kx[(i-1)],ky[(i-1)],kz[(i-1)]);
         CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
         CkExit();
	}
      }
    }//endif

//===================================================================================
// Read the state into the rundescriptor puppy dog
	
    if(!config.doublePack){
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("The rundescriptor needs some love for the non-double pack\n"); 
       CkPrintf("It is not consistent with new FFT logic due to input data order\n");
       CkPrintf("If the data is just reordered all should be well, %s\n",fromFile);
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkExit();
    }//endif

//===================================================================================
// A little output to the screen!

    *nplane_ret    = nplane;
    *nline_tot_ret = nline_tot;
    *nx_ret        = nx;
    *ny_ret        = ny;
    *nz_ret        = nz;

#ifdef _CP_UTIL_VERBOSE_
     CkPrintf("Done reading state from file: %s\n",fromFile);
#endif

//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================

void  readStateInfo(int &nPacked,int &minx, int &maxx, int &nx, int &ny, int &nz,
                    const char *fromFile, int ibinary_opt) {

//===================================================================================

#ifdef _CP_UTIL_VERBOSE_
     CkPrintf("Reading state info from file: %s\n",fromFile);
#endif

     if(ibinary_opt < 0 || ibinary_opt > 1){
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("Bad binary option\n",ibinary_opt); CkExit();
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkExit();
     }//endif

     int nktot;
     int nplane0;
     int n=1;
     if(ibinary_opt==0){

        FILE *fp=fopen(fromFile,"r");
         if(fp==NULL){
            CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
            CkPrintf("Can't open state file :%s", fromFile);
            CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
            CkExit();
         }//endif
         if(4!=fscanf(fp,"%d%d%d%d",&nPacked,&nx,&ny,&nz)){
            CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
            CkPrintf("Can't parse size line of file %s\n", fromFile);
            CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
            CkExit();
         }
	 nktot=0;
	 nplane0=0;
         for(int pNo=0;pNo<nPacked;pNo++) {
           double re,im; int x,y,z;
           if(5!=fscanf(fp,"%lg%lg%d%d%d",&re,&im,&x,&y,&z)){
             CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
             CkPrintf("Can't parse packed state location");
             CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
             CkExit();
           }
           if(pNo==0){minx=x; maxx=x;}
           if(x<minx){minx=x;}
           if(x>maxx){maxx=x;}
	   if(x==0){nplane0++;}
	   nktot++;
	   if(x==0 && y==0 && z==0 && config.doublePack)break;
         }//endfor
        fclose(fp);

     }else{

        FILE *fp=fopen(fromFile,"rb");
         if (fp==NULL){
           CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
           CkPrintf("Can't open state file :%s", fromFile);
           CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
           CkExit();
         }
         fread(&(nPacked),sizeof(int),n,fp);
         fread(&(nx),sizeof(int),n,fp);
         fread(&(ny),sizeof(int),n,fp);
         fread(&(nz),sizeof(int),n,fp);
	 nktot=0;
	 nplane0=0;
         for(int pNo=0;pNo<nPacked;pNo++) {
           double re,im; int x,y,z;
           fread(&(re),sizeof(double),n,fp);
           fread(&(im),sizeof(double),n,fp);
           fread(&(x),sizeof(int),n,fp);
           fread(&(y),sizeof(int),n,fp);
           fread(&(z),sizeof(int),n,fp);
           if(pNo==0){minx=x; maxx=x;}
           if(x<minx){minx=x;}
           if(x>maxx){maxx=x;}
	   if(x==0){nplane0++;}
	   nktot++;
	   if(x==0 && y==0 && z==0 && config.doublePack)break;
         }//endfor
        fclose(fp);

     }//endif

     if(minx<0){minx+=nx;}
     n    = minx; 
     minx = maxx; 
     maxx = n;

     // a few extra g-vectors are needed for plane0 than necessary

     if(config.doublePack){
       nPacked=nktot+nplane0-1;
     }//endif

//----------------------------------------------------------------------------------
   }//end routine
//===================================================================================



//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::print() {
//===================================================================================
        CkPrintf("\n");
        ckout     << "gSpacePPC: " << gSpacePPC << endl
                  << "realSpacePPC: " << realSpacePPC << endl
                  << "rhoGPPC: " << rhoGPPC << endl
                  << "dataPath: " << dataPath << endl
                  << "maxIter: " << maxIter << endl
                  << "sGrainSize: " << sGrainSize << endl
                  << "gSpaceNumChunks: " << gSpaceNumChunks << endl
                  << "rhoGHelpers: " << rhoGHelpers << endl
                  << "nstates: " << nstates << endl
                  << "useCommlib: " << useCommlib << endl
                  << "useGMulticast: " << useGMulticast << endl
                  << "useGReduction: " << useGReduction << endl
                  << "multicastDelayMS: " << multicastDelayMS << endl
                  << "numMulticastMsgs: " << numMulticastMsgs << endl
                  << "numPartialReduction: " << numPartialReduction << endl
                  << "useCommlibMulticast: " << useCommlibMulticast << endl
                  << "reductionDelayMS: " << reductionDelayMS << endl
                  << "checkForces: " << checkForces << endl
		  << "doublePack: " << doublePack << endl
		  << "inPlaceFFT: " << inPlaceFFT << endl
		  << "doublePack: " << doublePack << endl
                  << "low_x_size: " << low_x_size <<endl
                  << "high_x_size: " << high_x_size <<endl
	          << "conserveMemory: " << conserveMemory << endl
	          << "lbpaircalc: " << lbpaircalc << endl
   	          << "lbgspace: " << lbgspace << endl
	          << "pesPerState: " << pesPerState << endl
	          << "RpesPerState: " << RpesPerState << endl
	          << "GpesPerState: " << GpesPerState << endl
	          << "localSF: " << localSF << endl
	          << "delayCompStruct: " << delayCompStruct << endl
	          << "gspacesum: " << gspacesum << endl
	          << "parlambda: " << parlambda << endl
	          << "numSfGrps: " << numSfGrps << endl
	          << "numSfDups: " << numSfDups << endl
   	          << "sfpriority: " << sfpriority << endl
	          << "rsfftpriority: " << rsfftpriority << endl
	          << "gsfftpriority: " << gsfftpriority << endl
	          << "rsifftpriority: " << rsifftpriority << endl
	          << "lambdapriority: " << lambdapriority << endl
	          << "psipriority: " << psipriority << endl
	          << "rhorpriority: " << rhorpriority << endl
		  << "rhogpriority: " << rhogpriority << endl;

        CkPrintf("\n");

//----------------------------------------------------------------------------------
   }//end routine
//===================================================================================




//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::readConfig(const char* fileName, Config &config,
			int nstates_in, int nkf1, int nkf2, int nkf3, int maxIter_in,
			int ibinary_opt,int natm_nl)
//===================================================================================
    { // begin routine
//===================================================================================
// Initialize parameters

    config.nstates       = nstates_in;
    config.maxIter       = maxIter_in;
    config.checkForces   = 0;
    config.atomIndex     = 0;
    config.stateIndex    = 0;
    config.planeIndex    = 0;
    config.xIndex        = 0;
    config.yIndex        = 0;
    config.displacement  = 1e-6;
    config.delayComputeZ = 0;
    config.pesPerState   = 1; //Partition the states according to the planes
    config.RpesPerState  = 0; 
    config.GpesPerState  = 0; 
    config.doublePack    = 1;
    config.inPlaceFFT	 = 1;
    config.conserveMemory= 0;
    config.prioFFTMsg    = 0; 
    config.localSF       = 0;
    config.delayCompStruct=0;
    config.lbpaircalc    = 0;
    config.lbgspace      = 0;
    config.fftuseCommlib = 0;
    config.gspacesum     = 0;
    config.parlambda     = 0;
    config.numSfGrps     = 1;
    config.numSfDups     = 1;
    config.gSpacePPC     = 1;
    config.realSpacePPC  = 1;
    config.rhoGPPC       = 1;
    config.gSpaceNumChunks = 1;
    config.sfpriority=    10000000;
    config.rsfftpriority= 1000000;
    config.gsfftpriority= 1000000;
    config.rsifftpriority=100000000;
    config.gsifftpriority=200000000;
    config.lambdapriority=300000000;
    config.psipriority=   400000000;
    config.rhorpriority=  2000000;
    config.rhogpriority=  2000000;// unused
    config.priority=10; //unused?

//===================================================================================
// Read parameters

    CkPrintf("   Opening cpaimd config file : %s\n",fileName);
    ifstream configFile(fileName, ios::in);
    if (configFile.fail()) {
      CkAbort("Bad config file, trouble opening\n");
      CkExit();
    }
    char parameterName[MAX_CHAR_ARRAY_LENGTH];
    char parameterValue[MAX_CHAR_ARRAY_LENGTH];
    while (true) {
	configFile >> parameterName >> parameterValue;
	if (configFile.eof())
	    break;
        config.numSet++;
        if (!strcmp(parameterName, "gSpacePPC"))
            config.gSpacePPC = atoi(parameterValue);
        else if (!strcmp(parameterName, "realSpacePPC"))
            config.realSpacePPC = atoi(parameterValue);
        else if (!strcmp(parameterName, "rhoGPPC"))
            config.rhoGPPC = atoi(parameterValue);
        else if (!strcmp(parameterName, "dataPath"))
            strcpy(config.dataPath, parameterValue);
        else if (!strcmp(parameterName, "sGrainSize"))
            config.sGrainSize = atoi(parameterValue);
        else if (!strcmp(parameterName, "gSpaceNumChunks"))
            config.gSpaceNumChunks = atoi(parameterValue);
        else if (!strcmp(parameterName, "useCommlib"))
            config.useCommlib = atoi(parameterValue);
        else if (!strcmp(parameterName, "useGMulticast"))
            config.useGMulticast = atoi(parameterValue);
        else if (!strcmp(parameterName, "useGReduction"))
            config.useGReduction = atoi(parameterValue);
        else if (!strcmp(parameterName, "numMulticastMsgs"))
            config.numMulticastMsgs = atoi(parameterValue);
        else if (!strcmp(parameterName, "multicastDelayMS"))
            config.multicastDelayMS = atoi(parameterValue);
        else if (!strcmp(parameterName, "numPartialReduction"))
            config.numPartialReduction = atoi(parameterValue);
        else if (!strcmp(parameterName, "reductionDelayMS"))
            config.reductionDelayMS = atoi(parameterValue);
        else if (!strcmp(parameterName, "useCommlibMulticast"))
            config.useCommlibMulticast = atoi(parameterValue);
        else if (!strcmp(parameterName, "rhoGHelpers"))
            config.rhoGHelpers = atoi(parameterValue);
        else if (!strcmp(parameterName, "delayComputeZ"))
            config.delayComputeZ = atoi(parameterValue);
        else if (!strcmp(parameterName, "delayCompStruct"))
            config.delayCompStruct = atoi(parameterValue);
        else if (!strcmp(parameterName, "pesPerState"))
            config.pesPerState = atoi(parameterValue);
        else if (!strcmp(parameterName, "RpesPerState"))
            config.RpesPerState = atoi(parameterValue);
        else if (!strcmp(parameterName, "GpesPerState"))
            config.GpesPerState = atoi(parameterValue);
        //1 or 0 to enable debugging of force computation
        else if (!strcmp(parameterName, "checkForces"))
            config.checkForces = atoi(parameterValue);
        //which atom to displace
        else if (!strcmp(parameterName, "atomIndex"))
            config.atomIndex = atoi(parameterValue);
        //which state to use
        else if (!strcmp(parameterName, "stateIndex"))
            config.stateIndex = atoi(parameterValue);
        //which plane to displace
        else if (!strcmp(parameterName, "planeIndex"))
            config.planeIndex = atoi(parameterValue);
        else if (!strcmp(parameterName, "xIndex"))
            config.xIndex = atoi(parameterValue);
        else if (!strcmp(parameterName, "yIndex"))
            config.yIndex = atoi(parameterValue);
	else if (!strcmp(parameterName, "doublePack"))
            config.doublePack = atoi(parameterValue);
	else if (!strcmp(parameterName, "inPlaceFFT"))
            config.inPlaceFFT = atoi(parameterValue);
	else if (!strcmp(parameterName, "low_x_size"))
	    config.low_x_size = atoi(parameterValue);
	else if (!strcmp(parameterName, "high_x_size"))
	    config.high_x_size = atoi(parameterValue);
	else if (!strcmp(parameterName, "conserveMemory"))
	    config.conserveMemory = atoi(parameterValue);
	else if (!strcmp(parameterName, "prioFFTMsg"))
	    config.prioFFTMsg = atoi(parameterValue);
	else if (!strcmp(parameterName, "localSF"))
	    config.localSF = atoi(parameterValue);
	else if (!strcmp(parameterName, "lbgspace"))
	    config.lbgspace = atoi(parameterValue);
	else if (!strcmp(parameterName, "lbpaircalc"))
	    config.lbpaircalc = atoi(parameterValue);
        else if (!strcmp(parameterName, "fftuseCommlib"))
            config.fftuseCommlib = atoi(parameterValue);
        else if (!strcmp(parameterName, "gspacesum"))
            config.gspacesum = atoi(parameterValue);
        else if (!strcmp(parameterName, "parlambda"))
            config.parlambda = atoi(parameterValue);
        else if (!strcmp(parameterName, "numSfGrps"))
            config.numSfGrps = atoi(parameterValue);
        else if (!strcmp(parameterName, "numSfDups"))
            config.numSfDups = atoi(parameterValue);
        else if (!strcmp(parameterName, "sfpriority"))
            config.sfpriority = atoi(parameterValue);
        else if (!strcmp(parameterName, "rsfftpriority"))
            config.rsfftpriority = atoi(parameterValue);
        else if (!strcmp(parameterName, "gsfftpriority"))
            config.gsfftpriority = atoi(parameterValue);
        else if (!strcmp(parameterName, "rsifftpriority"))
            config.rsifftpriority = atoi(parameterValue);
        else if (!strcmp(parameterName, "gsifftpriority"))
            config.gsifftpriority = atoi(parameterValue);
        else if (!strcmp(parameterName, "lambdapriority"))
            config.lambdapriority = atoi(parameterValue);
        else if (!strcmp(parameterName, "psipriority"))
            config.psipriority = atoi(parameterValue);
        else if (!strcmp(parameterName, "rhorpriority"))
            config.rhorpriority = atoi(parameterValue);
        else if (!strcmp(parameterName, "rhogpriority"))
            config.rhogpriority = atoi(parameterValue);
        else {
            config.numSet --;
            ckout << "Unknown parameter: " << parameterName << endl;
        }
    }
    configFile.close();
    CkPrintf("   Closing cpaimd config file : %s\n",fileName);

    if(config.pesPerState>0 && config.RpesPerState <1){
      config.RpesPerState=config.pesPerState;
    }//endif
    if(config.pesPerState>0 && config.GpesPerState <1){
      config.GpesPerState=config.pesPerState;
    }//endif

//===================================================================================
// Consistency Checks on the input

//-----------------------------------------------------------------------------------
// Set FFT and g-space size


    config.numFFTPoints = nkf1 * nkf2 * nkf3;

    char fname[1024];
    int sizex,sizey,sizez,nPacked,minx,maxx;
    sprintf (fname, "%s/state1.out", config.dataPath);
    CkPrintf("   Opening state file : %s\n",fname);
    readStateInfo(nPacked,minx,maxx,sizex,sizey,sizez,fname,ibinary_opt);
    CkPrintf("   Closing state file : %s\n",fname);
    config.low_x_size  = minx+1;
    config.high_x_size = maxx-1;
    config.numData     = nPacked;

    if(sizex!=nkf1 || sizey!=nkf2 || sizez !=nkf3){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Incorrect FFT size in state files.\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    if(nkf1!=nkf2 || nkf1!=nkf3){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Only Cubic boxes for now\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    if (sizex % config.gSpacePPC != 0) {
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("x dimension should be divisible by gSpacePPC\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }
    if (sizey % config.realSpacePPC != 0) {
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("y dimension should be divisible by realSpacePPC\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }
    if (sizez % config.rhoGPPC != 0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("z dimension should be divisible by rhoGPPC\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }

    if (sizey % config.rhoGHelpers != 0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("y dimension should be divisible by rhoGHelpers\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }

//-----------------------------------------------------------------------------------
// Parameter values that are broken or must be within a certain range

    if (nstates_in % config.sGrainSize != 0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("number of states should be divisible by S matrix grain-size\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }

    if (nstates_in / config.numMulticastMsgs <= 0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Problem in the configuration of number of mcast msgs");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }

    if(config.doublePack!= 1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Non-double Pack code is broken\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    if(config.inPlaceFFT!=1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Non-place in FFT code is broken and not useful\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    if(config.useCommlibMulticast!=1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("No commlib, no work. Sameer is happy!\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    if(config.gSpacePPC!=1 || config.rhoGPPC!=1 || config.realSpacePPC!=1 || 
       config.gSpaceNumChunks!=1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("The PPC and NumChunks have to be unity or the code is horribly\n");
      CkPrintf("broken!\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    if(config.numSfGrps<1 || config.numSfGrps> natm_nl){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("The number of sf atm groups must be >=1 < natm_nl\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    if(config.numSfDups<1||config.numSfDups>config.nstates){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("The number of sf dup groups must be >=1 < num states\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================




//============================================================================
// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void create_line_decomp_descriptor(CPcharmParaInfo *sim)


//============================================================================
     { //begin rotunie
//============================================================================
// Set the file name and data points

    char fname[1000];
    sprintf(fname, "%s/state%d.out", config.dataPath,1);
    int numData     = config.numData;
    int ibinary_opt = sim->ibinary_opt;
    int sizeY = sim->sizeY;
    int sizeZ = sim->sizeZ;
    int doublePack = config.doublePack;

//============================================================================
// Get the complex data, Psi(g) and the run descriptor (z-lines in g-space)

    complex *complexPoints = new complex[numData];
    CkVec<RunDescriptor> runDescriptorVec;
    int nline_tot;
    readStateIntoRuns(numData,complexPoints,runDescriptorVec,fname,ibinary_opt,
                      &nline_tot,&(sim->nplane_x));
    int nplane=sim->nplane_x;

    if(config.low_x_size != nplane && config.doublePack){
       CkPrintf("Mismatch in allowed gspace chare arrays\n");
       CkExit();
    }//endif

//============================================================================
// Parse the run descriptor into integer vectors for use with decomp function
// There are two rundescriptors per line

    int *istrt_line   = new int [nline_tot];
    int *iend_line    = new int [nline_tot];
    int *npts_line    = new int [nline_tot];
    int *kx_line      = new int [nline_tot];
    int *ky_line      = new int [nline_tot];

    int nnn=0;
    for(int i=0,j=0;i<nline_tot;i++,j+=2){
      RunDescriptor *Desi  = &runDescriptorVec[j];
      RunDescriptor *Desi1 = &runDescriptorVec[(j+1)];
      kx_line[i]   = Desi->x;
      ky_line[i]   = Desi->y;
      npts_line[i] = Desi->length + Desi1->length;
      istrt_line[i]= nnn;
      nnn         += npts_line[i];
      iend_line[i] = nnn;
    }//endfor

//============================================================================
// Create the line decomposition and a sorted run descriptor
// There are two rundescriptors per line : Noah's arc sort

    int *istrt_lgrp   = new int [sizeX];
    int *iend_lgrp    = new int [sizeX];
    int *npts_lgrp    = new int [sizeX];
    int *nline_lgrp   = new int [sizeX];
    int *kx_str_lgrp  = new int [sizeX];
    int *kx_end_lgrp  = new int [sizeX];
    int *ky_str_lgrp  = new int [sizeX];
    int *ky_end_lgrp  = new int [sizeX];
    int nktot         = numData;

    ParaGrpParse::get_plane_line_prms(nktot,nplane,nline_tot,npts_line,kx_line,ky_line,
                      istrt_lgrp,iend_lgrp,npts_lgrp,nline_lgrp,
		      kx_str_lgrp,kx_end_lgrp,ky_str_lgrp,ky_end_lgrp);

   

    int nlines_max=0;
    for(int i=0;i<nplane;i++){nlines_max=MAX(nlines_max,nline_lgrp[i]);}
    int **index_tran_upack = cmall_int_mat(0,sizeX,0,nlines_max,"util.C");
   
    int yspace = sizeX;
    if(doublePack){yspace=sizeX/2+1;}
    for(int igrp=0;igrp<nplane;igrp++){
      for(int i=istrt_lgrp[igrp],j=0;i<iend_lgrp[igrp];i++,j++){
        index_tran_upack[igrp][j] = kx_line[i] + ky_line[i]*yspace;
      }//endfor
    }//endfor

    CkVec<RunDescriptor> *sortedRunDescriptors;
    sortedRunDescriptors = new CkVec<RunDescriptor> [sizeX];
    for(int igrp = 0; igrp < nplane; igrp++){
      for(int i=istrt_lgrp[igrp];i<iend_lgrp[igrp];i++){
 	 int j  = 2*i;
 	 int j1 = 2*i+1;
         sortedRunDescriptors[igrp].push_back(runDescriptorVec[j]);
         sortedRunDescriptors[igrp].push_back(runDescriptorVec[j1]);
      }//endfor
    }//endfor

    for(int igrp = nplane; igrp < sizeX; igrp++){
      sortedRunDescriptors[igrp].length() = 0;
    }//endfor

  for(int x = 0; x < nplane; x ++) {
      int runsToBeSent = sortedRunDescriptors[x].size();
      int numPoints    = 0;
      for (int j = 0; j < sortedRunDescriptors[x].size(); j++){
        numPoints += sortedRunDescriptors[x][j].length;
      }//endfor
  }

//============================================================================
// Pack up the stuff, clean up the memory and exit

    sim->npts_per_plane       = npts_lgrp;
    sim->index_tran_upack     = index_tran_upack;
    sim->nlines_max           = nlines_max;
    sim->nlines_per_plane     = nline_lgrp;
    sim->sortedRunDescriptors = sortedRunDescriptors;
    sim->npts_tot             = numData;
    sim->nlines_tot           = nline_tot;

    delete [] complexPoints;

    delete [] istrt_line;
    delete [] iend_line;
    delete [] npts_line;
    delete [] kx_line;
    delete [] ky_line;

    delete [] istrt_lgrp;
    delete [] iend_lgrp;
    delete [] kx_str_lgrp;
    delete [] kx_end_lgrp;
    delete [] ky_str_lgrp;
    delete [] ky_end_lgrp;

//============================================================================
  }//end routine
//============================================================================
