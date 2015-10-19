
//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================

void readState(int ibinary_opt, char *fromFile) 

  //===================================================================================
{//begin routine
  //===================================================================================
  // A little screen output for the fans

  CkPrintf("Reading state file: %s for chare %d.\n",fromFile,thisIndex);

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
      CkPrintf("Inconsistent number of coefficients in %s\n",localFile);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
  }

  //----------------------------------------------------------------------------------
}//end routine
//===================================================================================

