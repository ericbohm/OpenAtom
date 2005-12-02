//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
/** \file util.C
 *
 */
//===================================================================================

#define CHARM_ON
#include "../../src_piny_physics_v1.0/include/class_defs/piny_constants.h"
#include "util.h"
#include "para_grp_parse.h"
#include "CPcharmParaInfo.h"
#include "../../src_piny_physics_v1.0/friend_lib/proto_friend_lib_entry.h"
#include "../../src_mathlib/mathlib.h"
#include "../../src_piny_physics_v1.0/include/class_defs/allclass_gen.h"
#include "../../src_piny_physics_v1.0/include/class_defs/allclass_cp.h"
extern Config config;
extern int sizeX;
#if CMK_PROJECTIONS_USE_ZLIB
#include "zlib.h"
#endif

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif
//===================================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================

void make_rho_runs(CPcharmParaInfo *sim){

//===================================================================================

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP *cp                     = CP::get();

#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

//===================================================================================
// set up the rho kvectors

   int *kx;
   int *ky;
   int *kz;
   int nline_tot; 
   int nPacked;
   int sizeX     = sim->sizeX;
   int sizeY     = sim->sizeY;
   int sizeZ     = sim->sizeZ;
   double *hmati = gencell->hmati;
   double ecut4  = 8.0*cpcoeffs_info->ecut; // convert to Ryd extra factor of 2.0

   get_rho_kvectors(ecut4,hmati,&kx,&ky,&kz,&nline_tot,&nPacked,sizeX,sizeY,sizeZ);

//===================================================================================
// Reorder the kvectors to produce better balance for the lines : 

    int *kx_ind     = new int[nline_tot];
    int *kx_line    = new int[nline_tot];
    int *ky_line    = new int[nline_tot];
    int *istrt_line = new int [nline_tot];
    int *iend_line  = new int [nline_tot];
    int *npts_line  = new int [nline_tot];

    int nplane_x=0;
    int ic = 0;
    istrt_line[0] = 0;
    for(int i = 1;i<nPacked;i++){
      int iii = abs(kx[i]);
      nplane_x = (iii > nplane_x ? iii : nplane_x);
      if(kx[i]!=kx[(i-1)] || ky[i]!=ky[(i-1)]){
        iend_line[ic] = i;
        npts_line[ic] = iend_line[ic]-istrt_line[ic];
        ic++;
        istrt_line[ic] = i;
      }//endfor
    }//endfor
    iend_line[ic] = nPacked;
    npts_line[ic] = iend_line[ic]-istrt_line[ic];
    ic++;
    if(ic!=nline_tot){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Toasty Line Flip-lines!\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
    nplane_x      += 1;

    double temp    = ((double)nplane_x)*config.gExpandFactRho;
    int nchareRhoG = (int)temp;
    int nx             = sim->sizeX;
    int ny             = sim->sizeY;
    int nz             = sim->sizeZ;
    int *kxt        = new int[nPacked];
    int *kyt        = new int[nPacked];
    int *kzt        = new int[nPacked];
    memcpy(kxt,kx,(nPacked*sizeof(int)));
    memcpy(kyt,ky,(nPacked*sizeof(int)));
    memcpy(kzt,kz,(nPacked*sizeof(int)));

    int jc      = 0;
    int lc      = 0;
    for(int i=0;i<nchareRhoG; i++){
      for(int j=i;j<nline_tot;j+=nchareRhoG){
        for(int lt=istrt_line[j],l=jc;lt<iend_line[j];lt++,l++){
          kx[l]    = kxt[lt];
          ky[l]    = kyt[lt];
          kz[l]    = kzt[lt];
	}//endfor
        jc+=npts_line[j];
        lc++;
      }//endfor
    }//endfor
    if(jc!=nPacked)  {
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Toasty Line Flip-pts!\n");  
      CkExit();
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    }//endif
    if(lc!=nline_tot){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Toasty Line Flip-lines!\n");
      CkExit();
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    }//endif

    delete [] kxt;
    delete [] kyt;
    delete [] kzt;

    ic            = 0;
    istrt_line[0] = 0;
    kx_line[0]    = kx[0];
    ky_line[0]    = ky[0];
    if(ky_line[0]<0){ky_line[0]+=ny;}
    if(kx_line[0]<0){kx_line[0]+=nx;}
    for(int i = 1;i<nPacked;i++){
      if(kx[i]!=kx[(i-1)] || ky[i]!=ky[(i-1)]){
        iend_line[ic] = i;
        npts_line[ic] = iend_line[ic]-istrt_line[ic];
        ic++;
        istrt_line[ic] = i;
        kx_line[ic] = kx[i];
        ky_line[ic] = ky[i];
        if(ky_line[ic]<0){ky_line[ic]+=ny;}
        if(kx_line[ic]<0){kx_line[ic]+=nx;}
      }//endfor
    }//endfor
    iend_line[ic] = nPacked;
    npts_line[ic] = iend_line[ic]-istrt_line[ic];
    ic++;
    if(ic!=nline_tot){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Toasty Line Flip-lines.b!\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

//===================================================================================
// Create runs using reorder k-vectors

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
    CkVec<RunDescriptor> runs;

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
        runs.push_back(RunDescriptor(curr_x,curr_y,curr_z,run_length_sum,run_length,1));
        nrun_tot      +=1;
        run_length_sum += run_length;
        curr_x          = x;
        curr_y          = y;
        curr_z          = z;
        tmpz            = z;
        run_length      = 1;
        if(kz[pNo]==0 && kz[(pNo-1)]>=0){
          runs.push_back(RunDescriptor(curr_x,curr_y,curr_z,run_length_sum,0,1));
          nrun_tot      +=1;
	}//endif
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
       CkPrintf("The rundescriptor didn't assign all pts to rho runs \n"); 
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkExit();
    }//endif

    if((nrun_tot %2)!=0 || nrun_tot != 2*nline_tot){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("The rho rundescriptor MUST find an even number of half-lines\n");
      CkPrintf("The rho rundescriptor MUST find the correct number of lines\n");
      CkPrintf("%d %d %d %d\n",nrun_tot,nrun_tot/2,nline_tot,
                                  nline_tot_now);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
    for(int i=0;i<nline_tot;i+=2){
      RunDescriptor *Desi  = &runs[i];
      RunDescriptor *Desi1 = &runs[(i+1)];
      //      CkPrintf("i %d kx %d kx1 %d ky %d ky1 %d\n",i, Desi->x, Desi1->x, Desi->y, Desi1->y);
      if( (Desi->x != Desi1->x) || (Desi->y != Desi1->y) ){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("The rho rundescriptor MUST pair up the half-lines\n");
        CkPrintf("i %d kx %d kx1 %d ky %d ky1 %d\n",i, Desi->x, Desi1->x, Desi->y, Desi1->y);
	CkPrintf("or you will not be a happy camper :\n");
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit(); 
      }//endfor 
    }//endfor

//===================================================================================
// Decompose lines to balance points

    int *istrt_lgrp   = new int [sizeX];
    int *iend_lgrp    = new int [sizeX];
    int *npts_lgrp    = new int [sizeX];
    int *nline_lgrp   = new int [sizeX];

    ParaGrpParse::get_chareG_line_prms(nPacked,nchareRhoG,nline_tot,npts_line,
                               istrt_lgrp,iend_lgrp,npts_lgrp,nline_lgrp,false);

//============================================================================
// Create the line decomposition and a sorted run descriptor
// There are two rundescriptors per line : Noah's arc sort

    int nlines_max=0;
    for(int i=0;i<nchareRhoG;i++){nlines_max=MAX(nlines_max,nline_lgrp[i]);}
    int **index_tran_upack_rho = cmall_int_mat(0,sizeX,0,nlines_max,"util.C");
   
    int yspace=sizeX/2+1;

    for(int igrp=0;igrp<nchareRhoG;igrp++){
      for(int i=istrt_lgrp[igrp],j=0;i<iend_lgrp[igrp];i++,j++){
        index_tran_upack_rho[igrp][j] = kx_line[i] + ky_line[i]*yspace;
      }//endfor
    }//endfor

    CkVec<RunDescriptor> *RhosortedRunDescriptors;
    RhosortedRunDescriptors = new CkVec<RunDescriptor> [sizeX];
    for(int igrp = 0; igrp < nchareRhoG; igrp++){
      for(int i=istrt_lgrp[igrp];i<iend_lgrp[igrp];i++){
 	 int j  = 2*i;
 	 int j1 = 2*i+1;
         RhosortedRunDescriptors[igrp].push_back(runs[j]);
         RhosortedRunDescriptors[igrp].push_back(runs[j1]);
      }//endfor
    }//endfor

    for(int igrp = nchareRhoG; igrp < sizeX; igrp++){
      RhosortedRunDescriptors[igrp].length() = 0;
    }//endfor

  for(int x = 0; x < nchareRhoG; x ++) {
      int runsToBeSent = RhosortedRunDescriptors[x].size();
      int numPoints    = 0;
      for (int j = 0; j < RhosortedRunDescriptors[x].size(); j++){
        numPoints += RhosortedRunDescriptors[x][j].length;
      }//endfor
  }//endfor

//============================================================================
// variables that could be used for mapping but aren't
  double *pts_per_chare = new double[sizeX];
  double *lines_per_chare = new double[sizeX];
  for(int i=0;i<nchareRhoG;i++)
    {
      pts_per_chare[i]=(double) npts_lgrp[i];
      lines_per_chare[i]=(double) nline_lgrp[i];
    }
//============================================================================
// Pack up the stuff

    sim->nplane_rho_x            = nplane_x;
    sim->nchareRhoG              = nchareRhoG;
    sim->npts_per_chareRhoG      = npts_lgrp;
    sim->index_tran_upack_rho    = index_tran_upack_rho;
    sim->nlines_max_rho          = nlines_max;
    sim->nlines_per_chareRhoG    = nline_lgrp;
    sim->RhosortedRunDescriptors = RhosortedRunDescriptors;
    sim->npts_tot_rho            = nPacked;
    sim->nlines_tot_rho          = nline_tot;
    sim->lines_per_chareRhoG     = lines_per_chare;
    sim->pts_per_chareRhoG       = pts_per_chare;
//============================================================================
// Clean up the memory

    delete [] kx_ind;
    delete [] kx_line;
    delete [] ky_line;
    delete [] istrt_line;
    delete [] iend_line;
    delete [] npts_line;
    delete [] kx;
    delete [] ky;
    delete [] kz;
    delete [] istrt_lgrp;
    delete [] iend_lgrp;

//============================================================================
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void get_rho_kvectors(double ecut4, double *hmati, int **kx_ret, int **ky_ret, 
                      int **kz_ret, int *nline_tot_ret, int *nPacked_ret,
                      int sizeX, int sizeY, int sizeZ)

//============================================================================
  {//begin routine
//============================================================================
// count the k-vectors

  double tpi  = 2.0*M_PI;
  int ka_max  = sizeX/2-1;
  int kb_max  = sizeY/2-1;
  int kc_max  = sizeZ/2-1;

  int nPacked = 0;
  for(int ka=0;ka<=ka_max;ka++){
    for(int kb=-kb_max;kb<=kb_max;kb++){
      for(int kc=-kc_max;kc<=kc_max;kc++){
        double gx = tpi*(ka*hmati[1] + kb*hmati[2] + kc*hmati[3]);
        double gy = tpi*(ka*hmati[4] + kb*hmati[5] + kc*hmati[6]);
        double gz = tpi*(ka*hmati[7] + kb*hmati[8] + kc*hmati[9]);
        double g2 = gx*gx+gy*gy+gz*gz;
        if(g2<=ecut4){nPacked++;}
      }//endfor
    }//endfor
  }//endfor

//============================================================================
// fill the k-vectors

  int *kx = new int[nPacked];
  int *ky = new int[nPacked];
  int *kz = new int[nPacked];

  int ic = 0;
  for(int ka=0;ka<=ka_max;ka++){
    for(int kb=-kb_max;kb<=kb_max;kb++){
     for(int kc=-kc_max;kc<=kc_max;kc++){
       double gx = tpi*(ka*hmati[1] + kb*hmati[2] + kc*hmati[3]);
       double gy = tpi*(ka*hmati[4] + kb*hmati[5] + kc*hmati[6]);
       double gz = tpi*(ka*hmati[7] + kb*hmati[8] + kc*hmati[9]);
       double g2 = gx*gx+gy*gy+gz*gz;
       if(g2<=ecut4){
         kx[ic]=ka;
         ky[ic]=kb;
         kz[ic]=kc;
         ic++;
       }//endif
     }//endfor
   }//endfor
 }//endfor

//============================================================================
// Count the lines

  int nline_tot = 1;
  for(int i=1; i<nPacked;i++){
    if(kx[i]!=kx[(i-1)] || ky[i]!=ky[i-1]){nline_tot++;}
  }//endfor
  
//============================================================================
// Set return values

   *kx_ret        = kx;
   *ky_ret        = ky;
   *kz_ret        = kz;
   *nline_tot_ret = nline_tot;
   *nPacked_ret   = nPacked;

//============================================================================
   }//end routine
//============================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================

void readStateIntoRuns(int nPacked, complex *arrCP, CkVec<RunDescriptor> &runs, 
                       const char *fromFile,int ibinary_opt,
                       int *nline_tot_ret,int *nplane_ret,
                       int *istrt_lgrp,int *iend_lgrp,
                       int *npts_lgrp,int *nline_lgrp,
                       int **kx_line_ret, int **ky_line_ret,
                       int **kx_ret,int **ky_ret, int **kz_ret, int iget_decomp)

//===================================================================================
    {//begin routine
//===================================================================================
// A little screen output for the fans

#ifdef _CP_DEBUG_UTIL_VERBOSE_
    CkPrintf("Reading state from file: %s\n",fromFile);
#endif
    if(ibinary_opt < 0 || ibinary_opt > 3){
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
              nplane_ret, kx, ky, kz, &nx, &ny, &nz,
              istrt_lgrp,iend_lgrp,npts_lgrp,nline_lgrp,iget_decomp,0);
    int nplane    = (*nplane_ret);
    int nline_tot = (*nline_tot_ret);
    int nchareG   = config.nchareG;

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
        runs.push_back(RunDescriptor(curr_x,curr_y,curr_z,run_length_sum,run_length,1));
        nrun_tot      +=1;
        run_length_sum += run_length;
        curr_x          = x;
        curr_y          = y;
        curr_z          = z;
        tmpz            = z;
        run_length      = 1;
        if(kz[pNo]==0 && kz[(pNo-1)]>=0){
          runs.push_back(RunDescriptor(curr_x,curr_y,curr_z,run_length_sum,0,1));
          nrun_tot      +=1;
	}//endif
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

//===================================================================================

    int *kx_line  = new int [nline_tot];
    int *ky_line  = new int [nline_tot];
    int ic        = 0;
    kx_line[0]    = kx[0];
    ky_line[0]    = ky[0];
    if(kx_line[0]<0){kx_line[0]+=nx;}
    if(ky_line[0]<0){ky_line[0]+=ny;}

    for(int i = 1;i<nPacked;i++){
      if(kx[i]!=kx[(i-1)] || ky[i]!=ky[(i-1)]){
        ic++;
        kx_line[ic] = kx[i];
        ky_line[ic] = ky[i];
        if(ky_line[ic]<0){ky_line[ic]+=ny;}
        if(kx_line[ic]<0){kx_line[ic]+=nx;}
      }//endfor
    }//endfor
    ic++;

    if(ic!=nline_tot){
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("Incorrect number of lines : %s\n",fromFile);
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkExit();
    }//endif

    *kx_line_ret  = kx_line;
    *ky_line_ret  = ky_line;

//===================================================================================

    *kx_ret =  kx;
    *ky_ret =  ky;
    *kz_ret =  kz;

//===================================================================================
// A little output to the screen!

#ifdef _CP_DEBUG_UTIL_VERBOSE_
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
               int *nx_ret, int *ny_ret, int *nz_ret,
               int *istrt_lgrp,int *iend_lgrp,int *npts_lgrp,int *nline_lgrp,
               int iget_decomp,int iget_vstate) 

//===================================================================================
   {//begin routine
//===================================================================================
// A little screen output for the fans

    int nchareG = config.nchareG;
    char stuff[25];
    if(iget_vstate==0){
      strcpy(stuff,"coef");
    }else{
      strcpy(stuff,"velocity");
    }
#ifdef _CP_DEBUG_UTIL_VERBOSE_
      CkPrintf("Reading %s state file: %s\n",stuff,fromFile);
#endif
    if(ibinary_opt < 0 || ibinary_opt > 3){
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
            CkPrintf("Can't open %s state file %s\n",stuff,fromFile);
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
              CkPrintf("Can't parse packed %s state location %s\n",stuff,fromFile);
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
    }
#ifdef ZLIB_H
    else if(ibinary_opt==2){
      char bigenough[1000];  //we know our lines are shorter than this
	char localFile[1000]; // fromFile is const
	strcpy(localFile,fromFile);
	strcat(localFile,".gz");
	gzFile zfp=gzopen(localFile,"rb");
         if (zfp==NULL){
            CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
            CkPrintf("Can't open %s state file %s\n",stuff,localFile);
            CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
            CkExit();
         }//endif
         int nPackedLoc;
	 if(gzgets(zfp,bigenough,1000)!=Z_NULL)
	   {
	     if(4!=sscanf(bigenough,"%d%d%d%d",&nPackedLoc,&nx,&ny,&nz)){
	       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	       CkPrintf("Can't parse size line of file %s\n",localFile);
	       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	       CkExit();
	     }//endif
	   }
	 else
	   {
             CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
             CkPrintf("Can't parse size line of file %s\n",localFile);
             CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
             CkExit();
	   }
         for(int pNo=0;pNo<nPacked;pNo++) {
           int x,y,z;
           double re,im;
	   if(gzgets(zfp,bigenough,1000)!=Z_NULL)
	     {
	       if(5!=sscanf(bigenough,"%lg%lg%d%d%d",&re,&im,&x,&y,&z)){
		 CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		 CkPrintf("Can't parse packed %s state location %s\n",stuff,localFile);
		 CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
		 CkExit();
	       }//endif
	     }
	   else
	     {
	       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	       CkPrintf("Can't parse packed %s state location %s\n",stuff,localFile);
	       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	       CkExit();
	     }
           arrCP[pNo] = complex(re, im);
           kx[pNo]    = x;
           ky[pNo]    = y;
           kz[pNo]    = z;
           nktot++;
           if(config.doublePack && x==0 && y==0 && z==0){break;}
	}//endfor
       gzclose(zfp);
    }
    else if (ibinary_opt==3)
      {
	char localFile[1000]; // fromFile is const
	strcpy(localFile,fromFile);
	strcat(localFile,".gz");
	gzFile zfp=gzopen(localFile,"rb");
         if (zfp==NULL){
            CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
            CkPrintf("Can't open %s state file %s\n",stuff,localFile);
            CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
            CkExit();
         }
         int nPackedLoc;
         int n=1;
         gzread(zfp,&(nPackedLoc),sizeof(int));
         gzread(zfp,&(nx),sizeof(int));
         gzread(zfp,&(ny),sizeof(int));
         gzread(zfp,&(nz),sizeof(int));
         for(int pNo=0;pNo<nPacked;pNo++) {
           int x,y,z;
           double re,im;
           gzread(zfp,&(re),sizeof(double));
           gzread(zfp,&(im),sizeof(double));
           gzread(zfp,&(x),sizeof(int));
           gzread(zfp,&(y),sizeof(int));
           gzread(zfp,&(z),sizeof(int));
           arrCP[pNo] = complex(re, im);
           kx[pNo]    = x;
           ky[pNo]    = y;
           kz[pNo]    = z;
           nktot++;
           if(config.doublePack && x==0 && y==0 && z==0){break;}
	 }
	 gzclose(zfp);
    
      }
#endif
    else{
       FILE *fp=fopen(fromFile,"rb");
         if (fp==NULL){
            CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
            CkPrintf("Can't open %s state file %s\n",stuff,fromFile);
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
      }//endif
      if(kx[i]==kx[(i-1)]){
       if(ky[i]<ky[(i-1)]){
         CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
         CkPrintf("Bad y-flip in readState() %s\n",fromFile); 
         CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
         CkExit();
       }//endif
      }//endif
      if(kx[i]==kx[(i-1)] && ky[i]==ky[(i-1)]){
        if(kz[i]!=(kz[(i-1)]+1)){
         CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
         CkPrintf("Bad y-flip in readState() %s\n",fromFile); 
         CkPrintf("  %d %d %d\n",kx[i],ky[i],kz[i]);
         CkPrintf("  %d %d %d\n",kx[(i-1)],ky[(i-1)],kz[(i-1)]);
         CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
         CkExit();
	}//endif
      }//endif
    }//endfor : pts in g-space

    if(config.nchareG>nline_tot){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Dude, too much juice on the chunks. Chill on gExpandFact\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_Error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
    
//===================================================================================
// Reorder the data to produce better balance for the lines : 

    int *kx_ind     = new int[nline_tot];
    int *kx_line    = new int[nline_tot];
    int *ky_line    = new int[nline_tot];
    int *istrt_line = new int [nline_tot];
    int *iend_line  = new int [nline_tot];
    int *npts_line  = new int [nline_tot];
    int ic = 0;
    istrt_line[0] = 0;
    for(int i = 1;i<nPacked;i++){
      if(kx[i]!=kx[(i-1)] || ky[i]!=ky[(i-1)]){
        iend_line[ic] = i;
        npts_line[ic] = iend_line[ic]-istrt_line[ic];
        ic++;
        istrt_line[ic] = i;
      }//endfor
    }//endfor
    iend_line[ic] = nPacked;
    npts_line[ic] = iend_line[ic]-istrt_line[ic];
    ic++;
    if(ic!=nline_tot){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Toasty Line Flip-lines!\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    complex *arrCPt = new complex[nPacked];
    int *kxt        = new int[nPacked];
    int *kyt        = new int[nPacked];
    int *kzt        = new int[nPacked];
    memcpy(arrCPt,arrCP,(nPacked*sizeof(complex)));
    memcpy(kxt,kx,(nPacked*sizeof(int)));
    memcpy(kyt,ky,(nPacked*sizeof(int)));
    memcpy(kzt,kz,(nPacked*sizeof(int)));

    int jc      = 0;
    int lc      = 0;
    for(int i=0;i<nchareG; i++){
      for(int j=i;j<nline_tot;j+=nchareG){
        for(int lt=istrt_line[j],l=jc;lt<iend_line[j];lt++,l++){
          kx[l]    = kxt[lt];
          ky[l]    = kyt[lt];
          kz[l]    = kzt[lt];
          arrCP[l] = arrCPt[lt];
	}//endfor
        jc+=npts_line[j];
        lc++;
      }//endfor
    }//endfor
    if(jc!=nPacked)  {
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Toasty Line Flip-pts!\n");  
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
    if(lc!=nline_tot){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Toasty Line Flip-lines!\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    ic            = 0;
    istrt_line[0] = 0;
    int kxmax     = kx[0];
    kx_line[0]    = kx[0];
    ky_line[0]    = ky[0];
    for(int i = 1;i<nPacked;i++){
      kxmax = MAX(kx[i],kxmax);
      if(kx[i]!=kx[(i-1)] || ky[i]!=ky[(i-1)]){
        iend_line[ic] = i;
        npts_line[ic] = iend_line[ic]-istrt_line[ic];
        ic++;
        istrt_line[ic] = i;
        kx_line[ic] = kx[i];
        ky_line[ic] = ky[i];
      }//endfor
    }//endfor
    iend_line[ic] = nPacked;
    npts_line[ic] = iend_line[ic]-istrt_line[ic];
    ic++;
    if(ic!=nline_tot){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Toasty Line Flip-lines.b!\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

//===================================================================================
// Decompose if necessary : note istrt_lgrp, iend_lgrp only define when iget_decomp==1

    if(iget_decomp==1){
      if(istrt_lgrp==NULL || iend_lgrp == NULL){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("Toasty Line Flip memory!\n");
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit();
      }//endif
      ParaGrpParse::get_chareG_line_prms(nPacked,nchareG,nline_tot,npts_line,
                               istrt_lgrp,iend_lgrp,npts_lgrp,nline_lgrp,true);
    }//endif

//===================================================================================
// For each decomposd chunk : sort on kx
//      note istrt_lgrp, iend_lgrp only define when iget_decomp==1

    int mal_size = MAX(nline_tot,(kxmax+1));
    int *kx_tmp = new int[nline_tot];
    int *ky_tmp = new int[nline_tot];
    int *k_tmp  = new int[mal_size];
    memcpy(arrCPt,arrCP,(nPacked*sizeof(complex)));
    memcpy(kxt,kx,(nPacked*sizeof(int)));
    memcpy(kyt,ky,(nPacked*sizeof(int)));
    memcpy(kzt,kz,(nPacked*sizeof(int)));

    int loff = 0;
    int joff = 0;
    for(int i=0;i<nchareG;i++){
      for(int l=0;l<nline_lgrp[i];l++){
        kx_tmp[l] = kx_line[(l+loff)];
        ky_tmp[l] = ky_line[(l+loff)];
        kx_ind[l] = l;
      }//endfor
      if(nline_lgrp[i]>1){sort_kxky(nline_lgrp[i],kx_tmp,ky_tmp,kx_ind,k_tmp,ny);}
      for(int l=0;l<nline_lgrp[i];l++){
        int istrt = istrt_line[(kx_ind[l]+loff)];
        int iend  = iend_line[(kx_ind[l]+loff)];
        for(int j=istrt,jk=joff;j<iend;j++,jk++){
          arrCP[jk] = arrCPt[j];
          kx[jk]    = kxt[j];
          ky[jk]    = kyt[j];
          kz[jk]    = kzt[j];
	}//endfor
        joff += npts_line[(kx_ind[l]+loff)];
      }//endfor
      loff += nline_lgrp[i];
    }//endfor
    if(joff!=nPacked)  {
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Toasty Line Flip-pts.2!\n");  
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
    if(loff!=nline_tot){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Toasty Line Flip-lines.2!\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

//===================================================================================
// Debug output

#ifdef _CP_DEBUG_LINE_
    double norm = 0;
    for(int i=1;i<nPacked;i++){
       double wght_now = 2.0;
       if(kx[i]==0 && ky[i]<0){wght_now=0.0;}
       if(kx[i]==0 && ky[i]==0 && kz[i]<0){wght_now=0.0;}
       if(kx[i]==0 && ky[i]==0 && kz[i]==0){wght_now=1.0;}
       norm += (wght_now)*arrCP[i].getMagSqr();
    }//endif
    double normt = 0;
    for(int i=1;i<nPacked;i++){
       double wght_now = 2.0;
       if(kxt[i]==0 && kyt[i]<0){wght_now=0.0;}
       if(kxt[i]==0 && kyt[i]==0 && kzt[i]<0){wght_now=0.0;}
       if(kxt[i]==0 && kyt[i]==0 && kzt[i]==0){wght_now=1.0;}
       normt += (wght_now)*arrCPt[i].getMagSqr();
    }//endif
    CkPrintf("state : %g %g\n",norm,normt);
#endif

//===================================================================================
// Fix !double pack

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

    delete [] arrCPt;
    delete [] kxt;
    delete [] kyt;
    delete [] kzt;
    delete [] istrt_line;
    delete [] iend_line;
    delete [] npts_line;
    delete [] kx_line;
    delete [] ky_line;
    delete [] kx_tmp;
    delete [] ky_tmp;
    delete [] k_tmp;
    delete [] kx_ind;

#ifdef _CP_DEBUG_UTIL_VERBOSE_
     CkPrintf("Done reading %s state from file: %s\n",stuff,fromFile);
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

#ifdef _CP_DEBUG_UTIL_VERBOSE_
  CkPrintf("Reading state info from file: %s\n",fromFile);
#endif

  //===================================================================================

  if(ibinary_opt < 0 || ibinary_opt > 3){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Bad binary option\n",ibinary_opt);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  int nktot;
  int nplane0;
  int n=1;
  if(ibinary_opt==0)
    {

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
    }
#ifdef ZLIB_H
  else if(ibinary_opt==2){
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
    int nPackedLoc;
    if(gzgets(zfp,bigenough,1000)!=Z_NULL)
      {
	if(4!=sscanf(bigenough,"%d%d%d%d",&nPacked,&nx,&ny,&nz)){
	  CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	  CkPrintf("Can't parse size line of file %s\n", localFile);
	  CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	  CkExit();
	}
      }
    else
      {
	CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	CkPrintf("Can't parse size line of file %s\n", localFile);
	CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	CkExit();
      }
    nktot=0;
    nplane0=0;
    for(int pNo=0;pNo<nPacked;pNo++) {
      double re,im; int x,y,z;
      if(gzgets(zfp,bigenough,1000)!=Z_NULL)
	{
	  if(5!=sscanf(bigenough,"%lg%lg%d%d%d",&re,&im,&x,&y,&z)){
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
	}//endif
      else
	{
	  CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	  CkPrintf("Can't parse size line of file %s\n", localFile);
	  CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	  CkExit();
	}
    }//endfor
    gzclose(zfp);
  }
  else if(ibinary_opt==3){
    char localFile[1000]; // fromFile is const
    strcpy(localFile,fromFile);
    strcat(localFile,".gz");
    gzFile zfp=gzopen(localFile,"rb");
      if (zfp==NULL){
	CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	CkPrintf("Can't open state file :%s", localFile);
	CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	CkExit();
      }//endif
      gzread(zfp,&(nPacked),sizeof(int));
      gzread(zfp,&(nx),sizeof(int));
      gzread(zfp,&(ny),sizeof(int));
      gzread(zfp,&(nz),sizeof(int));
      nktot=0;
      nplane0=0;
      for(int pNo=0;pNo<nPacked;pNo++) {
	double re,im; int x,y,z;
	gzread(zfp,&(re),sizeof(double));
	gzread(zfp,&(im),sizeof(double));
	gzread(zfp,&(x),sizeof(int));
	gzread(zfp,&(y),sizeof(int));
	gzread(zfp,&(z),sizeof(int));
	if(pNo==0){minx=x; maxx=x;}
	if(x<minx){minx=x;}
	if(x>maxx){maxx=x;}
	if(x==0){nplane0++;}
	nktot++;
	if(x==0 && y==0 && z==0 && config.doublePack)break;
      }//endfor
      gzclose(zfp);
    }
#endif
    else{
      FILE *fp=fopen(fromFile,"rb");
      if (fp==NULL){
	CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	CkPrintf("Can't open state file :%s", fromFile);
	CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	CkExit();
      }//endif
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
void Config::print(char *fname_in) {
//===================================================================================

   char fname[1024];
   sprintf(fname,"%s.out",fname_in);
   CkPrintf("   Writing all parameters both specified and default, to %s\n",fname);

   FILE *fp = fopen(fname,"w");
     fprintf(fp,"dataPath: %s\n",dataPath);
     fprintf(fp,"gSpacePPC: %d\n",gSpacePPC);
     fprintf(fp,"realSpacePPC: %d\n",realSpacePPC);
     fprintf(fp,"rhoGPPC: %d\n",rhoGPPC);
     fprintf(fp,"maxIter: %d\n",maxIter);
     fprintf(fp,"sGrainSize: %d\n",sGrainSize);
     fprintf(fp,"gSpaceNumChunks: %d\n",gSpaceNumChunks);
     fprintf(fp,"rhoGHelpers: %d\n",rhoGHelpers);
     fprintf(fp,"nstates: %d\n",nstates);
     fprintf(fp,"useCommlib: %d\n",useCommlib);
     fprintf(fp,"useGMulticast: %d\n",useGMulticast);
     fprintf(fp,"useGReduction: %d\n",useGReduction);
     fprintf(fp,"multicastDelayMS: %d\n",multicastDelayMS);
     fprintf(fp,"numMulticastMsgs: %d\n",numMulticastMsgs);
     fprintf(fp,"numPartialReduction: %d\n",numPartialReduction);
     fprintf(fp,"useCommlibMulticast: %d\n",useCommlibMulticast);
     fprintf(fp,"reductionDelayMS: %d\n",reductionDelayMS);
     fprintf(fp,"checkForces: %d\n",checkForces);
     fprintf(fp,"doublePack: %d\n",doublePack);
     fprintf(fp,"inPlaceFFT: %d\n",inPlaceFFT);
     fprintf(fp,"doublePack: %d\n",doublePack);
     fprintf(fp,"low_x_size: %d\n",low_x_size);
     fprintf(fp,"high_x_size: %d\n",high_x_size);
     fprintf(fp,"conserveMemory: %d\n",conserveMemory);
     fprintf(fp,"lbpaircalc: %d\n",lbpaircalc);
     fprintf(fp,"lbgspace: %d\n",lbgspace);
     fprintf(fp,"pesPerState: %d\n",pesPerState);
     fprintf(fp,"RpesPerState: %d\n",RpesPerState);
     fprintf(fp,"GpesPerState: %d\n",GpesPerState);
     fprintf(fp,"localSF: %d\n",localSF);
     fprintf(fp,"delayCompStruct: %d\n",delayCompStruct);
     fprintf(fp,"gspacesum: %d\n",gspacesum);
     fprintf(fp,"numSfGrps: %d\n",numSfGrps);
     fprintf(fp,"numSfDups: %d\n",numSfDups);
     fprintf(fp,"sfpriority: %d\n",sfpriority);
     fprintf(fp,"rsfftpriority: %d\n",rsfftpriority);
     fprintf(fp,"gsfftpriority: %d\n",gsfftpriority);
     fprintf(fp,"rsifftpriority: %d\n",rsifftpriority);
     fprintf(fp,"lambdapriority: %d\n",lambdapriority);
     fprintf(fp,"psipriority: %d\n",psipriority);
     fprintf(fp,"rhogpriority: %d\n",rhogpriority);
     fprintf(fp,"fftprogresssplit: %d\n",fftprogresssplit);
     fprintf(fp,"stateOutputOn: %d\n",stateOutputOn);
     fprintf(fp,"gExpandFact: %g\n",gExpandFact);
     fprintf(fp,"gExpandFactRho: %g\n",gExpandFactRho);
     fprintf(fp,"toleranceInterval: %d\n",toleranceInterval);
   fclose(fp);


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

    config.rhoGHelpers     = 1;
    config.nstates         = nstates_in;
    config.maxIter         = maxIter_in;
    config.checkForces     = 0;
    config.atomIndex       = 0;
    config.stateIndex      = 0;
    config.planeIndex      = 0;
    config.xIndex          = 0;
    config.yIndex          = 0;
    config.displacement    = 1e-6;
    config.delayComputeZ   = 0;
    config.pesPerState     = 1;    //Partition the states g-chares amongst the procs
    config.RpesPerState    = 0; 
    config.GpesPerState    = 0; 
    config.doublePack      = 1;
    config.inPlaceFFT	   = 1;
    config.conserveMemory  = 0;
    config.prioFFTMsg      = 0; 
    config.localSF         = 0;
    config.delayCompStruct  =0;
    config.lbpaircalc      = 0;
    config.lbgspace        = 0;
    config.fftuseCommlib   = 0;
    config.gspacesum       = 0;
    config.numSfGrps       = 1;
    config.numSfDups       = 1;
    config.gSpacePPC       = 1;
    config.realSpacePPC    = 1;
    config.rhoGPPC         = 1;
    config.gSpaceNumChunks = 1;
    config.sfpriority      = 10000000;
    config.rsfftpriority   = 1000000;
    config.gsfftpriority   = 1000000;
    config.rsifftpriority  = 100000000;
    config.gsifftpriority  = 200000000;
    config.lambdapriority  = 300000000;
    config.psipriority     = 400000000;
    config.rhorpriority    = 2000000;
    config.rhogpriority    = 2000000;   // unused?
    config.priority        = 10;        // unused?
    config.gExpandFact     = 1.0;
    config.gExpandFactRho  = 1.0;
    config.fftprogresssplit= 20;
    config.stateOutputOn  =  0;
    config.toleranceInterval=5;
//===================================================================================
// Read parameters

    CkPrintf("   Opening cpaimd config file : %s\n",fileName);
    ifstream configFile(fileName, ios::in);

    if (configFile.fail()) {
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkAbort("Bad config file, trouble opening\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

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
        else if (!strcmp(parameterName, "fftprogresssplit"))
            config.fftprogresssplit = atoi(parameterValue);
        else if (!strcmp(parameterName, "stateOutputOn"))
            config.stateOutputOn = atoi(parameterValue);
        else if (!strcmp(parameterName, "toleranceInterval"))
            config.toleranceInterval = atoi(parameterValue);
        else if (!strcmp(parameterName, "parlambda"))
            CkPrintf("      Warning : parlambda is compulsory. It can't be set.\n");
        else if (!strcmp(parameterName, "gExpandFact")){
               sscanf(parameterValue,"%lg",&(config.gExpandFact));
               }
        else if (!strcmp(parameterName, "gExpandFactRho")){
               sscanf(parameterValue,"%lg",&(config.gExpandFactRho));
               }
        else {
            config.numSet --;
            CkPrintf("@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@\n");
            ckout << "Unknown parameter: " << parameterName << endl;
            CkPrintf("@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@\n");
            CkExit();
        }//endif
    }//end while reading
    configFile.close();
    CkPrintf("   Closing cpaimd config file : %s\n\n",fileName);

    if(config.pesPerState>0 && config.RpesPerState <1){
      config.RpesPerState=config.pesPerState;
    }//endif

    if(config.pesPerState>0 && config.GpesPerState <1){
      config.GpesPerState=config.pesPerState;
    }//endif

//===================================================================================
// Set FFT and g-space size

    char fname[1024];
    int sizex,sizey,sizez,nPacked,minx,maxx;
    sprintf (fname, "%s/state1.out", config.dataPath);
    CkPrintf("   Opening state file : %s\n",fname);
    readStateInfo(nPacked,minx,maxx,sizex,sizey,sizez,fname,ibinary_opt);
    CkPrintf("   Closing state file : %s\n\n",fname);

    config.numFFTPoints = nkf1 * nkf2 * nkf3;
    config.low_x_size   = minx+1;
    config.high_x_size  = maxx-1;
    config.numData      = nPacked;
    int nplane_x        = minx+1;
    double temp         = (config.gExpandFact)*((double)nplane_x);
    int nchareG         = ((int)temp);
//    nchareG             = MIN(nchareG,sizex);
    config.nchareG      = nchareG;
    int nplane_x_rho        = 2*minx+1;
    double temp_rho         = (config.gExpandFactRho)*((double)nplane_x_rho);
    int nchareRhoG         = ((int)temp_rho);
    config.nchareRhoG      = nchareRhoG;

//===================================================================================
// Consistency Checks on the input

    if(config.gExpandFact<1.0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Chare array expansion factor out of range\n");
      CkPrintf("This probably could work but I'd check first.\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    if(nchareG<nplane_x){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Too few g-space chares %d %d\n",nplane_x,nchareG);
      CkPrintf("This probably could work but I'd check first.\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    if(nchareRhoG<nplane_x_rho){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Too few rhog-space chares %d %d\n",nplane_x_rho,nchareRhoG);
      CkPrintf("This probably could work but I'd check first.\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    if(config.pesPerState>config.nchareG){
      CkPrintf("Warning : pesPerState > %d %g %d\n",config.nchareG,
                 config.gExpandFact,sizex);
    }//endif

    if(config.nchareG>sizex){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Error: nchareG> sizex : reduce gExpandFact\n");
      CkPrintf("Memory allocations needs love in pups and elsewhere!\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    if(config.nchareRhoG>sizex){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Error: nchareRhoG> sizex : reduce gExpandFactRho\n");
      CkPrintf("Memory allocations needs love in pups and elsewhere!\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

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
    }//endif

    if (sizey % config.realSpacePPC != 0) {
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("y dimension should be divisible by realSpacePPC\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    if (sizez % config.rhoGPPC != 0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("z dimension should be divisible by rhoGPPC\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    if (sizey % config.rhoGHelpers != 0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("y dimension should be divisible by rhoGHelpers\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

//-----------------------------------------------------------------------------------
// Parameter values that are broken or must be within a certain range

    if (nstates_in % config.sGrainSize != 0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("number of states should be divisible by S matrix grain-size\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    if (nstates_in / config.numMulticastMsgs <= 0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Problem in the configuration of number of mcast msgs");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

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

    if(config.stateOutputOn<0 ||config.stateOutputOn>1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("The state output flag must be 1(on) or 0 (off) \n");
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

    char fname[1000]; sprintf(fname,"%s/state%d.out",config.dataPath,1);
    int numData        = config.numData;
    int ibinary_opt    = sim->ibinary_opt;
    int sizeY          = sim->sizeY;
    int sizeZ          = sim->sizeZ;
    int doublePack     = config.doublePack;
    double gExpandFact = config.gExpandFact;

//============================================================================
// Get the complex data, Psi(g) and the run descriptor (z-lines in g-space)

    complex *complexPoints = new complex[numData];
    CkVec<RunDescriptor> runDescriptorVec;
    int nline_tot;
    int *istrt_lgrp   = new int [sizeX];
    int *iend_lgrp    = new int [sizeX];
    int *npts_lgrp    = new int [sizeX];
    int *nline_lgrp   = new int [sizeX];
    int *kx_line      = NULL;
    int *ky_line      = NULL;
    int *kx           = NULL;
    int *ky           = NULL;
    int *kz           = NULL;

    readStateIntoRuns(numData,complexPoints,runDescriptorVec,fname,ibinary_opt,
                      &nline_tot,&(sim->nplane_x),istrt_lgrp,iend_lgrp,
                      npts_lgrp,nline_lgrp,&kx_line,&ky_line,&kx,&ky,&kz,1);
    int nplane  = sim->nplane_x;
    int nchareG = sim->nchareG;

    if(config.low_x_size != nplane && config.doublePack){
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("Mismatch in allowed gspace planes %d %d\n",config.low_x_size,nplane);
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkExit();
    }//endif

    double temp  = ((double)nplane)*gExpandFact;
    int mychareG = (int)temp;
    if(mychareG!=nchareG || mychareG!=config.nchareG){
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("Mismatch in allowed gspace chare arrays %d %d %d\n",
             mychareG,nchareG,config.nchareG);
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkExit();
    }//endif

//============================================================================
// Create the line decomposition and a sorted run descriptor
// There are two rundescriptors per line : Noah's arc sort

    int nlines_max=0;
    for(int i=0;i<nchareG;i++){nlines_max=MAX(nlines_max,nline_lgrp[i]);}
    int **index_tran_upack = cmall_int_mat(0,sizeX,0,nlines_max,"util.C");
   
    int yspace = sizeX;
    if(doublePack){yspace=sizeX/2+1;}
    for(int igrp=0;igrp<nchareG;igrp++){
      for(int i=istrt_lgrp[igrp],j=0;i<iend_lgrp[igrp];i++,j++){
        index_tran_upack[igrp][j] = kx_line[i] + ky_line[i]*yspace;
      }//endfor
    }//endfor

    CkVec<RunDescriptor> *sortedRunDescriptors;
    sortedRunDescriptors = new CkVec<RunDescriptor> [sizeX];
    for(int igrp = 0; igrp < nchareG; igrp++){
      for(int i=istrt_lgrp[igrp];i<iend_lgrp[igrp];i++){
 	 int j  = 2*i;
 	 int j1 = 2*i+1;
         sortedRunDescriptors[igrp].push_back(runDescriptorVec[j]);
         sortedRunDescriptors[igrp].push_back(runDescriptorVec[j1]);
      }//endfor
    }//endfor

    for(int igrp = nchareG; igrp < sizeX; igrp++){
      sortedRunDescriptors[igrp].length() = 0;
    }//endfor

    int *index_output_off = new int[nchareG];
    int numPoints    = 0;
    for(int x = 0; x < nchareG; x ++) {
      index_output_off[x] = numPoints;
      int nnn = 0;
      int runsToBeSent = sortedRunDescriptors[x].size();
      for (int j = 0; j < runsToBeSent; j++){
        numPoints += sortedRunDescriptors[x][j].length;
        nnn       += sortedRunDescriptors[x][j].length;
      }//endfor
      if(nnn != npts_lgrp[x]){
         CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
         CkPrintf("Incorrect number of points in gspace chare\n");
         CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
         CkExit();
      }//endif
    }//endfor
 
    if(numPoints!=numData){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Incorrect number of total g-space points\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

//============================================================================
// Figure out the unique/redundent communication :
//   Send out the unique : Receive and over write the redudant.

    // find the unique and redundent g's on each chare
    int *num_uni = new int [nchareG];
    int *num_red = new int [nchareG];
    int nk0_max=0;
    for(int i = 0; i < nchareG; i ++) {
      num_uni[i]=0;
      num_red[i]=0;
      for(int j=0;j<npts_lgrp[i];j++){
        int iii = j+index_output_off[i];
        if(kx[iii]==0 && ky[iii]>0){num_uni[i]++;}
        if(kx[iii]==0 && ky[iii]<0){num_red[i]++;}
        if(kx[iii]==0 && ky[iii]==0 && kz[iii]>=0){num_uni[i]++;}
        if(kx[iii]==0 && ky[iii]==0 && kz[iii]<0){num_red[i]++;}
      }//endif
      nk0_max=MAX(num_uni[i],nk0_max);
      nk0_max=MAX(num_red[i],nk0_max);
    }//endif

    // Find where my unique guys go and make a list so I can send them.
    // Make a list of where the unique guys arrive so I can receive them.
    RedundantCommPkg *RCommPkg = new RedundantCommPkg [nchareG]; 
    for(int i=0;i<nchareG;i++){RCommPkg[i].Init(nk0_max,nchareG);}

    for(int i=0;i<nchareG;i++){
      int  *num_send = RCommPkg[i].num_send;
      int **lst_send = RCommPkg[i].lst_send;
      for(int j=0;j<nchareG;j++){
        int  *num_recv = RCommPkg[j].num_recv;
        int **lst_recv = RCommPkg[j].lst_recv;
        for(int ip=0;ip<num_uni[i];ip++){
          int iii = ip+index_output_off[i] + num_red[i];
          for(int jp=0;jp<num_red[j];jp++){
            int jjj = jp+index_output_off[j];
            if(ky[iii]==-ky[jjj] && kz[iii]==-kz[jjj]){
              lst_send[j][num_send[j]] = ip+num_red[i];
              lst_recv[i][num_recv[i]] = jp;
              num_send[j]++;
              num_recv[i]++;
	    }//endif
  	  }//endfor
        }//endfor
      }//endfor
    }//endfor

    for(int i=0;i<nchareG;i++){
      int  *num_send   = RCommPkg[i].num_send;
      int  *num_recv   = RCommPkg[i].num_recv;
      int num_recv_tot = 0;
      int num_send_tot = 0;
      for(int j=0;j<nchareG;j++){ 
        if(num_send[j]>0){num_send_tot++;}
        if(num_recv[j]>0){num_recv_tot++;}
      }//endif
      RCommPkg[i].num_recv_tot = num_recv_tot;
      RCommPkg[i].num_send_tot = num_send_tot;
    }//endfor

    for(int i=0;i<nchareG;i++){
      int  *num_send = RCommPkg[i].num_send;
      int **lst_send = RCommPkg[i].lst_send;
      for(int j=0;j<nchareG;j++){
        int  *num_recv = RCommPkg[j].num_recv;
        int **lst_recv = RCommPkg[j].lst_recv;
        if(num_send[j]!=num_recv[i]){
          CkPrintf("Bad Num of Redundent g-vectors %d %d\n",i,j);
          CkExit();
	}//endif
        for(int ip=0;ip<num_send[j];ip++){
          int iii = index_output_off[i]+lst_send[j][ip];
          int jjj = index_output_off[j]+lst_recv[i][ip];
          if(kx[iii]!=0 || kx[jjj]!=0 || ky[iii]!=-ky[jjj] || kz[iii]!=-kz[jjj] ||
             ky[iii]<0 || (kz[iii]<0 && ky[iii]==0)){
            CkPrintf("Bad Num of Redundent g-vectors %d %d %d %d %d\n",iii,jjj,ip,i,j);
            CkPrintf("%d %d %d : %d %d %d\n",kx[iii],ky[iii],kz[iii],
                        		     kx[jjj],ky[jjj],kz[jjj]);
            CkExit();
   	  }//endif
	}//endfor
      }//endfor
    }//endfor

//============================================================================
// Pack up the stuff, clean up the memory and exit

    sim->npts_per_chareG      = npts_lgrp;
    sim->index_tran_upack     = index_tran_upack;
    sim->nlines_max           = nlines_max;
    sim->nlines_per_chareG    = nline_lgrp;
    sim->sortedRunDescriptors = sortedRunDescriptors;
    sim->npts_tot             = numData;
    sim->nlines_tot           = nline_tot;
    sim->index_output_off     = index_output_off;
    sim->RCommPkg             = RCommPkg;

    delete [] istrt_lgrp;
    delete [] iend_lgrp;
    delete [] complexPoints;
    delete [] kx_line;
    delete [] ky_line;
    delete [] kx;
    delete [] ky;
    delete [] kz;
    delete [] num_uni;
    delete [] num_red;

//============================================================================
  }//end routine
//============================================================================


//=============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=============================================================================
void writeStateFile(int ncoef,complex *psi,complex *vpsi,
                    int *k_x,int *k_y,int *k_z,int cp_min_opt,
                    int sizeX,int sizeY,int sizeZ,char *psiName,char *vpsiName,
                    int ibinary_write_opt)
//=============================================================================
  { //begin rotunie
//=============================================================================

  int *index = new int [ncoef];
  int *ktemp = new int [ncoef];
  int istrt=0;
  sort_psi_output(ncoef,k_x,k_y,k_z,index,ktemp,&istrt);
  int ncoef_true=ncoef-istrt;

   if(ibinary_write_opt==0){
     FILE *fp  = fopen(psiName,"w");
     fprintf(fp,"%d %d %d %d\n",ncoef_true,sizeX,sizeY,sizeZ);
     for(int i=istrt+1;i<ncoef;i++){
       fprintf(fp,"%g %g %d %d %d \n",psi[index[i]].re,psi[index[i]].im,
	       k_x[index[i]],k_y[index[i]],k_z[index[i]]);
     }//endfor
     int i = istrt;
     fprintf(fp,"%g %g %d %d %d \n",psi[index[i]].re,psi[index[i]].im,
	     k_x[index[i]],k_y[index[i]],k_z[index[i]]);
     fclose(fp);
   }
#ifdef ZLIB_H
   else if(ibinary_write_opt==2){
     strcat(psiName,".gz");
     gzFile zfp  = gzopen(psiName,"w");
     gzprintf(zfp,"%d %d %d %d\n",ncoef_true,sizeX,sizeY,sizeZ);
     for(int i=istrt+1;i<ncoef;i++){
       gzprintf(zfp,"%g %g %d %d %d \n",psi[index[i]].re,psi[index[i]].im,
	       k_x[index[i]],k_y[index[i]],k_z[index[i]]);
     }//endfor
     int i = istrt;
     gzprintf(zfp,"%g %g %d %d %d \n",psi[index[i]].re,psi[index[i]].im,
	     k_x[index[i]],k_y[index[i]],k_z[index[i]]);
     gzclose(zfp);
   }
   else if(ibinary_write_opt==3){
     strcat(psiName,".gz");
     gzFile zfp  = gzopen(psiName,"w");
     int n=1;
     gzwrite(zfp,&ncoef_true,sizeof(int));
     gzwrite(zfp,&sizeX,sizeof(int));
     gzwrite(zfp,&sizeY,sizeof(int));
     gzwrite(zfp,&sizeZ,sizeof(int));
     for(int i=istrt+1;i<ncoef;i++){
       gzwrite(zfp,&psi[index[i]].re,sizeof(double));
       gzwrite(zfp,&psi[index[i]].im,sizeof(double));
       gzwrite(zfp,&k_x[index[i]],sizeof(int));
       gzwrite(zfp,&k_y[index[i]],sizeof(int));
       gzwrite(zfp,&k_z[index[i]],sizeof(int));
     }//endfor
     int i = istrt;
     gzwrite(zfp,&psi[index[i]].re,sizeof(double));
     gzwrite(zfp,&psi[index[i]].im,sizeof(double));
     gzwrite(zfp,&k_x[index[i]],sizeof(int));
     gzwrite(zfp,&k_y[index[i]],sizeof(int));
     gzwrite(zfp,&k_z[index[i]],sizeof(int));
     gzclose(zfp);
   }
#endif
   else{
     FILE *fp  = fopen(psiName,"w");
     int n=1;
     fwrite(&ncoef_true,sizeof(int),n,fp);
     fwrite(&sizeX,sizeof(int),n,fp);
     fwrite(&sizeY,sizeof(int),n,fp);
     fwrite(&sizeZ,sizeof(int),n,fp);
     for(int i=istrt+1;i<ncoef;i++){
       fwrite(&psi[index[i]].re,sizeof(double),n,fp);
       fwrite(&psi[index[i]].im,sizeof(double),n,fp);
       fwrite(&k_x[index[i]],sizeof(int),n,fp);
       fwrite(&k_y[index[i]],sizeof(int),n,fp);
       fwrite(&k_z[index[i]],sizeof(int),n,fp);
     }//endfor
     int i = istrt;
     fwrite(&psi[index[i]].re,sizeof(double),n,fp);
     fwrite(&psi[index[i]].im,sizeof(double),n,fp);
     fwrite(&k_x[index[i]],sizeof(int),n,fp);
     fwrite(&k_y[index[i]],sizeof(int),n,fp);
     fwrite(&k_z[index[i]],sizeof(int),n,fp);
     fclose(fp);
   }//endif


  if(cp_min_opt==0){

    if(ibinary_write_opt==0){
      FILE *fp  = fopen(vpsiName,"w");
      fprintf(fp,"%d %d %d %d\n",ncoef_true,sizeX,sizeY,sizeZ);
      for(int i=istrt+1;i<ncoef;i++){
        fprintf(fp,"%g %g %d %d %d \n",vpsi[index[i]].re,vpsi[index[i]].im,
                   k_x[index[i]],k_y[index[i]],k_z[index[i]]);
      }//endfor
      int i = istrt;
      fprintf(fp,"%g %g %d %d %d \n",vpsi[index[i]].re,vpsi[index[i]].im,
	      k_x[index[i]],k_y[index[i]],k_z[index[i]]);
      fclose(fp);
    }
#ifdef ZLIB_H
   else if(ibinary_write_opt==2){
     strcat(vpsiName,".gz");
     gzFile zfp=gzopen(vpsiName,"w");
      gzprintf(zfp,"%d %d %d %d\n",ncoef_true,sizeX,sizeY,sizeZ);
      for(int i=istrt+1;i<ncoef;i++){
        gzprintf(zfp,"%g %g %d %d %d \n",vpsi[index[i]].re,vpsi[index[i]].im,
                   k_x[index[i]],k_y[index[i]],k_z[index[i]]);
      }//endfor
      int i = istrt;
      gzprintf(zfp,"%g %g %d %d %d \n",vpsi[index[i]].re,vpsi[index[i]].im,
		k_x[index[i]],k_y[index[i]],k_z[index[i]]);
      gzclose(zfp);

   }
   else if(ibinary_write_opt==3){
     strcat(vpsiName,".gz");
     gzFile zfp=gzopen(vpsiName,"w");
     int n=1;
     gzwrite(zfp,&ncoef_true,sizeof(int));
     gzwrite(zfp,&sizeX,sizeof(int));
     gzwrite(zfp,&sizeY,sizeof(int));
     gzwrite(zfp,&sizeZ,sizeof(int));
     for(int i=istrt+1;i<ncoef;i++){
       gzwrite(zfp,&vpsi[index[i]].re,sizeof(double));
       gzwrite(zfp,&vpsi[index[i]].im,sizeof(double));
       gzwrite(zfp,&k_x[index[i]],sizeof(int));
       gzwrite(zfp,&k_y[index[i]],sizeof(int));
       gzwrite(zfp,&k_z[index[i]],sizeof(int));
     }//endfor
     int i = istrt;
     gzwrite(zfp,&vpsi[index[i]].re,sizeof(double));
     gzwrite(zfp,&vpsi[index[i]].im,sizeof(double));
     gzwrite(zfp,&k_x[index[i]],sizeof(int));
     gzwrite(zfp,&k_y[index[i]],sizeof(int));
     gzwrite(zfp,&k_z[index[i]],sizeof(int));
     gzclose(zfp);
   }//endif
#endif
else{
      FILE *fp  = fopen(vpsiName,"w");
      int n=1;
      fwrite(&ncoef_true,sizeof(int),n,fp);
      fwrite(&sizeX,sizeof(int),n,fp);
      fwrite(&sizeY,sizeof(int),n,fp);
      fwrite(&sizeZ,sizeof(int),n,fp);
      for(int i=istrt+1;i<ncoef;i++){
        fwrite(&vpsi[index[i]].re,sizeof(double),n,fp);
        fwrite(&vpsi[index[i]].im,sizeof(double),n,fp);
        fwrite(&k_x[index[i]],sizeof(int),n,fp);
        fwrite(&k_y[index[i]],sizeof(int),n,fp);
        fwrite(&k_z[index[i]],sizeof(int),n,fp);
      }//endfor
        int i = istrt;
        fwrite(&vpsi[index[i]].re,sizeof(double),n,fp);
        fwrite(&vpsi[index[i]].im,sizeof(double),n,fp);
        fwrite(&k_x[index[i]],sizeof(int),n,fp);
        fwrite(&k_y[index[i]],sizeof(int),n,fp);
        fwrite(&k_z[index[i]],sizeof(int),n,fp);
	fclose(fp);
    }//endif
  }//endif

  delete [] index;
  delete [] ktemp;

//============================================================================
  }//end routine
//============================================================================


//=============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=============================================================================
void sort_psi_output(int n,int *kx,int *ky,int *kz,int *index,int *ktemp,int *istrt_ret){
//=============================================================================
// Find max and min

  int kxmax = kx[0]; int kxmin=kx[0];
  int kymax = ky[0]; int kymin=ky[0];
  int kzmax = kz[0]; int kzmin=kz[0];
  for(int i=0;i<n;i++){
    kxmax = MAX(kxmax,kx[i]); kxmin = MIN(kxmin,kx[i]);
    kymax = MAX(kymax,ky[i]); kymin = MIN(kymin,ky[i]);
    kzmax = MAX(kzmax,kz[i]); kzmin = MIN(kzmin,kz[i]);
  }//endfor
  int nx = kxmax-kxmin+1;
  int ny = kymax-kymin+1;
  int nz = kzmax-kzmin+1;

//=============================================================================
// Create a sortable 1d array and sort


  for(int i=0;i<n;i++){index[i]=i;}
  for(int i=0;i<n;i++){ktemp[i]=(kz[i]-kzmin)+(ky[i]-kymin)*nz+(kx[i]-kxmin)*nz*ny;}
  sort_commence(n,ktemp,index);

//=============================================================================
// Find g=0 term

  int istrt=0;
  for(int i=0;i<n;i++){
    if(kx[index[i]]==0&&ky[index[i]]==0&&kz[index[i]]==0){
      istrt=i; break;
    }//endif
  }//endfor
  (*istrt_ret) = istrt;
   
//-----------------------------------------------------------------------------
  }// end routine
//=============================================================================


//=============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=============================================================================
void sort_kxky(int n,int *kx,int *ky,int *index,int *ktemp,int sizeY){
//=============================================================================
// Sort on kx

  for(int i=0;i<n;i++){index[i]=i;}
  sort_commence(n,kx,index);
  for(int i=0;i<n;i++){ktemp[i]=ky[i];}
  for(int i=0;i<n;i++){ky[i]=ktemp[index[i]];}

  if(kx[0]<0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Double pack only in sort kxky\n");  
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

//=============================================================================
// Sort on ky for each kx

  int kmin = kx[0];
  int kmax = kx[(n-1)];
  for(int i=kmin;i<=kmax;i++){ktemp[i]=0;}
  for(int i=0;i<n;i++){ktemp[kx[i]]++;}

  int ioff=0;
  for(int i=kmin;i<=kmax;i++){
    if(ktemp[i]>1){sort_commence(ktemp[i],&ky[ioff],&index[ioff]);}
    ioff += ktemp[i];
  }//endfor
   
//============================================================================
  }//end routine
//============================================================================

//=============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=============================================================================
void sort_kxky_old(int n,int *kx,int *ky,int *index,int *kyt){
//=============================================================================

 int nk0=0;
 for(int i=0;i<n;i++){index[i]=i;}
 for(int i=0;i<n;i++){
   if(kx[i]==0){
     int itmp   = index[nk0];
     index[nk0] = index[i];
     index[i]   = itmp;
     nk0++;
   }//endif
 }//endfor

 for(int i=0;i<nk0;i++){kyt[i]=ky[index[i]];}
 sort_commence(nk0,kyt,index);

//============================================================================
  }//end routine
//============================================================================
