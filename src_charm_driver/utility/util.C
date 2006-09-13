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
#ifdef CMK_PROJECTIONS_USE_ZLIB
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
   int rhoGHelpers = config.rhoGHelpers;
   int sizeXEext = sim->ngrid_eext_a;
   int sizeX     = sim->sizeX;
   int sizeY     = sim->sizeY;
   int sizeZ     = sim->sizeZ;
   double *hmati = gencell->hmati;
   double ecut4  = 8.0*cpcoeffs_info->ecut; // convert to Ryd extra factor of 2.0

   get_rho_kvectors(ecut4,hmati,&kx,&ky,&kz,&nline_tot,&nPacked,sizeX,sizeY,sizeZ);

//===================================================================================
// Reorder the kvectors to produce better balance for the lines : 

    int *kx_ind      = new int[nline_tot];
    int *kx_line     = new int[nline_tot];
    int *ky_line     = new int[nline_tot];
    int *kx_line_ext = new int[nline_tot];
    int *ky_line_ext = new int[nline_tot];
    int *istrt_line  = new int [nline_tot];
    int *iend_line   = new int [nline_tot];
    int *npts_line   = new int [nline_tot];

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

    double temp     = ((double)nplane_x)*config.gExpandFactRho;
    int nchareRhoG  = (int)temp;
    int nx          = sim->sizeX;
    int ny          = sim->sizeY;
    int nz          = sim->sizeZ;
    int nx_ext      = sim->ngrid_eext_a;
    int ny_ext      = sim->ngrid_eext_b;
    int *kxt        = new int[nPacked];
    int *kyt        = new int[nPacked];
    int *kzt        = new int[nPacked];
    memcpy(kxt,kx,(nPacked*sizeof(int)));
    memcpy(kyt,ky,(nPacked*sizeof(int)));
    memcpy(kzt,kz,(nPacked*sizeof(int)));
    int nsplit = (3*nplane_x)/2;


    int jc      = 0;
    int lc      = 0;
    for(int i=0;i<nsplit; i++){
      for(int j=i;j<nline_tot;j+=nsplit){
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
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
    if(lc!=nline_tot){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Toasty Line Flip-lines!\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    delete [] kxt;
    delete [] kyt;
    delete [] kzt;

    ic              = 0;
    istrt_line[0]   = 0;
    kx_line[ic]     = kx[ic];
    ky_line[ic]     = ky[ic];
    kx_line_ext[ic] = kx[ic];
    ky_line_ext[ic] = ky[ic];
    if(kx_line[ic]    <0){kx_line[ic]    +=nx;}
    if(ky_line[ic]    <0){ky_line[ic]    +=ny;}
    if(kx_line_ext[ic]<0){kx_line_ext[ic]+=nx_ext;}
    if(ky_line_ext[ic]<0){ky_line_ext[ic]+=ny_ext;}
    for(int i = 1;i<nPacked;i++){
      if(kx[i]!=kx[(i-1)] || ky[i]!=ky[(i-1)]){
        iend_line[ic]   = i;
        npts_line[ic]   = iend_line[ic]-istrt_line[ic];
        ic++;
        istrt_line[ic]  = i;
        kx_line[ic]     = kx[i];
        ky_line[ic]     = ky[i];
        kx_line_ext[ic] = kx[i];
        ky_line_ext[ic] = ky[i];
        if(kx_line[ic]    <0){kx_line[ic]    +=nx;}
        if(ky_line[ic]    <0){ky_line[ic]    +=ny;}
        if(kx_line_ext[ic]<0){kx_line_ext[ic]+=nx_ext;}
        if(ky_line_ext[ic]<0){ky_line_ext[ic]+=ny_ext;}
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
        runs.push_back(RunDescriptor(curr_x,curr_y,curr_z,run_length_sum,run_length,1,nz));
        nrun_tot      +=1;
        run_length_sum += run_length;
        curr_x          = x;
        curr_y          = y;
        curr_z          = z;
        tmpz            = z;
        run_length      = 1;
        if(kz[pNo]==0 && kz[(pNo-1)]>=0){ // test the new line to see if it is of length 1
          runs.push_back(RunDescriptor(curr_x,curr_y,curr_z,run_length_sum,0,1,nz));
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
    runs.push_back(RunDescriptor(curr_x,curr_y,curr_z,run_length_sum,run_length,1,nz));
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
     //CkPrintf("i %d kx %d kx1 %d ky %d ky1 %d\n",i,Desi->x,Desi1->x,Desi->y,Desi1->y);
      if( (Desi->x != Desi1->x) || (Desi->y != Desi1->y) ){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("The rho rundescriptor MUST pair up the half-lines\n");
        CkPrintf("i %d kx %d kx1 %d ky %d ky1 %d\n",i,Desi->x,Desi1->x,Desi->y,Desi1->y);
	CkPrintf("or you will not be a happy camper :\n");
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit(); 
      }//endfor 
      if(Desi1->z!=0 || Desi1->length == 0){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("The rho rundescriptor MUST have 2nd z ==0 and len>0\n");
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit(); 
      }
      if(Desi->z==0 && Desi->length != 0){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("The rho rundescriptor with 1st z == 0 must have 0 lngth\n");
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit(); 
      }
    }//endfor

//===================================================================================
// Decompose lines to balance points

    int nchareRhoGEext   = nchareRhoG*rhoGHelpers;
    int *istrt_lgrp      = new int [nchareRhoG];
    int *iend_lgrp       = new int [nchareRhoG];
    int *npts_lgrp       = new int [nchareRhoG];
    int *nline_lgrp      = new int [nchareRhoG];
    int *nline_lgrp_eext = new int [nchareRhoGEext];

    ParaGrpParse::get_chareG_line_prms(nPacked,nchareRhoG,nline_tot,npts_line,
                               istrt_lgrp,iend_lgrp,npts_lgrp,nline_lgrp,false);

//============================================================================
// Create the line decomposition and a sorted run descriptor
// There are two rundescriptors per line : Noah's arc sort

    int nlines_min=nline_lgrp[0];
    int nlines_max=0;
    for(int i=0;i<nchareRhoG;i++){
      nlines_max=MAX(nlines_max,nline_lgrp[i]);
      nlines_min=MIN(nlines_min,nline_lgrp[i]);
    }//endfor
    int **index_tran_upack_rho  = cmall_int_mat(0,nchareRhoG,0,nlines_max,"util.C");
    int **index_tran_upack_eext = cmall_int_mat(0,nchareRhoGEext,0,nlines_max,"util.C");
   
    int yspace=sizeX/2+1;
    for(int igrp=0;igrp<nchareRhoG;igrp++){
      for(int i=istrt_lgrp[igrp],j=0;i<iend_lgrp[igrp];i++,j++){
        index_tran_upack_rho[igrp][j] = kx_line[i] + ky_line[i]*yspace;
      }//endfor
    }//endfor

    int yspaceEext=sizeXEext/2+1;
    for(int igrp=0,jgrp=0;igrp<nchareRhoG;igrp++){
      int nlTot = nline_lgrp[igrp];
      int istrt = istrt_lgrp[igrp];
      for(int k=0;k<rhoGHelpers;k++,jgrp++){
        int kstrt,kend,nl;
        getSplitDecomp(&kstrt,&kend,&nl,nlTot,rhoGHelpers,k);
        kstrt += istrt;
        nline_lgrp_eext[jgrp] = nl;
        for(int j=0,i=kstrt;j<nl;j++,i++){
          index_tran_upack_eext[jgrp][j] = kx_line_ext[i] + ky_line_ext[i]*yspaceEext;
        }//endfor
      }//endfor
    }//endfor

    CkVec<RunDescriptor> *RhosortedRunDescriptors;
    RhosortedRunDescriptors = new CkVec<RunDescriptor> [nchareRhoG];
    for(int igrp = 0; igrp < nchareRhoG; igrp++){
      for(int i=istrt_lgrp[igrp];i<iend_lgrp[igrp];i++){
 	 int j  = 2*i;
 	 int j1 = 2*i+1;
         RhosortedRunDescriptors[igrp].push_back(runs[j]);
         RhosortedRunDescriptors[igrp].push_back(runs[j1]);
      }//endfor
    }//endfor

  for(int x = 0; x < nchareRhoG; x ++) {
      int runsToBeSent = RhosortedRunDescriptors[x].size();
      int numPoints    = 0;
      for (int j = 0; j < RhosortedRunDescriptors[x].size(); j++){
        numPoints += RhosortedRunDescriptors[x][j].length;
      }//endfor
  }//endfor

//============================================================================
// variables that could be used for mapping but aren't yet.

  double *pts_per_chare = new double[nchareRhoG];
  double *lines_per_chare = new double[nchareRhoG];
  for(int i=0;i<nchareRhoG;i++){
      pts_per_chare[i]  =(double) npts_lgrp[i];
      lines_per_chare[i]=(double) nline_lgrp[i];
  }//endfor

//============================================================================
// Check for rhoghelper size

  if(nlines_min==0){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("No lines in a RhoG collection. Your RhoG decomp stinks.\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit(); 
  }//endif

  if(nlines_min<config.rhoGHelpers){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("RhoGHelper parameter, %d, must be <= minimum number\n",config.rhoGHelpers);
     CkPrintf("of lines in any RhoG chare array element, %d.\n",nlines_min);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit(); 
  }//endif

//============================================================================
// Pack up the stuff

    sim->nplane_rho_x            = nplane_x;
    sim->nchareRhoG              = nchareRhoG;
    sim->nchareRhoGEext          = nchareRhoGEext;
    sim->npts_per_chareRhoG      = npts_lgrp;
    sim->index_tran_upack_rho    = index_tran_upack_rho;
    sim->index_tran_upack_eext   = index_tran_upack_eext;
    sim->nlines_max_rho          = nlines_max;
    sim->nlines_per_chareRhoG    = nline_lgrp;
    sim->nlines_per_chareRhoGEext= nline_lgrp_eext;
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
    delete [] kx_line_ext;
    delete [] ky_line_ext;
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
        runs.push_back(RunDescriptor(curr_x,curr_y,curr_z,run_length_sum,run_length,1,nz));
        nrun_tot      +=1;
        run_length_sum += run_length;
        curr_x          = x;
        curr_y          = y;
        curr_z          = z;
        tmpz            = z;
        run_length      = 1;
        if(kz[pNo]==0 && kz[(pNo-1)]>=0){// test the new line to see if it of length 0
          runs.push_back(RunDescriptor(curr_x,curr_y,curr_z,run_length_sum,0,1,nz));
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
    runs.push_back(RunDescriptor(curr_x,curr_y,curr_z,run_length_sum,run_length,1,nz));
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
      if(Desi->z==0 && Desi->length != 0){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("The rundescriptor with 1st z == 0 must have 0 lngth\n");
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit(); 
      }
      if(Desi1->z!=0 || Desi1->length == 0){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("The rundescriptor MUST have 2nd z == 0\n");
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit(); 
      }
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

#ifndef CMK_PROJECTIONS_USE_ZLIB
    if(ibinary_opt>1)
      {
	CkPrintf("Attempt to use ZLIB Failed! Please review compilation\n");
	CkPrintf("Macro cmk-projections-use-zlib  is %d \n", CMK_PROJECTIONS_USE_ZLIB);
	CkExit();
      }
#endif
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
#ifdef CMK_PROJECTIONS_USE_ZLIB
    else if(ibinary_opt==2){
      //      CkPrintf("Using ZLIB to load ascii states\n");
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
	//	CkPrintf("Using ZLIB to load binary states\n");
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
#ifndef CMK_PROJECTIONS_USE_ZLIB
    if(ibinary_opt>1)
      {
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	CkPrintf("Attempt to use ZLIB Failed! Please review compilation\n");
	CkPrintf("Macro cmk-projections-use-zlib  is %d \n", CMK_PROJECTIONS_USE_ZLIB);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	CkExit();
      }
#endif

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
#ifdef CMK_PROJECTIONS_USE_ZLIB
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
     fprintf(fp,"stateOutputOn: %d\n",stateOutputOn);
     fprintf(fp,"atmOutputOn: %d\n",atmOutputOn);
     fprintf(fp,"sGrainSize: %d\n",sGrainSize);
     fprintf(fp,"orthoStride: %d\n",orthoStride);
     fprintf(fp,"orthoGrainSize: %d\n",orthoGrainSize);
     fprintf(fp,"launchNLeesFromRho: %d\n",launchNLeesFromRho);
     fprintf(fp,"useBWBarrier: %d\n",useBWBarrier);
     fprintf(fp,"useOrthoDirect: %d\n",useOrthoDirect);
     fprintf(fp,"useOrthoSection: %d\n",useOrthoSection);
     fprintf(fp,"useOrthoSectionRed: %d\n",useOrthoSectionRed);
     fprintf(fp,"useOrthoHelpers: %d\n",useOrthoHelpers);
     fprintf(fp,"lambdaGrainSize: %d\n",lambdaGrainSize);
     fprintf(fp,"pesPerState: %d\n",pesPerState);
     fprintf(fp,"gExpandFact: %g\n",gExpandFact);
     fprintf(fp,"gExpandFactRho: %g\n",gExpandFactRho);
     fprintf(fp,"rhoGHelpers: %d\n",rhoGHelpers);
     fprintf(fp,"numSfGrps: %d\n",numSfGrps);
     fprintf(fp,"numSfDups: %d\n",numSfDups);
     fprintf(fp,"toleranceInterval: %d\n",toleranceInterval);
     fprintf(fp,"localAtomBarrier: %d\n",localAtomBarrier);
     fprintf(fp,"localEnergyBarrier: %d\n",localEnergyBarrier);
     fprintf(fp,"lbpaircalc: %d\n",lbpaircalc);
     fprintf(fp,"lbgspace: %d\n",lbgspace);
     fprintf(fp,"lbdensity: %d\n",lbdensity);
     fprintf(fp,"useCommlib: %d\n",useCommlib);
     fprintf(fp,"usePairEtoM: %d\n",usePairEtoM);
     fprintf(fp,"usePairDirectSend: %d\n",usePairDirectSend);
     fprintf(fp,"PCCollectTiles: %d\n",PCCollectTiles);
     fprintf(fp,"PCdelayBWSend: %d\n",PCdelayBWSend);
     fprintf(fp,"PCstreamBWout: %d\n",PCstreamBWout);
     fprintf(fp,"PCstreamFWblock: %d\n",PCstreamFWblock);
     fprintf(fp,"useGHartInsRhoRP: %d\n", useGHartInsRhoRP);
     fprintf(fp,"useGIns0RhoRP: %d\n", useGIns0RhoRP);
     fprintf(fp,"useGIns1RhoRP: %d\n", useGIns1RhoRP);
     fprintf(fp,"useGIns2RhoRP: %d\n", useGIns2RhoRP);
     fprintf(fp,"useGIns3RhoRP: %d\n", useGIns3RhoRP);
     fprintf(fp,"useGByrdInsRhoRBP: %d\n", useGByrdInsRhoRBP);
     fprintf(fp,"useRInsRhoGP: %d\n", useRInsRhoGP);
     fprintf(fp,"useRInsIGXRhoGP: %d\n", useRInsIGXRhoGP);
     fprintf(fp,"useRInsIGYRhoGP: %d\n", useRInsIGYRhoGP);
     fprintf(fp,"useRInsIGZRhoGP: %d\n", useRInsIGZRhoGP);

     fprintf(fp,"useGssInsRealP: %d\n", useGssInsRealP);
     fprintf(fp,"useGssInsRealP: %d\n", useGssInsRealPP);
     fprintf(fp,"useMssInsGP: %d\n", useMssInsGP);
     fprintf(fp,"useMssInsGPP: %d\n", useMssInsGPP);
     fprintf(fp,"useGHartInsRHart %d\n",useGHartInsRHart);
     fprintf(fp,"useRHartInsGHart %d\n",useRHartInsGHart);

     fprintf(fp,"useGMulticast: %d\n",useGMulticast);
     fprintf(fp,"useCommlibMulticast: %d\n",useCommlibMulticast);
     fprintf(fp,"numMulticastMsgs: %d\n",numMulticastMsgs);
     fprintf(fp,"PCSpanFactor: %d\n",PCSpanFactor);
     fprintf(fp,"OrthoRedSpanFactor: %d\n",OrthoRedSpanFactor);
     fprintf(fp,"OrthoMcastSpanFactor: %d\n",OrthoMcastSpanFactor);

     fprintf(fp,"fftprogresssplit: %d\n",fftprogresssplit);
     fprintf(fp,"fftprogresssplitReal: %d\n",fftprogresssplitReal);
     fprintf(fp,"RpesPerState: %d\n",RpesPerState);
     fprintf(fp,"GpesPerState: %d\n",GpesPerState);

     fprintf(fp,"prioFFTMsg %d\n",prioFFTMsg);
     fprintf(fp,"rsfftpriority: %d\n",rsfftpriority);
     fprintf(fp,"gsfftpriority: %d\n",gsfftpriority);
     fprintf(fp,"rsifftpriority: %d\n",rsifftpriority);
     fprintf(fp,"gsifftpriority: %d\n",gsifftpriority);
     fprintf(fp,"rhorpriority: %d\n",rhorpriority);
     fprintf(fp,"rhogpriority: %d\n",rhogpriority);

     fprintf(fp,"sfpriority: %d\n",sfpriority);
     fprintf(fp,"lambdapriority: %d\n",lambdapriority);
     fprintf(fp,"psipriority: %d\n",psipriority);

     fprintf(fp,"prioNLFFTMsg %d\n",prioNLFFTMsg);
     fprintf(fp,"prioEextFFTMsg %d\n",prioEextFFTMsg);
     fprintf(fp,"gsNLfftpriority: %d\n",gsNLfftpriority);
     fprintf(fp,"rsNLfftpriority: %d\n",rsNLfftpriority);
     fprintf(fp,"rhorHartpriority %d\n", rhorHartpriority);
     fprintf(fp,"rhogHartpriority %d\n", rhogHartpriority);

     fprintf(fp,"doublePack: %d\n", doublePack);
     fprintf(fp, "conserveMemory: %d\n", conserveMemory);
     fprintf(fp, "Gstates_per_pe: %d\n", Gstates_per_pe);
     fprintf(fp, "scalc_per_plane: %d\n", scalc_per_plane);
     fprintf(fp, "Rstates_per_pe: %d\n", Rstates_per_pe);
     fprintf(fp,"numChunks: %d\n",numChunks);
     fprintf(fp,"gSpaceSum: %d\n",gSpaceSum);
     fprintf(fp,"phantomSym: %d\n",phantomSym);
     fprintf(fp,"prioBW: %d\n",prioBW);
     fprintf(fp,"numChunksSym: %d\n",numChunksSym);
     fprintf(fp,"numChunksAsym: %d\n",numChunksAsym);
     fprintf(fp,"gStreamPeriod: %g\n",gStreamPeriod);
     fprintf(fp,"rStreamPeriod: %g\n",rStreamPeriod);
     fprintf(fp,"gBucketSize: %d\n",gBucketSize);
     fprintf(fp,"rBucketSize: %d\n",rBucketSize);
     fprintf(fp,"useCuboidMap: %d\n",useCuboidMap);
     fprintf(fp,"useCuboidMapRS: %d\n",useCuboidMapRS);
     fprintf(fp,"useCentroidMap: %d\n",useCentroidMap);
     fprintf(fp,"useCentroidMapRho: %d\n",useCentroidMapRho);

     //     fprintf(fp,"nstates: %d\n",nstates);
     //     fprintf(fp,"nchareG %d\n",nchareG);
     //     fprintf(fp,"nchareRhoG %d\n",nchareRhoG);
     //     fprintf(fp,"maxIter: %d\n",maxIter);
     //     fprintf(fp,"low_x_size: %d\n",low_x_size);
     //     fprintf(fp,"high_x_size: %d\n",high_x_size);
     //     fprintf(fp,"numData %d\n",numData);
     //     fprintf(fp,"numFFTPoints %d\n",numFFTPoints);
   fclose(fp);


//----------------------------------------------------------------------------------
   }//end routine
//===================================================================================




//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::readConfig(const char* fileName, Config &config,
			int nstates_in, int nkf1, int nkf2, int nkf3, int maxIter_in,
			int ibinary_opt,int natm_nl, int ees_nonlocal_on)
//===================================================================================
    { // begin routine
//===================================================================================
// Initialize parameters

    config.nstates         = nstates_in;
    config.maxIter         = maxIter_in;

    config.sGrainSize           = nstates_in;
    config.orthoGrainSize           =     config.sGrainSize;
    config.orthoStride           = 0;
    config.useBWBarrier           = 0;
    config.launchNLeesFromRho           = 0;
    config.useOrthoSection       = 0;
    config.useOrthoDirect       = 0;
    config.useOrthoSectionRed    = 0;
    config.lambdaGrainSize       =     config.sGrainSize;
    config.useOrthoHelpers       = 0;
    config.rhoGHelpers          = 1;
    config.pesPerState          = 1;
    config.RpesPerState         = 0; 
    config.GpesPerState         = 0; 
    config.doublePack           = 1;
    config.conserveMemory       = 0;
    config.numMulticastMsgs     = 10;
    config.PCSpanFactor         = 2;
    config.OrthoRedSpanFactor   = 2;
    config.OrthoMcastSpanFactor = 16;
    config.useGMulticast        = 0;
    config.useCommlibMulticast  = 1;
    config.useCommlib           = 1;
    config.usePairEtoM           = 0;
    config.usePairDirectSend     = 0;
    config.useCuboidMap           = 0;
    config.useCuboidMapRS           = 0;
    config.useCentroidMap        = 0;
    config.useCentroidMapRho        = 0;
    config.PCCollectTiles       = 1;
    config.PCstreamBWout       = 0;
    config.PCstreamFWblock       = 0;
    config.PCdelayBWSend       = 1;

    // Density FFT comlib flags
    config.useGHartInsRhoRP	= config.useCommlib;
    config.useGIns0RhoRP	= config.useCommlib;
    config.useGIns1RhoRP	= config.useCommlib;
    config.useGIns2RhoRP	= config.useCommlib;
    config.useGIns3RhoRP	= config.useCommlib;
    config.useGByrdInsRhoRBP	= config.useCommlib;
    config.useRInsRhoGP		= config.useCommlib;
    config.useRInsIGXRhoGP	= config.useCommlib;
    config.useRInsIGYRhoGP	= config.useCommlib;
    config.useRInsIGZRhoGP	= config.useCommlib;

    // state real and state g FFT comlib flags
    config.useGssInsRealP	= config.useCommlib;
    config.useMssInsGP		= config.useCommlib;

    // Ees methods for NLPP and EESEext FFT comblib flags
    config.useGssInsRealPP	= config.useCommlib;
    config.useMssInsGPP		= config.useCommlib;
    config.useGHartInsRHart	= config.useCommlib;
    config.useRHartInsGHart	= config.useCommlib;

    config.lbpaircalc           = 0;
    config.lbgspace             = 0;
    config.lbdensity            = 0;
    config.numSfGrps            = 1;
    config.numSfDups            = 1;

    // density and state fft prios
    config.prioFFTMsg           = 1; 
    config.rsfftpriority        = 1000000;
    config.gsfftpriority        = 1000000;
    config.rsifftpriority       = 100000000;
    config.gsifftpriority       = 200000000;
    config.rhorpriority         = 2000000;
    config.rhogpriority         = 2000000; 

    // PC and SF prios
    config.sfpriority           = 10000000;
    config.lambdapriority       = 300000000;
    config.psipriority          = 400000000;

    // ees method prios
    config.prioNLFFTMsg         = 1; 
    config.prioEextFFTMsg       = 1; 
    config.rsNLfftpriority      = 2300000;
    config.gsNLfftpriority      = 2500000;
    config.rhorHartpriority     = 2000000;
    config.rhogHartpriority     = 2000000;

    config.gExpandFact          = 1.0;
    config.gExpandFactRho       = 1.0;
    config.fftprogresssplit     = 20;
    config.fftprogresssplitReal = 5;
    config.stateOutputOn        = 0;
    config.atmOutputOn          = 0;
    config.localAtomBarrier     = 1;
    config.localEnergyBarrier   = 1;
    config.toleranceInterval    = 1;
    config.Gstates_per_pe	= config.nstates;
    config.scalc_per_plane	= 1;
    config.Rstates_per_pe	= config.nstates;
    config.gStreamPeriod         = 2.0;
    config.rStreamPeriod         = 2.0;
    config.gBucketSize         = 5;
    config.rBucketSize         = 5;
    config.numChunks            = 1;
    config.gSpaceSum             = 0;
    config.phantomSym             = 0;
    config.numChunksSym            = 1;
    config.numChunksAsym            = 1;
    strcpy(config.dataPath,"./");

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
        if (!strcmp(parameterName, "dataPath"))
            strcpy(config.dataPath, parameterValue);
        else if (!strcmp(parameterName, "sGrainSize"))
            config.sGrainSize = atoi(parameterValue);
        else if (!strcmp(parameterName, "orthoGrainSize"))
            config.orthoGrainSize = atoi(parameterValue);
        else if (!strcmp(parameterName, "orthoStride"))
            config.orthoStride = atoi(parameterValue);
        else if (!strcmp(parameterName, "useBWBarrier"))
            config.useBWBarrier = atoi(parameterValue);
        else if (!strcmp(parameterName, "launchNLeesFromRho"))
            config.launchNLeesFromRho = atoi(parameterValue);
        else if (!strcmp(parameterName, "useOrthoDirect"))
            config.useOrthoDirect = atoi(parameterValue);
        else if (!strcmp(parameterName, "useOrthoHelpers"))
            config.useOrthoHelpers = atoi(parameterValue);
        else if (!strcmp(parameterName, "useOrthoSection"))
            config.useOrthoSection = atoi(parameterValue);
        else if (!strcmp(parameterName, "useOrthoSectionRed"))
            config.useOrthoSectionRed = atoi(parameterValue);
        else if (!strcmp(parameterName, "lambdaGrainSize"))
            config.lambdaGrainSize = atoi(parameterValue);
        else if (!strcmp(parameterName, "useCommlib"))
            config.useCommlib = atoi(parameterValue);
        else if (!strcmp(parameterName, "useCuboidMap"))
            config.useCuboidMap = atoi(parameterValue);
        else if (!strcmp(parameterName, "useCuboidMapRS"))
            config.useCuboidMapRS = atoi(parameterValue);
        else if (!strcmp(parameterName, "useCentroidMap"))
            config.useCentroidMap = atoi(parameterValue);
        else if (!strcmp(parameterName, "useCentroidMapRho"))
            config.useCentroidMapRho = atoi(parameterValue);
        else if (!strcmp(parameterName, "usePairEtoM"))
            config.usePairEtoM = atoi(parameterValue);
        else if (!strcmp(parameterName, "usePairDirectSend"))
            config.usePairDirectSend = atoi(parameterValue);
        else if (!strcmp(parameterName, "PCCollectTiles"))
            config.PCCollectTiles = atoi(parameterValue);
        else if (!strcmp(parameterName, "PCdelayBWSend"))
            config.PCdelayBWSend = atoi(parameterValue);
        else if (!strcmp(parameterName, "PCstreamBWout"))
            config.PCstreamBWout = atoi(parameterValue);
        else if (!strcmp(parameterName, "PCstreamFWblock"))
            config.PCstreamFWblock = atoi(parameterValue);
        else if (!strcmp(parameterName, "useGHartInsRhoRP"))
            config.useGHartInsRhoRP = atoi(parameterValue);
        else if (!strcmp(parameterName, "useGHartInsRHart"))
            config.useGHartInsRHart = atoi(parameterValue);
        else if (!strcmp(parameterName, "useRHartInsGHart"))
            config.useRHartInsGHart = atoi(parameterValue);
        else if (!strcmp(parameterName, "useGIns0RhoRP"))
            config.useGIns0RhoRP = atoi(parameterValue);
        else if (!strcmp(parameterName, "useGIns1RhoRP"))
            config.useGIns1RhoRP = atoi(parameterValue);
        else if (!strcmp(parameterName, "useGIns2RhoRP"))
            config.useGIns2RhoRP = atoi(parameterValue);
        else if (!strcmp(parameterName, "useGIns3RhoRP"))
            config.useGIns3RhoRP = atoi(parameterValue);
        else if (!strcmp(parameterName, "useGByrdInsRhoRBP"))
            config.useGByrdInsRhoRBP = atoi(parameterValue);
        else if (!strcmp(parameterName, "useRInsRhoGP"))
            config.useRInsRhoGP = atoi(parameterValue);
        else if (!strcmp(parameterName, "useRInsIGXRhoGP"))
            config.useRInsIGXRhoGP = atoi(parameterValue);
        else if (!strcmp(parameterName, "useRInsIGYRhoGP"))
            config.useRInsIGYRhoGP = atoi(parameterValue);
        else if (!strcmp(parameterName, "useRInsIGZRhoGP"))
            config.useRInsIGZRhoGP = atoi(parameterValue);
        else if (!strcmp(parameterName, "useGssInsRealP"))
            config.useGssInsRealP = atoi(parameterValue);
        else if (!strcmp(parameterName, "useGssInsRealPP"))
            config.useGssInsRealPP = atoi(parameterValue);
        else if (!strcmp(parameterName, "useMssInsGP"))
            config.useMssInsGP = atoi(parameterValue);
        else if (!strcmp(parameterName, "useMssInsGPP"))
            config.useMssInsGPP = atoi(parameterValue);
        else if (!strcmp(parameterName, "useGMulticast"))
            config.useGMulticast = atoi(parameterValue);
        else if (!strcmp(parameterName, "numMulticastMsgs"))
            config.numMulticastMsgs = atoi(parameterValue);
        else if (!strcmp(parameterName, "PCSpanFactor"))
            config.PCSpanFactor = atoi(parameterValue);
        else if (!strcmp(parameterName, "OrthoRedSpanFactor"))
            config.OrthoRedSpanFactor = atoi(parameterValue);
        else if (!strcmp(parameterName, "OrthoMcastSpanFactor"))
            config.OrthoMcastSpanFactor = atoi(parameterValue);
        else if (!strcmp(parameterName, "useCommlibMulticast"))
            config.useCommlibMulticast = atoi(parameterValue);
        else if (!strcmp(parameterName, "rhoGHelpers"))
            config.rhoGHelpers = atoi(parameterValue);
        else if (!strcmp(parameterName, "pesPerState"))
            config.pesPerState = atoi(parameterValue);
        else if (!strcmp(parameterName, "RpesPerState"))
            config.RpesPerState = atoi(parameterValue);
        else if (!strcmp(parameterName, "GpesPerState"))
            config.GpesPerState = atoi(parameterValue);
	else if (!strcmp(parameterName, "doublePack"))
            config.doublePack = atoi(parameterValue);
	else if (!strcmp(parameterName, "conserveMemory"))
	    config.conserveMemory = atoi(parameterValue);
	else if (!strcmp(parameterName, "lbgspace"))
	    config.lbgspace = atoi(parameterValue);
	else if (!strcmp(parameterName, "phantomSym"))
	    config.phantomSym = atoi(parameterValue);
	else if (!strcmp(parameterName, "lbpaircalc"))
	    config.lbpaircalc = atoi(parameterValue);
	else if (!strcmp(parameterName, "lbdensity"))
	    config.lbdensity = atoi(parameterValue);
        else if (!strcmp(parameterName, "numSfGrps"))
            config.numSfGrps = atoi(parameterValue);
        else if (!strcmp(parameterName, "numSfDups"))
            config.numSfDups = atoi(parameterValue);

	else if (!strcmp(parameterName, "prioFFTMsg"))
	    config.prioFFTMsg = atoi(parameterValue);
        else if (!strcmp(parameterName, "rsfftpriority"))
            config.rsfftpriority = atoi(parameterValue);
        else if (!strcmp(parameterName, "gsfftpriority"))
            config.gsfftpriority = atoi(parameterValue);
        else if (!strcmp(parameterName, "rsifftpriority"))
            config.rsifftpriority = atoi(parameterValue);
        else if (!strcmp(parameterName, "gsifftpriority"))
            config.gsifftpriority = atoi(parameterValue);
        else if (!strcmp(parameterName, "rhorpriority"))
            config.rhorpriority = atoi(parameterValue);
        else if (!strcmp(parameterName, "rhogpriority"))
            config.rhogpriority = atoi(parameterValue);

        else if (!strcmp(parameterName, "sfpriority"))
            config.sfpriority = atoi(parameterValue);
        else if (!strcmp(parameterName, "lambdapriority"))
            config.lambdapriority = atoi(parameterValue);
        else if (!strcmp(parameterName, "psipriority"))
            config.psipriority = atoi(parameterValue);

	else if (!strcmp(parameterName, "prioNLFFTMsg"))
	    config.prioNLFFTMsg = atoi(parameterValue);
	else if (!strcmp(parameterName, "prioEextFFTMsg"))
	    config.prioEextFFTMsg = atoi(parameterValue);
        else if (!strcmp(parameterName, "rhorHartpriority"))
            config.rhorHartpriority = atoi(parameterValue);
        else if (!strcmp(parameterName, "rhogHartpriority"))
            config.rhogHartpriority = atoi(parameterValue);
        else if (!strcmp(parameterName, "rsNLfftpriority"))
            config.rsNLfftpriority = atoi(parameterValue);
        else if (!strcmp(parameterName, "gsNLfftpriority"))
            config.gsNLfftpriority = atoi(parameterValue);


        else if (!strcmp(parameterName, "fftprogresssplit"))
            config.fftprogresssplit = atoi(parameterValue);
        else if (!strcmp(parameterName, "fftprogresssplitReal"))
            config.fftprogresssplitReal = atoi(parameterValue);
        else if (!strcmp(parameterName, "localAtomBarrier"))
            config.localAtomBarrier = atoi(parameterValue);
        else if (!strcmp(parameterName, "localEnergyBarrier"))
            config.localEnergyBarrier = atoi(parameterValue);
        else if (!strcmp(parameterName, "atmOutputOn"))
            config.atmOutputOn = atoi(parameterValue);
        else if (!strcmp(parameterName, "stateOutputOn"))
            config.stateOutputOn = atoi(parameterValue);
        else if (!strcmp(parameterName, "toleranceInterval"))
            config.toleranceInterval = atoi(parameterValue);
	else if (!strcmp(parameterName, "Gstates_per_pe"))
            config.Gstates_per_pe = atoi(parameterValue);
        else if (!strcmp(parameterName, "Rstates_per_pe"))
            config.Rstates_per_pe = atoi(parameterValue);
	else if (!strcmp(parameterName, "gSpaceSum"))
            config.gSpaceSum = atoi(parameterValue);
	else if (!strcmp(parameterName, "numChunks"))
            config.numChunks = atoi(parameterValue);
	else if (!strcmp(parameterName, "numChunksSym"))
            config.numChunksSym = atoi(parameterValue);
	else if (!strcmp(parameterName, "numChunksAsym"))
            config.numChunksAsym = atoi(parameterValue);
	else if (!strcmp(parameterName, "prioBW"))
            config.prioBW = atoi(parameterValue);
        else if (!strcmp(parameterName, "gExpandFact")){
               sscanf(parameterValue,"%lg",&(config.gExpandFact));
               }
        else if (!strcmp(parameterName, "gExpandFactRho")){
               sscanf(parameterValue,"%lg",&(config.gExpandFactRho));
               }
        else if (!strcmp(parameterName, "gStreamPeriod")){
               sscanf(parameterValue,"%lg",&(config.gStreamPeriod));
               }
        else if (!strcmp(parameterName, "rStreamPeriod")){
               sscanf(parameterValue,"%lg",&(config.rStreamPeriod));
               }
	else if (!strcmp(parameterName, "gBucketSize"))
            config.gBucketSize = atoi(parameterValue);
	else if (!strcmp(parameterName, "rBucketSize"))
            config.rBucketSize = atoi(parameterValue);
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

      
//===================================================================================
// Set FFT and g-space size from state file

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
    int nplane_x_rho        = 2*minx+1;

//===================================================================================
//come up with sane values if the user didn't bother to try
    if(config.numChunks>config.numChunksSym)
      config.numChunksSym=config.numChunks;

    if(config.numChunks>config.numChunksAsym)
      config.numChunksAsym=config.numChunks;

    if(config.lambdaGrainSize==config.nstates && config.orthoGrainSize!=config.nstates)
      config.lambdaGrainSize=config.orthoGrainSize;

    config.guesstimateParms(natm_nl, sizez);

    if(config.pesPerState>0 && config.RpesPerState <1){
      config.RpesPerState=config.pesPerState;
    }//endif

    if(config.pesPerState>0 && config.GpesPerState <1){
      config.GpesPerState=config.pesPerState;
    }//endif

//===================================================================================
// Set FFT and g-space size based on config post guesstimate

    double temp         = (config.gExpandFact)*((double)nplane_x);
    int nchareG         = ((int)temp);
    config.nchareG      = nchareG;
    double temp_rho         = (config.gExpandFactRho)*((double)nplane_x_rho);
    int nchareRhoG         = ((int)temp_rho);
    config.nchareRhoG      = nchareRhoG;
    config.scalc_per_plane = (config.nstates/config.sGrainSize)*(config.nstates/config.sGrainSize);

//===================================================================================
// Check the parameter ranges 
    rangeExit(config.launchNLeesFromRho,"launchNLeesFromRho",1);
    rangeExit(config.prioFFTMsg,"prioFFTMsg",1);
    rangeExit(config.stateOutputOn,"stateOutputOn",1);
    rangeExit(config.atmOutputOn,"atmOutputOn",1);
    rangeExit(config.localAtomBarrier,"localAtomBarrier",1);
    rangeExit(config.localEnergyBarrier,"localEnergyBarrier",1);
    rangeExit(config.lbpaircalc,"lbpaircalc",1);
    rangeExit(config.lbgspace,"lbgspace",1);
    rangeExit(config.lbdensity,"lbgspace",1);
    rangeExit(config.useGMulticast,"useGMulticast",1);
    rangeExit(config.useCommlibMulticast,"useCommlibMulticast",1);
    rangeExit(config.useCommlib,"useCommlib",1);
    rangeExit(config.usePairEtoM,"usePairEtoM",1);
    rangeExit(config.usePairDirectSend,"usePairDirectSend",1);
    rangeExit(config.PCCollectTiles,"PCCollectTiles",1);
    rangeExit(config.PCdelayBWSend,"PCdelayBWSend",1);    
    rangeExit(config.PCstreamBWout,"PCstreamBWout",1);
    rangeExit(config.doublePack,"doublePack",1);
    rangeExit(config.conserveMemory,"conserveMemory",1);
    rangeExit(config.fftprogresssplit,"fftprogresssplit",0);
    rangeExit(config.fftprogresssplitReal,"fftprogresssplitReal",0);
    rangeExit(config.rhoGHelpers,"rhoGHelpers",0);
    rangeExit(config.numMulticastMsgs,"numMulticastMsgs",0);
    rangeExit(config.PCSpanFactor,"numMulticastMsgs",0);
    rangeExit(config.pesPerState,"pesPerState;",0);
    rangeExit(config.GpesPerState,"GpesPerState;",0);
    rangeExit(config.RpesPerState,"RpesPerState;",0);
    rangeExit(config.toleranceInterval,"toleranceInterval;",0);
    rangeExit(config.numChunks,"numChunks;",0);
    rangeExit(config.gSpaceSum,"gSpaceSum;",1);
    rangeExit(config.gSpaceSum,"phantomSym;",1);
    rangeExit(config.useCuboidMap,"useCuboidMap;",1);
    rangeExit(config.useCuboidMapRS,"useCuboidMapRS;",1);
    rangeExit(config.useCentroidMap,"useCentroidMap;",1);
    rangeExit(config.useCentroidMapRho,"useCentroidMapRho;",1);
    rangeExit(config.numChunksAsym,"numChunksAsym;",0);
    rangeExit(config.numChunksSym,"numChunksSym;",0);

//===================================================================================
// Consistency Checks on the input
#ifndef CMK_VERSION_BLUEGENE
    if(config.useCuboidMap || config.useCuboidMapRS)
      {
	CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	CkPrintf("useCuboidMap requires CMK_VERSION_BLUEGENE\n");
	CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	CkExit();
      }
#endif
    if(config.gSpaceSum && !config.usePairDirectSend){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("gSpaceSum requires usePairDirectSend\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }

    if(config.gExpandFact<1.0){
      CkPrintf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      CkPrintf("Chare array expansion factor out of range %g\n",config.gExpandFact);
      CkPrintf("This probably could work but I'd check first.\n");
      CkPrintf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      if(config.gExpandFact<=0.0){
       CkExit();
      }
    }//endif


    if(config.gExpandFactRho<1.0){
      CkPrintf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      CkPrintf("RhoChare array expansion factor out of range %g\n",config.gExpandFactRho);
      CkPrintf("This probably could work but I'd check first.\n");
      CkPrintf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      if(config.gExpandFactRho<=0.0){
       CkExit();
      }
    }//endif

    if(nchareG<nplane_x){
      CkPrintf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      CkPrintf("Too few g-space chares %d %d\n",nplane_x,nchareG);
      CkPrintf("This probably could work but I'd check first.\n");
      CkPrintf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    }//endif

    if(nchareRhoG<nplane_x_rho){
      CkPrintf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      CkPrintf("Too few rhog-space chares %d %d\n",nplane_x_rho,nchareRhoG);
      CkPrintf("This probably could work but I'd check first.\n");
      CkPrintf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    }//endif

    if(config.pesPerState>config.nchareG){
      CkPrintf("Warning : pesPerState > %d %g %d\n",config.nchareG,
                 config.gExpandFact,sizex);
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

    if(config.scalc_per_plane<1.0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Scalc per plane cannot be less then 1\n");
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

    if (config.sGrainSize %config.orthoGrainSize != 0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("S matrix grain-size should be divisible by orthoGrainSize\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    if (config.sGrainSize %config.lambdaGrainSize != 0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("S matrix grain-size should be divisible by lambdaGrainSize\n");
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

    if(config.usePairEtoM==1 && config.useCommlib!=1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("EachToMany pairCalc requires Commlib!\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }
    /*
    if(config.useCommlibMulticast!=1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("No commlibMulticast no work. Sameer is happy!\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
    */
    if(config.useCommlibMulticast+config.useGMulticast!=1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("You can't use both the g and commlib multicast\n");
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

    if(config.Gstates_per_pe<1 || config.Gstates_per_pe>config.nstates){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("The number of Gstates per pe must be >=1 < num states\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
    if(config.phantomSym && !config.gSpaceSum)
      {
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Current implementation of phantomSym requires gSpaceSum\n");
      CkPrintf("The price of midnight hacking sessions.\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
      }
    if(config.Rstates_per_pe<1 || config.Rstates_per_pe>config.nstates){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("The number of Rstates per pe must be >=1 < num states\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================


//============================================================================
// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Config::rangeExit(int param, char *name, int iopt){
//============================================================================
  switch(iopt){
   case 0: 
      if(param<1){
          CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          CkPrintf("The parameter %s must be >0 not %d \n",name,param);
          CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          CkExit();
      }//endif
      break;
   case 1: 
      if(param<0 || param>1){
          CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          CkPrintf("The parameter %s must be 1(on) or 0 (off) \n",name);
          CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          CkExit();
      }//endif
      break;
  }//end switch

//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================


//=============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=============================================================================
/**
 * Guesses decent values for configuration parameters based on user
 * set values, system size, and the number of processes available.
 *
 * For large pe runs, scheme is thus: 
 *    1. find power2 nchareG. (typically 32 or 64)
 *    2. determine numchunks and grainsize
 *      a. find numPes()/nchareG.
 *      b. set sGrainsize to nstates, set numchunks to 1
 *      c. set numgrains to (nstates/sGrainSize)^2
 *      
 *    3. if phantoms are in use, numchunks is same for asym and sym
 *    4. set rstates_per_pe and gstates_per_pe as high as they will
 *       go for numPes and this system (constrained by nstates,
 *       nplanes, nchareg
 *    5. if numPes is large enable gSpaceSum, cuboid mapping,
 *       centroid mapping, useDirectSend and probably a few other things
 *
 * Supported parameters: sGrainSize, gExpandFact, gExpandFactRho,
 * fftprogresssplit, fftprogresssplitReal, numSfGrps, numSfDups,
 * pesPerState, rhoGHelpers, numMulticastMsgs

 */
//=============================================================================
void Config::guesstimateParms(int natm_nl, int sizez){
//=============================================================================
    // If the user hasn't set a value in the config, try to come up with
    // something suitable for the system size and number of processes
    // based on what they have set.
  
    // Else accept user set value

    // should come up with an ifdef macro approach so we don't
    // have the extra function call in non BG/L case
//=============================================================================

#ifndef CMK_VERSION_BLUEGENE
    fftprogresssplit=1000;
    fftprogresssplitReal=1000;
#else
    // just use the initial value or user override
#endif  
    int numPes=CkNumPes();  //cache this
    int sqrtpes=(int) sqrt((double)numPes);
    int sqrtstates=(int) sqrt((double)nstates);
    if(gExpandFact==1.0) 
    { 
      //gives us more gspace chares
      //only worth expanding if we have enough pes
      // suitable range is 1.0 to 2.0

      if(numPes>low_x_size)
	{
	  int i=1;
	  double mypow=1;
	  while((mypow=pow(2.0, (double)i)) <= low_x_size)
	    {
	      i++;
	    }
	  //	  CkPrintf("i is %d from low_x_size %d\n");
	  gExpandFact= mypow / (double) low_x_size;
	  nchareG=(int)( gExpandFact * (double) low_x_size);
	}
    }

    if(numChunks==1&&numChunksSym==1&&numChunksAsym==1)
      {
	//	if(numPes>nstates)
	//	  sGrainSize/2;
	int numGrains = nstates/sGrainSize;
	numGrains*=numGrains;
	numChunks=numPes/(nchareG*numGrains);
	if(numChunks<1)
	  numChunks=1;
	while(numChunks>16)
	  {
	    sGrainSize=sGrainSize/2;
	    numGrains = nstates/sGrainSize;
	    numGrains*=numGrains;
	    numChunks=numPes/(nchareG*numGrains);
	  }
	//	if(numPes>=nstates)
	//	  {
	    phantomSym=1;
	    numChunksSym=numChunks;
	    numChunksAsym=numChunks;
	    //	  }
	
      }
    
    if(numPes!=1 && Gstates_per_pe==nstates)
    {

      if(numPes<=128)
	Gstates_per_pe=nstates/4;
      else if(numPes>128 && numPes<=512)
	Gstates_per_pe=nstates/16;
      else 
	Gstates_per_pe= nchareG*nstates/numPes;
      if (Gstates_per_pe==0)
	Gstates_per_pe=1;
    }

    if(numPes!=1 && Rstates_per_pe==nstates)
    {

      if(numPes<=128)
	Rstates_per_pe=nstates/4;
      else 
	Rstates_per_pe= sizez*nstates/numPes;
      if (Rstates_per_pe==0)
	Rstates_per_pe=1;
			 
    }

    if((sGrainSize%orthoGrainSize !=0)|| (sGrainSize==orthoGrainSize))
    {
      // this is lame and should be replace with something which finds
      // an even mod of any sGrainSize
      orthoGrainSize=32;
      if(orthoGrainSize<32)
	orthoGrainSize=32;
    }

    if((sGrainSize%lambdaGrainSize !=0)|| (sGrainSize==lambdaGrainSize))
    {
      // this is lame and should be replace with something which finds
      // an even mod of any (non prime) sGrainSize
      lambdaGrainSize=orthoGrainSize;
      //      if(lambdaGrainSize<32)
      //	lambdaGrainSize=32;

    }


    if(gExpandFactRho==1.0) 
    { //gives us more gspace chares
      //only worth expanding if we have enough pes
      // suitable range is 1.0 to 2.0
	if(numPes>low_x_size*4)
	{
	    gExpandFactRho+=fabs((double) (numPes/2-low_x_size*4)/ (double)( numPes));
	  
	}
	else if(numPes>low_x_size)
	{
	    gExpandFactRho+=(double) sqrtpes/ (double)( numPes);
	}

    }
    if(numSfGrps==1)
    {// number of groups to chop atom calc into
	// needs to grow with size and number of PEs
	// range 1->natm_nl 
	int atmstates=natm_nl*nstates;
	if(numPes<atmstates)
	{
	    double ratio=(double)numPes/(double)atmstates;
	    numSfGrps=(int) (ratio*(double)natm_nl);
	}
	else //there is only so far we can go
	{ 
	    numSfGrps=natm_nl-1;
	}
	if(numSfGrps==0)
	  numSfGrps=1;
    }
    if(numSfDups==1)
    {
	// numbers of duplicate caches to create needs to grow with
	// number of atom chunks mapped to multiple PEs 
	// which is map dependant... yikes. can't set this here.
	// Real range is 1 -> nstates
	int atmstates=natm_nl*nstates;
	if(numPes<atmstates)
	{
	    double ratio=(double)numPes/(double)atmstates;
	    numSfDups=(int) (ratio*(double)nstates);
	}
	else //there is only so far we can go
	{ 
	    numSfDups=nstates-1;
	}
	if(numSfDups==0)
	  numSfDups=1;

    }

    // ranges from 1 to numPes/numChareRhoG
    int temp_rho         = (int) (gExpandFactRho*2.0*((double) low_x_size + 1.0));
    if(rhoGHelpers==1)
    {
	if(numPes>temp_rho)
	{
	    rhoGHelpers=numPes/temp_rho;
	}
    }
//============================================================================
}//end routine
//============================================================================

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
    int nchareG        = sim->nchareG;
    int sizeXNL        = sim->ngrid_nloc_a;
    int nxNL           = sim->ngrid_nloc_a;
    int nyNL           = sim->ngrid_nloc_b;
    int doublePack     = config.doublePack;
    double gExpandFact = config.gExpandFact;

//============================================================================
// Get the complex data, Psi(g) and the run descriptor (z-lines in g-space)

    complex *complexPoints = new complex[numData];
    CkVec<RunDescriptor> runDescriptorVec;
    int nline_tot;
    int *istrt_lgrp   = new int [nchareG];
    int *iend_lgrp    = new int [nchareG];
    int *npts_lgrp    = new int [nchareG];
    int *nline_lgrp   = new int [nchareG];
    int *kx_line      = NULL;
    int *ky_line      = NULL;
    int *kx           = NULL;
    int *ky           = NULL;
    int *kz           = NULL;

    readStateIntoRuns(numData,complexPoints,runDescriptorVec,fname,ibinary_opt,
                      &nline_tot,&(sim->nplane_x),istrt_lgrp,iend_lgrp,
                      npts_lgrp,nline_lgrp,&kx_line,&ky_line,&kx,&ky,&kz,1);
    int nplane  = sim->nplane_x;

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
// setup the non-local stuff

    int *kx_lineNL  = new int[nline_tot];
    int *ky_lineNL  = new int[nline_tot];
   
    int ic = 0;
    kx_lineNL[ic]    = kx[ic];
    ky_lineNL[ic]    = ky[ic];
    if(kx_lineNL[ic]<0){kx_lineNL[ic]+=nxNL;}
    if(ky_lineNL[ic]<0){ky_lineNL[ic]+=nyNL;}
    for(int i = 1;i<numData;i++){
      if(kx[i]!=kx[(i-1)] || ky[i]!=ky[(i-1)]){
        ic++;
        kx_lineNL[ic] = kx[i];
        ky_lineNL[ic] = ky[i];
        if(kx_lineNL[ic]<0){kx_lineNL[ic]+=nxNL;}
        if(ky_lineNL[ic]<0){ky_lineNL[ic]+=nyNL;}
      }//endfor
    }//endfor
    ic++;

    if(ic!=nline_tot){
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("Incorrect number of lines in util.C\n");
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkExit();
    }//endif

//============================================================================
// Create the line decomposition and a sorted run descriptor
// There are two rundescriptors per line : Noah's arc sort

    int nlines_max=0;
    for(int i=0;i<nchareG;i++){nlines_max=MAX(nlines_max,nline_lgrp[i]);}
    int **index_tran_upack   = cmall_int_mat(0,nchareG,0,nlines_max,"util.C");
    int **index_tran_upackNL = cmall_int_mat(0,nchareG,0,nlines_max,"util.C");
   
    int yspace   = sizeX;
    int yspaceNL = sizeXNL;
    if(doublePack){yspace  =sizeX/2+1;}
    if(doublePack){yspaceNL=sizeXNL/2+1;}
    for(int igrp=0;igrp<nchareG;igrp++){
      for(int i=istrt_lgrp[igrp],j=0;i<iend_lgrp[igrp];i++,j++){
        index_tran_upack[igrp][j] = kx_line[i] + ky_line[i]*yspace;
        index_tran_upackNL[igrp][j] = kx_lineNL[i] + ky_lineNL[i]*yspaceNL;
      }//endfor
    }//endfor

    CkVec<RunDescriptor> *sortedRunDescriptors;
    sortedRunDescriptors = new CkVec<RunDescriptor> [nchareG];
    for(int igrp = 0; igrp < nchareG; igrp++){
      for(int i=istrt_lgrp[igrp];i<iend_lgrp[igrp];i++){
 	 int j  = 2*i;
 	 int j1 = 2*i+1;
         sortedRunDescriptors[igrp].push_back(runDescriptorVec[j]);
         sortedRunDescriptors[igrp].push_back(runDescriptorVec[j1]);
      }//endfor
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
    sim->index_tran_upackNL   = index_tran_upackNL;
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
    delete [] kx_lineNL;
    delete [] ky_lineNL;
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
#ifdef CMK_PROJECTIONS_USE_ZLIB
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
#ifdef CMK_PROJECTIONS_USE_ZLIB
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



//============================================================================
// Initialization Function
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void getSplitDecomp(int *istrt_ret,int *iend_ret,int *n_ret,
                    int ntot, int ndiv,int idiv) 
//============================================================================
  {//begin routine
//============================================================================

   if(idiv>=ndiv || ntot< ndiv){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Incorrect input to RhoGHart collection creator.\n");
     CkPrintf("idiv %d ndiv %d, ntot %d.\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
   }//endif

   int n     = (ntot/ndiv);
   int r     = (ntot%ndiv);

   int istrt = n*idiv;
   if(idiv>=r){istrt += r;}
   if(idiv<r) {istrt += idiv;}
   if(idiv<r) {n++;}
   int iend  = n+istrt;

   if(n==0){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("No lines in a RhoGHart collection!!\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
   }//endif
  
   (*n_ret)     = n;
   (*istrt_ret) = istrt;
   (*iend_ret)  = iend;

//---------------------------------------------------------------------------
  }//end routine
//============================================================================
