//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
/** \file util.C
 *
 */
//===================================================================================

#define CHARM_ON
#include "src_piny_physics_v1.0/include/class_defs/piny_constants.h"
#include "util.h"
#include "para_grp_parse.h"
#include "main/CPcharmParaInfoGrp.h"
#include "src_piny_physics_v1.0/friend_lib/proto_friend_lib_entry.h"
#include "src_mathlib/mathlib.h"
#include "src_piny_physics_v1.0/include/class_defs/allclass_gen.h"
#include "src_piny_physics_v1.0/include/class_defs/allclass_cp.h"
#include "src_piny_physics_v1.0/include/class_defs/PINY_INIT/PhysicsParamTrans.h"
#include <assert.h>
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

#include "src_piny_physics_v1.0/include/class_defs/allclass_strip_gen.h"
#include "src_piny_physics_v1.0/include/class_defs/allclass_strip_cp.h"

//===================================================================================
// set up the rho kvectors

   int *kx;
   int *ky;
   int *kz;
   int nline_tot; 
   int nPacked;
   int rhoRsubplanes = config.rhoRsubplanes;
   int rhoGHelpers   = config.rhoGHelpers;

   int sizeXEext     = sim->ngrid_eext_a;
   int sizeYEext     = sim->ngrid_eext_b;
   int sizeZEext     = sim->ngrid_eext_c;
   int sizeX         = sim->sizeX;
   int sizeY         = sim->sizeY;
   int sizeZ         = sim->sizeZ;

   double *hmati     = gencell->hmati;
   double ecut4      = 8.0*cpcoeffs_info->ecut; // convert to Ryd extra factor of 2.0

   int ka_max = 2*(sim->kx_max);
   int kb_max = 2*(sim->ky_max);
   int kc_max = 2*(sim->kz_max);

   get_rho_kvectors(ecut4,hmati,&kx,&ky,&kz,&nline_tot,&nPacked,ka_max,kb_max,kc_max);

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
      CkPrintf("Toasty Rho Line Flip-lines!\n");
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
    int nz_ext      = sim->ngrid_eext_c;

    int *kxt        = new int[nPacked];
    int *kyt        = new int[nPacked];
    int *kzt        = new int[nPacked];
    CmiMemcpy(kxt,kx,(nPacked*sizeof(int)));
    CmiMemcpy(kyt,ky,(nPacked*sizeof(int)));
    CmiMemcpy(kzt,kz,(nPacked*sizeof(int)));

    int *mapl      = new int[nline_tot];
    for(int i=0;i<nline_tot;i++){mapl[i]=i;}

    if(config.rhoLineOrder==0){
       int nsplit = (3*nplane_x)/2;
       int jj=0;
       for(int i=0;i<nsplit; i++){
       for(int j=i;j<nline_tot;j+=nsplit){
         mapl[jj]=j; jj++;
       }}//endfor
    }//end switch
    if(config.rhoLineOrder==1){
      long seed      = 174571;
      for(int j=0;j<4;j++){
      for(int i=0;i<nline_tot;i++){
        double stuff = altRandom(&seed);
        int index    = (int)( ((double) nline_tot)*stuff );
        index        = MIN(index,nline_tot-1);
        int itemp    = mapl[i];
        mapl[i]      = mapl[index];
        mapl[index]  = itemp;     
      }}//endfor
    }//endif

    int jc      = 0;
    for(int i=0;i<nline_tot;i++){
      int j = mapl[i];
      for(int lt=istrt_line[j],l=jc;lt<iend_line[j];lt++,l++){
        kx[l]    = kxt[lt];
        ky[l]    = kyt[lt];
        kz[l]    = kzt[lt];
      }//endfor
      jc+=npts_line[j];
    }//endfor
    if(jc!=nPacked)  {
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Toasty Rho Line Flip-pts!\n");  
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    delete [] mapl;
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
      CkPrintf("Toasty Line Rho Flip-lines.b!\n");
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
 
    if(kz[0]>=0){ // add neg half-line of size 0
       runs.push_back(RunDescriptor(curr_x,curr_y,0,run_length_sum,0,1,nz));
       nrun_tot      +=1;
    }//endif
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
        nrun_tot       +=1;
        run_length_sum += run_length;
        if((curr_x!=x || curr_y!=y) && kz[pNo-1]<0){ // add pos half-line of size 0
          runs.push_back(RunDescriptor(curr_x,curr_y,0,run_length_sum,0,1,nz));
          nrun_tot      +=1;
	}//endif
        if((curr_x!=x || curr_y!=y) && kz[pNo]>=0){ // add neg half-line of size 0
          runs.push_back(RunDescriptor(x,y,0,run_length_sum,0,1,nz));
          nrun_tot      +=1;
	}//endif
        curr_x          = x;
        curr_y          = y;
        curr_z          = z;
        tmpz            = z;
        run_length      = 1;
      }//endif
      if(kx[pNo]!=kx[(pNo-1)] || ky[pNo]!=ky[(pNo-1)] ){
        nline_tot_now++;
        if( (nrun_tot-1)/2 != nline_tot_now-1){
          CkPrintf("\n");
          CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	  CkPrintf("Broken Rho Run Desscriptor : %d %d %d : %d %d %d: %d %d\n",
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
    if(kz[nPacked-1]<0){ // add pos half-line of size 0
       runs.push_back(RunDescriptor(curr_x,curr_y,0,run_length_sum,0,1,nz));
       nrun_tot      +=1;
    }//endif


    if(run_length_sum!=nPacked){
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("The rho rundescriptor didn't assign all pts to rho runs \n"); 
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
      }//endif
      if(Desi->z==0 && Desi->length != 0){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("The rho rundescriptor with 1st z == 0 must have 0 lngth\n");
        CkPrintf("%d %d %d : %d %d %d\n",Desi->x,Desi->y,Desi->z,
		                         Desi1->x,Desi1->y,Desi1->z);
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit(); 
      }//endif
      if(Desi1->z<0 || (Desi->z<=nz/2 && Desi->length != 0)){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("The rho rundescriptor MUST have 1st z <=0 and 2nd >=0\n");
        CkPrintf("%d %d %d : %d %d %d\n",Desi->x,Desi->y,Desi->z,
		                         Desi1->x,Desi1->y,Desi1->z);
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit(); 
      }//endif
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

    int nlines_max_eext = 0;
    int yspaceEext=sizeXEext/2+1;
    for(int igrp=0,jgrp=0;igrp<nchareRhoG;igrp++){
      int nlTot = nline_lgrp[igrp];
      int istrt = istrt_lgrp[igrp];
      for(int k=0;k<rhoGHelpers;k++,jgrp++){
        int kstrt,kend,nl;
        getSplitDecomp(&kstrt,&kend,&nl,nlTot,rhoGHelpers,k);
        kstrt += istrt;
        nline_lgrp_eext[jgrp] = nl;
        nlines_max_eext = MAX(nlines_max_eext,nl);
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

  if(rhoRsubplanes > nplane_x){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("RhoGsubplanes parameter, %d, must be <= number\n",rhoRsubplanes);
     CkPrintf("number of gx values which span 0<=gx <=%d.\n",nplane_x);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit(); 
  }//endif

//============================================================================
// Rejigger the stuff for subplane parallelization : RHOG

    int ***index_tran_upack_rho_y  = NULL;
    int ***index_tran_pack_rho_y   = NULL;
    int **nline_send_rho_y         = NULL;
    int ***index_tran_upack_eext_y = NULL;
    int ***index_tran_upack_eext_ys= NULL;
    int ***index_tran_pack_eext_y  = NULL;
    int ***index_tran_pack_eext_ys = NULL;
    int **nline_send_eext_y        = NULL;
    int *numSubGx                  = NULL;
    int **listSubGx                = NULL;
    int listSubFlag                = 0;
    int ngxSubMax                  = nplane_x/rhoRsubplanes+1;
    if(rhoRsubplanes==1){ngxSubMax=0;}

    if(rhoRsubplanes>1){
      CkPrintf("\n");
      PRINT_LINE_STAR;
      CkPrintf("Creating the subPlane maps\n");
      PRINT_LINE_DASH;
    }//endif

    if(rhoRsubplanes>1){//subplanes on
     //----------------------------------------------------------------
     // for normal rho work
      index_tran_upack_rho_y  = cmall_itens3(0,nchareRhoG,0,rhoRsubplanes,
                                             0,nlines_max,"util.C");
      index_tran_pack_rho_y   = cmall_itens3(0,nchareRhoG,0,rhoRsubplanes,
                                             0,nlines_max,"util.C");
      nline_send_rho_y        = cmall_int_mat(0,nchareRhoG,0,rhoRsubplanes,"util.C");
     //----------------------------------------------------------------
     // for eext-ees rho work
      index_tran_upack_eext_y  = cmall_itens3(0,nchareRhoGEext,0,rhoRsubplanes,
                                              0,nlines_max_eext,"util.C");
      index_tran_upack_eext_ys = cmall_itens3(0,nchareRhoGEext,0,rhoRsubplanes,
                                              0,nlines_max_eext,"util.C");
      index_tran_pack_eext_y   = cmall_itens3(0,nchareRhoGEext,0,rhoRsubplanes,
                                              0,nlines_max_eext,"util.C");
      index_tran_pack_eext_ys  = cmall_itens3(0,nchareRhoGEext,0,rhoRsubplanes,
                                              0,nlines_max_eext,"util.C");
      nline_send_eext_y        = cmall_int_mat(0,nchareRhoGEext,0,rhoRsubplanes,"util.C");
     //----------------------------------------------------------------
     // Balanced SubPlane decomp
      numSubGx                 = new int [rhoRsubplanes];
      listSubGx                = cmall_int_mat(0,rhoRsubplanes,0,ngxSubMax,"util.C");
      listSubFlag              = 0;

     //----------------------------------------------------------------
     // Group the Gx in the subplanes so as to minimize # of message sent from RtoG
      int *listGx    = new int [nplane_x];
      int *mapGrpGx  = new int [nplane_x];
      int *mapMemGx  = new int [nplane_x];

      int div = (nplane_x/rhoRsubplanes);
      int rem = (nplane_x % rhoRsubplanes);
      for(int igrp=0;igrp<rhoRsubplanes;igrp++){
        int max        = (igrp < rem ? 1 : 0);
        numSubGx[igrp] = max+div;
      }//endfor
      create_subPlane_decomp(nplane_x,listGx,mapGrpGx,nchareRhoGEext,numSubGx,
                             nline_lgrp_eext,kx_line,nline_send_eext_y,rhoRsubplanes);
      if(config.rhoSubPlaneBalance==1){
	for(int i=0;i<nplane_x;i++){mapGrpGx[listGx[i]] = i;}
        create_gx_decomp(nline_tot,nplane_x,kx_line,mapGrpGx,rhoRsubplanes,numSubGx);
      }//endfif
      
      listSubFlag=0;
      for(int i=0;i<nplane_x;i++){
        if(listGx[i]!=i){listSubFlag==1;}
      }//endif
      if(listSubFlag==0){
        CkPrintf("This is a straight run through gx on the subplanes\n");
      }//endif
      int iii = 0;
      for(int igrp=0;igrp<rhoRsubplanes;igrp++){
        int num = numSubGx[igrp];
        if(num>1){sort_me(num,&listGx[iii]);}  //order the gx you have
        CkPrintf("subplane[%d] : %d : Gx { ",igrp,num);
        for(int jc=0,ic=iii;ic<iii+num;ic++,jc++){
          listSubGx[igrp][jc]  = listGx[ic];
          mapGrpGx[listGx[ic]] = igrp;
          mapMemGx[listGx[ic]] = jc;
          CkPrintf("%d ",listGx[ic]);
        }//endfor
        CkPrintf("}\n");
        iii += num;
      }//endfor

      //----------------------------------------------------------------
      // RhoR(gx,gy,z) parallelized by gx(subPlane) and z
      // RhoG(gx,gy,z) parallelized by collections of {gx,gy}
      // pack and upack indicies for rhoG(gx,gy,z) <-> rhoR(gx,gy,z)
      for(int igrp=0,i=0;igrp<nchareRhoG;igrp++){
      for(int ic=0;ic<rhoRsubplanes;ic++){
        nline_send_rho_y[igrp][ic]=0;
      }}//endfor

      for(int igrp=0;igrp<nchareRhoG;igrp++){
        for(int i=istrt_lgrp[igrp],j=0;i<iend_lgrp[igrp];i++,j++){
  	      int ic = mapGrpGx[kx_line[i]];       //subPlane index where kx is located
              int ip = mapMemGx[kx_line[i]];       //which kx in the subplane this is
              int jc = nline_send_rho_y[igrp][ic]; //cnt pts sent from subplane to g-chare
              nline_send_rho_y[igrp][ic]++;
              index_tran_pack_rho_y[igrp][ic][jc]  = j*sizeZ;
              index_tran_upack_rho_y[igrp][ic][jc] = ip*sizeY +  ky_line[i];
        }//endfor
      }//endfor

      int nsend_min=10000000;
      int nsend_max=0;
      for(int igrp=0;igrp<nchareRhoG;igrp++){
	for(int ic=0;ic<rhoRsubplanes;ic++){
	  if(nline_send_rho_y[igrp][ic]>0){
	    nsend_min= MIN(nsend_min, nline_send_rho_y[igrp][ic]);
  	    nsend_max= MAX(nsend_max, nline_send_rho_y[igrp][ic]);
	  }//endif
	}//endif
      }//endfor
      CkPrintf("Msg size Imbalance : rhoG <-> rhoR min %d max %d\n",
                nsend_min, nsend_max);

      //----------------------------------------------------------------
      // EextR(gx,gy,z) parallelized by gx(subPlane) and z : bigger Z than rho
      // EextG(gx,gy,z) parallelized by collections of {gx,gy} (subdivided)
      // pack and upack indicies for EextG(gx,gy,z) <-> EextR(gx,gy,z)
      for(int igrp=0,i=0;igrp<nchareRhoGEext;igrp++){
      for(int ic=0;ic<rhoRsubplanes;ic++){
        nline_send_eext_y[igrp][ic]=0;
      }}//endfor

      for(int igrp=0,i=0;igrp<nchareRhoGEext;igrp++){
        for(int j=0;j<nline_lgrp_eext[igrp];j++,i++){
  	      int ic = mapGrpGx[kx_line[i]];        //subPlane index
              int ip = mapMemGx[kx_line[i]];        //which kx in the subplane
              int jc = nline_send_eext_y[igrp][ic]; // cnt num pts sent
              nline_send_eext_y[igrp][ic]++;
              index_tran_pack_eext_y[igrp][ic][jc]  = j*sizeZEext;
              index_tran_pack_eext_ys[igrp][ic][jc] = j*sizeZ;
              index_tran_upack_eext_y[igrp][ic][jc] = ip*sizeYEext +  ky_line_ext[i];
              index_tran_upack_eext_ys[igrp][ic][jc]= ip*sizeY     +  ky_line[i];
        }//endfor
      }//endfor

      nsend_min=10000000;
      nsend_max=0;
      for(int igrp=0;igrp<nchareRhoGEext;igrp++){
	for(int ic=0;ic<rhoRsubplanes;ic++){
           if(nline_send_eext_y[igrp][ic]>0){
	    nsend_min= MIN(nsend_min, nline_send_eext_y[igrp][ic]);
	    nsend_max= MAX(nsend_max, nline_send_eext_y[igrp][ic]);
           }//endif
        }//endfor
      }//endfor
      CkPrintf("Msg size Imbalance : rhoGhart <-> rhoRhart min %d max %d\n",
                nsend_min, nsend_max);

      delete [] listGx;
      delete [] mapGrpGx;
      delete [] mapMemGx;

    }//endif : subplanes are on.

//============================================================================
// Pack up the stuff

    config.nchareRhoG            = nchareRhoG;

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

    sim->rhoRsubplanes           = config.rhoRsubplanes;
    sim->nlines_max_eext         = nlines_max_eext;
    sim->nline_send_eext_y       = nline_send_eext_y;
    sim->nline_send_rho_y        = nline_send_rho_y;
    sim->index_tran_upack_rho_y  = index_tran_upack_rho_y;
    sim->index_tran_upack_eext_y = index_tran_upack_eext_y;
    sim->index_tran_upack_eext_ys= index_tran_upack_eext_ys;
    sim->index_tran_pack_rho_y   = index_tran_pack_rho_y;
    sim->index_tran_pack_eext_y  = index_tran_pack_eext_y;
    sim->index_tran_pack_eext_ys = index_tran_pack_eext_ys;

    sim->ngxSubMax               = ngxSubMax;  // max number of gx values in any grp
    sim->numSubGx                = numSubGx;   // number of gx values in each grp
    sim->listSubGx               = listSubGx;  // gx values in grp
    sim->listSubFlag             = listSubFlag;
      
//=================================================================================
// Analyze the Send and Receives

    if(rhoRsubplanes>1){
	int Rhart_max=0;
	int Rhart_min=10000000;
        double total=0;
	for( int j=0; j< nchareRhoGEext;j++){
	    int recvCountFromRHartExt = 0;
	    for(int i=0;i<rhoRsubplanes;i++){
		if(sim->nline_send_eext_y[j][i]>0)
		  recvCountFromRHartExt++;
            }//endfor
	    recvCountFromRHartExt*=sizeZEext;
            total += recvCountFromRHartExt;
	    Rhart_max=MAX(Rhart_max,recvCountFromRHartExt);
	    Rhart_min=MIN(Rhart_min,recvCountFromRHartExt);
        }//endfor
        total /= (double)nchareRhoGEext;
	CkPrintf("GHart recv %d min msg %d max msg avg %g from RHart\n",
                  Rhart_min, Rhart_max,total);
    }//endif

    if(rhoRsubplanes>1){
	int Rho_max=0;
	int Rho_min=10000000;
        double total=0;
	for( int j=0; j< nchareRhoGEext;j++){
	    int recvCountFromRho = 0;
	    for(int i=0;i<rhoRsubplanes;i++){
		if(sim->nline_send_eext_y[j][i]>0)
		  recvCountFromRho++;
	    }//endfor
	    recvCountFromRho*=sizeZ;
            total += recvCountFromRho;
	    Rho_max=MAX(Rho_max,recvCountFromRho);
	    Rho_min=MIN(Rho_min,recvCountFromRho);
	}//endfor
        total /= (double)nchareRhoGEext;
	CkPrintf("GHart recv %d min msg %d max msg avg %g from RRho\n",
                    Rho_min, Rho_max,total);
    }//endif

    if(rhoRsubplanes>1){
	int RRho_max=0;
	int RRho_min=10000000;
        double total=0;
	for( int j=0; j< nchareRhoG;j++){
	    int recvCountFromRRho = 0;
	    for(int i=0;i<rhoRsubplanes;i++){
		if(sim->nline_send_rho_y[j][i]>0)
		  recvCountFromRRho++;
	    }//endfor
	    recvCountFromRRho*=sizeZ;
            total += recvCountFromRRho;
	    RRho_max=MAX(RRho_max,recvCountFromRRho);
	    RRho_min=MIN(RRho_min,recvCountFromRRho);
	  }//endfor
        total /= (double)nchareRhoG;
	CkPrintf("GRho recv %d min msg %d max msg avg %g from RRho\n",
                  RRho_min,RRho_max,total);
    }//endif

    if(rhoRsubplanes>1){
	int GRho_max=0;
	int GRho_min=10000000;
	for( int j=0; j< rhoRsubplanes;j++){
	    int recvCountFromGRho = 0;
	    for(int i=0;i<nchareRhoG;i++){
		if(sim->nline_send_rho_y[i][j]>0)
		  recvCountFromGRho++;
	    }//endfor
	    GRho_max=MAX(GRho_max,recvCountFromGRho);
	    GRho_min=MIN(GRho_min,recvCountFromGRho);
	}//endfor
	CkPrintf("RRho recv %d min msg %d max msg from RhoG\n",GRho_min, GRho_max);
    }//endif

    if(rhoRsubplanes>1){
	int GHart_max=0;
	int GHart_min=10000000;
	for( int j=0; j< rhoRsubplanes;j++){
	    int recvCountFromGHartExt = 0;
	    for(int i=0;i<nchareRhoGEext;i++){
		if(sim->nline_send_eext_y[i][j]>0)
		  recvCountFromGHartExt++;
	    }//endfor
	    GHart_max=MAX(GHart_max,recvCountFromGHartExt);
	    GHart_min=MIN(GHart_min,recvCountFromGHartExt);
	}//endfor
	CkPrintf("RRho recv %d min msg %d max msg from GHart\n",GHart_min, GHart_max);
    }//endif

    if(rhoRsubplanes>1){
	int GHart_max=0;
	int GHart_min=10000000;
	for( int j=0; j< rhoRsubplanes;j++){
	    int recvCountFromGHartExt = 0;
	    for(int i=0;i<nchareRhoGEext;i++){
		if(sim->nline_send_eext_y[i][j]>0)
		  recvCountFromGHartExt++;
	    }//endfor
	    GHart_max=MAX(GHart_max,recvCountFromGHartExt);
	    GHart_min=MIN(GHart_min,recvCountFromGHartExt);
	}//endfor
	CkPrintf("RHart recv %d min msg %d max msg from GHart\n",GHart_min, GHart_max);
    }//endif

    if(rhoRsubplanes>1){
      PRINT_LINE_DASH;
      CkPrintf("Completed subPlane map creation.\n");
      PRINT_LINE_STAR; CkPrintf("\n\n");
    }//endif
 
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
                      int ka_max, int kb_max, int kc_max)

//============================================================================
  {//begin routine
//============================================================================
// count the k-vectors : preserve nice symmetry for non-cubic boxes even
//                       though you need more kvectors   
 
  int iii;
  double tpi  = 2.0*M_PI;

  int nPacked = 0;
  int nPacked_np = 0;
  for(int ka=0;ka<=ka_max;ka++){
    for(int kb=-kb_max;kb<=kb_max;kb++){
      for(int kc=-kc_max;kc<=kc_max;kc++){
        double gx = tpi*(ka*hmati[1] + kb*hmati[2] + kc*hmati[3]);
        double gy = tpi*(ka*hmati[4] + kb*hmati[5] + kc*hmati[6]);
        double gz = tpi*(ka*hmati[7] + kb*hmati[8] + kc*hmati[9]);
        double g2 = gx*gx+gy*gy+gz*gz;
        if(g2<=ecut4){nPacked++;}
        if(g2<=ecut4){
          nPacked_np++;
          if(ka==0 && kb<0){nPacked_np--;}
          if(ka==0 && kb==0 && kc<=0){nPacked_np--;}
        }
      }//endfor:kc
    }//endfor:kb
  }//endfor:ka

  CkPrintf("Rho kvectors perfect sphere %d : %d %d %d\n",nPacked_np,ka_max,kb_max,kc_max);

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
	}/*endif*/
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

void readStateIntoRuns(int nPacked, int ncoef,complex *arrCP, CkVec<RunDescriptor> &runs, 
                       const char *fromFile,int ibinary_opt,
                       int *nline_tot_ret,int *nplane_ret,
                       int *istrt_lgrp,int *iend_lgrp,
                       int *npts_lgrp,int *nline_lgrp,
                       int **kx_line_ret, int **ky_line_ret,
                       int **kx_ret,int **ky_ret, int **kz_ret, int iget_decomp,
                       int gen_wave,int nx_in, int ny_in , int nz_in)

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

    int doublePack = config.doublePack;
    int nx,ny,nz;
    int *kx= new int [nPacked];
    int *ky= new int [nPacked];
    int *kz= new int [nPacked];

    if(gen_wave==0){
       readState(nPacked,arrCP,fromFile,ibinary_opt,nline_tot_ret, 
                 nplane_ret,kx,ky,kz,&nx,&ny,&nz,
                 istrt_lgrp,iend_lgrp,npts_lgrp,nline_lgrp,iget_decomp,0);
    }else{
      nx = nx_in;
      ny = ny_in;
      nz = nz_in;
      kx -= 1;  ky -= 1; kz -=1;
      PhysicsParamTransfer::fetch_state_kvecs(kx,ky,kz,ncoef,doublePack);
      kx += 1;  ky += 1; kz +=1;
      processState(nPacked,ncoef,arrCP,fromFile,ibinary_opt,nline_tot_ret,
                   nplane_ret,kx,ky,kz,
                   istrt_lgrp,iend_lgrp,npts_lgrp,nline_lgrp,iget_decomp,0,0,ny); 
    }//endif

    int nplane    = (*nplane_ret);
    int nline_tot = (*nline_tot_ret);
    int nchareG   = config.nchareG;

    if(nx!=nx_in || ny!=ny_in || nz!=nz_in){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Bad State fft size Input\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

//===================================================================================
// Read the state into the rundescriptor puppy dog

#ifdef _INPUT_FOR_KPTS_MESSES_UP_
    if(!config.doublePack){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_warning_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("The rundescriptor needs some love for the non-double pack\n"); 
      CkPrintf("It is not consistent with new FFT logic due to input data order\n");
      CkPrintf("If the data is just reordered all should be well, %s\n",fromFile);
      CkPrintf("Raz is on the job.\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_warning_@@@@@@@@@@@@@@@@@@@@\n");
    }//endif
#endif

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

    if(kz[0]>=0){ // add neg half-line of size 0
      runs.push_back(RunDescriptor(curr_x,curr_y,curr_z,run_length_sum,0,1,nz));
      nrun_tot      +=1;
    }//endif
 
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
        if((curr_x!=x || curr_y!=y) && kz[pNo-1]<0){ // add pos half-line of size 0
          runs.push_back(RunDescriptor(curr_x,curr_y,0,run_length_sum,0,1,nz));
          nrun_tot      +=1;
	}//endif
        if((curr_x!=x || curr_y!=y) && kz[pNo]>=0){ // add neg half-line of size 0
          runs.push_back(RunDescriptor(x,y,0,run_length_sum,0,1,nz));
          nrun_tot      +=1;

	}//endif
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
	  CkPrintf("Broken State Run Desscriptor : %d %d %d : %d %d %d: %d %d\n",
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
    if(kz[nPacked-1]<0){ // add pos half-line of size 0
       runs.push_back(RunDescriptor(curr_x,curr_y,0,run_length_sum,0,1,nz));
       nrun_tot      +=1;
    }//endif

    if(run_length_sum!=nPacked){
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("The state rundescriptor didn't assign all pts to runs %s\n",fromFile); 
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkExit();
    }//endif

    if((nrun_tot %2)!=0 || nrun_tot != 2*nline_tot){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("The state rundescriptor MUST find an even number of half-lines\n");
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
        CkPrintf("The state rundescriptor MUST pair up the half-lines\n");
        CkPrintf("or you will not be a happy camper : %s\n",fromFile);
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit();
      }//endif
      if(Desi->z==0 && Desi->length != 0){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("The state rundescriptor with 1st z == 0 must have 0 lngth\n");
        CkPrintf("%d %d %d : %d %d %d\n",Desi->x,Desi->y,Desi->z,
		                         Desi1->x,Desi1->y,Desi1->z);
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit(); 
      }//endif
      if(Desi1->z<0 || (Desi->z<=nz/2 && Desi->length != 0)){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("The state rundescriptor MUST have 1st z <=0 and 2nd >=0\n");
        CkPrintf("%d %d %d : %d %d %d\n",Desi->x,Desi->y,Desi->z,
		                         Desi1->x,Desi1->y,Desi1->z);
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit(); 
      }//endif
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
       CkPrintf("Incorrect number of (state) lines : %s\n",fromFile);
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
         assert(fread(&(nPackedLoc),sizeof(int),n,fp)>0);
         assert(fread(&(nx),sizeof(int),n,fp));
         assert(fread(&(ny),sizeof(int),n,fp));
         assert(fread(&(nz),sizeof(int),n,fp));
         for(int pNo=0;pNo<nPacked;pNo++) {
           int x,y,z;
           double re,im;
           assert(fread(&(re),sizeof(double),n,fp));
           assert(fread(&(im),sizeof(double),n,fp));
           assert(fread(&(x),sizeof(int),n,fp));
           assert(fread(&(y),sizeof(int),n,fp));
           assert(fread(&(z),sizeof(int),n,fp));
           arrCP[pNo] = complex(re, im);
           kx[pNo]    = x;
           ky[pNo]    = y;
           kz[pNo]    = z;
           nktot++;
           if(config.doublePack && x==0 && y==0 && z==0){break;}
	 }//endfor
       fclose(fp);
    }//endif

#ifdef _CP_DEBUG_UTIL_VERBOSE_
     CkPrintf("Done reading %s state from file: %s\n",stuff,fromFile);
#endif

//===================================================================================

    processState(nPacked,nktot,arrCP,fromFile,ibinary_opt,nline_tot_ret,nplane_ret,kx,ky,kz, 
                 istrt_lgrp,iend_lgrp,npts_lgrp,nline_lgrp,iget_decomp,iget_vstate,1,ny); 

    *nx_ret = nx;
    *ny_ret = ny;
    *nz_ret = nz;

//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================

void processState(int nPacked, int nktot, complex *arrCP, const char *fromFile,int ibinary_opt,
  	          int *nline_tot_ret,int *nplane_ret, int *kx, int *ky, int *kz, 
                  int *istrt_lgrp,int *iend_lgrp,int *npts_lgrp,int *nline_lgrp,
                  int iget_decomp,int iget_vstate, int iopt,int ny) {

//===================================================================================
// Set up some variables 

#ifdef _CP_DEBUG_UTIL_VERBOSE_
     CkPrintf("Proceesing %s state %s\n",stuff,fromFile);
#endif

    int nchareG = config.nchareG;

    char stuff[25];
    if(iget_vstate==0){
      strcpy(stuff,"coef");
    }else{
      strcpy(stuff,"velocity");
    }/*endif*/

//===================================================================================
// Add the bottom half of plane zero because the code likes to have it.

    if(config.doublePack){
       int n_ret;
       ParaGrpParse::flip_data_set(nktot,&n_ret,kx,ky,kz,arrCP,iopt);
       if(n_ret!=nPacked){
          CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          CkPrintf("Bad num pts in readState() %s: %d %d \n",nktot,n_ret,fromFile); 
          CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          CkExit();
       }//endif
    }//endif

   // We are hoping here that the sensible input for !doublePack works corrrectly

//===================================================================================
// Process the state data

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

    complex *arrCPt;
    if(iopt==1){
      arrCPt = new complex[nPacked];
      CmiMemcpy(arrCPt,arrCP,(nPacked*sizeof(complex)));
    }//endif

    int *kxt        = new int[nPacked];
    int *kyt        = new int[nPacked];
    int *kzt        = new int[nPacked];
    CmiMemcpy(kxt,kx,(nPacked*sizeof(int)));
    CmiMemcpy(kyt,ky,(nPacked*sizeof(int)));
    CmiMemcpy(kzt,kz,(nPacked*sizeof(int)));

    int jc      = 0;
    int lc      = 0;
    for(int i=0;i<nchareG; i++){
      for(int j=i;j<nline_tot;j+=nchareG){
        for(int lt=istrt_line[j],l=jc;lt<iend_line[j];lt++,l++){
          kx[l]    = kxt[lt];
          ky[l]    = kyt[lt];
          kz[l]    = kzt[lt];
          if(iopt==1){arrCP[l] = arrCPt[lt];}
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
    if(iopt==1){CmiMemcpy(arrCPt,arrCP,(nPacked*sizeof(complex)));}
    CmiMemcpy(kxt,kx,(nPacked*sizeof(int)));
    CmiMemcpy(kyt,ky,(nPacked*sizeof(int)));
    CmiMemcpy(kzt,kz,(nPacked*sizeof(int)));

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
          if(iopt==1){arrCP[jk] = arrCPt[j];}
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
    if(iopt==1 && doublePack==1){
      double norm = 0;
      for(int i=1;i<nPacked;i++){
        double wght_now = 2.0;
        if(kx[i]==0 && ky[i]<0){wght_now=0.0;}
        if(kx[i]==0 && ky[i]==0 && kz[i]<0){wght_now=0.0;}
        if(kx[i]==0 && ky[i]==0 && kz[i]==0){wght_now=1.0;}
        norm += (wght_now)*arrCP[i].getMagSqr();
      }//endfor
      double normt = 0;
      for(int i=1;i<nPacked;i++){
        double wght_now = 2.0;
        if(kxt[i]==0 && kyt[i]<0){wght_now=0.0;}
        if(kxt[i]==0 && kyt[i]==0 && kzt[i]<0){wght_now=0.0;}
        if(kxt[i]==0 && kyt[i]==0 && kzt[i]==0){wght_now=1.0;}
        normt += (wght_now)*arrCPt[i].getMagSqr();
      }//endif
      CkPrintf("state : %g %g\n",norm,normt);
    }//endif

    if(iopt==1 && doublePack==0){
      double norm = 0;
      for(int i=1;i<nPacked;i++){
        norm += arrCP[i].getMagSqr();
      }//endfor
      double normt = 0;
      for(int i=1;i<nPacked;i++){
        normt += arrCPt[i].getMagSqr();
      }//endif
      CkPrintf("state : %g %g\n",norm,normt);
    }//endif
#endif

//===================================================================================
// Fix !double pack

#ifdef _INPUT_FOR_KPTS_MESSES_UP_
    if(!config.doublePack){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_warning_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("The rundescriptor needs some love for the non-double pack\n"); 
      CkPrintf("It is not consistent with new FFT logic due to input data order\n");
      CkPrintf("If the data is just reordered all should be well, %s\n",fromFile);
      CkPrintf("Raz is on the job.\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_warning_@@@@@@@@@@@@@@@@@@@@\n");
    }//endif
#endif

//===================================================================================
// A little output to the screen!

    *nplane_ret    = nplane;
    *nline_tot_ret = nline_tot;

    if(iopt==1){delete [] arrCPt;}
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
     CkPrintf("Done proceesing %s state %s\n",stuff,fromFile);
#endif

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
    sprintf(fname,"%s.s0.k0.b0.t0/state%d.out",config.dataPath,1);

    int numData        = config.numData;
    int ibinary_opt    = sim->ibinary_opt;
    int sizeY          = sim->sizeY;
    int sizeZ          = sim->sizeZ;
    int nchareG        = sim->nchareG;
    int sizeXNL        = sim->ngrid_nloc_a;
    int nxNL           = sim->ngrid_nloc_a;
    int nyNL           = sim->ngrid_nloc_b;
    int gen_wave       = sim->gen_wave; 
    int ncoef          = sim->ncoef;
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

    readStateIntoRuns(numData,ncoef,complexPoints,runDescriptorVec,fname,ibinary_opt,
                      &nline_tot,&(sim->nplane_x),istrt_lgrp,iend_lgrp,
                      npts_lgrp,nline_lgrp,&kx_line,&ky_line,&kx,&ky,&kz,1,gen_wave,
                      sizeX,sizeY,sizeZ);
    int nplane  = sim->nplane_x;

    if(config.nGplane_x != nplane && config.doublePack){
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("Mismatch in allowed gspace planes %d %d\n",config.nGplane_x,nplane);
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

    int kx_min = kx[0];
    int kx_max = kx[0];
    int ky_min = ky[0];
    int ky_max = ky[0];
    int kz_min = kz[0];
    int kz_max = kz[0];
    for(int i = 1;i<numData;i++){
      kx_min = (kx_min <= kx[i] ? kx_min : kx[i]);
      ky_min = (ky_min <= ky[i] ? ky_min : ky[i]);
      kz_min = (kz_min <= kz[i] ? kz_min : kz[i]);
      kx_max = (kx_max >= kx[i] ? kx_max : kx[i]);
      ky_max = (ky_max >= ky[i] ? ky_max : ky[i]);
      kz_max = (kz_max >= kz[i] ? kz_max : kz[i]);
    }//endfor
    kx_max = (kx_max >= -kx_min ? kx_max : -kx_min);
    ky_max = (ky_max >= -ky_min ? ky_max : -ky_min);
    kz_max = (kz_max >= -kz_min ? kz_max : -kz_min);

    sim->kx_max = kx_max;
    sim->ky_max = ky_max;
    sim->kz_max = kz_max;

    if(2*kx_max >= sizeX/2 || 2*ky_max >= sizeY/2 || 2*kz_max >= sizeZ/2){
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("Bad FFT sizes : (%d,%d), (%d,%d), (%d,%d)\n",
                 2*kx_max,sizeX/2,2*ky_max,sizeY/2,2*kz_max,sizeZ/2);
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
        index_tran_upack[igrp][j]   = kx_line[i]   + ky_line[i]*yspace;
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
      if(doublePack==1){
        for(int j=0;j<npts_lgrp[i];j++){
          int iii = j+index_output_off[i];
          if(kx[iii]==0 && ky[iii]>0){num_uni[i]++;}
          if(kx[iii]==0 && ky[iii]<0){num_red[i]++;}
          if(kx[iii]==0 && ky[iii]==0 && kz[iii]>=0){num_uni[i]++;}
          if(kx[iii]==0 && ky[iii]==0 && kz[iii]<0){num_red[i]++;}
        }//endfor
        nk0_max=MAX(num_uni[i],nk0_max);
        nk0_max=MAX(num_red[i],nk0_max);
      }//endif
    }//endif

    // Find where my unique guys go and make a list so I can send them.
    // Make a list of where the unique guys arrive so I can receive them.

    RedundantCommPkg *RCommPkg = new RedundantCommPkg [nchareG]; 
    for(int i=0;i<nchareG;i++){RCommPkg[i].Init(nk0_max,nchareG);} // mallocs and zeros 

    if(doublePack==1){

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
      CkPrintf("send recv stuff for %d: %d %d \n",i,num_send_tot,num_recv_tot);
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

    }//endif :: Do redundancy when doublepack is off

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
                    int ibinary_write_opt,int iteration, int istate)
//=============================================================================
  { //begin rotunie
//=============================================================================

  int *index = new int [ncoef];
  int *ktemp = new int [ncoef];
  int istrt=0;
  sort_psi_output(ncoef,k_x,k_y,k_z,index,ktemp,&istrt);
  int ncoef_true=ncoef-istrt;

  if(istate==1){
    char fname[1000];
    sprintf(fname,"%s.s0.k0.b0.t0/TimeStamp",config.dataPathOut);
    FILE *fp = fopen(fname,"w");
     fprintf(fp,"time step = %d\n",iteration);
    fclose(fp);
  }//endif

   if(ibinary_write_opt==0){
     FILE *fp  = fopen(psiName,"w");
     assert(fp!=NULL);
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
#if CMK_PROJECTIONS_USE_ZLIB
   else if(ibinary_write_opt==2){
     strcat(psiName,".gz");
     gzFile zfp  = gzopen(psiName,"w");
     assert(zfp!=NULL);
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
     assert(zfp!=NULL);
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
     assert(fp!=NULL);
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
      assert(fp!=NULL);
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
#if CMK_PROJECTIONS_USE_ZLIB
   else if(ibinary_write_opt==2){
     strcat(vpsiName,".gz");
     gzFile zfp=gzopen(vpsiName,"w");
     assert(zfp!=NULL);
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
     assert(zfp!=NULL);
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
      assert(fp!=NULL);
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

//=============================================================================
// Sort on ky for each kx

  int kmin = kx[0];
  int kmax = kx[(n-1)];
  for(int i=0;i<=kmax-kmin;i++){ktemp[i]=0;}
  for(int i=0;i<n;i++){ktemp[kx[i]-kmin]++;}

  int ioff=0;
  for(int i=0;i<=kmax-kmin;i++){
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



//============================================================================
// Create some decompositions and find the best one.
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void create_subPlane_decomp(int nplane_x,int *listGx,int *mapGrpGx,
                            int nchareRhoGEext,int *numSubGx,
                            int *nline_lgrp_eext,int *kx_line,
                            int **nline_send_eext_y, int rhoRsubplanes){
//============================================================================
// try out some several decomps and find the best

  int *list1 = new int [nplane_x];
  long seed  = 1145;

  int score_max = 0;
  for(int j=0;j<nplane_x;j++){list1[j]=j;}
  for(int ntry=1;ntry<100;ntry++){

   //------------------------------------------------------
   // Mix up the kx
    if(ntry>1){
      double stuff = altRandom(&seed);
      int index    = (int) (((double)nplane_x)*stuff);
      index        = MIN(index,nplane_x-1);
      stuff        = altRandom(&seed);
      int jndex    = (int) (((double)nplane_x)*stuff);
      jndex        = MIN(jndex,nplane_x-1);
      int itemp    = list1[index];
      list1[index] = list1[jndex];
      list1[jndex] = itemp;
    }//endif

   //------------------------------------------------------
   // Create the decomp
    int iii = 0;
    for(int igrp=0;igrp<rhoRsubplanes;igrp++){
      int num = numSubGx[igrp];
      for(int jc=0,ic=iii;ic<iii+num;ic++,jc++){
        mapGrpGx[list1[ic]] = igrp;
      }//endfor
      iii += num;
    }//endfor

   //------------------------------------------------------
   // Score the decomp
    int score;
    score_subPlane_decomp(nchareRhoGEext,rhoRsubplanes,nline_lgrp_eext,
                          mapGrpGx,kx_line,nline_send_eext_y,&score);

   //------------------------------------------------------
   // Keep the best decomp
    //    CkPrintf("try %d score %d score_max %d : ",ntry,score,score_max);
    //    for(int i=0;i<nplane_x;i++){CkPrintf("%d ",list1[i]);} CkPrintf("\n");
    if(score<score_max || ntry==1){
      score_max = score;
      for(int i=0;i<nplane_x;i++){listGx[i]=list1[i];}
    }else{
      for(int i=0;i<nplane_x;i++){list1[i]=listGx[i];}
    }//endif

  }//endfor

  delete [] list1;

//============================================================================
 }//end routine
//============================================================================


//============================================================================
// Score a decomposition
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void score_subPlane_decomp(int nchareRhoGEext,int rhoRsubplanes, int *nline_lgrp_eext,
                           int *mapGrpGx, int *kx_line,int **nline_send_eext_y, 
                           int *Rhart_max){
//============================================================================
// Count the sends

  for(int igrp=0,i=0;igrp<nchareRhoGEext;igrp++){
   for(int ic=0;ic<rhoRsubplanes;ic++){
     nline_send_eext_y[igrp][ic]=0;
  }}//endfor

  for(int igrp=0,i=0;igrp<nchareRhoGEext;igrp++){
    for(int j=0;j<nline_lgrp_eext[igrp];j++,i++){
      int ic = mapGrpGx[kx_line[i]];   //subPlane index
      nline_send_eext_y[igrp][ic]++;
    }//endfor
  }//endfor

//============================================================================
// Find the total number of send/recvs

  Rhart_max[0] = 0;
  for( int j=0; j< nchareRhoGEext;j++){
    int recvCountFromRHartExt = 0;
    for(int i=0;i<rhoRsubplanes;i++){
      if(nline_send_eext_y[j][i]>0){recvCountFromRHartExt++;}
    }//endfor
    Rhart_max[0]+=recvCountFromRHartExt;
  }//endfor

//============================================================================
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void create_gx_decomp(int nktot, int nline, int *kx_line, int *mapl, 
                      int nchare,int *nline_grp){
//============================================================================
// count the pts in the lines

  int *npts = new int [nline];
  for(int i =0;i<nline;i++){npts[i]=0;}
  for(int i =0;i<nktot;i++){npts[mapl[kx_line[i]]]++;} // ky pts per kx line

//============================================================================
// Find the best division

  double dev_min  = 0.0;
  int nbal_min    = 0;
  int ibal_min    = 0;
  int ifirst      = 0;
  int mmm         = (nktot/nchare);
  for(int ibal=0;ibal<=mmm-1;ibal++){
    int nbal  = ibal;
    int ntarg = (nktot/nchare);
    if(ntarg > nbal){ntarg -= nbal;}
    int nmax  = 0;
    int nmin  = nktot;
    int nnow  = 0;
    int ic    = 0;
    for(int i=0;i<nline;i++){
      nnow += npts[i];
      if( (nnow>=ntarg) && (ic<(nchare-1)) ){
        ic+=1;
        nmin = MIN(nmin,nnow);
        nmax = MAX(nmax,nnow);
        nnow = 0;
      }//endif
    }//endfor
    nmin = MIN(nmin,nnow);
    nmax = MAX(nmax,nnow);
    double dev = 100.0*((double)(nmax-nmin))/((double)MAX(nmin,1));
    if(ic==nchare-1){
     if(dev<dev_min || ifirst==0){
       ifirst   = 1;
       dev_min  = dev;
       nbal_min = nbal;
       ibal_min = ibal;
     }//endif
    }//endif
  }//endfor

//==========================================================================
// Store the good decomposition

  int ntarg = (nktot/nchare);
  if(ntarg > nbal_min){ntarg = ntarg-nbal_min;}

  int ic   = 0;
  int nnow = 0;
  for(int i=0;i<nchare;i++){nline_grp[i]=0;}
  for(int i=0;i<nline;i++){
     nline_grp[ic]++;
     nnow += npts[i];
     if( (nnow>=ntarg) && (ic<(nchare-1)) ){
       ic  +=1;
       nnow = 0;
    }//endif
  }//endfor

  delete [] npts;

//============================================================================
  }//end routine
//============================================================================
