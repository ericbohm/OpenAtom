//=============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=============================================================================
#ifndef _parainfo_h_
#define _parainfo_h_

#include "../include/RunDescriptor.h"
#include "../src_piny_physics_v1.0/friend_lib/proto_friend_lib_entry.h"

//=============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=============================================================================
class RedundantCommPkg {
  public :
      int nk0_max;
      int nchareG;
      int num_recv_tot; //how many chares do I recieve from
      int num_send_tot; //how many chares do I send to
      int *num_send;
      int **lst_send;
      int *num_recv;
      int **lst_recv;

  RedundantCommPkg(){
     nk0_max=0; nchareG=0; 
     num_send = NULL;
     lst_send = NULL;
     num_recv = NULL;
     lst_recv = NULL;
  }

  void Init(int nk0_max_,int nchareG_){
    nk0_max  = nk0_max_;
    nchareG  = nchareG_;
    num_send = new int[nchareG];
    num_recv = new int[nchareG];
    for(int i=0;i<nchareG;i++){num_send[i] = 0;}
    for(int i=0;i<nchareG;i++){num_recv[i] = 0;}

    lst_send = new int *[nchareG];
    lst_recv = new int *[nchareG];
    for(int i=0;i<nchareG;i++){
      lst_send[i]=new int[nk0_max];
      lst_recv[i]=new int[nk0_max];
    }//endfor
    for(int i=0;i<nchareG;i++){
    for(int j=0;j<nk0_max;j++){
      lst_send[i][j] = 0;
      lst_recv[i][j] = 0;
    }}
  }

  void Init(RedundantCommPkg *R){
    nk0_max      = R->nk0_max;
    nchareG      = R->nchareG;
    num_recv_tot = R->num_recv_tot;
    num_send_tot = R->num_send_tot;
    Init(nk0_max,nchareG);
    for(int i=0;i<nchareG;i++){
      num_send[i]=R->num_send[i];
      num_recv[i]=R->num_recv[i];
    }//endfor
    for(int i=0;i<nchareG;i++){
      for(int j=0;j<num_send[i];j++){
        lst_send[i][j]=R->lst_send[i][j];
      }//endfor
      for(int j=0;j<num_recv[i];j++){
        lst_recv[i][j]=R->lst_recv[i][j];
      }//endfor
    }//endfor
  }

 ~RedundantCommPkg(){
    delete [] num_send;
    delete [] num_recv;
    for(int j=0;j<nchareG;j++){
      delete [] lst_send[j];
      delete [] lst_recv[j];
    }//endfor
    delete [] lst_send;
    delete [] lst_recv;
  }

  void pup(PUP::er &p){
    p|nk0_max;        p|nchareG;
    p|num_recv_tot;   p|num_send_tot;
    int *lst_send_all = new int [nk0_max*nchareG];
    int *lst_recv_all = new int [nk0_max*nchareG];
    if(p.isUnpacking()) {
      Init(nk0_max,nchareG);
    }else{
      int iii=0;
      for(int i=0;i<nchareG;i++){
      for(int j=0;j<nk0_max;j++){
        lst_send_all[iii]=lst_send[i][j];
        lst_recv_all[iii]=lst_recv[i][j];
        iii++;
      }}
    }//endif
    p(num_send,nchareG);
    p(num_recv,nchareG);
    p(lst_send_all,nk0_max*nchareG);
    p(lst_recv_all,nk0_max*nchareG);
    if(p.isUnpacking()) {
      int iii=0;
      for(int i=0;i<nchareG;i++){
      for(int j=0;j<nk0_max;j++){
        lst_send[i][j] = lst_send_all[iii];
        lst_recv[i][j] = lst_recv_all[iii];
        iii++;
      }}
    }//endif unpacking
    delete [] lst_send_all;
    delete [] lst_recv_all;
  }

//----------------------------------------------------------------------------
   };// end class
//=============================================================================


//=============================================================================
PUPmarshall(RedundantCommPkg);
//=============================================================================

//=============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=============================================================================
class CPcharmParaInfo {
//=============================================================================
 public:

   double vol,dt,tol_norb,tol_cp_min;
   int ndump_frq;
   int istart_typ_cp;
   int cp_grad_corr_on;
   int cp_opt; 
   int cp_std;
   int cp_wave;
   int cp_min_opt;
   int cp_min_cg;
   int cp_min_std;
   int sizeX, sizeY, sizeZ;  
   int nstates;
   int ntime;
   int ibinary_opt;
   int ibinary_write_opt;
   int nplane_x;   // # of non-zero planes of gx
   int nchareG;    // # of collections of lines in g-space
   int nlines_tot; // total number of lines in g-space
   int nlines_max; // maximum number of lines in any collection
   int npts_tot;   // total number of pts in g-space
   int natm_tot;
   int natm_nl;
   int numSfGrps;
   int natm_nl_grp_max;
   double *lines_per_chareG;  // lines in each collection (g-space)
   double *pts_per_chareG;    // pts in each collection   (g-space)
   CkVec<RunDescriptor> *sortedRunDescriptors; // description of collection
   int *nlines_per_chareG;   // lines in each collection (g-space)
   int *npts_per_chareG;     // pts in each collection   (g-space)
   int *index_output_off;    // offset to unpack output fore each chareg
   int **index_tran_upack;   // unpack indicies for each chareG


   int nlines_tot_rho; // total number of lines in rhog-space
   int nlines_max_rho; // maximum number of lines in any rho collection
   int npts_tot_rho;   // total number of pts in rhog-space

   int nplane_rho_x;   // # of non-zero planes of rho gx
   int nchareRhoG;    // # of collections of lines in rhog-space
   double *lines_per_chareRhoG;  // lines in each collection (rhog-space)
   double *pts_per_chareRhoG;    // pts in each collection   (rhog-space)
   int *nlines_per_chareRhoG;   // lines in each collection (rhog-space)
   int *npts_per_chareRhoG;     // pts in each collection   (rhog-space)
   int **index_tran_upack_rho;   // unpack indicies for each chareRhoG
   CkVec<RunDescriptor> *RhosortedRunDescriptors; // description of collection
   RedundantCommPkg *RCommPkg;  // communication of redundant elements
//=============================================================================


//=============================================================================
   CPcharmParaInfo(CPcharmParaInfo &s){
//=============================================================================
#ifdef _CP_DEBUG_PARAINFO_VERBOSE_
     CkPrintf("CPcharmParaInfo copy constructor\n");
#endif
     vol         = s.vol;
     tol_norb    = s.tol_norb;
     tol_cp_min  = s.tol_cp_min;
     dt          = s.dt;
     ndump_frq   = s.ndump_frq;
     istart_typ_cp = s.istart_typ_cp;
     cp_grad_corr_on = s.cp_grad_corr_on;
     cp_opt      = s.cp_opt; 
     cp_std      = s.cp_std;
     cp_wave     = s.cp_wave;
     cp_min_opt  = s.cp_min_opt;
     cp_min_cg   = s.cp_min_cg;
     cp_min_std  = s.cp_min_std;
     sizeX       = s.sizeX;
     sizeY       = s.sizeY;
     sizeZ       = s.sizeZ;
     nstates     = s.nstates;
     ntime       = s.ntime;
     ibinary_opt = s.ibinary_opt;
     ibinary_write_opt = s.ibinary_write_opt;
     nplane_x    = s.nplane_x;
     nchareG     = s.nchareG;
     nlines_tot  = s.nlines_tot;
     npts_tot    = s.npts_tot;
     nlines_max  = s.nlines_max;
     natm_tot    = s.natm_tot;
     natm_nl     = s.natm_nl;
     numSfGrps   = s.numSfGrps;
     natm_nl_grp_max = s.natm_nl_grp_max;
     if(nplane_x==0 || nplane_x > sizeX){
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("Error in CPcharmParaInfo constructor\n");
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkExit();
     }//endif
     if(nchareG==0 || nchareG > sizeX || nchareG < nplane_x){
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("Error in CPcharmParaInfo constructor\n");
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkExit();
     }//endif
     nlines_tot_rho  = s.nlines_tot_rho;
     npts_tot_rho    = s.npts_tot_rho;
     nlines_max_rho  = s.nlines_max_rho;
     nplane_rho_x    = s.nplane_rho_x;
     nchareRhoG     = s.nchareRhoG;
     npts_per_chareRhoG   = new int[nchareRhoG];
     nlines_per_chareRhoG = new int[nchareRhoG];
     lines_per_chareRhoG  = new double[nchareRhoG];
     pts_per_chareRhoG    = new double[nchareRhoG];
     if(nplane_rho_x==0 || nplane_rho_x > sizeX){
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("Error in CPcharmParaInfo constructor\n");
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkExit();
     }//endif
     if(nchareRhoG==0 || nchareRhoG > sizeX || nchareRhoG < nplane_rho_x){
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("Error in CPcharmParaInfo constructor\n");
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkExit();
     }//endif
     for(int i=0;i<nchareRhoG;i++){
       nlines_per_chareRhoG[i] = s.nlines_per_chareRhoG[i];
       lines_per_chareRhoG[i]  = s.lines_per_chareRhoG[i];
       pts_per_chareRhoG[i]    = s.pts_per_chareRhoG[i];
       npts_per_chareRhoG[i]   = s.npts_per_chareRhoG[i];
     }//endfor
     index_tran_upack_rho = cmall_int_mat(0,nchareRhoG,0,nlines_max_rho,"cpcharmparainfo.h");
     for(int i=0;i<nchareRhoG;i++){
      for(int j=0;j<nlines_per_chareRhoG[i];j++){
        index_tran_upack_rho[i][j] = s.index_tran_upack_rho[i][j];
      }//endfor
     }//endfor

     sortedRunDescriptors = new CkVec<RunDescriptor> [sizeX];
     for(int i=0;i<nchareG;i++){
       for(int j=0;j<s.sortedRunDescriptors[i].size();j++){
          sortedRunDescriptors[i].push_back(s.sortedRunDescriptors[i][j]);
       }//endfor
     }//endfor

     RhosortedRunDescriptors = new CkVec<RunDescriptor> [nchareRhoG];
     for(int i=0;i<nchareRhoG;i++){
       for(int j=0;j<s.RhosortedRunDescriptors[i].size();j++){
          RhosortedRunDescriptors[i].push_back(s.RhosortedRunDescriptors[i][j]);
       }//endfor
     }//endfor
     npts_per_chareG   = new int[nchareG];
     index_output_off  = new int[nchareG];
     nlines_per_chareG = new int[nchareG];
     lines_per_chareG  = new double[nchareG];
     pts_per_chareG    = new double[nchareG];
     for(int i=0;i<nchareG;i++){
       nlines_per_chareG[i] = s.nlines_per_chareG[i];
       lines_per_chareG[i]  = s.lines_per_chareG[i];
       pts_per_chareG[i]    = s.pts_per_chareG[i];
       npts_per_chareG[i]   = s.npts_per_chareG[i];
       index_output_off[i]  = s.index_output_off[i];
     }//endfor
     index_tran_upack = cmall_int_mat(0,nchareG,0,nlines_max,"cpcharmparainfo.h");
     for(int i=0;i<nchareG;i++){
      for(int j=0;j<nlines_per_chareG[i];j++){
        index_tran_upack[i][j] = s.index_tran_upack[i][j];
      }//endfor
     }//endfor

     RCommPkg = new RedundantCommPkg [nchareG]; 
     for(int i=0;i<nchareG;i++){RCommPkg[i].Init(&s.RCommPkg[i]);}

     LBTurnInstrumentOff();
   }//end constructor
//=============================================================================

//=============================================================================
   CPcharmParaInfo() {
//=============================================================================
#ifdef _CP_DEBUG_PARAINFO_VERBOSE_
     CkPrintf("CPcharmParaInfo constructor\n");
#endif
       lines_per_chareG=NULL; 
       pts_per_chareG=NULL;
       nlines_per_chareG=NULL; 
       npts_per_chareG=NULL;
       index_output_off=NULL;

       lines_per_chareRhoG=NULL; 
       pts_per_chareRhoG=NULL;
       nlines_per_chareRhoG=NULL; 
       npts_per_chareRhoG=NULL;
       RCommPkg  = NULL;
   }
//=============================================================================

//=============================================================================
  ~CPcharmParaInfo() {
      delete []  lines_per_chareG;  lines_per_chareG = NULL;
      delete [] nlines_per_chareG; nlines_per_chareG = NULL;
      delete []  pts_per_chareG;    pts_per_chareG   = NULL;
      delete [] npts_per_chareG;   npts_per_chareG   = NULL;
      delete [] index_output_off;   index_output_off   = NULL;

      delete []  lines_per_chareRhoG;  lines_per_chareRhoG = NULL;
      delete [] nlines_per_chareRhoG; nlines_per_chareRhoG = NULL;
      delete []  pts_per_chareRhoG;    pts_per_chareRhoG   = NULL;
      delete [] npts_per_chareRhoG;   npts_per_chareRhoG   = NULL;
      delete [] RCommPkg;             RCommPkg= NULL;
      /*      for(int i=0;i<nchareG;i++)
	delete [] index_tran_upack[i];
      delete [] index_tran_upack;
      for(int i=0;i<nchareRhoG;i++)
	delete [] index_tran_upack_rho[i];
      delete [] index_tran_upack_rho;
      */

  }//end destructor
//=============================================================================

//=============================================================================
  void pup(PUP::er &p){
//=============================================================================
#ifdef _CP_DEBUG_PARAINFO_VERBOSE_
     CkPrintf("CPcharmParaInfo pup\n");
#endif
      p|vol;        p|dt;         p|tol_norb;   p|tol_cp_min;
      p|ndump_frq;  p|istart_typ_cp; p|cp_grad_corr_on;
      p|cp_opt;     p|cp_std;     p|cp_wave;
      p|cp_min_opt; p|cp_min_std; p|cp_min_cg;
      p|sizeX;      p|sizeY;      p|sizeZ;  
      p|nplane_x;   p|nchareG;    p|natm_tot;    p|natm_nl;
      p|nstates;    p|ntime;      p|ibinary_opt; p|ibinary_write_opt;
      p|natm_tot;   p|natm_nl;    p|numSfGrps;
      p|natm_nl_grp_max;  p|nlines_tot; p|npts_tot;
      p|nlines_max;
      p|nplane_rho_x; p|nchareRhoG; 
      p|nlines_max_rho; p|nlines_tot_rho; p|npts_tot_rho;
      if(p.isUnpacking()) {
        lines_per_chareRhoG  = new double[nchareRhoG];
        pts_per_chareRhoG    = new double[nchareRhoG];
        nlines_per_chareRhoG = new int[nchareRhoG];
        npts_per_chareRhoG   = new int[nchareRhoG];
#ifdef _CP_DEBUG_PARAINFO_VERBOSE_
	CkPrintf("nchareRhoG %d nlines_max_rho %d sizeX %d\n",
             nchareRhoG, nlines_max_rho, sizeX);
#endif
        index_tran_upack_rho = cmall_int_mat(0,nchareRhoG,0,nlines_max_rho,
                                             "cpcharmparainfo.h");
        RhosortedRunDescriptors = new CkVec<RunDescriptor> [nchareRhoG];
      }//enddif : unpacking
#ifdef _CP_DEBUG_PARAINFO_VERBOSE_
      CkPrintf("CPcharmParaInfo pup 2 \n");
#endif
      p(lines_per_chareRhoG, nchareRhoG);
      p(nlines_per_chareRhoG, nchareRhoG);
      p(pts_per_chareRhoG, nchareRhoG);
      p(npts_per_chareRhoG, nchareRhoG);
      for(int igrp=0;igrp<nchareRhoG;igrp++){
	p|RhosortedRunDescriptors[igrp];
      }
      for(int igrp=0;igrp<nchareRhoG;igrp++){
	p(index_tran_upack_rho[igrp],nlines_per_chareRhoG[igrp]);
      }
#ifdef _CP_DEBUG_PARAINFO_VERBOSE_
      CkPrintf("CPcharmParaInfo pup 3 \n");
#endif
      if(p.isUnpacking()) {
        lines_per_chareG  = new double[nchareG];
        pts_per_chareG    = new double[nchareG];
        nlines_per_chareG = new int[nchareG];
        npts_per_chareG   = new int[nchareG];
        index_output_off  = new int[nchareG];
        index_tran_upack = cmall_int_mat(0,nchareG,0,nlines_max,"cpcharmparainfo.h");
        sortedRunDescriptors = new CkVec<RunDescriptor> [nchareG];
      }//endif : unpacking
#ifdef _CP_DEBUG_PARAINFO_VERBOSE_
      CkPrintf("CPcharmParaInfo pup 4 \n");
#endif
      p(lines_per_chareG,nchareG);
      p(nlines_per_chareG,nchareG);
      p(pts_per_chareG,nchareG);
      p(npts_per_chareG,nchareG);
      p(index_output_off,nchareG);
      for(int igrp=0;igrp<nchareG;igrp++){
	p|sortedRunDescriptors[igrp];
      }//endfor
      for(int igrp=0;igrp<nchareG;igrp++){
	p(index_tran_upack[igrp],nlines_per_chareG[igrp]);
      }//endfor
#ifdef _CP_DEBUG_PARAINFO_VERBOSE_
      CkPrintf("CPcharmParaInfo pup 5 \n");
#endif
      if(p.isUnpacking()){
	  if(sizeX<0|sizeY<0|sizeZ<0|RhosortedRunDescriptors[0].size()<=0){
	    int *p=NULL;
	    *p=12;  //segfault now please
#ifdef _CP_DEBUG_PARAINFO_VERBOSE_
	  }else{
	    CkPrintf("unpacked RhosortedRunDescriptors[0].size()=%d\n",
              RhosortedRunDescriptors[0].size());
#endif
	  }//endif
      }//endif
      if(p.isUnpacking()){
         RCommPkg = new RedundantCommPkg [nchareG]; 
      }//endif
      for(int igrp=0;igrp<nchareG;igrp++){
        RCommPkg[igrp].pup(p);
      }//endfor
   };
//----------------------------------------------------------------------------
   }; // end class
//=============================================================================


//=============================================================================
PUPmarshall(CPcharmParaInfo);
//=============================================================================


//=============================================================================
#endif
