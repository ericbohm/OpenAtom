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
   if(num_send!=NULL)
     delete [] num_send;
   if(num_recv!=NULL)
    delete [] num_recv;
   if(lst_send!=NULL)
     {
       for(int j=0;j<nchareG;j++){
	 if(lst_send[j]!=NULL)
	   delete [] lst_send[j];
	 if(lst_recv[j]!=NULL)
	   delete [] lst_recv[j];
       }//endfor
       delete [] lst_send;
       delete [] lst_recv;
     }
  }

  void pup(PUP::er &p){
    p|nk0_max;        p|nchareG;
    p|num_recv_tot;   p|num_send_tot;
    if(p.isUnpacking()) {
      Init(nk0_max,nchareG);
    }
    PUParray(p,num_send,nchareG);
    PUParray(p,num_recv,nchareG);
    for(int i=0;i<nchareG;i++)
      PUParray(p,lst_send[i],nk0_max);
    for(int i=0;i<nchareG;i++)
      PUParray(p,lst_recv[i],nk0_max);
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
   double tol_cp_dyn;
   int fftopt;
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
   int rhoRsubplanes;
   int ees_eext_on;         //Opt: EES option on for external energy
   int ees_nloc_on;         //Opt: EES option on for non-local energy
   int ngrid_nloc_a, ngrid_nloc_b, ngrid_nloc_c; 
   int ngrid_eext_a, ngrid_eext_b, ngrid_eext_c; 
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
   int nlIters;        // number of non-local iterations
   int nmem_zmat_max;  // total size of zmatrix memory
   int nmem_zmat_tot;  // total size of zmatrix memory
   int *ioff_zmat;     // offset into zmatrix memory
   int *nmem_zmat;     // zmatrix memory for each iteration
   int natm_nl;
   int natm_typ;
   int numSfGrps;
   int natm_nl_grp_max;
   double *lines_per_chareG;  // lines in each collection (g-space)
   double *pts_per_chareG;    // pts in each collection   (g-space)
   CkVec<RunDescriptor> *sortedRunDescriptors; // description of collection
   int *nlines_per_chareG;   // lines in each collection (g-space)
   int *npts_per_chareG;     // pts in each collection   (g-space)
   int *index_output_off;    // offset to unpack output fore each chareg
   int **index_tran_upack;   // unpack indicies for each chareG
   int **index_tran_upackNL; // unpack indicies for each chareG


   int nlines_tot_rho; // total number of lines in rhog-space
   int nlines_max_rho; // maximum number of lines in any rho collection
   int npts_tot_rho;   // total number of pts in rhog-space
   int nlines_max_eext;// rhoghelpers subdivided g-collections

   int nplane_rho_x;   // # of non-zero planes of rho gx
   int nchareRhoG;    // # of collections of lines in rhog-space
   int nchareRhoGEext;  // # of collections of lines in rhog-space for eext
   double *lines_per_chareRhoG;  // lines in each collection (rhog-space)
   double *pts_per_chareRhoG;    // pts in each collection   (rhog-space)

   int *nlines_per_chareRhoG;   // lines in each collection (rhog-space)
   int *nlines_per_chareRhoGEext;  // lines in each collection (rhog-space)
   int *npts_per_chareRhoG;     // pts in each collection   (rhog-space)
   int **index_tran_upack_rho;   // unpack indicies for each chareRhoG
   int **index_tran_upack_eext;   // unpack indicies for each chareRhoG

   int **nline_send_eext_y;       // g-collections send/recv R to G in Rsubplane decomp
   int **nline_send_rho_y;        // g-collections send/recv R to G in Rsubplane decomp
   int ***index_tran_upack_rho_y; // fancy indices to pack/upack collections
   int ***index_tran_upack_eext_y;// in Rsubplane decomp.
   int ***index_tran_upack_eext_ys;// in Rsubplane decomp.
   int ***index_tran_pack_rho_y;
   int ***index_tran_pack_eext_y;
   int ***index_tran_pack_eext_ys;

   CkVec<RunDescriptor> *RhosortedRunDescriptors; // description of collection
   RedundantCommPkg *RCommPkg;  // communication of redundant elements

//=============================================================================

//=============================================================================
   CPcharmParaInfo(CPcharmParaInfo &s){
//=============================================================================
#ifdef _CP_DEBUG_PARAINFO_VERBOSE_
    CkPrintf("CPcharmParaInfo copy constructor\n");
#endif
     vol          = s.vol;
     tol_norb     = s.tol_norb;
     tol_cp_min   = s.tol_cp_min;
     tol_cp_dyn   = s.tol_cp_dyn;
     dt           = s.dt;
     fftopt       = s.fftopt;
     ndump_frq    = s.ndump_frq;
     istart_typ_cp= s.istart_typ_cp;
     cp_grad_corr_on = s.cp_grad_corr_on;
     cp_opt       = s.cp_opt; 
     cp_std       = s.cp_std;
     cp_wave      = s.cp_wave;
     cp_min_opt   = s.cp_min_opt;
     cp_min_cg    = s.cp_min_cg;
     cp_min_std   = s.cp_min_std;
     rhoRsubplanes= s.rhoRsubplanes;
     sizeX        = s.sizeX;
     sizeY        = s.sizeY;
     sizeZ        = s.sizeZ;
     nstates      = s.nstates;
     ees_eext_on  = s.ees_eext_on;
     ees_nloc_on  = s.ees_nloc_on; 
     ngrid_nloc_a = s.ngrid_nloc_a; 
     ngrid_nloc_b = s.ngrid_nloc_b; 
     ngrid_nloc_c = s.ngrid_nloc_c; 
     ngrid_eext_a = s.ngrid_eext_a; 
     ngrid_eext_b = s.ngrid_eext_b; 
     ngrid_eext_c = s.ngrid_eext_c; 
     ntime        = s.ntime;
     ibinary_opt  = s.ibinary_opt;
     ibinary_write_opt = s.ibinary_write_opt;
     nplane_x     = s.nplane_x;
     nchareG      = s.nchareG;
     nlines_tot   = s.nlines_tot;
     npts_tot     = s.npts_tot;
     nlines_max   = s.nlines_max;
     natm_tot     = s.natm_tot;
     natm_typ     = s.natm_typ;
     natm_nl      = s.natm_nl;
     numSfGrps    = s.numSfGrps;
     natm_nl_grp_max = s.natm_nl_grp_max;
     if(nplane_x==0 || nplane_x > sizeX){
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("Error in CPcharmParaInfo constructor\n");
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkExit();
     }//endif
     if(nchareG==0){
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("Error in CPcharmParaInfo constructor\n");
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkExit();
     }//endif
     nlines_tot_rho  = s.nlines_tot_rho;
     npts_tot_rho    = s.npts_tot_rho;
     nlines_max_rho  = s.nlines_max_rho;
     nlines_max_eext = s.nlines_max_eext;
     nplane_rho_x    = s.nplane_rho_x;
     nchareRhoG      = s.nchareRhoG;
     nchareRhoGEext  = s.nchareRhoGEext;
     npts_per_chareRhoG   = new int[nchareRhoG];
     nlines_per_chareRhoG = new int[nchareRhoG];
     nlines_per_chareRhoGEext = new int[nchareRhoGEext];
     lines_per_chareRhoG  = new double[nchareRhoG];
     pts_per_chareRhoG    = new double[nchareRhoG];
     if(nplane_rho_x==0 || nplane_rho_x > sizeX){
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("Error in CPcharmParaInfo constructor\n");
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkExit();
     }//endif
     if(nchareRhoG==0){
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
     index_tran_upack_rho = cmall_int_mat(0,nchareRhoG,0,nlines_max_rho,
                                          "cpcharmparainfo.h");
     for(int i=0;i<nchareRhoG;i++){
      for(int j=0;j<nlines_per_chareRhoG[i];j++){
        index_tran_upack_rho[i][j] = s.index_tran_upack_rho[i][j];
      }//endfor
     }//endfor

     for(int i=0;i<nchareRhoGEext;i++){
       nlines_per_chareRhoGEext[i] = s.nlines_per_chareRhoGEext[i];
     }//endfor
     index_tran_upack_eext=cmall_int_mat(0,nchareRhoGEext,0,nlines_max_rho,
                                         "cpcharmparainfo.h");
     for(int i=0;i<nchareRhoGEext;i++){
      for(int j=0;j<nlines_per_chareRhoGEext[i];j++){
        index_tran_upack_eext[i][j] = s.index_tran_upack_eext[i][j];
      }//endfor
     }//endfor

     sortedRunDescriptors = new CkVec<RunDescriptor> [nchareG];
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

     index_tran_upackNL = cmall_int_mat(0,nchareG,0,nlines_max,"cpcharmparainfo.h");
     for(int i=0;i<nchareG;i++){
      for(int j=0;j<nlines_per_chareG[i];j++){
        index_tran_upackNL[i][j] = s.index_tran_upackNL[i][j];
      }//endfor
     }//endfor

     RCommPkg = new RedundantCommPkg [nchareG]; 
     for(int i=0;i<nchareG;i++){RCommPkg[i].Init(&s.RCommPkg[i]);}

     nlIters   = s.nlIters;
     ioff_zmat = new int [nlIters];
     nmem_zmat = new int [nlIters];;
     for(int i =0;i<nlIters;i++){
       ioff_zmat[i] = s.ioff_zmat[i];
       nmem_zmat[i] = s.nmem_zmat[i];
     }//endif

     if(rhoRsubplanes>1){
       index_tran_upack_rho_y  = cmall_itens3(0,nchareRhoG,0,rhoRsubplanes,
                                              0,nlines_max_rho,"util.C");
       index_tran_pack_rho_y   = cmall_itens3(0,nchareRhoG,0,rhoRsubplanes,
                                              0,nlines_max_rho,"util.C");
       nline_send_rho_y        = cmall_int_mat(0,nchareRhoG,0,rhoRsubplanes,"util.C");
       index_tran_upack_eext_y = cmall_itens3(0,nchareRhoGEext,0,rhoRsubplanes,
                                              0,nlines_max_eext,"util.C");
       index_tran_upack_eext_ys= cmall_itens3(0,nchareRhoGEext,0,rhoRsubplanes,
                                              0,nlines_max_eext,"util.C");
       index_tran_pack_eext_y  = cmall_itens3(0,nchareRhoGEext,0,rhoRsubplanes,
                                              0,nlines_max_eext,"util.C");
       index_tran_pack_eext_ys = cmall_itens3(0,nchareRhoGEext,0,rhoRsubplanes,
                                              0,nlines_max_eext,"util.C");
       nline_send_eext_y       = cmall_int_mat(0,nchareRhoGEext,0,rhoRsubplanes,"util.C");
       for(int igrp=0,i=0;igrp<nchareRhoG;igrp++){
          for(int ic=0;ic<rhoRsubplanes;ic++){
            nline_send_rho_y[igrp][ic] = s.nline_send_rho_y[igrp][ic];
            for(int jc=0;jc<nline_send_rho_y[igrp][ic];jc++){
              index_tran_upack_rho_y[igrp][ic][jc] = s.index_tran_upack_rho_y[igrp][ic][jc];
              index_tran_pack_rho_y[igrp][ic][jc] = s.index_tran_pack_rho_y[igrp][ic][jc];
	    }//endfor
	  }//endfor
       }//endfor
       for(int igrp=0,i=0;igrp<nchareRhoGEext;igrp++){
        for(int ic=0;ic<rhoRsubplanes;ic++){
         nline_send_eext_y[igrp][ic] = s.nline_send_eext_y[igrp][ic];
         for(int jc=0;jc<nline_send_eext_y[igrp][ic];jc++){
          index_tran_upack_eext_y[igrp][ic][jc]= s.index_tran_upack_eext_y[igrp][ic][jc];
          index_tran_upack_eext_ys[igrp][ic][jc]= s.index_tran_upack_eext_ys[igrp][ic][jc];
          index_tran_pack_eext_y[igrp][ic][jc] = s.index_tran_pack_eext_y[igrp][ic][jc];
          index_tran_pack_eext_ys[igrp][ic][jc] = s.index_tran_pack_eext_ys[igrp][ic][jc];
         }//endfor
	}//endfor
       }//endfor
     }//endif

     LBTurnInstrumentOff();
   }//end constructor
//=============================================================================

//=============================================================================
   CPcharmParaInfo() {
//=============================================================================
#ifdef _CP_DEBUG_PARAINFO_VERBOSE_
   CkPrintf("CPcharmParaInfo constructor\n");
#endif
       rhoRsubplanes = 1;
       lines_per_chareG=NULL; 
       pts_per_chareG=NULL;
       nlines_per_chareG=NULL; 
       npts_per_chareG=NULL;
       index_output_off=NULL;

       lines_per_chareRhoG=NULL; 
       pts_per_chareRhoG=NULL;
       nlines_per_chareRhoG=NULL; 
       nlines_per_chareRhoGEext=NULL; 
       npts_per_chareRhoG=NULL;
       RCommPkg  = NULL;

       ioff_zmat = NULL;
       nmem_zmat = NULL;

       index_tran_upack_rho_y   = NULL;
       index_tran_pack_rho_y    = NULL;
       nline_send_rho_y         = NULL;
       index_tran_upack_eext_y  = NULL;
       index_tran_upack_eext_ys = NULL;
       index_tran_pack_eext_y   = NULL;
       index_tran_pack_eext_ys  = NULL;
       nline_send_eext_y        = NULL;

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
      delete [] nlines_per_chareRhoGEext; nlines_per_chareRhoGEext = NULL;
      delete []  pts_per_chareRhoG;    pts_per_chareRhoG   = NULL;
      delete [] npts_per_chareRhoG;   npts_per_chareRhoG   = NULL;
      delete [] RCommPkg;             RCommPkg= NULL;
      delete [] ioff_zmat;
      delete [] nmem_zmat;

      cfree_int_mat(index_tran_upack,0,nchareG,0,nlines_max);
      cfree_int_mat(index_tran_upackNL,0,nchareG,0,nlines_max);
      cfree_int_mat(index_tran_upack_rho,0,nchareRhoG,0,nlines_max_rho);
      cfree_int_mat(index_tran_upack_eext,0,nchareRhoGEext,0,nlines_max_rho);

     if(rhoRsubplanes>1){
      cfree_itens3(index_tran_upack_rho_y,0,nchareRhoG,0,rhoRsubplanes,0,nlines_max_rho);
      cfree_itens3(index_tran_pack_rho_y,0,nchareRhoG,0,rhoRsubplanes,0,nlines_max_rho);
      cfree_itens3(index_tran_upack_eext_y,0,nchareRhoGEext,0,rhoRsubplanes,
		                           0,nlines_max_eext);
      cfree_itens3(index_tran_upack_eext_ys,0,nchareRhoGEext,0,rhoRsubplanes,
		                           0,nlines_max_eext);
      cfree_itens3(index_tran_pack_eext_y,0,nchareRhoGEext,0,rhoRsubplanes,
  		                          0,nlines_max_eext);
      cfree_itens3(index_tran_pack_eext_ys,0,nchareRhoGEext,0,rhoRsubplanes,
  		                          0,nlines_max_eext);
      cfree_int_mat(nline_send_rho_y,0,nchareRhoG,0,rhoRsubplanes);
      cfree_int_mat(nline_send_eext_y,0,nchareRhoGEext,0,rhoRsubplanes);
     }//endif
      
  }//end destructor
//=============================================================================

//=============================================================================
  void pup(PUP::er &p){
//=============================================================================
#ifdef _CP_DEBUG_PARAINFO_VERBOSE_
     CkPrintf("CPcharmParaInfo pup\n");
#endif
      p|vol;        p|dt;         p|tol_norb;   p|tol_cp_min; p|tol_cp_dyn;
      p|fftopt;
      p|ndump_frq;  p|istart_typ_cp; p|cp_grad_corr_on;
      p|cp_opt;     p|cp_std;     p|cp_wave;
      p|cp_min_opt; p|cp_min_std; p|cp_min_cg; p|rhoRsubplanes;
      p|sizeX;      p|sizeY;      p|sizeZ;  
      p|ees_eext_on;    p|ees_nloc_on;
      p|ngrid_nloc_a;  p|ngrid_nloc_b;   p|ngrid_nloc_c;
      p|ngrid_eext_a;  p|ngrid_eext_b;   p|ngrid_eext_c;
      p|nplane_x;   p|nchareG;    p|natm_tot;    p|natm_nl;
      p|nstates;    p|ntime;      p|ibinary_opt; p|ibinary_write_opt;
      p|natm_typ;   p|natm_nl;    p|numSfGrps; p|nlIters; 
      p|nmem_zmat_tot; p|nmem_zmat_max;
      p|natm_nl_grp_max;  p|nlines_tot; p|npts_tot;
      p|nlines_max; p|nlines_max_eext;
      p|nplane_rho_x; p|nchareRhoG; p|nchareRhoGEext; 
      p|nlines_max_rho; p|nlines_tot_rho; p|npts_tot_rho;
      if(p.isUnpacking()) {
        lines_per_chareRhoG      = new double[nchareRhoG];
        pts_per_chareRhoG        = new double[nchareRhoG];
        nlines_per_chareRhoG     = new int[nchareRhoG];
        nlines_per_chareRhoGEext = new int[nchareRhoGEext];
        npts_per_chareRhoG       = new int[nchareRhoG];
#ifdef _CP_DEBUG_PARAINFO_VERBOSE_
	CkPrintf("nchareRhoG %d nlines_max_rho %d sizeX %d\n",
             nchareRhoG, nlines_max_rho, sizeX);
#endif
        index_tran_upack_rho = cmall_int_mat(0,nchareRhoG,0,nlines_max_rho,
                                             "cpcharmparainfo.h");
        index_tran_upack_eext= cmall_int_mat(0,nchareRhoGEext,0,nlines_max_rho,
                                             "cpcharmparainfo.h");
        RhosortedRunDescriptors = new CkVec<RunDescriptor> [nchareRhoG];
        if(rhoRsubplanes>1){
          index_tran_upack_rho_y  = cmall_itens3(0,nchareRhoG,0,rhoRsubplanes,
                                                 0,nlines_max_rho,"cpcharmparainfo.h");
          index_tran_pack_rho_y   = cmall_itens3(0,nchareRhoG,0,rhoRsubplanes,
                                                 0,nlines_max_rho,"cpcharmparainfo.h");
          nline_send_rho_y        = cmall_int_mat(0,nchareRhoG,0,rhoRsubplanes,
                                                  "cpcharmparainfo.h");
          index_tran_upack_eext_y = cmall_itens3(0,nchareRhoGEext,0,rhoRsubplanes,
                                                 0,nlines_max_eext,"cpcharmparainfo.h");
          index_tran_upack_eext_ys= cmall_itens3(0,nchareRhoGEext,0,rhoRsubplanes,
                                                 0,nlines_max_eext,"cpcharmparainfo.h");
          index_tran_pack_eext_y  = cmall_itens3(0,nchareRhoGEext,0,rhoRsubplanes,
                                                 0,nlines_max_eext,"cpcharmparainfo.h");
          index_tran_pack_eext_ys = cmall_itens3(0,nchareRhoGEext,0,rhoRsubplanes,
                                                 0,nlines_max_eext,"cpcharmparainfo.h");
          nline_send_eext_y       = cmall_int_mat(0,nchareRhoGEext,0,rhoRsubplanes,
                                                  "cpcharmparainfo.h");
        }//endif
      }//enddif : unpacking
#ifdef _CP_DEBUG_PARAINFO_VERBOSE_
      CkPrintf("CPcharmParaInfo pup 2 \n");
#endif
      PUParray(p,lines_per_chareRhoG, nchareRhoG);
      PUParray(p,nlines_per_chareRhoG, nchareRhoG);
      PUParray(p,nlines_per_chareRhoGEext, nchareRhoGEext);
      PUParray(p,pts_per_chareRhoG, nchareRhoG);
      PUParray(p,npts_per_chareRhoG, nchareRhoG);
      for(int igrp=0;igrp<nchareRhoG;igrp++){
	p|RhosortedRunDescriptors[igrp];
      }
      for(int igrp=0;igrp<nchareRhoG;igrp++){
	PUParray(p,index_tran_upack_rho[igrp],nlines_per_chareRhoG[igrp]);
      }
      for(int igrp=0;igrp<nchareRhoGEext;igrp++){
	PUParray(p,index_tran_upack_eext[igrp],nlines_per_chareRhoGEext[igrp]);
      }
      if(rhoRsubplanes>1){
        for(int igrp=0;igrp<nchareRhoG;igrp++){
          PUParray(p,nline_send_rho_y[igrp],rhoRsubplanes);
	}
        for(int igrp=0;igrp<nchareRhoGEext;igrp++){
          PUParray(p,nline_send_eext_y[igrp],rhoRsubplanes);
	}
        for(int igrp=0,i=0;igrp<nchareRhoG;igrp++){
          for(int ic=0;ic<rhoRsubplanes;ic++){
            PUParray(p,index_tran_pack_rho_y[igrp][ic],nline_send_rho_y[igrp][ic]);
            PUParray(p,index_tran_upack_rho_y[igrp][ic],nline_send_rho_y[igrp][ic]);
	  }
	}
        for(int igrp=0,i=0;igrp<nchareRhoGEext;igrp++){
          for(int ic=0;ic<rhoRsubplanes;ic++){
            PUParray(p,index_tran_pack_eext_y[igrp][ic],nline_send_eext_y[igrp][ic]);
            PUParray(p,index_tran_pack_eext_ys[igrp][ic],nline_send_eext_y[igrp][ic]);
            PUParray(p,index_tran_upack_eext_y[igrp][ic],nline_send_eext_y[igrp][ic]);
            PUParray(p,index_tran_upack_eext_ys[igrp][ic],nline_send_eext_y[igrp][ic]);
	  }
	}
      }//endif : subplanes are in use.
#ifdef _CP_DEBUG_PARAINFO_VERBOSE_
      CkPrintf("CPcharmParaInfo pup 3 \n");
#endif
      if(p.isUnpacking()) {
        ioff_zmat         = new int [nlIters];
        nmem_zmat         = new int [nlIters];;
        lines_per_chareG  = new double[nchareG];
        pts_per_chareG    = new double[nchareG];
        nlines_per_chareG = new int[nchareG];
        npts_per_chareG   = new int[nchareG];
        index_output_off  = new int[nchareG];
        index_tran_upack  = cmall_int_mat(0,nchareG,0,nlines_max,"cpcharmparainfo.h");
        index_tran_upackNL= cmall_int_mat(0,nchareG,0,nlines_max,"cpcharmparainfo.h");
        sortedRunDescriptors = new CkVec<RunDescriptor> [nchareG];
      }//endif : unpacking
#ifdef _CP_DEBUG_PARAINFO_VERBOSE_
      CkPrintf("CPcharmParaInfo pup 4 \n");
#endif
      p(ioff_zmat,nlIters);
      p(nmem_zmat,nlIters);
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
      for(int igrp=0;igrp<nchareG;igrp++){
	p(index_tran_upackNL[igrp],nlines_per_chareG[igrp]);
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
      PUParray(p,RCommPkg,nchareG);

#ifdef _CP_DEBUG_PARAINFO_VERBOSE_
     CkPrintf("end CPcharmParaInfo pup\n");
#endif

   };
//----------------------------------------------------------------------------
   }; // end class
//=============================================================================


//=============================================================================
PUPmarshall(CPcharmParaInfo);
//=============================================================================


//=============================================================================
#endif
