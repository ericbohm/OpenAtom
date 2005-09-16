#ifndef _parainfo_h_
#define _parainfo_h_

#include "../include/RunDescriptor.h"
#include "../src_piny_physics_v1.0/friend_lib/proto_friend_lib_entry.h"

class CPcharmParaInfo {
 public:

   double vol;
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
   int **index_tran_upack;   // unpack indicies for each chareG

   //------------------------------------
   CPcharmParaInfo(CPcharmParaInfo &s){
     vol         = s.vol;
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
       CkPrintf("Error in CPcharmParaInfo constructor\n");
       CkExit();
     }//endif
     if(nchareG==0 || nchareG > sizeX || nchareG < nplane_x){
       CkPrintf("Error in CPcharmParaInfo constructor\n");
       CkExit();
     }//endif
     sortedRunDescriptors = new CkVec<RunDescriptor> [sizeX];
     for(int i=0;i<nchareG;i++){
       for(int j=0;j<s.sortedRunDescriptors[i].size();j++){
          sortedRunDescriptors[i].push_back(s.sortedRunDescriptors[i][j]);
       }//endfor
     }//endfor
     npts_per_chareG   = new int[sizeX];
     nlines_per_chareG = new int[sizeX];
     lines_per_chareG  = new double[sizeX];
     pts_per_chareG    = new double[sizeX];
     for(int i=0;i<nchareG;i++){
       nlines_per_chareG[i] = s.nlines_per_chareG[i];
       lines_per_chareG[i]  = s.lines_per_chareG[i];
       pts_per_chareG[i]    = s.pts_per_chareG[i];
       npts_per_chareG[i]   = s.npts_per_chareG[i];
     }//endfor
     for(int i=nchareG;i<sizeX;i++){
       pts_per_chareG[i]    = 0.0;
       lines_per_chareG[i]  = 0.0;
       npts_per_chareG[i]    = 0;
       nlines_per_chareG[i]  = 0;
     }//endfor
     index_tran_upack = cmall_int_mat(0,nchareG,0,nlines_max,"cpcharmparainfo.h");
     for(int i=0;i<nchareG;i++){
      for(int j=0;j<nlines_per_chareG[i];j++){
        index_tran_upack[i][j] = s.index_tran_upack[i][j];
      }//endfor
     }//endfor
     LBTurnInstrumentOff();
   }//end constructor

   //------------------------------------
   CPcharmParaInfo() {
       lines_per_chareG=NULL; 
       pts_per_chareG=NULL;
       nlines_per_chareG=NULL; 
       npts_per_chareG=NULL;
   }

   //------------------------------------
  ~CPcharmParaInfo() {
      delete []  lines_per_chareG;  lines_per_chareG = NULL;
      delete [] nlines_per_chareG; nlines_per_chareG = NULL;
      delete []  pts_per_chareG;    pts_per_chareG   = NULL;
      delete [] npts_per_chareG;   npts_per_chareG   = NULL;
   }

   //------------------------------------
   void pup(PUP::er &p){
      p|vol; 
      p|cp_opt;     p|cp_std;     p|cp_wave;
      p|cp_min_opt; p|cp_min_std; p|cp_min_cg;
      p|sizeX;      p|sizeY;      p|sizeZ;  
      p|nplane_x;   p|nchareG;   p|natm_tot;   p|natm_nl;
      p|nstates;    p|ntime;      p|ibinary_opt;
      p|natm_tot;   p|natm_nl;    p|numSfGrps;
      p|natm_nl_grp_max;   p|nlines_tot; p|npts_tot;
      p|nlines_max;

      int *itemp = new int[nchareG*nlines_max];
      memset(itemp,0,nchareG*nlines_max*sizeof(int));
      RunDescriptor *Desi  = new RunDescriptor[(2*nlines_tot)];

      if(p.isUnpacking()) {
        lines_per_chareG  = new double[sizeX];
        pts_per_chareG    = new double[sizeX];
        nlines_per_chareG = new int[sizeX];
        npts_per_chareG   = new int[sizeX];
        index_tran_upack = cmall_int_mat(0,nchareG,0,nlines_max,"cpcharmparainfo.h");
        sortedRunDescriptors = new CkVec<RunDescriptor> [sizeX];
      }else{
        int iii = 0;
        int jjj = 0;
        for(int igrp=0;igrp<nchareG;igrp++){
          for(int i=0;i<nlines_per_chareG[igrp];i++){
            itemp[iii] = index_tran_upack[igrp][i];
            iii++;
          }//endfor
          for(int i=0;i<2*nlines_per_chareG[igrp];i++){
            Desi[jjj]  = sortedRunDescriptors[igrp][i];
            jjj++;
	  }//endfor
        }//endfor
      }/*endif*/
      p(lines_per_chareG,nchareG);
      p(nlines_per_chareG,nchareG);
      p(pts_per_chareG,nchareG);
      p(npts_per_chareG,nchareG);
      p(npts_per_chareG,nchareG);
      for(int i=0;i<2*nlines_tot;i++){Desi[i].pup(p);}
      p(itemp,nlines_max*nchareG);

      if(p.isUnpacking()) {
        int iii = 0;
        int jjj = 0;
        for(int igrp=0;igrp<nchareG;igrp++){
          for(int i=0;i<nlines_per_chareG[igrp];i++){
            index_tran_upack[igrp][i] = itemp[iii];
            iii++;
          }//endfor
          for(int i=0;i<2*nlines_per_chareG[igrp];i++){
            sortedRunDescriptors[igrp].push_back(Desi[jjj]);
            jjj++;
	  }//endfor
        }//endfor
        for(int i=nchareG;i<sizeX;i++){
          pts_per_chareG[i]    = 0.0;
          lines_per_chareG[i]  = 0.0;
	  nlines_per_chareG[i] = 0;
          npts_per_chareG[i]   = 0;
        }//endfor

      }//endif : unpacking

      delete [] itemp;
      delete [] Desi;
   };

};

PUPmarshall(CPcharmParaInfo);
#endif
