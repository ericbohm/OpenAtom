#ifndef _eesDataClass_h_
#define _eesDataClass_h_

#include "../include/RunDescriptor.h"

//============================================================================
// Data class
//============================================================================
class RPPDATA {
 public:
   int index;
   int ngrid_a,ngrid_b,ngrid_c;   // fft grid size
   int n_interp;                  // ees interpolation order
   int natm;                      // non-local atms
   int *plane_index;              // lth: natm : atm is on plane? 1/0
   int **igrid;                   // lth: natm*ninterp^2 : where is atm
   double **mn;                   // lth: natm*ninterp^2 : b-spline
   double **dmn_x,**dmn_y,**dmn_z;// lth: natm*ninterp^2 : nabla (b-spline)
   RPPDATA(){};   
  ~RPPDATA(){};
   void init(int );
};
//============================================================================


//============================================================================
// Data class
//============================================================================
class RHORHARTDATA {
 public:
   int index;
   int ngrid_a,ngrid_b,ngrid_c;   // fft grid size
   int n_interp;                  // ees interpolation order
   int natm;                      // non-local atms
   int *plane_index;              // lth: natm : atm is on plane? 1/0
   int **igrid;                   // lth: natm*ninterp^2 : where is atm
   double **mn;                   // lth: natm*ninterp^2 : b-spline
   double **dmn_x,**dmn_y,**dmn_z;// lth: natm*ninterp^2 : nabla (b-spline)
   RHORHARTDATA(){};   
  ~RHORHARTDATA(){};   
   void init(int );
};
//============================================================================


//============================================================================
// Data class
//============================================================================
class GPPDATA {
 public:
   int index;
   int ngrid_a,ngrid_b,ngrid_c;   //fft grid size
   int n_interp;                  //interpolation order
   int natm;                      //number of non-local atoms
   int ncoef;                     //number of g-pts in this collection
   double *b_re, *b_im;           //lth: ncoef : ees g-space scaling factor
   double *h_gspl;                //lth: ncoef : pseudo-spline lookup parameters 
   int *ind_gspl; 
   GPPDATA(){};   
  ~GPPDATA(){};   
   void init(int ,int , int *,int *,int *);
};
//============================================================================

//============================================================================
// Data class
//============================================================================
class RHOGHARTDATA {
 public:
   int index;
   int ngrid_a,ngrid_b,ngrid_c;   //fft grid size
   int n_interp;                  //interpolation order
   int natm;                      //number of non-local atoms
   int ncoef;                     //number of g-pts in this collection
   double *b_re, *b_im;           //lth: ncoef : ees g-space scaling factor
   double *h_gspl;int *ind_gspl;  //lth: ncoef : pseudo-spline lookup parameters 
   RHOGHARTDATA(){};   
  ~RHOGHARTDATA(){};   
   void init(int ,int , int *,int *,int *);
};
//============================================================================

//============================================================================
// Data class : Minimal copies of big redundant data sets
//============================================================================
class GCHAREPKG {
 public:
   int ihave_g000;
   int ind_g000;
   int ihave_kx0;
   int nkx0;   
   int nkx0_uni;
   int nkx0_red;
   int nkx0_zero;
   int kx0_strt;
   int kx0_end; 
   GCHAREPKG(){};   
  ~GCHAREPKG(){};   
};
//============================================================================


//============================================================================
// Data class : Minimal copies of big redundant data sets
//============================================================================
class GSPDATA {
 public:
   int index;                     //plane index
   int ngrid_a,ngrid_b,ngrid_c;   //fft grid size
   int ncoef;                     //number of g-pts in this collection
   int numLines;
   int numRuns;                   //2x the number of z-lines in collection
   int *ka, *kb, *kc;             //lth: ncoef : g-space
   double *g,*g2;
   RunDescriptor *runs;           //lth:numruns : k's in fft order
   GSPDATA(){};   
  ~GSPDATA(){};   
   void init(int);
   GCHAREPKG gCharePkg;
};
//============================================================================

//============================================================================
// Data class : Minimal copies of big redundant data sets
//============================================================================
class RHOGDATA {
 public:
   int index;                     //plane index
   int ngrid_a,ngrid_b,ngrid_c;   //fft grid size
   int ncoef;                     //number of g-pts in this collection
   int numRuns;                   //2x the number of z-lines in collection
   int *ka, *kb, *kc;             //lth: ncoef : g-space
   RunDescriptor *runs;           //lth:numruns : k's in fft order
   RHOGDATA(){};   
  ~RHOGDATA(){};   
   void init(int ,int , int, int, int, int);
};
//============================================================================
#endif
