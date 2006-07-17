#ifndef _eesDataClass_h_
#define _eesDataClass_h_

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

#endif
