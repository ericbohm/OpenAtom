#ifndef _PMEGSPACEDATA_H_
#define _PMEGSPACEDATA_H_

#include "leanMD.h"
#include "BasePMEGSpaceData.h"

#include "../include/pentium_par/standard_include.h"

class PMEGSpaceData : public BasePMEGSpaceData {
  public:
    PMEGSpaceData() {}
    ~PMEGSpaceData ();
    
    void initialize (int sy_, int nPlanes_, double alp_ewd_, 
                     double vol_, double* hmati_, int n_interp_);

    void compute_recip_energy (double* hmati, StepOutput& out);
    
    complex* getQGrid ();
    
    void pup(PUP::er& p);
    
    static void calcGridSize (double* _hmati, double _alpha_ewald,
                              int& _ngrid_a, int& _ngrid_b, int& _ngrid_c);
    
  private:
    double   alp_ewd;
    double   vol;         // simulation box volume
    double   ecut_recip;
  
    int      n_interp;
    
    int      ngrid_a;
    int      ngrid_b;
    int      ngrid_c;
    
    int      ka_max;
    int      kb_max;
    int      kc_max;
    
    int      nktot;
    
    int      qgridSize;
    
    int      sy;       // start plane index
    int      nPlanes;  // [sy, sy+nPlanes-1] planes belong to me
    
    complex* qgrid_gspace;  
    int*     ka;         // size [nktot]
    int*     kb;         // size [nktot]
    int*     kc;         // size [nktot]
  
    double*  bweight_tot;// size [nktot]
    
    double*  aj;         // size [interpOrder]
    double*  rn;         // size [interpOrder]
    double*  rn1;        // size [interpOrder]
    
    void count_recip (double* hmati);
    void fill_recip (double* hmati);
    void set_pme_wght (int pme_b_opt, double* bfact_r, double* bfact_i);
    void get_bspline_wght1d(int ngrid, double *mn_k, double *uk, double *bden_r,
                            double *bden_i, double *bweight);
    void clear ();
};

#endif
