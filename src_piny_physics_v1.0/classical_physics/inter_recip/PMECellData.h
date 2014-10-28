#ifndef _PMECELLDATA_H_
#define _PMECELLDATA_H_

#include "leanMD.h"
#include "BasePMECellData.h"

#include "../include/pentium_par/standard_include.h"

class PMECellData : public BasePMECellData {
  //----------------
  public:

    PMECellData();
    ~PMECellData();

    void initialize (int nAtoms, int nLen, int interpOrder,
        int gridSizeX, int gridSizeY, int gridSizeZ);

    void reinitialize (int nAtoms, int nLen,int interpOrder,
        int gridSizeX, int gridSizeY, int gridSizeZ);

    void clear ();

    void pup (PUP::er &p);

    // some of arrays can be precomputed
    void precompute ();

    void computeMN (Vector* _pos, double* _hmati);

    void computeRSpaceGrid (Vector* _pos, double* _q, double* _hmati);

    double* getRSpaceGrid ();

    void calcForces (Vector* _force, double* _q, double* _hmati);

  private:
    int natm;
    int n_interp;
    int ngrid_a,ngrid_b,ngrid_c;
    int nlen;

    int *iatemp; // size [natm]
    int *ibtemp; // size [natm]
    int *ictemp; // size [natm]

    int **igrid_a; // size [n_interp*nlen]
    int **igrid_b; // size [n_interp*nlen]
    int **igrid_c; // size [n_interp*nlen]

    int **igrid_now; // size [n_interp*nlen]

    double *aj; // size [n_interp]
    double *rn; // size [n_interp]
    double *rn1; // size [n_interp]

    double *frac_a; // size [natm]
    double *frac_b; // size [natm]
    double *frac_c; // size [natm]

    double *qgrid_rspace; // size [ngrid_a*ngrid_b*ngrid_c]
    double **qgrid_now; // size [n_interp*nlen]
    double **ua; // size [n_interp*nlen]
    double **ub; // size [n_interp*nlen]
    double **uc; // size [n_interp*nlen]

    double **mn_a; // size [n_interp*nlen]
    double **mn_b; // size [n_interp*nlen]
    double **mn_c; // size [n_interp*nlen]

    double **dmn_a; // size [n_interp*nlen]
    double **dmn_b; // size [n_interp*nlen]
    double **dmn_c; // size [n_interp*nlen]
}; // PMECellData

#endif // _PMECELLDATA_H_
