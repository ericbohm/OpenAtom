#ifndef SYSINFO_H
#define SYSINFO_H

#include <cstdio>

class SYSINFO{
    
  public:

    int nspin;       // number of spin
    int nkpt;        // number of k points in each spin channel
    int nstate;      // total number of states
    int ntotkpt;     // nspin*nkpt
    int nocc;        // total number of valence (filled) state
    int nunocc;      // total number of conduction (empty) state
    
    double a1[3], a2[3], a3[3];  // lattice vectors (atomic unit)
    double b1[3], b2[3], b3[3];  // reciprocal lattice vectors
    double shift[3];             // shift vector for Vcoulb(q+G) when q=0 and G=0
    
    double vol;      // volume of the simulation box (atomic unit)
    double alat;     // lattice constant (atomic unit)
    
    double *kwt;
    double **kcart;  // cartesian coordinates
    double **qcart;  // cartesian coordinates
    double **kvec;   // crystal coordinates (reciprocal basis)
    double **qvec;   // crystal coordinates (reciprocal basis)
    
    int *npwk;       // number of planewaves at each k point
    int nfftDen[3];  // number of dense fft grid (for density)

    // member function declaration
    inline void get_ntotkpt(int _nspin, int _nkpt){
    
        ntotkpt = _nspin*_nkpt;

        printf("--------------------------------------\n");
        printf("Number of spin channels: %d\n", _nspin);
        printf("Number of k points in each spin channel: %d\n", _nkpt);
        printf("Number of total k points %d\n", ntotkpt);
        printf(" \n");
    }
        
};




#endif
