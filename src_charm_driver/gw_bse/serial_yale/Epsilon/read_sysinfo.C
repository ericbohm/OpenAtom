// read system information
#include "class_defs/sysinfo.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "read_sysinfo.h"
#include "util.h"

using namespace std;

// function definitions

void read_sysinfo(SYSINFO &sys){
    
    ifstream fp("sysinfo.dat");       // system information file from quantum espresso
    
    if(!fp){                            // if fails to open, exit program
        cerr << "Failed to open file.\n";
        exit(EXIT_FAILURE);
    }
    
    if(fp){
        
        // read system information
        
        // lattice constant 
        fp >> sys.alat;


        // lattice vectors
	// atomic unit
        fp >> sys.a1[0] >> sys.a1[1] >> sys.a1[2];
        fp >> sys.a2[0] >> sys.a2[1] >> sys.a2[2];
        fp >> sys.a3[0] >> sys.a3[1] >> sys.a3[2];

        // reciprocal lattice vectors 
	// this has unit of  1/(2pi/alat) a.u.^-1
        fp >> sys.b1[0] >> sys.b1[1] >> sys.b1[2];
        fp >> sys.b2[0] >> sys.b2[1] >> sys.b2[2];
        fp >> sys.b3[0] >> sys.b3[1] >> sys.b3[2];
        
        // nspin, nkpt, nstate 
        fp >> sys.nspin >> sys.nkpt >> sys.nstate;
        
        
        // number of planewaves at each k points 
	sys.npwk = new int [sys.nkpt];
        sys.kwt = new double [sys.nkpt];
        sys.kvec = new double *[sys.nkpt];
	sys.kcart = new double *[sys.nkpt];
        for (int i=0; i < sys.nkpt ; i++) {
            sys.kcart[i] = new double [3];
        }
      
        // number of planwvaves in each k point
        for (int i=0; i<sys.nkpt; i++) {
           fp >> sys.npwk[i];
        }
        
        /* kvectors */
        // in cartesian unit (2pi/alat factor is not multiplied in quantum espresso)
        for (int i=0; i<sys.nkpt; i++) {
            fp >> sys.kcart[i][0] >> sys.kcart[i][1] >> sys.kcart[i][2] >> sys.kwt[i];
        }

	// dense fft grid size
	fp >> sys.nfftDen[0] >> sys.nfftDen[1] >> sys.nfftDen[2];

    }
    
    fp.close();
    
// for debugging
sys.nkpt = 1;    
    // calculate total number of k points in all spin channel
    sys.get_ntotkpt(sys.nspin, sys.nkpt);
    
    // to calculate q vectors, it is better handle k vectors with crystal coordinates
    // kvec = b^-1 kcart  (where b is the reciprocal lattice 3x3 matrix) 
    cartesian_to_crystal(sys);

    // calculate q vectors
    calc_qvec(sys);
    
    // calculate volume
    calc_vol(sys);

    
}




void cartesian_to_crystal(SYSINFO &sys){
    
    sys.kvec = new double *[sys.nkpt];
    for (int i=0; i<sys.nkpt; i++) {
        sys.kvec[i] = new double [3];
    }
    
    double b[3][3] = {{sys.b1[0],sys.b2[0],sys.b3[0]},{sys.b1[1],sys.b2[1],sys.b3[1]},{sys.b1[2],sys.b2[2],sys.b3[2]}};
    double invb[3][3];
    inverse( b, invb );
    
    
    for (int ik=0; ik<sys.nkpt; ik++) {
        matvec3( invb, sys.kcart[ik], sys.kvec[ik] );
    }
    
}




void calc_qvec(SYSINFO &sys){
    
    sys.qvec = new double *[sys.nkpt];
    for (int i=0; i <sys.nkpt; i++){
        sys.qvec[i] = new double [3];
    }
    
    for (int i=0; i<sys.nkpt; i++){
        for (int j=0; j<3; j++) {
            sys.qvec[i][j] = sys.kvec[i][j]-sys.kvec[0][j];
        }
        
    }
    
}



/* calculate system volme */
void calc_vol(SYSINFO &sys){

    double a[3][3];
    double m[3][3];
    double det, del;
        
    for (int i; i<3; i++){
        m[0][i] = sys.a1[i];
        m[1][i] = sys.a2[i];
        m[2][i] = sys.a3[i];
    }

    /* compute matrix of cofactors */
    a[0][0] =  m[1][1]*m[2][2] - m[1][2]*m[2][1];
    a[1][0] = -m[1][0]*m[2][2] + m[1][2]*m[2][0];
    a[2][0] =  m[1][0]*m[2][1] - m[1][1]*m[2][0];
    a[0][1] = -m[0][1]*m[2][2] + m[0][2]*m[2][1];
    a[1][1] =  m[0][0]*m[2][2] - m[0][2]*m[2][0];
    a[2][1] = -m[0][0]*m[2][1] + m[0][1]*m[2][0];
    a[0][2] =  m[0][1]*m[1][2] - m[0][2]*m[1][1];
    a[1][2] = -m[0][0]*m[1][2] + m[0][2]*m[1][0];
    a[2][2] =  m[0][0]*m[1][1] - m[0][1]*m[1][0];

    det = m[0][0]*a[0][0] + m[0][1]*a[1][0] + m[0][2]*a[2][0];

    del = 1.0d-5;

    if (abs(det) < del ){
        Die("cannot compute volume of the cell!!");
    }

    sys.vol = det;

    printf("-----------------------------------------\n");
    printf("volume of the simulation cell: %f\n",sys.vol);
    printf(" \n");
        
}


