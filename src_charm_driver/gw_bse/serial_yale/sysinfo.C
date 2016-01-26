// member functions for sysinfo class
// it does read input files and calculate q vectors and volume of the simulation cell
#include "class_defs/sysinfo.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "mtrxop_3x3.h"
#include "util.h"

using namespace std;

void SYSINFO::read_sysinfo(){
    
    ifstream fp("sysinfo.dat");       // system information file from quantum espresso
    
    if(!fp){                            // if fails to open, exit program
        cerr << "Failed to open file.\n";
        exit(EXIT_FAILURE);
    }
    
    if(fp){
        
        // read system information
        
        // lattice constant 
        fp >> alat;


        // lattice vectors
	// atomic unit
        fp >> a1[0] >> a1[1] >> a1[2];
        fp >> a2[0] >> a2[1] >> a2[2];
        fp >> a3[0] >> a3[1] >> a3[2];

        // reciprocal lattice vectors 
	// this has unit of  1/(2pi/alat) a.u.^-1
        fp >> b1[0] >> b1[1] >> b1[2];
        fp >> b2[0] >> b2[1] >> b2[2];
        fp >> b3[0] >> b3[1] >> b3[2];
        
        // nspin, nkpt, nstate 
        fp >> nspin >> nkpt >> nstate;
        
        
        // number of planewaves at each k points 
	npwk = new int [nkpt];
        kwt = new double [nkpt];
        kvec = new double *[nkpt];
	kcart = new double *[nkpt];
        for (int i=0; i < nkpt ; i++) {
            kcart[i] = new double [3];
        }
      
        // number of planwvaves in each k point
        for (int i=0; i<nkpt; i++) {
           fp >> npwk[i];
        }
        
        /* kvectors */
        // in cartesian unit (2pi/alat factor is not multiplied in quantum espresso)
        for (int i=0; i<nkpt; i++) {
            fp >> kcart[i][0] >> kcart[i][1] >> kcart[i][2] >> kwt[i];
        }

	// dense fft grid size
	fp >> nfftDen[0] >> nfftDen[1] >> nfftDen[2];

    }
    
    fp.close();
}


// to calculate q vectors, it is better handle k vectors with crystal coordinates
// kvec = b^-1 kcart  (where b is the reciprocal lattice 3x3 matrix)              
void SYSINFO::cartesian_to_crystal(){
    
    kvec = new double *[nkpt];
    for (int i=0; i<nkpt; i++) {
        kvec[i] = new double [3];
    }
    
    double b[3][3] = {{b1[0],b2[0],b3[0]},{b1[1],b2[1],b3[1]},{b1[2],b2[2],b3[2]}};
    double invb[3][3];
    inverse( b, invb );
    
    
    for (int ik=0; ik<nkpt; ik++) {
        matvec3( invb, kcart[ik], kvec[ik] );
    }
    
}


void SYSINFO::calc_qvec(){
    
    qvec = new double *[nkpt];
    for (int i=0; i <nkpt; i++){
        qvec[i] = new double [3];
    }
    
    for (int i=0; i<nkpt; i++){
        for (int j=0; j<3; j++) {
            qvec[i][j] = kvec[i][j] - kvec[0][j];
        }
        
    }

}



/* calculate system volme */
void SYSINFO::calc_vol(){

    double a[3][3];
    double m[3][3];
    double det, del;
        
    for (int i; i<3; i++){
        m[0][i] = a1[i];
        m[1][i] = a2[i];
        m[2][i] = a3[i];
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

    del = 1.0e-5;

    if (abs(det) < del ){
        Die("cannot compute volume of the cell!!");
    }

    vol = det;

    printf("-----------------------------------------\n");
    printf("volume of the simulation cell: %f\n",vol);
    printf("-----------------------------------------\n");
        
}


