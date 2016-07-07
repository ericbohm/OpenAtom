// define class for 2-dimensional array
#ifndef MATRIX_H
#define MATRIX_H

#include "../include/ckcomplex.h"
#include "../util.h"
#include "gspace.h"

using namespace std;

// complex variable matrix (rank2 array)

class CMATRIX{

  private:
     
  public:
    int thisIndex;  // index if array
    int nrow, ncol; // number of rows and colums
    int ndata;      // number of data ( = nrow * ncol )
    complex *m;
    bool T; // if the matrix is transposed or not
   
    // constructor  
    CMATRIX(){};
    CMATRIX(int _nrow, int _ncol){
        nrow = _nrow;
        ncol = _ncol;
        ndata = _nrow * _ncol;
	
        // allocate memory for matrix m
        m = new complex [ndata];

        //  initialization  
        for (int i=0; i<ndata; i++)
            m[i] = (0, 0);
    }

    // access function
    inline complex& get(int i, int j){
        T = false;
	if(!T)
	    return m[i*ncol+j];
	if(T)
	    return m[i+j*nrow];
    }

    // in place transpose of the matrix
    void transpose();

    // in place conjugate tranpose of the matrix
    void ctranspose();

    //printint matrix elements
    void printRow(int, char*);
    void printColumn(int, char*);
    void printRowEps(int, char*, GSPACE*);
    void printDiag(char*);
    void printIJdiff(int,char*);
    void printallG(char*,GSPACE*);

};




#endif
