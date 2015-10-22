#include <cstdlib>
#include "class_defs/matrix.h"

// member functions

// in place transpose
void CMATRIX::transpose(){

    if(nrow!=ncol){
	mymessage("The program has tried to transpose a matrix, but the matrix is not squre. I am not sure if this is what you want.");
    }
    else{
	for (int i=0; i<nrow; i++){
	    for (int j=0; j<i; j++){
		complex tmp = m[i*ncol+j];
		m[i*ncol+j] = m[j*ncol+i];
		m[j*ncol+i] = tmp;
	    }
	}
    }	
}	

void CMATRIX::ctranspose(){
    
    if(nrow!=ncol){
	mymessage("The program has tried to conjugate transpose a matrix, but the matrix is not squre. I am not sure if this is what you want.");
    }
    else{
	for (int i=0; i<nrow; i++){
	    for (int j=0; j<i; j++){
		complex tmp = m[i*ncol+j];
		m[i*ncol+j] = m[j*ncol+i].conj();
		m[j*ncol+i] = tmp.conj();
	    }
	}
    }	
}	

