#include <cstdlib>
#include "class_defs/matrix.h"
#include "print_util.h"

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



// for debugging purpose... 
// print row for P matrix 
void CMATRIX::printRow(int row, char* fileName){
    FILE* fp = fopen(fileName, "w");
    printf("print row %d in CMATRIX\n", row);

    int rowGidx[3];
    PdataIndexG(row,rowGidx);
    
    // loop over column
    fprintf(fp,"row G index   column G index\n");
    for (int ic=0; ic<ncol; ic++){
	complex data;
	data = get(row, ic);
	int colGidx[3];
	PdataIndexG(ic,colGidx);
	fprintf(fp, "(%d, %d, %d)  (%d, %d, %d)   %lg  %lg \n", rowGidx[0], rowGidx[1], rowGidx[2],
		colGidx[0], colGidx[1], colGidx[2], data.re, data.im);
    }
    fclose(fp);
    
}

// print row for epsilon matrix
void CMATRIX::printRowEps(int row, char* fileName, GSPACE* geps){
    FILE* fp = fopen(fileName, "w");
    printf("print row %d in CMATRIX\n", row);

    int rowGidx[3];
    rowGidx[0] = geps->ig[row];
    rowGidx[1] = geps->jg[row];
    rowGidx[2] = geps->kg[row];
	
    
    // loop over column
    fprintf(fp,"row G index   column G index\n");
    for (int ic=0; ic<ncol; ic++){
	complex data;
	data = get(row, ic);
	int colGidx[3];
	colGidx[0] = geps->ig[ic];
        colGidx[1] = geps->jg[ic];
        colGidx[2] = geps->kg[ic];
	fprintf(fp, "(%d, %d, %d)  (%d, %d, %d)   %lg  %lg \n", rowGidx[0], rowGidx[1], rowGidx[2],
		colGidx[0], colGidx[1], colGidx[2], data.re, data.im);
    }
    fclose(fp);
    
}







// print column 
void CMATRIX::printColumn(int column, char* fileName){
    FILE* fp = fopen(fileName, "w");
    printf("print column %d in CMATRIX\n", column);

    int colGidx[3];
    PdataIndexG(column,colGidx);
    
    // loop over column
    fprintf(fp,"row G index   column G index\n");
    for (int ir=0; ir<nrow; ir++){
	complex data;
	data = get(ir, column);
	int rowGidx[3];
	PdataIndexG(ir,rowGidx);
	fprintf(fp, "(%d, %d, %d)  (%d, %d, %d)   %lg  %lg \n", rowGidx[0], rowGidx[1], rowGidx[2],
		colGidx[0], colGidx[1], colGidx[2], data.re, data.im);
    }
    fclose(fp);
}
