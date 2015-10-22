// do matrix inversion using iterative methods
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "include/ckcomplex.h"
#include "class_defs/matrix.h"
#include "class_defs/usrinput.h"
#include "iter_invmtrx.h"
#include "util.h"

using namespace std;

void iter_invmtrx(CMATRIX* eps, USRINPUT usrin, int N){

    // allocate matrix.
    // myGEMM takes 1D array, so it has to be complex *, not complex ** //
    complex *X, *X1, *M1, *A;
    A = new complex[N*N]; // input matrix (eps->m)
    X = new complex[N*N];
    X1 = new complex[N*N];
    M1 = new complex[N*N];
  
    // calculate matrix X  : row major
    // X = m * m^T
    // the purpose of this matrix is find parameter alpha
    // we will destruct X once we find alpha
    for (int i=0; i<N; i++){
	for (int j=0; j<N; j++){
	    X[i*N+j] = 0;
	    for (int k=0; k<N; k++){
		X[i*N+j] += eps->get(i,k)*eps->get(j,k).conj();
	    }
	}
    }
    
    // find alpha 
    double R, alpha;
    alpha = 0;
    for(int i=0; i<N; i++){
        R=0;
        for(int j=0; j<N; j++){
            R += abs(X[i*N+j]);
        }
        if ( R > alpha ) {
            alpha = R;
        }
    }
    alpha = 1/alpha;


    // X matrix :: colum major (to call lapack routine)
    // initial guess. X is now different from previous
    // X = alpha * m^T   (warning!! column major)
    for (int j=0; j<N; j++) {
        for (int i=0; i<N; i++) {
            X[i+j*N] = alpha * eps->get(j,i).conj();
        }
    }

    // assign A = eps->m (warning!! A is column major)
    for (int i=0; i<N; i++){
	for (int j=0; j<N; j++){
	    A[i+j*N] = eps->get(i,j);
	}
    }
    

    // matrix multiplication

    const int maxiter = usrin.iter_maxiter; // maximum number of iteration (user provided)
    double convg = usrin.iter_convg;        // convergence (user provided)
    
    complex compl_two(2.0, 0.0);   

    double resid;
    int niter;                   // number of iteration performed
    char transformT = 'C';       // conjugate transpose (Hermitian conjugate)
    char transform = 'N';        // Do nothing for this matrix
    double Lalpha = double(1.0); // Scale A*B by this scalar factor
    double Lbeta = double(0.0);  // Scale initial valuce of C by this factor

    for(int iter=0; iter < maxiter; iter++){

#ifdef USE_LAPACK
	// Call GEMM for matrix multiplication
	// Step 1: M1 = 2I - A * X
	//     Step 1-(a) : M1 = -A * X
	Lalpha = double(-1.0);
	myGEMM(&transform, &transform, &N, &N, &N, &Lalpha, A, &N, X, &N, &Lbeta, M1, &N);
	//     Step 1-(b) : M1 = 2I-M1 
	for(int i=0; i<N; i++){
	    M1[i*N+i] += compl_two;
	}

	// Step 2: X1 = X * M1
	Lalpha = double(1.0);
	myGEMM(&transform, &transform, &N, &N, &N, &Lalpha, X, &N, M1, &N, &Lbeta, X1, &N);

#else

	// fill this
	// Step 1: M1 = 2I - A * X


	// Step 2: X1 = X * M1
	

#endif


        // check if it can stop.
        // compare X vs X1

        resid = convergence_check(X1, X, N);        

        // Step 3: X = X1
        for(int i=0; i<N*N; i++){
            X[i] = X1[i];
        }//end step 3
        
        
        if (resid <= convg){
            niter = iter+1;
            break;
	}
        
    }//end iteration

    delete[] A;
    delete[] X1;
    delete[] M1;

    
    // epsilon becomes epsilon inverse!
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            eps->get(i,j) = X[i+j*N];
        }
    }//end step 3

    delete[] X;
    
    // print out the number of iteration performed
    cout << " ------------------------------------ " << endl;
    cout << " number of iteration performed : " << niter << endl; 
    cout << " ------------------------------------ " << endl;
    
}


double convergence_check(complex *M1, complex *M0, int N){

    // find the largest element of abs(M0_i,j - M1_i,j)    
    double Rmax=0;  // the largest element
    double tmp;

    for(int i=0; i<N*N; i++){
        tmp = abs( M0[i] - M1[i] );
        if( tmp > Rmax ){ Rmax = tmp; }
    }
    return Rmax;
}
