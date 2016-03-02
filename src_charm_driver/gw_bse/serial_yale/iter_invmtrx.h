// LAPACK linking

// Define the GEMM macros that paircalc will use to invoke the appropriate matrix multiplys

#include "include/ckcomplex.h"
#include "include/mylapack.h"



double convergence_check(complex *, complex *, int);
