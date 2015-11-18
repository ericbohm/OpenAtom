#ifndef _diagonalizer_h_
#define _diagonalizer_h_
//#define INTEROP 1

//#include "preprocessor.h"
//#include "mpi.h"

template <typename T>
struct diagData_t {
  int numStatesOA;
  int orthoGrainSizeOA;
  int numOrthosPerDimOA;
  int totalOrthosOA;
  //int blockDim;
  //int dim;
  //bool elpaInit;
  //int myRow;
  //int myCol;
  //MPI_Fint mpiCommRows;
  //MPI_Fint mpiCommCols;

  //void *msg;
  T*** globalmatrix;
  T* plambda;
  T* rlambda;
  int pelements;
  //T *lambda;
  //T *lambda_evec;
  //T *lambda_eval;
};

//template <typename T>
//  void diagonalize(void *data);

#endif
