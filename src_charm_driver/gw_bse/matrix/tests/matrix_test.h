#ifndef TESTMATRIX_H
#define TESTMATRIX_H

#include "matrix.h"
#include "matrix_test.decl.h"

#define ZERO    0
#define RANDOM  1
#define LOCAL   2
#define GLOBAL  3

class MatrixTester : public CBase_MatrixTester {
  MatrixTester_SDAG_CODE;
  private:
    MatrixConfig config1D, config2D;
    CkCallback callback;
    CProxy_Matrix mat1, mat2, mat3;
  public: 
    MatrixTester(CkArgMsg* msg);
};

class TestMatrix : public CBase_TestMatrix {
  private:
    void setInitial(int);
  public:
    TestMatrix(MatrixConfig,int);
    TestMatrix(int,int,int,int,int);
};

#endif
