#ifndef TESTMATRIX_H
#define TESTMATRIX_H

#include "matrix_test.decl.h"

#include "matrix.h"

class MatrixTester : public CBase_MatrixTester {
  MatrixTester_SDAG_CODE;
  private:
    int mat_rows, mat_cols, tile_rows, tile_cols, chare_rows, chare_cols;
    CProxy_Matrix mat1, mat2;
  public: 
    MatrixTester(CkArgMsg* msg);
};

class TestMatrix : public CBase_TestMatrix {
  private:
    void setInitial(int type = 0);
  public:
    TestMatrix(int,int,int,int,CkCallback);
};

#endif
