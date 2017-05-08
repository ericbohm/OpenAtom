#include "matrix_test.h"

MatrixTester::MatrixTester(CkArgMsg* msg) {
  thisProxy.test();
}

TestMatrix::TestMatrix(MatrixConfig config, int t) : CBase_TestMatrix(config) {
  setInitial(t);
}

TestMatrix::TestMatrix(int mr, int mc, int tr, int tc, int t) : CBase_TestMatrix(mr, mc, tr, tc) {
  setInitial(t);
}

void TestMatrix::setInitial(int type) {
  for (int r = 0; r < config.tile_rows; r++) {
    for (int c = 0; c < config.tile_cols; c++) {
      switch (type) {
        case ZERO:
          break;
        case RANDOM:
          data[r*config.tile_cols + c].re = (double)rand();
          data[r*config.tile_cols + c].im = (double)rand();
          break;
        case LOCAL:
          data[r*config.tile_cols + c].re = r;
          data[r*config.tile_cols + c].im = c;
          break;
        case GLOBAL:
          data[r*config.tile_cols + c].re = start_row + r;
          data[r*config.tile_cols + c].im = start_col + c;
          break;
        default:
          CkAbort("Bad initialization type\n");
          break;
      }
    }
  }
}

#include "matrix_test.def.h"
