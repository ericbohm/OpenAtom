#ifndef MATRIX_H
#define MATRIX_H

struct MatrixConfig {
  unsigned int mat_rows, mat_cols;
  unsigned int tile_rows, tile_cols;

  inline unsigned int chareRows() const { return mat_rows/tile_rows; }
  inline unsigned int chareCols() const { return mat_cols/tile_cols; }
};

MatrixConfig convertTo1D(const MatrixConfig& in, int tile_rows);

#include "matrix.decl.h"

struct MatrixHandle {
  MatrixConfig config;
  CProxy_Matrix proxy;
};
PUPbytes(MatrixConfig);
PUPbytes(MatrixHandle);

// API
void matrixDestroy(CProxy_Matrix matrix);
void matrixCopy(CProxy_Matrix src, CProxy_Matrix dest, CkCallback cb);
void matrixCompare(CProxy_Matrix src, CProxy_Matrix dest, CkCallback cb);
void matrixRead(CProxy_Matrix dest, std::string prefix, CkCallback cb);
void matrixWrite(CProxy_Matrix src, std::string prefix, CkCallback cb);
void matrixVerify(CProxy_Matrix src, std::string prefix, CkCallback cb);

struct DataMessage : public CMessage_DataMessage {
  unsigned int row, col, num_rows, num_cols;
  complex* data;
};

class Matrix : public CBase_Matrix {
  Matrix_SDAG_CODE;
  protected:
    MatrixConfig config;
    unsigned int start_row, start_col;
    complex* data;

    bool busy;
    unsigned int data_received, total_data;
    void (Matrix::*dataHandler)(DataMessage*);

    void initialize();

    void packMsg(DataMessage* msg);
    void unpackMsg(DataMessage* msg);
    void compareMsg(DataMessage* msg);

  public:
    Matrix() {}
    Matrix(MatrixConfig config);
    Matrix(int mr, int mc, int tr, int tc);

    // Entry methods
    void copy(CProxy_Matrix src, CkCallback cb);
    void compare(CProxy_Matrix src, CkCallback cb);
    void read(std::string prefix, CkCallback cb);
    void write(std::string prefix, CkCallback cb);
    void verify(std::string prefix, CkCallback cb);
    void sendData(CProxy_Matrix dest, int dest_rows, int dest_cols);
};

#endif
