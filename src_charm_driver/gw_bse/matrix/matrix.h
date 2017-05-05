#ifndef MATRIX_H
#define MATRIX_H

#include "matrix.decl.h"

// API
void matrixCopy(CProxy_Matrix src, CProxy_Matrix dest, CkCallback cb);
void matrixCompare(CProxy_Matrix src, CProxy_Matrix dest, CkCallback cb);

struct DataMessage : public CMessage_DataMessage {
  unsigned int row, col, num_rows, num_cols;
  complex* data;
};

// TODO: Instead of passing around all kinds of ints
struct MatrixConfig {};

// TODO: Maybe make a wrapper class around these, that wraps proxy and config
// that way the user doesn't have to know about the other matrix config when
// doing things like copying. It's error prone.

class Matrix : public CBase_Matrix {
  Matrix_SDAG_CODE;
  protected:
    unsigned int matrix_rows, matrix_cols;
    unsigned int num_rows, num_cols;
    unsigned int start_row, start_col;
    complex* data;

    bool busy;
    unsigned int data_received, total_data;
    void (Matrix::*dataHandler)(DataMessage*);

    void initialize(int mr, int mc, int tr, int tc, CkCallback cb);

    void sendTo(CProxy_Matrix dest, int dest_rows, int dest_cols);

    void packMsg(DataMessage* msg);
    void unpackMsg(DataMessage* msg);
    void compareMsg(DataMessage* msg);

  public:
    Matrix() {}
    Matrix(int mr, int mc, int tr, int tc, CkCallback cb);
    void startCopy(CProxy_Matrix src, CkCallback cb);
    void copy(CProxy_Matrix dest, int dest_rows, int dest_cols);
    void startCompare(CProxy_Matrix src, CkCallback cb);
    void compare(CProxy_Matrix dest, int dest_rows, int dest_cols);

    void read(std::string prefix, CkCallback cb);
    void write(std::string prefix, CkCallback cb);
};




#endif
