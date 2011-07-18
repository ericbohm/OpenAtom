#include "ckcomplex.h"

#ifndef MATRIX_2_FILE_H
#define MATRIX_2_FILE_H

#include "debug_flags.h"

//@{
/// Matrix read/write utils
void dumpMatrix(const char *infilename,  double *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric);
void dumpMatrix(const char *infilename, complex *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric);
void loadMatrix(const char *infilename,  double *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric);
void loadMatrix(const char *infilename, complex *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric);
void dumpMatrix2DDouble(const char *infilename, double **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric);
void loadMatrix2DDouble(const char *infilename, double **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric);
void dumpMatrix2DInt(const char *infilename, int **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric);
void loadMatrix2DInt(const char *infilename, int **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric);
//@}

#endif // MATRIX_2_FILE_H

