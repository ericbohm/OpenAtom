#ifndef CLA_Matrix_H
#define CLA_Matrix_H

#include "CLA_Matrix.decl.h"
#include "ckmulticast.h"
#include <assert.h>

/* Class definitions */

/* The CLA_Matrix class should not be used directy by the user. See comments
 * below regarding CLA_Matrix_interface. */
class CLA_Matrix : public CBase_CLA_Matrix{
  friend class CLA_Matrix_interface;

  public:
    CLA_Matrix(){}
    CLA_Matrix(CkMigrateMessage *m){}

    /* For 2D algorihtm */
    CLA_Matrix(int M, int K, int N, int m, int k, int n, int part,
     CProxy_CLA_Matrix other1, CProxy_CLA_Matrix other2,
     CProxy_ArrayElement bound, CkCallback ready);
    void receiveA(CLA_Matrix_msg *m);
    void receiveB(CLA_Matrix_msg *m);
  private:
    void multiply(double alpha, double beta, double *data,
     void (*fptr) (void *), void *usr_data);
    void multiply();

    /* shared */
    int M, K, N, m, k, n, um, uk, un;
    int M_chunks, K_chunks, N_chunks;
    int part;
    int algorithm;
    CProxy_CLA_Matrix other1; // For A, B. For B, A. For C, A.
    CProxy_CLA_Matrix other2; // For A, C. For B, C. For C, B.
    CProxy_ArrayElement bound;
    void (*fcb) (void *obj);
    void *user_data;
    double alpha, beta;

    /* For 2D algorithm */
    CProxySection_CLA_Matrix commGroup; // used by A and B
    double *tmpA, *tmpB, *dest; // used by C
    int row_count, col_count; // used by C
    bool got_start; // used by C
};

/* This class below and the make_multiplier function below are the only way in
 * which a user of the library should interact with the libary. Users should
 * never explicitly create char CLA_Matrix object or call their entry methods.
 */

class CLA_Matrix_interface {
  public:
    CLA_Matrix_interface(){}
    inline void multiply(double alpha, double beta, double *data,
     void (*fptr) (void *), void *usr_data, int x, int y){
      assert(p(x, y).ckLocal() != NULL);
      p(x, y).ckLocal()->multiply(alpha, beta, data, fptr, usr_data);
    }
    void pup(PUP::er &per){ per | p; }
  private:
    CProxy_CLA_Matrix p;
    inline void setProxy(CProxy_CLA_Matrix pp){ p = pp; }

  friend int make_multiplier(CLA_Matrix_interface *A, CLA_Matrix_interface *B,
   CLA_Matrix_interface *C, CProxy_ArrayElement bindA,
   CProxy_ArrayElement bindB, CProxy_ArrayElement bindC,
   int M, int K, int N, int m, int k, int n, CkCallback cbA, CkCallback cbB,
   CkCallback cbC, int algorithm);
};

class CLA_Matrix_msg : public CMessage_CLA_Matrix_msg, public CkMcastBaseMsg {
  public:
    CLA_Matrix_msg(double *data, int d1, int d2, int fromX, int fromY);
    double *data;
    int d1, d2;
    int fromX, fromY;
};

/* Function below creates the necessary interfaces so that
 * C = beta * C + alpha * A * B
 * can be computed. X will be bound to bindX for
 * X in {A, B, C}. A is M x K, B is K x N, C is M x N. M, K, and N are
 * decomposed into chunks of size m, k, and n, respectively. Along a given
 * dimension, the array must be of size ceil(1.0 * Y / y) (Y in {M, K, N},
 * y in {m, k, n}). If y does not divide Y, the last element must have only
 * Y % y elements along the given dimension. When matrix X is ready, it does a
 * callback to cbX. The algorithm used to multiply the matrices is determined
 * by the value of 'algorithm', which can take the values defined below.
 *
 * Return value: Zero is returned upon success. A negative value is returned
 * if an error occurs.
 *  ERR_INVALID_ALG: an invalid algorithm was selected
 *  ERR_INVALID_DIM: invalid dimensions were given
 */

/* Valid values for 'algorithm' are given below. As new ones are added,
 * they should be given incremental numbers. MM_ALG_MIN should not be changed,
 * and MM_ALG_MAX should have the value of the greatest defined algorithm.
 * MM_ALG_MIN MM_ALG_MAX are used to validate user input, so they should always
 * be updated as new algorithms are added.
 */
#define MM_ALG_MIN	1
#define MM_ALG_2D	1
#define MM_ALG_3D	2
#define MM_ALG_MAX	2

int make_multiplier(CLA_Matrix_interface *A, CLA_Matrix_interface *B,
 CLA_Matrix_interface *C, CProxy_ArrayElement bindA,
 CProxy_ArrayElement bindB, CProxy_ArrayElement bindC,
 int M, int K, int N, int m, int k, int n, CkCallback cbA, CkCallback cbB,
 CkCallback cbC, int algorithm);

/* Error codes */
#define ERR_INVALID_ALG		-1
#define ERR_INVALID_DIM		-2

#endif
