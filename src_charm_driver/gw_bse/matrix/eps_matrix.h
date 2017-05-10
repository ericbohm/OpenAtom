
#include "matrix.h"
#include "eps_matrix.decl.h"

#include "mylapack.h"
#include "CLA_Matrix.h"
#include "ckcomplex.h"

class EpsMatrix : public CBase_EpsMatrix {
  private:
    unsigned L; // Number of occupied psis
    int* nfft; // number of fft grids in each direction
    unsigned qindex;

    unsigned data_received;
    double total_time;
    CLA_Matrix_interface matrix;

  public:
    EpsMatrix();
    EpsMatrix(MatrixConfig config);
    EpsMatrix(CkMigrateMessage* msg) {}

    void checkReady();
    void createTranspose(CProxy_EpsMatrix other, bool todo);
    void receiveTranspose(std::vector<complex> incoming);
    void calc_vcoulb();
    void calc_Eps(Phase3Message* msg);
    void receiveFs(Phase3Message* msg);
    void multiply(double alpha, double beta);
    void round_done(void);
    void findAlpha(void);
    void screenedExchange();
    void bareExchange();
    void scalar_multiply(double alpha);
    void convergence_check(CProxy_EpsMatrix cmp_proxy);
    void add_compl_two();
    void setI(CLA_Matrix_interface mat, bool clean);
    void receiveConvCheck(std::vector<complex> incoming);
    static void done_cb(void *obj){
     ((EpsMatrix*) obj)->round_done();
    }
};
