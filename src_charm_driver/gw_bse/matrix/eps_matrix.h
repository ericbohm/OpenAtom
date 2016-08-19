
#include "eps_matrix.decl.h"

#include "CLA_Matrix.h"
#include "ckcomplex.h"

#define ZGEMM zgemm
#define ZGERC zgerc

extern "C" {
  void ZGERC (int*, int*, complex*, complex*, int*, complex*, int*, complex*, int*);
  void ZGEMM (char*, char*, int*, int*, int*, complex*, complex*, int*, complex*, int*, complex*, complex*, int*);
}

class EpsMatrix2D : public CBase_EpsMatrix2D {
  EpsMatrix2D_SDAG_CODE
  public:
    EpsMatrix2D();
    EpsMatrix2D(CkMigrateMessage* msg) {}
    EpsMatrix2D(CLA_Matrix_interface mat);

    void checkReady();
    void setSize(int size);
    void applyFs();
    void sendTo1D();
    void receiveChunk(Phase2Message*);
    void createTranspose();
    void receiveTranspose(std::vector<complex> incoming);
    void sendTo(CProxy_EpsMatrix2D receiver_proxy);
    void receiveData(std::vector<complex> incoming);
    void calc_vcoulb();
    void calc_Eps(Phase3Message* msg);
    void receiveFs(complex eps_data[140*140]);
    void reportPTime();
    void print_res();
    void multiply(double alpha, double beta);
    void round_done(void);
    void findAlpha(void);
    long double max_fn(int size);
    void scalar_multiply(double alpha);
    void convergence_check(CProxy_EpsMatrix2D cmp_proxy);
    void add_compl_two();
    void setI(CLA_Matrix_interface mat, bool clean);
    void receiveConvCheck(std::vector<complex> incoming);
    static void done_cb(void *obj){
     ((EpsMatrix2D*) obj)->round_done();                                                                     
    } 
    
    virtual void pup(PUP::er &p){
      CBase_EpsMatrix2D::pup(p);
      p | num_rows;
      p | num_cols;
      p | matrix;
      if(p.isUnpacking())
        data = new complex[num_rows * num_cols];
      PUParray(p, data, num_cols * num_rows);
    }
    CLA_Matrix_interface matrix;
  private:
    unsigned L; // Number of occupied psis
    unsigned matrix_dimension; // Size of the entire matrix
    unsigned num_rows, num_cols, start_row, start_col; // The shape of our data
    unsigned trans_count, num_chares; // SDAG variables
    complex* data;
    int* nfft; // number of fft grids in each direction
    unsigned qindex;
    unsigned local_mtx_size_1d_y;
    int receive_counter;
    complex total[144];
    void kqIndex(unsigned, unsigned&, int*);
    complex* umklapp_factor;
    void getUmklappFactor(complex*, int[3]);

    double total_time;
};
extern /* readonly */ CProxy_EpsMatrix2D pmatrix2D_bproxy;
extern /* readonly */ CProxy_EpsMatrix2D pmatrix2D_cproxy;

extern /* readonly */ CProxy_EpsMatrix2D pmatrix2D_aaproxy;
extern /* readonly */ CProxy_EpsMatrix2D pmatrix2D_bbproxy;
extern /* readonly */ CProxy_EpsMatrix2D pmatrix2D_ccproxy;
