/* 
   header file for epsilon program
*/
#include "class_defs/sysinfo.h"
#include "class_defs/usrinput.h"
#include "class_defs/states.h"
#include "class_defs/matrix.h"
#include "class_defs/gspace.h"
#include <string>
#include "util.h"
#include "iter_invmtrx.h"
#include "get_fftsize.h"
#include "fft_routines.h"
#include "calc_P.h"
#include "print_util.h"

using namespace std;

// function declaration here. If functions are called in the main function, 
// the function declaraions is in this header file
void check_inputs(USRINPUT, SYSINFO);
void read_states(STATES *);
void do_fft_states(STATES**, STATES**, SYSINFO, int[3]);
void get_k_plus_q_index(int, int, int &, SYSINFO, int (&)[3]);
void calc_Pmtrx_Rspace(CMATRIX*, SYSINFO, STATES*, STATES*, int [3], int [3]);
void Pmtrx_R_to_G(int, CMATRIX*, int [3], double);
void get_geps(GSPACE*, USRINPUT, SYSINFO, int, int[3], bool*);
void calc_Epsmat(double*, int, SYSINFO, GSPACE*, CMATRIX*, CMATRIX*, int, bool*);
void iter_invmtrx(CMATRIX*, USRINPUT, int);

