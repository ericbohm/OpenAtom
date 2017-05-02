#include "class_defs/sysinfo.h"
#include "class_defs/matrix.h"
#include "class_defs/states.h"
#include "class_defs/interpolator.h"
#include "class_defs/laplace.h"
#include "include/ckcomplex.h"
#include "fftw3.h"
#include "util.h"
#include "my_fftw.h"

void update_Pmtrx(complex *, double, complex *, double, int, SYSINFO, CMATRIX*);
void modify_state_Uproc(complex *, int [3], int [3], SYSINFO);
void CalcPmtrxRspaceInterpolation(CMATRIX*, INTERPOLATOR*, STATES*, int [3], int [3], SYSINFO);
void CalcPmtrxLaplace(CMATRIX*, LAPLACE&, STATES*, STATES*, int [3], int [3], SYSINFO);
