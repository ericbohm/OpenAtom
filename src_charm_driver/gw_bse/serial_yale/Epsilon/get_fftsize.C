/* Set size of the FFT */

#include "get_fftsize.h"
#include <cstdlib>
#include <cmath>
#include <iostream>

void wfn_fftsize(USRINPUT usrin, SYSINFO sys, int (&nfft)[3]){

    double frac = usrin.wfnFFTsize;

    nfft[0] = int(sys.nfftDen[0] * frac);
    nfft[1] = int(sys.nfftDen[1] * frac);
    nfft[2] = int(sys.nfftDen[2] * frac);

}
