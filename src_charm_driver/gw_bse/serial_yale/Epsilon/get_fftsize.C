/* Set size of the FFT */

#include "get_fftsize.h"
#include <cstdlib>
#include <cmath>
#include <iostream>

// at this moment, this function takes the size of the dense fft grid from QEspresso output
// and just multiply some number (which is given by user) to them, and make it the fftgrid for wavefunction
// TODO maybe we should change this routine
void wfn_fftsize(USRINPUT usrin, SYSINFO sys, int (&nfft)[3]){

    double frac = usrin.wfnFFTsize;

    nfft[0] = int(sys.nfftDen[0] * frac);
    nfft[1] = int(sys.nfftDen[1] * frac);
    nfft[2] = int(sys.nfftDen[2] * frac);

}
