#include <cstdlib>
#include <iostream>
#include "class_defs/gspace.h"
#include "class_defs/usrinput.h"
#include "class_defs/sysinfo.h"
#include "constant.h"
#include "fft_routines.h"
#include "util.h"

using namespace std;



void get_geps(GSPACE *geps, USRINPUT usrin, SYSINFO sys, int iq, int nfft[3], bool *accept){
    
    int ndata = nfft[0]*nfft[1]*nfft[2];
    
    int *gx, *gy, *gz;
    
    gx = new int [ndata];
    gy = new int [ndata];
    gz = new int [ndata];
    
    fftidx_to_gidx( gx, gy, gz, nfft);

    
    double gxtmp, gytmp, gztmp;
    double vtmp[3];
    double Ekin;
    int count = 0;
    
    for (int i=0; i<ndata; i++) {
        gxtmp = gx[i] + sys.qvec[iq][0];
        gytmp = gy[i] + sys.qvec[iq][1];
        gztmp = gz[i] + sys.qvec[iq][2];
        /* transfer to cartesian unit to calculate energy */
        Ekin = 0;
        for (int j=0; j<3; j++) {
            vtmp[j] = gxtmp*sys.b1[j] + gytmp*sys.b2[j] + gztmp*sys.b3[j];
            vtmp[j] *= 2*PI/sys.alat;
            Ekin += 0.5 * vtmp[j] * vtmp[j];
        }
        
        if (Ekin <= usrin.EpsCut) {
            accept[i] = true;
            count += 1;
        }
        else{
            accept[i] = false;
        }
        
    }

    // set values
    geps->ng = count;
    geps->ig = new int [count];
    geps->jg = new int [count];
    geps->kg = new int [count];
    
    int j=0;
    
    for (int i=0; i<ndata; i++) {
        if (accept[i]) {
            geps->ig[j] = gx[i];
            geps->jg[j] = gy[i];
            geps->kg[j] = gz[i];
            j += 1;
        }
    }
    
    if ( j!= count ) {
        cout << " Oops. Error when reducing gspace!!!" << endl;
    }
    
    delete[] gx;
    delete[] gy;
    delete[] gz;
}
