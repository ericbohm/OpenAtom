/* this routine calculates coulomb potential in G space,
   i.e., vcoulb = 4*PI/(Nk*vol)/(q+G)^2  (q and G should be in cartesian coordinates (1/a.u.), so we transform it from crystal coordinates to the cartesian coordinates by multipying reciprocal lattice vectors) 
   the total g vectors are passed into this routine
*/

#include <cstdlib>
#include "class_defs/sysinfo.h"
#include "class_defs/gspace.h"
#include "constant.h"

void calc_vcoulb(int iq, double *vcoulb, GSPACE *geps, SYSINFO sys) {
    
    double gx, gy, gz;
    double gq[3];
    const double fact = 4*PI/sys.vol/sys.nkpt;
    
    for (int i=0; i<geps->ng; i++) {
 
 
        if (iq==0){
            gx = geps->ig[i] + sys.shift[0];
            gy = geps->jg[i] + sys.shift[1];
            gz = geps->kg[i] + sys.shift[2];
        }
        else{
            gx = geps->ig[i] + sys.qvec[iq][0];
            gy = geps->jg[i] + sys.qvec[iq][1];
            gz = geps->kg[i] + sys.qvec[iq][2];
        }
        
        vcoulb[i] = 0;
        for (int j=0; j<3; j++) {
            gq[j] =  gx*sys.b1[j] + gy*sys.b2[j] + gz*sys.b3[j];
            gq[j] *= 2*PI/sys.alat;

            vcoulb[i] += gq[j]*gq[j];
        }
        vcoulb[i] = 1/vcoulb[i];
        vcoulb[i] *= fact;
    }
    
}
