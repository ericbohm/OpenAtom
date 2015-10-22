#include "class_defs/sysinfo.h"
#include <cmath>
// this routine finds k+q vector index
// ik is current k vector index 
// iq is current q vector index
// ikq is k+q vector index

using namespace std;
void get_k_plus_q_index(int iq, int ik, int &ikq, SYSINFO sys, int (&uklpp)[3]){

    // temporary space to save k/q/k+q vectors
    double this_k[3], k_plus_q[3], k_plus_q_orig[3];
    
    for (int i=0; i<3; i++) {
        this_k[i] = sys.kvec[ik][i]; // k vector
        k_plus_q[i] = sys.qvec[iq][i]; // q vector
        // calculate k+q vector 
        k_plus_q[i] += this_k[i]; // k+q vector
        k_plus_q_orig[i] = k_plus_q[i]; // save it for Umklapp scattering
	// if k+q[i] is not between 0 and 1, adjust k+q
        if ( k_plus_q[i] >= 1 ) {
            k_plus_q[i] -= 1;
        }
        else if( k_plus_q[i] < 0 ){
            k_plus_q[i] += 1;
        }
    }
    
    double diff;
    
    // find k+q vector index 
    for (int kk=0; kk<sys.nkpt; kk++) {
        diff = 0;
	//this_k is now a difference between kvec[kk] (trial kvec) and kvec[ik]+qvec[iq]
        for (int i=0; i<3; i++) {
            this_k[i] = abs( sys.kvec[kk][i] - k_plus_q[i] );
            diff += this_k[i];
        }
        if (diff == 0) {
            ikq = kk;  //k+q index is found
        }
    }
    // save umklapp scattering information
    for (int i=0; i<3; i++) {
        uklpp[i] = int( k_plus_q_orig[i] - k_plus_q[i] );
    }
}

