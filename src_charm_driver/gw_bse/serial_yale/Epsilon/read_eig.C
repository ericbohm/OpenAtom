/*    This subroutine reads eigenvalues from a data file "Eigenvalue" in each wavefunction directory 
      This subroutine takes psi[is][ik] pointer, then save eigenvalue and occupancy at psi[is][ik][ib]->eig and psi[is][ik][ib]->occ
*/


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include "class_defs/states.h"
#include "util.h"

using namespace std;

void read_eig(STATES **psi, int nstate){
    
    int ispin = psi[0]->ispin;
    int ikpt = psi[0]->ikpt;
    bool shift = psi[0]->shifted;
    
    /* open Eigenvalue file */
    stringstream ss;
    stringstream sk;
    ss << ispin;
    sk << ikpt;

    string fname;
    
    if(!shift){
        fname = "./Spin." + ss.str() + "_Kpt." + sk.str() + "_Bead.0_Temper.0/Eigenvalue";
    }
    if(shift){
        fname = "./Spin." + ss.str() + "_Kpt.0" + sk.str() + "_Bead.0_Temper.0/Eigenvalue";
    }

    
    ifstream fp;
    fp.open( fname.c_str() );
    
    if(!fp){
	Die("Failed to open Eigenvalue file");
    }
    
    if(fp){
	double dummy;
        for (int i=0; i < nstate; i++) {
            fp >> psi[i]->eig >> psi[i]->occ;
        }
    }
    fp.close();
}
