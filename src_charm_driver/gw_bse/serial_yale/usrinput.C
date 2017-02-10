/* 
 read input file from stdin and get params
*/
#include <iostream>
#include <fstream>                         // class to handle files
#include <cstdlib>                         // for exit()
#include "class_defs/usrinput.h"
#include "class_defs/sysinfo.h"
#include "util.h"

using namespace std;

void USRINPUT::read_usrinput(){
    
    ifstream fp ("epsilon.in");              //fp is a pointer for epsilon.in file
    if (!fp){
	char errmsg[] = "Failed to open input file (epsilon.in)!!!";
	Die(errmsg);
    }
    
    if (fp){
	fp >> nstate;
	fp >> nocc;
	fp >> nunocc;
	fp >> nkpt;
	fp >> nspin;
	
        fp >> EcutFFT;
        fp >> EpsCut;
        fp >> shift[0] >> shift[1] >> shift[2];

	fp >> iter_maxiter;
        fp >> iter_convg;
#ifdef USE_P_INTERPOLATION
        fp >> nPitp;
#endif
#ifdef USE_P_LAPLACE
        fp >> nPitp; // dummy variable in laplace method
        fp >> ptol;
#endif

    }
    
    fp.close();
}

