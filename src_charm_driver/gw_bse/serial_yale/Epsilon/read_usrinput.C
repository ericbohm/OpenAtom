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

void read_usrinput(USRINPUT &usrin, SYSINFO &sys){
    
    ifstream fp ("gw.input");              //fp is a pointer for gw.input file
    if (!fp){
	char errmsg[] = "Failed to open input file (gw.input)!!!";
	Die(errmsg);
    }
    
    if (fp){
	fp >> usrin.nstate;
	fp >> usrin.nocc;
	fp >> usrin.nunocc;
	fp >> usrin.nkpt;
	fp >> usrin.nspin;
	
        fp >> usrin.wfnFFTsize;
        fp >> usrin.EpsCut;
        fp >> usrin.shift[0] >> usrin.shift[1] >> usrin.shift[2];

	fp >> usrin.iter_maxiter;
        fp >> usrin.iter_convg;
    }
    
    fp.close();
    
    //check if the input file variables matche with the data file variables
    if ( usrin.nstate != sys.nstate ){
        char errmsg[] = "nstate in gw.input is different from nstate in sysinfo.dat file";
        Die(errmsg);
    }
        
    if ( usrin.nspin != sys.nspin ){
        char errmsg[] = "nspin in gw.input is different from nspin in sysinfo.dat file";
        Die(errmsg);
    }
    
    if ( usrin.nkpt != sys.nkpt ){
        char errmsg[] = "nkpt in gw.input is different from nkpt in sysinfo.dat file";
        Die(errmsg);
    }

    // save number of occupied and unoccupied states to sys object
    sys.nocc = usrin.nocc;
    sys.nunocc = usrin.nunocc;
}


// note: you cannot jump to the next line automatically as can be done in fortran
// I choose the simplest option. No commenting inside of the input file
