#include "class_defs/sysinfo.h"
#include "class_defs/usrinput.h"
#include "util.h"

void check_inputs(USRINPUT usrin, SYSINFO sys){   
    //check if the input file variables matche with the data file variables
    if ( usrin.nstate != sys.nstate ){
        char errmsg[] = "nstate in epsilon.in is different from nstate in sysinfo.dat file";
        Die(errmsg);
    }
        
    if ( usrin.nspin != sys.nspin ){
        char errmsg[] = "nspin in epsilon.in is different from nspin in sysinfo.dat file";
        Die(errmsg);
    }
    
    if ( usrin.nkpt != sys.nkpt ){
        char errmsg[] = "nkpt in epsilon.in is different from nkpt in sysinfo.dat file";
        Die(errmsg);
    }

}

