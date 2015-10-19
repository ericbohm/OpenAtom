//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                         GW_BSE:                                           
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

//num_pup :: 


/* class definition header files here
*/
#include "class_gwbseopts.h"
#include "class_gw_epsilon.h"
#include "class_gw_sigma.h"
#include "class_gw_parallel.h"


//==========================================================================

#ifndef _GWBSE_
#define _GWBSE_

#ifdef CHARM_ON
// declare an instance of the class, called readonly_name, as extern or global
class GWBSE; extern GWBSE readonly_gwbse; 
#endif

class GWBSE {

  public:
    // classes that you defined in ../class_defs/GWBSE/xxx.h will be written here 
    GWBSEOPTS         gwbseopts; // GWBSE options. read from user input files
    GW_EPSILON        gw_epsilon; // variables for epsilon.
    GW_SIGMA          gw_sigma;   // variables for sigma
    int num_pup;

    // Parallel configuration
    GW_PARALLEL       gw_parallel; // variables for gw parallelization

    //-------------------------------------------------------------------------
    GWBSE(){num_pup=0;};
    ~GWBSE(){};

#ifdef CHARM_ON
    static GWBSE *get(){
      return &readonly_gwbse;  // return the pointer of the global instance
    }
#endif

    //-------------------------------------------------------------------------
    void state_class_out(){
      PRINTF("\n");
      PRINT_LINE_STAR;
      PRINTF("GWBSE state_class_out %d\n",CKMYPE());
   
      gwbseopts.state_class_out();
      gw_epsilon.state_class_out();
      gw_sigma.state_class_out();
      gw_parallel.state_class_out();

      PRINT_LINE_STAR;
      PRINTF("\n");
    }// end routine

    //-------------------------------------------------------------------------
#ifdef PUP_ON
    void pup(PUP::er &p){

      gwbseopts.pup(p);
      gw_epsilon.pup(p);
      gw_sigma.pup(p);
      gw_parallel.pup(p);


      num_pup++;
#ifdef _DEBUG_PUP_
      if(num_pup==1){
	state_class_out();//every processor write this
        PRINTF("\n");
        PRINT_LINE_STAR;
        PRINTF("GWBSE PUP COMPLETED %d\n",CKMYPE());
        PRINT_LINE_STAR;
        PRINTF("\n");
      }//endif
#endif
    }
#endif

    //-------------------------------------------------------------------------
}; // end class GWBSE
//==========================================================================

#ifdef PUP_ON
PUPmarshall(GWBSE);
#endif

#endif

//==========================================================================



