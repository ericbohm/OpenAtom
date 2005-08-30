//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                         PI_MD:                                           
//             The future of simulation technology                          
//             ------------------------------------                         
//                   Structures: typedefs_gen.h                             
//                                                                          
// The include file with the typedefs of all the option structures          
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#include "../class_defs/GEN/class_genstatepoint.h"
#include "../class_defs/GEN/class_gencell.h"
#include "../class_defs/GEN/class_gentimeinfo.h"
#include "../class_defs/GEN/class_genptens.h"
#include "../class_defs/GEN/class_genmdstat_avg.h"
#include "../class_defs/GEN/class_gencpstat_avg.h"
#include "../class_defs/GEN/class_genpistat_avg.h"
#include "../class_defs/GEN/class_genewald.h"
#include "../class_defs/GEN/class_gensimopts.h"
#include "../class_defs/GEN/class_genminopts.h"
#include "../class_defs/GEN/class_genensopts.h"
#include "../class_defs/GEN/class_genfilenames.h"

//==========================================================================

#ifndef _GENERAL_DATA_
#define _GENERAL_DATA_

#ifdef CHARM_ON
// declare an instance of the class, called readonly_name, as extern or global
class GENERAL_DATA; extern GENERAL_DATA readonly_general_data;  
#endif

class GENERAL_DATA {

  public :
   int num_pup;
   double        tot_memory;
   GENSTATEPOINT genstatepoint;
   GENCELL       gencell;
   GENTIMEINFO   gentimeinfo;
   GENPTENS      genptens;
   GENMDSTAT_AVG genmdstat_avg;
   GENCPSTAT_AVG gencpstat_avg;
   GENPISTAT_AVG genpistat_avg;
   GENEWALD      genewald;
   GENSIMOPTS    gensimopts;
   GENMINOPTS    genminopts;
   GENENSOPTS    genensopts;
   GENFILENAMES  genfilenames;

//--------------------------------------------------------------------------

   GENERAL_DATA(){num_pup=0;tot_memory=0;};
  ~GENERAL_DATA(){};

#ifdef CHARM_ON
  static GENERAL_DATA *get(){
    return &readonly_general_data;  // return the pointer of the global instance
  }
#endif

//--------------------------------------------------------------------------
  void state_class_out(){
   PRINTF("\n");
   PRINT_LINE_STAR;
   PRINTF("GENERAL_DATA state_class_out\n");
   genstatepoint.state_class_out();
   gencell.state_class_out();
   gentimeinfo.state_class_out();
   genptens.state_class_out();
   genmdstat_avg.state_class_out();
   gencpstat_avg.state_class_out();
   genpistat_avg.state_class_out();
   genewald.state_class_out();
   gensimopts.state_class_out();
   genminopts.state_class_out();
   genensopts.state_class_out();
   genfilenames.state_class_out();
   PRINT_LINE_STAR;
   PRINTF("\n");
  }// end routine

//--------------------------------------------------------------------------
#ifdef PUP_ON
  void pup(PUP::er &p){
   p | tot_memory;
   p | num_pup;
   genstatepoint.pup(p);
   gencell.pup(p);
   gentimeinfo.pup(p);
   genptens.pup(p);
   genmdstat_avg.pup(p);
   gencpstat_avg.pup(p);
   genpistat_avg.pup(p);
   genewald.pup(p);
   gensimopts.pup(p);
   genminopts.pup(p);
   genensopts.pup(p);
   genfilenames.pup(p);
   num_pup++;
   if(num_pup==2){
     PUP_PRINTF("\n");
     PUP_PRINT_LINE_STAR;
     PUP_PRINTF("GEN PUP COMPLETED\n");
     PUP_PRINT_LINE_STAR;
     PUP_PRINTF("\n");
   }//endif
  }
#endif

}; // end class general_data
//--------------------------------------------------------------------------

#ifdef PUP_ON
PUPmarshall(GENERAL_DATA);
#endif

#endif

//==========================================================================
