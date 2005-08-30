/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                      class_mdnbr_list.h                                  */
/*                                                                          */
/*    Class definition for intermolecular neibhor list methods              */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#ifndef _MDNBR_LIST_
#define _MDNBR_LIST_

#include "class_mdverlist.h"
#include "class_mdlnklist.h"
#include "class_mdbrnch_root.h"

class MDNBR_LIST {
 public:
  int natm_tot;                /* Num: Number of atoms               */
  int iver;                    /* Opt: Verlet list option             */ 
  int ilnk;                    /* Opt: Lnk list option                */
  int nolst;                   /* Opt: No list option                 */     
  int brnch_root_list_opt;     /* Num: Option to shave brnch interactions */
  double *x0,*y0,*z0;          /* Lst: Pos of the atms when the list
                                       is calculated;        
                                  Lth: natm_tot                       */
//=====================================================================================
// Declare instances of subclasses 

  MDVERLIST mdverlist;
  MDLNKLIST mdlnklist;
  MDBRNCH_ROOT mdbrnch_root;

//=====================================================================================
// Default constructor/destructor

   MDNBR_LIST(){
    natm_tot = 0;            
    iver     = 0;                
    ilnk     = 0;                
    nolst    = 0;               
    brnch_root_list_opt = 0; 

   }
  ~MDNBR_LIST(){}
//=====================================================================================
// Pack/Unpack

#ifdef PUP_ON
  void pup(PUP::er &p){

    // PUP ints
    p | natm_tot;
    p | iver;
    p | ilnk;
    p | nolst;
    p | brnch_root_list_opt;

    // PUP arrays

#ifdef ANYLIST_IMPLEMENTED
    if(iver == 1){
      pup1d_dbl(p,&x0,natm_tot);
      pup1d_dbl(p,&y0,natm_tot);
      pup1d_dbl(p,&z0,natm_tot);
    } /* endif */
#endif
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif
  } /* end pack unpack */
#endif

//=====================================================================================
// Print out state of class

  void state_class_out(){
    
     char fileName [255];
     sprintf (fileName, "%d_mdnbr_list.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");

     fprintf(fp,"mdnbr_list:  natm_tot %d\n",natm_tot);
     fprintf(fp,"mdnbr_list:  iver %d\n",iver);
     fprintf(fp,"mdnbr_list:  ilnk %d\n",ilnk);
     fprintf(fp,"mdnbr_list:  nolst %d\n",nolst);
     fprintf(fp,"mdnbr_list:  brnch_root_list_opt %d\n",brnch_root_list_opt);

     fclose(fp);
  }/* end member function */

 
}; /* end class definitions */

#ifdef PUP_ON
PUPmarshall(MDNBR_LIST);
#endif

#endif





