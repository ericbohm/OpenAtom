/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

/** \file configure.h
 *
 */

#ifndef _Configure_GWBSE_
#define _Configure_GWBSE_

#include "dictionary.h"

//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
class Config {
//===================================================================================
 public:

  //==================================
  // input or derived parameters
  //----------------------------------


  //==================================
  // Class Functions
  //----------------------------------
   Config(){};
  ~Config(){};
   void readConfig(char*, GWBSE*);
   void readStateInfo(int &, int &, int &, int &, int &, int &,const char *, int,int,int,int,int); // no need?
   void simpleRangeCheck(); // check the input that user provides sensible
   void rangeExit(int, const char *, int); // if the range is not right, exit
   void Finale(int, int, int, int, int, int); // finish the code

   void set_config_dict_fun    (int *, DICT_WORD **); // functional keywords


   void set_config_dict_gen_GW            (int *, DICT_WORD **);
   void set_config_dict_GW_epsilon        (int *, DICT_WORD **);
   void set_config_dict_GW_sigma          (int *, DICT_WORD **);
   void set_config_dict_GW_file        (int *, DICT_WORD **);
   void set_config_dict_GW_charm_input        (int *, DICT_WORD **);

   void set_config_params_gen_GW  (DICT_WORD *, char *, char *, GWBSEOPTS *);
   void set_config_params_GW_epsilon (DICT_WORD *, char *, char *, GW_EPSILON *);
   void set_config_params_GW_sigma (DICT_WORD *, char *, char *, GW_SIGMA *);
   void set_config_params_GW_file (DICT_WORD *, char *, char *, GWBSEOPTS *);
   //   void guesstimateParmsConfig ( ); // need this
   void write_cpaimd_config    (FILE *, DICT_WORD *, int, char *);
//   void load_cpaimd_config (DICT_WORD *, int, PINY_NAME *, PINY_NAME *, int, int *);


   void read_lattice(GWBSEOPTS *);
   void read_klist(GWBSEOPTS *);


//============================================================================
  //==================================

//-----------------------------------------------------------------------------------
   }; // end class
//===================================================================================


//===================================================================================
#endif
