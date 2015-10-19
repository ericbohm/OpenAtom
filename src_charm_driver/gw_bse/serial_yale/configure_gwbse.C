//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                         GW-BSE:
//             The future of simulation technology                          
//             ------------------------------------                         
//                     Module: configure.c
//                                                                          
//                                                                          
// This subprogram reads in user simulation params and echoes them
// and the default parameters to a file                                     
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================


//#ifdef _DEBUG_GWINPUT_
//#include "standard_include_gwbse_dbg.h"
//#else 
#include "standard_include_gwbse.h"
//#endif

#include "allclass_gwbse.h"
#include "configure_gwbse.h"

// these are openatom dictionary parsers
#include "proto_friend_lib_entry.h"
#include "proto_handle_entry.h"

//===================================================================================

//===================================================================================
// dummy main program for scalar debugging
//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
#ifdef _DEBUG_GWINPUT_
int main();
int main(){
    GWBSE gwbse;
    char input_name[10000];
    Config gw_configure;
    
    sprintf (input_name,"minjung_test_file");
    gw_configure.readConfig (input_name,&gwbse);
    gwbse.state_class_out();
    
    return 1;
}
#endif
//===================================================================================





//===================================================================================
// read configuration file using openatom parser
//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::readConfig(char* input_name, GWBSE *gwbse)
  //===================================================================================
{//begin routine
  //===================================================================================

  int num_dict_fun; // number of meta-keywords
  int num_dict_gen_GW;
  int num_dict_GW_epsilon;
  int num_dict_GW_sigma;
  int num_dict_GW_file;
  int num_dict_GW_charm_input;
  DICT_WORD *dict_fun;
  DICT_WORD *dict_gen_GW;
  DICT_WORD *dict_GW_epsilon;
  DICT_WORD *dict_GW_sigma;
  DICT_WORD *dict_GW_file;
  DICT_WORD *dict_GW_charm_input;
  DICT_WORD word;            

  int nline;
  int nkey,nfun_key;
  char *fun_key,*fname;
  int ind_key;
  FILE *fp;
#include "allclass_strip_gwbse.h"

  fun_key = (char *)cmalloc(PINY_MAXWORD*sizeof(char),"Config::readConfig");
  fname   = (char *)cmalloc(1024*sizeof(char),"Config::readConfig");

  //===================================================================================
 
   //===================================================================================
  // Tell everyone you are busy 

  PRINTF("  =============================================================\n");
  PRINTF("  Reading GW-BSE input from file : %s\n",input_name);
  PRINTF("  -------------------------------------------------------------\n\n");

  if(PINY_MAXWORD!=MAXWORD || PINY_MAXLINE != MAXLINE){
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("   Incorrect word and line sizes\n");
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

  //===================================================================================
  //===================================================================================
  //===================================================================================
  // Set up the dictionaries   // fill the dictionary  with keywords and key arguments
  set_config_dict_fun( &num_dict_fun, &dict_fun );
  set_config_dict_gen_GW( &num_dict_gen_GW, &dict_gen_GW );
  set_config_dict_GW_epsilon( &num_dict_GW_epsilon, &dict_GW_epsilon );
  set_config_dict_GW_sigma( &num_dict_GW_sigma, &dict_GW_sigma );
  //set_config_dict_GW_file( &num_dict_GW_file, &dict_GW_file );
  //set_config_dict_GW_charm_input( &num_dict_GW_charm_input, &dict_GW_charm_input );

  //===================================================================================
  // Read the input file and fill the dictionaries with user input

  fp = cfopen((const char *) input_name,"r");

  nline = 1;  nfun_key=0;
  while(get_fun_key(fp,fun_key,&nline,&nfun_key,input_name)){ // finds the meta keywords
    get_fun_key_index(fun_key,num_dict_fun,dict_fun,nline,nfun_key,
        input_name,&ind_key);  //finds which one it is
    nkey  = 0;
    while(get_word(fp,&word,&nline,&nkey,nfun_key,input_name)){
      switch(ind_key){
        case 1 : put_word_dict(&word,dict_gen_GW,num_dict_gen_GW,fun_key,nline,
                     nkey,nfun_key,input_name);break;
        case 2 : put_word_dict(&word,dict_GW_epsilon,num_dict_GW_epsilon,fun_key,nline,
                     nkey,nfun_key,input_name);break;
        case 3 : put_word_dict(&word,dict_GW_sigma,num_dict_GW_sigma,fun_key,nline,
                     nkey,nfun_key,input_name);break;
        case 4 : put_word_dict(&word,dict_GW_file,num_dict_GW_file,fun_key,nline,
                     nkey,nfun_key,input_name);break;
        case 5 : put_word_dict(&word,dict_GW_charm_input,num_dict_GW_charm_input,fun_key,nline,
                     nkey,nfun_key,input_name);break;
      }//end switch
    }// end while 
  }//end while
// default values assigned
  fclose(fp);  

  //===================================================================================
  // Take the information out of the dictionary and put it in the class

  set_config_params_gen_GW      (dict_gen_GW,  dict_fun[1].keyword, input_name, gwbseopts);
  set_config_params_GW_epsilon  (dict_GW_epsilon, dict_fun[2].keyword, input_name, gw_epsilon);
  set_config_params_GW_sigma    (dict_GW_sigma, dict_fun[3].keyword, input_name, gw_sigma);
  //  set_config_params_GW_file     (dict_GW_file, dict_fun[4].keyword, input_name, gwbseopts);
//  set_config_params_GW_charm_input   (dict_GW_charm_input, dict_fun[5].keyword, input_name, iflag);
// my input values
//  simpleRangeCheck_gwbse(); // redundant checking   need to write this


  //===================================================================================
  // Improve user parameters and/or try to optimize unset parameters

  // guesstimateParmsConfig(sizez,dict_gen,dict_rho,dict_state,dict_pc,dict_nl,dict_map, 
  //    nchareRhoRHart, nplane_x_rho, natm_typ);

  //===================================================================================
  // Final consistency checks

//  Finale(nkf1,nkf2,nkf3,nplane_x,nplane_x_rho,cp_min_opt);

  //===================================================================================
  // Output your parameter choices to the screen

  // The cpaimd config output file is written in the current directory (and not where the input file is located)
  sprintf(fname,"%s.out",input_name);
  fp = cfopen((const char*) fname,"w");
  write_cpaimd_config(fp, dict_gen_GW,  num_dict_gen_GW,  dict_fun[1].keyword);
  write_cpaimd_config(fp, dict_GW_epsilon, num_dict_GW_epsilon, dict_fun[2].keyword);
  write_cpaimd_config(fp, dict_GW_sigma, num_dict_GW_sigma, dict_fun[3].keyword);
  //  write_cpaimd_config(fp, dict_GW_file, num_dict_GW_file, dict_fun[4].keyword);
//  write_cpaimd_config(fp,dict_GW_charm_input,   num_dict_GW_charm_input,   dict_fun[5].keyword);
  fclose(fp);
  //===================================================================================
  // Free memory : 

  cfree(fun_key,"Config::readCconfig");  
  cfree(fname  ,"Config::readCconfig");  

  cfree(&dict_fun[1], "Config::readCconfig");
  cfree(&dict_gen_GW[1], "Config::readCconfig");
  cfree(&dict_GW_epsilon[1], "Config::readCconfig");
  cfree(&dict_GW_sigma[1], "Config::readCconfig"); 
  //  cfree(&dict_GW_file[1]   ,"Config::readCconfig");
  //  cfree(&dict_GW_charm_input[1]   ,"Config::readCconfig");

  //===================================================================================
  // Tell Everyone you are done

  PRINTF("  -------------------------------------------------------------\n");
  PRINTF("  Completed reading GWBSE input from file : %s\n",input_name);
  PRINTF("  =============================================================\n\n");



  ///==================================================================================
  // read cell and k list
  read_lattice( gwbseopts );
  read_klist( gwbseopts );



  


  //----------------------------------------------------------------------------------
}//end routine
//===================================================================================












//===================================================================================
/* set_config_dict_fun sets metakeywords */
//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_fun  (int *num_dict  ,DICT_WORD **dict){
  //==================================================================================
  //  I) Malloc the dictionary                                              

  num_dict[0] = 5;   // how many functional keywords you have
  *dict = (DICT_WORD *)cmalloc(num_dict[0]*sizeof(DICT_WORD),"set_config_dict_fun")-1;

  //=================================================================================
  //  II) Initialize the user set option(did the user set the key word      

  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iuset    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iflag    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].key_type = 1;} //spec only once

  //=================================================================================
  // III) Set up the dictionary
  int ind;
  //------------------------------------------------------------------------------
  //  1)~gen_GW[ ]
  ind = 1;
  strcpy((*dict)[ind].error_mes," ");
  strcpy((*dict)[ind].keyword,"gen_GW");
  strcpy((*dict)[ind].keyarg," ");
  //------------------------------------------------------------------------------
  //  2)~GW_epsilon[ ]
  ind = 2;
  strcpy((*dict)[ind].error_mes," ");
  strcpy((*dict)[ind].keyword,"GW_epsilon");
  strcpy((*dict)[ind].keyarg," ");
  //------------------------------------------------------------------------------
  //  3)~GW_sigma[ ]
  ind = 3;
  strcpy((*dict)[ind].error_mes," ");
  strcpy((*dict)[ind].keyword,"GW_sigma");
  strcpy((*dict)[ind].keyarg," ");  
  //------------------------------------------------------------------------------
  //  4)~GW_filenames[ ]
  ind = 4;
  strcpy((*dict)[ind].error_mes," ");
  strcpy((*dict)[ind].keyword,"GW_filenames");
  strcpy((*dict)[ind].keyarg," ");
  //------------------------------------------------------------------------------
  //  5)~GW_charm_input[ ]
  ind = 5;
  strcpy((*dict)[ind].error_mes," ");
  strcpy((*dict)[ind].keyword,"GW_charm_input");
  strcpy((*dict)[ind].keyarg," ");
  //------------------------------------------------------------------------------
}//end routine
//===================================================================================






//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_gen_GW  (int *num_dict ,DICT_WORD **dict){
  //==================================================================================
  //  I) Malloc the dictionary                                              

  num_dict[0] = 8;
  *dict = (DICT_WORD *)cmalloc(num_dict[0]*sizeof(DICT_WORD),"set_dict_gen_GW")-1;

  //=================================================================================
  //  II) Initialize the user set option(did the user set the key word      

  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iuset    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iflag    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].key_type = 1;} //spec only once

  //=================================================================================
  // III) Set up the dictionary
  int ind;
  //-----------------------------------------------------------------------------
  //  1)\gwbse_opt{}
  ind =   1;
  strcpy((*dict)[ind].keyword,"gwbse_opt");
  strcpy((*dict)[ind].keyarg,"off");
  strcpy((*dict)[ind].error_mes,"on or off");
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  //  2)\num_tot_state{}
  ind =   2;   
  strcpy((*dict)[ind].keyword,"num_tot_state");
  strcpy((*dict)[ind].keyarg,"2");
  strcpy((*dict)[ind].error_mes,"an integer > 1 ");
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  //  3)\num_occ_state{}
  ind =   3;   
  strcpy((*dict)[ind].keyword,"num_occ_state");
  strcpy((*dict)[ind].keyarg,"1");
  strcpy((*dict)[ind].error_mes,"an integer > 0 ");
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  //  4)\num_unocc_state{}
  ind =   4;   
  strcpy((*dict)[ind].keyword,"num_unocc_state");
  strcpy((*dict)[ind].keyarg,"1");
  strcpy((*dict)[ind].error_mes,"an integer > 0");
  //-----------------------------------------------------------------------------  

  //-----------------------------------------------------------------------------
  //  5)\num_kpoint{}
  ind =   5;
  strcpy((*dict)[ind].keyword,"num_kpoint");
  strcpy((*dict)[ind].keyarg,"1");
  strcpy((*dict)[ind].error_mes,"an integer > 0");
  //----------------------------------------------------------------------------- 

  //-----------------------------------------------------------------------------
  //  6)\num_spin{}
  ind =   6;
  strcpy((*dict)[ind].keyword,"num_spin");
  strcpy((*dict)[ind].keyarg,"1");
  strcpy((*dict)[ind].error_mes,"an integer > 0");
  //----------------------------------------------------------------------------- 

  //-----------------------------------------------------------------------------
  //  7)\coulb_trunc_opt{}
  ind =   7;
  strcpy((*dict)[ind].keyword,"coulb_trunc_opt");
  strcpy((*dict)[ind].keyarg,"0");
  strcpy((*dict)[ind].error_mes,"0-no truncation, 1-wire, 2-sheet, 3-molecule");
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  //  8)\statefile_binary_opt{}
  ind =   8;
  strcpy((*dict)[ind].keyword,"statefile_binary_opt");
  strcpy((*dict)[ind].keyarg,"off");
  strcpy((*dict)[ind].error_mes,"off,off_gzip,on,on_gzip");

}//end routine
//===================================================================================





//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_GW_epsilon  (int *num_dict ,DICT_WORD **dict){
  //==================================================================================
  //  I) Malloc the dictionary                                              
  num_dict[0] = 2;
  *dict = (DICT_WORD *)cmalloc(num_dict[0]*sizeof(DICT_WORD),"set_dict_gen_GW")-1;

  //=================================================================================
  //  II) Initialize the user set option(did the user set the key word      

  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iuset    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iflag    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].key_type = 1;} //spec only once

  //=================================================================================
  // III) Set up the dictionary
  int ind;
  //-----------------------------------------------------------------------------
  //  1)\Ecuteps{}
  ind =   1;
  strcpy((*dict)[ind].keyword,"Ecuteps");
  strcpy((*dict)[ind].keyarg,"10");
  strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  //  2)\tol_iter_mtxinv{}
  ind =   2;   
  strcpy((*dict)[ind].keyword,"tol_iter_mtxinv");
  strcpy((*dict)[ind].keyarg,"0.001");
  strcpy((*dict)[ind].error_mes,"an integer > 1 ");
  //-----------------------------------------------------------------------------

  
}//end routine
//===================================================================================





//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_GW_sigma  (int *num_dict ,DICT_WORD **dict){
  //==================================================================================
  //  I) Malloc the dictionary                                              
  num_dict[0] = 5;
  *dict = (DICT_WORD *)cmalloc(num_dict[0]*sizeof(DICT_WORD),"set_dict_gen_GW")-1;

  //=================================================================================
  //  II) Initialize the user set option(did the user set the key word      

  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iuset    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iflag    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].key_type = 1;} //spec only once

  //=================================================================================
  // III) Set up the dictionary
  int ind;
  //-----------------------------------------------------------------------------
  //  1)\PP_num_mode{}
  ind =   1;
  strcpy((*dict)[ind].keyword,"PP_num_mode");
  strcpy((*dict)[ind].keyarg,"1");
  strcpy((*dict)[ind].error_mes,"an integer > 0");
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  //  2)\band_index_min{}
  ind =   2;   
  strcpy((*dict)[ind].keyword,"band_index_min");
  strcpy((*dict)[ind].keyarg,"1");
  strcpy((*dict)[ind].error_mes,"an integer > 0 ");
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  //  3)\band_index_max{}
  ind =   3;   
  strcpy((*dict)[ind].keyword,"band_index_max");
  strcpy((*dict)[ind].keyarg,"2");
  strcpy((*dict)[ind].error_mes,"an integer > 1 ");
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  //  4)\screened_coulomb_cutoff{}
  ind =   4;   
  strcpy((*dict)[ind].keyword,"screened_coulomb_cutoff");
  strcpy((*dict)[ind].keyarg,"1");
  strcpy((*dict)[ind].error_mes,"a number > 0 ");
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  //  5)\screened_coulomb_cutoff{}
  ind =   5;   
  strcpy((*dict)[ind].keyword,"bare_coulomb_cutoff");
  strcpy((*dict)[ind].keyarg,"1");
  strcpy((*dict)[ind].error_mes,"a number > 0 ");
  //-----------------------------------------------------------------------------
}//end routine
//===================================================================================








//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_GW_file  (int *num_dict ,DICT_WORD **dict){
  //==================================================================================
  //  I) Malloc the dictionary                                              
  num_dict[0] = 5;
  *dict = (DICT_WORD *)cmalloc(num_dict[0]*sizeof(DICT_WORD),"set_dict_gen_GW")-1;

  //=================================================================================
  //  II) Initialize the user set option(did the user set the key word      

  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iuset    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iflag    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].key_type = 1;} //spec only once

  //=================================================================================
  // III) Set up the dictionary
  int ind;
  //-----------------------------------------------------------------------------
  //  1)\read_kpt_cell{}
  ind =   1;
  strcpy((*dict)[ind].keyword,"read_kpt_cell");
  strcpy((*dict)[ind].keyarg,"");
  strcpy((*dict)[ind].error_mes,"a file name that reads kpoints and cell info");
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  //  2)\write_EpsMatInv{}
  ind =   2;   
  strcpy((*dict)[ind].keyword,"write_EpsMatInv");
  strcpy((*dict)[ind].keyarg,"off");
  strcpy((*dict)[ind].error_mes,"on or off");
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  //  3)\restart_EpsMatInv{}
  ind =   3;   
  strcpy((*dict)[ind].keyword,"restart_EpsMatInv");
  strcpy((*dict)[ind].keyarg,"EPSMATINV");
  strcpy((*dict)[ind].error_mes,"a file name that reads epsilon matrix inverse");
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  //  4)\read_density{}
  ind =   4;   
  strcpy((*dict)[ind].keyword,"read_density");
  strcpy((*dict)[ind].keyarg,"RHO");
  strcpy((*dict)[ind].error_mes,"a file name that reads density for GPP calculations");
  //-----------------------------------------------------------------------------
}//end routine
//===================================================================================







//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_GW_charm_input  (int *num_dict ,DICT_WORD **dict){
  //==================================================================================
  //  I) Malloc the dictionary                                              
  num_dict[0] = 5;
  *dict = (DICT_WORD *)cmalloc(num_dict[0]*sizeof(DICT_WORD),"set_dict_gen_GW")-1;

  //=================================================================================
  //  II) Initialize the user set option(did the user set the key word      

  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iuset    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iflag    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].key_type = 1;} //spec only once

  //=================================================================================
  // III) Set up the dictionary
  int ind;

}//end routine
//===================================================================================





//===================================================================================
/* set parameters to the class variables */
//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_params_gen_GW  (DICT_WORD *dict, char *fun_key, char *input_name, GWBSEOPTS *gwbseopts){
  //===================================================================================
  double real_arg;
  int int_arg;
  int ifound;
  int ind,ierr;

  //===================================================================================
  // Fill me with joy.

  //-----------------------------------------------------------------------------
  //  1)\gwbse_opt{}
  ind =   1;   
  parse_on_off(dict[ind].keyarg,&int_arg,&ierr);
  if (ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  gwbseopts->gwbse_on = (int_arg==1?true:false);  //if int_arg ==1, then true 

  //-----------------------------------------------------------------------------
  //  2)\num_tot_state{}
  ind =   2;   
  sscanf(dict[ind].keyarg,"%d",&int_arg);
  if (int_arg<2){keyarg_barf(dict,input_name,fun_key,ind);}  
  gwbseopts->nstate = int_arg;
  
  //-----------------------------------------------------------------------------
  //  3)\num_occ_state{}
  ind =   3;   
  sscanf(dict[ind].keyarg,"%d",&int_arg);
  if (int_arg<1){keyarg_barf(dict,input_name,fun_key,ind);}     
  gwbseopts->nocc = int_arg;

  //-----------------------------------------------------------------------------
  //  4)\num_unocc_state{}
  ind =   4;   
  sscanf(dict[ind].keyarg,"%d",&int_arg);
  if (int_arg<0){keyarg_barf(dict,input_name,fun_key,ind);}     
  gwbseopts->nunocc = int_arg;

  //-----------------------------------------------------------------------------
  //  5)\num_kpoint{}
  ind =   5;   
  sscanf(dict[ind].keyarg,"%d",&int_arg);
  if (int_arg<1){keyarg_barf(dict,input_name,fun_key,ind);}     
  gwbseopts->nkpt = int_arg;

  //-----------------------------------------------------------------------------
  //  6)\num_spin{}
  ind =   6;   
  sscanf(dict[ind].keyarg,"%d",&int_arg);
  if (int_arg<1){keyarg_barf(dict,input_name,fun_key,ind);}     
  gwbseopts->nspin = int_arg;

  //-----------------------------------------------------------------------------
  //  7)\coulb_trunc_opt{}
  ind =   7;   
  sscanf(dict[ind].keyarg,"%d",&int_arg);
  if (int_arg<0){keyarg_barf(dict,input_name,fun_key,ind);}     
  gwbseopts->coulb_trunc_opt = int_arg;

  //-----------------------------------------------------------------------------
  //  7)\statefile_binary_opt{}
  ind =   8;   

  ifound = 0;
  if(strcasecmp(dict[ind].keyarg,"off")==0)
    {gwbseopts->ibinary_opt = 0;ifound++;}
  if(strcasecmp(dict[ind].keyarg,"on")==0)
    {gwbseopts->ibinary_opt = 1;ifound++;}
  if(strcasecmp(dict[ind].keyarg,"off_gzip")==0)
    {gwbseopts->ibinary_opt = 2;ifound++;}
  if(strcasecmp(dict[ind].keyarg,"on_gzip")==0)
    {gwbseopts->ibinary_opt = 3;ifound++;}
  if(ifound == 0) keyarg_barf(dict,input_name,fun_key,ind);
  

  //----------------------------------------------------------------------------- 
  // Die if the number of states don't add up
  if(gwbseopts->nstate!=(gwbseopts->nunocc+gwbseopts->nocc)){
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("   total number of states %d is not equal to conduction %d + valence %d \n",
           gwbseopts->nstate, gwbseopts->nunocc, gwbseopts->nocc);
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }// end if

  //----------------------------------------------------------------------------- 
}// end routine
//================================================================================






//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_params_GW_epsilon  (DICT_WORD *dict, char *fun_key, char *input_name, GW_EPSILON *gw_epsilon){
  //===================================================================================
  double real_arg;
  int int_arg;

  int ind,ierr;

  //===================================================================================
  // Fill me with joy.

  //-----------------------------------------------------------------------------
  //  1)\Ecuteps{}
  ind =   1;   
  sscanf(dict[ind].keyarg,"%lg",&real_arg);
  if (real_arg<0){keyarg_barf(dict,input_name,fun_key,ind);}
  gw_epsilon->Ecuteps =real_arg; 

  //-----------------------------------------------------------------------------
  //  2)\tol_iter_mtxinv{}
  ind =   2;   
  sscanf(dict[ind].keyarg,"%lg",&real_arg);
  if (real_arg<0){keyarg_barf(dict,input_name,fun_key,ind);}  
  gw_epsilon->tol_iter = real_arg;
  //----------------------------------------------------------------------------- 
}// end routine
//================================================================================





//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_params_GW_sigma  (DICT_WORD *dict, char *fun_key, char *input_name, GW_SIGMA *gw_sigma){
  //===================================================================================
  double real_arg;
  int int_arg;

  int ind,ierr;

  //===================================================================================
  // Fill me with joy.

  //-----------------------------------------------------------------------------
  //  1)\PP_num_mode{}
  ind =   1;   
  sscanf(dict[ind].keyarg,"%d",&int_arg);
  if (int_arg<0){keyarg_barf(dict,input_name,fun_key,ind);}
  gw_sigma->PP_nmode = int_arg; 

  //-----------------------------------------------------------------------------
  //  2)\band_index_min{}
  ind =   2;   
  sscanf(dict[ind].keyarg,"%d",&int_arg);
  if (int_arg<0){keyarg_barf(dict,input_name,fun_key,ind);}  
  gw_sigma->band_index_min = int_arg;
  //-----------------------------------------------------------------------------
  //  3)\band_index_max{}
  ind =   3;   
  sscanf(dict[ind].keyarg,"%d",&int_arg);
  if (int_arg<1){keyarg_barf(dict,input_name,fun_key,ind);}
  gw_sigma->band_index_max = int_arg; 

  //-----------------------------------------------------------------------------
  //  4)\screened_coulomb_cutoff{}
  ind =   4;   
  sscanf(dict[ind].keyarg,"%lg",&real_arg);
  if (real_arg<0){keyarg_barf(dict,input_name,fun_key,ind);}  
  gw_sigma->screened_coulomb_cutoff = real_arg;

  //-----------------------------------------------------------------------------
  //  5)\bare_coulomb_cutoff{}
  ind =   5;   
  sscanf(dict[ind].keyarg,"%lg",&real_arg);
  if (real_arg<0){keyarg_barf(dict,input_name,fun_key,ind);}  
  gw_sigma->bare_coulomb_cutoff = real_arg;
  //----------------------------------------------------------------------------- 
}// end routine
//================================================================================




//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_params_GW_file  (DICT_WORD *dict, char *fun_key, char *input_name, GWBSEOPTS *gwbseopts){
  //===================================================================================
  double real_arg;
  int int_arg;

  int ind,ierr;

  //===================================================================================
  // Fill me with joy.

  //-----------------------------------------------------------------------------
  //  1)\write_EpsMatInv{}
  ind =   2;
  parse_on_off(dict[ind].keyarg,&int_arg,&ierr);
  if (ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  gwbseopts->write_epsmatinv = (int_arg==1?true:false);
 
  //-----------------------------------------------------------------------------
  //  2)\restart_EpsMatInv{}
  ind =   3;   
  sscanf(dict[ind].keyarg,"%s",gwbseopts->fileEpsMatInv);
  /* put here the restart option 
  if (int_arg<1){keyarg_barf(dict,input_name,fun_key,ind);}
  */
  
  //-----------------------------------------------------------------------------
  //  3)\read_density{}
  ind =   4;   
  sscanf(dict[ind].keyarg,"%s",gwbseopts->fileRho);
  /* put here sigma option to check the error
  if (real_arg<0){keyarg_barf(dict,input_name,fun_key,ind);}  
  */
  
  //----------------------------------------------------------------------------- 
}// end routine
//================================================================================










//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::readStateInfo(int &nPacked,int &minx, int &maxx, int &nx, int &ny, int &nz,
    const char *fromFile, int ibinary_opt,
    int nkf1, int nkf2, int nkf3,int ncoef) {
  //===================================================================================
  // Check for errors
#ifdef when_we_are_ready
  if(ibinary_opt < 0 || ibinary_opt > 3){
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("   Bad binary option\n",ibinary_opt);
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

#if !CMK_PROJECTIONS_USE_ZLIB
  if(ibinary_opt>1){
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("   Attempt to use ZLIB Failed! Please review compilation\n");
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif
#endif

  //===================================================================================
  // Read the file

  int nktot,nplane0,n=1;

  if(gen_wave==0){
    //---------------------------------------------------------------------------------
    // Ascii
    if(ibinary_opt==0){

      FILE *fp=fopen(fromFile,"r");
      if(fp==NULL){
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("   Can't open state file :%s", fromFile);
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        EXIT(1);
      }//endif
      if(4!=fscanf(fp,"%d%d%d%d",&nPacked,&nx,&ny,&nz)){
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("   Can't parse size line of file %s\n", fromFile);
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        EXIT(1);
      }
      nktot=0;
      nplane0=0;
      for(int pNo=0;pNo<nPacked;pNo++) {
        double re,im; int x,y,z;
        if(5!=fscanf(fp,"%lg%lg%d%d%d",&re,&im,&x,&y,&z)){
          PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          PRINTF("   Can't parse packed state location");
          PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          EXIT(1);
        }
        if(pNo==0){minx=x; maxx=x;}
        if(x<minx){minx=x;}
        if(x>maxx){maxx=x;}
        if(x==0){nplane0++;}
        nktot++;
        if(x==0 && y==0 && z==0 && doublePack)break;
      }//endfor
      fclose(fp);

    }//endif:: acii

#if CMK_PROJECTIONS_USE_ZLIB
    //---------------------------------------------------------------------------------
    // Zipped ascii
    if(ibinary_opt==2){

      char bigenough[1000];  //we know our lines are shorter than this
      char localFile[1000]; // fromFile is const
      strncpy(localFile,fromFile,1000);
      strncat(localFile,".gz",1000);
      gzFile zfp=gzopen(localFile,"rb");
      if (zfp==NULL){
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("   Can't open state file %s\n",localFile);
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        EXIT(1);
      }//endif
      int nPackedLoc;
      if(gzgets(zfp,bigenough,1000)!=Z_NULL){
        if(4!=sscanf(bigenough,"%d%d%d%d",&nPacked,&nx,&ny,&nz)){
          PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          PRINTF("   Can't parse size line of file %s\n", localFile);
          PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          EXIT(1);
        }
      }else{
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("   Can't parse size line of file %s\n", localFile);
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        EXIT(1);
      }//endif
      nktot=0;
      nplane0=0;
      for(int pNo=0;pNo<nPacked;pNo++) {
        double re,im; int x,y,z;
        if(gzgets(zfp,bigenough,1000)!=Z_NULL){
          if(5!=sscanf(bigenough,"%lg%lg%d%d%d",&re,&im,&x,&y,&z)){
            PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
            PRINTF("   Can't parse packed state location");
            PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
            EXIT(1);
          }
          if(pNo==0){minx=x; maxx=x;}
          if(x<minx){minx=x;}
          if(x>maxx){maxx=x;}
          if(x==0){nplane0++;}
          nktot++;
          if(x==0 && y==0 && z==0 && doublePack)break;
        }else{
          PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          PRINTF("   Can't parse size line of file %s\n", localFile);
          PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          EXIT(1);
        }//endif
      }//endfor
      gzclose(zfp);

    }//endif:: zipped ascii

    //---------------------------------------------------------------------------------
    // Zipped binary
    if(ibinary_opt==3){

      char localFile[1000]; // fromFile is const
      strncpy(localFile,fromFile,1000);
      strncat(localFile,".gz",1000);
      gzFile zfp=gzopen(localFile,"rb");
      if (zfp==NULL){
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("   Can't open state file :%s", localFile);
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        EXIT(1);
      }//endif
      gzread(zfp,&(nPacked),sizeof(int));
      gzread(zfp,&(nx),sizeof(int));
      gzread(zfp,&(ny),sizeof(int));
      gzread(zfp,&(nz),sizeof(int));
      nktot=0;
      nplane0=0;
      for(int pNo=0;pNo<nPacked;pNo++) {
        double re,im; int x,y,z;
        gzread(zfp,&(re),sizeof(double));
        gzread(zfp,&(im),sizeof(double));
        gzread(zfp,&(x),sizeof(int));
        gzread(zfp,&(y),sizeof(int));
        gzread(zfp,&(z),sizeof(int));
        if(pNo==0){minx=x; maxx=x;}
        if(x<minx){minx=x;}
        if(x>maxx){maxx=x;}
        if(x==0){nplane0++;}
        nktot++;
        if(x==0 && y==0 && z==0 && doublePack)break;
      }//endfor
      gzclose(zfp);

    }//endif :: zipped binary
#endif

    //---------------------------------------------------------------------------------
    // Binary
    if(ibinary_opt==1){

      FILE *fp=fopen(fromFile,"rb");
      if (fp==NULL){
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("   Can't open state file :%s", fromFile);
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        EXIT(1);
      }//endif
      if(fread(&(nPacked),sizeof(int),n,fp)){}
      if(fread(&(nx),sizeof(int),n,fp)){}
      if(fread(&(ny),sizeof(int),n,fp)){}
      if(fread(&(nz),sizeof(int),n,fp)){}
      nktot=0;
      nplane0=0;
      for(int pNo=0;pNo<nPacked;pNo++) {
        double re,im; int x,y,z;
        if(fread(&(re),sizeof(double),n,fp)){}
        if(fread(&(im),sizeof(double),n,fp)){}
        if(fread(&(x),sizeof(int),n,fp)){}
        if(fread(&(y),sizeof(int),n,fp)){}
        if(fread(&(z),sizeof(int),n,fp)){}
        if(pNo==0){minx=x; maxx=x;}
        if(x<minx){minx=x;}
        if(x>maxx){maxx=x;}
        if(x==0){nplane0++;}
        nktot++;
        if(x==0 && y==0 && z==0 && doublePack)break;
      }//endfor
      fclose(fp);

    }//endif::binary
  }//endif::gen_wave

  //===================================================================================
  // If we are generating the wave function from scratch

  if(gen_wave==1){
    nx    = nkf1;
    ny    = nkf2;
    nz    = nkf3;
    nktot = ncoef;

    int *ka       = (int *)cmalloc(nktot*sizeof(int),"parainfo")-1;
    int *kb       = (int *)cmalloc(nktot*sizeof(int),"parainfo")-1;
    int *kc       = (int *)cmalloc(nktot*sizeof(int),"parainfo")-1;

    PhysicsParamTransfer::fetch_state_kvecs(ka,kb,kc,nktot,doublePack);

    nplane0 = 0;
    minx    = ka[1]; 
    maxx    = ka[1];
    for(int i=1;i<=nktot;i++){
      if(ka[i]<minx){minx=ka[i];}
      if(ka[i]>maxx){maxx=ka[i];}
      if(ka[i]==0){nplane0++;}
    }//endfor

    cfree(&ka[1],"configures.C"); 
    cfree(&kb[1],"configures.C"); 
    cfree(&kc[1],"configures.C"); 

  }//endif

  //===================================================================================
  // Set a few parameters before you go home

  if(doublePack){nPacked=nktot+nplane0-1;}
  else {nPacked = nktot;}
//-------------------------------------------------
#endif //when_we_are_ready
  //----------------------------------------------------------------------------------
}//end routine
//===================================================================================


//===================================================================================
// Some simple range checking : needs some love
//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================

//============================================================================
// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Config::rangeExit(int param, const char *name, int iopt){
  //============================================================================

  switch(iopt){
    case 0: 
      if(param<1){
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("   The parameter %s must be >0 not %d \n",name,param);
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        EXIT(1);
      }//endif
      break;
    case 1: 
      if(param<0 || param>1){
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("   The parameter %s must be 1(on) or 0 (off) \n",name);
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        EXIT(1);
      }//endif
      break;
  }//end switch

  //----------------------------------------------------------------------------------
}//end routine
//===================================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void Config::write_cpaimd_config(FILE *fp,DICT_WORD *dict, int num_dict, char *fun_key){
  //========================================================================

  int iuset,itype;

  //========================================================================
  //     I) Write out meta key word 

  fprintf(fp,"================================================\n");
  fprintf(fp,"cccccccccccccccccccccccccccccccccccccccccccccccc\n");
  fprintf(fp,"================================================\n");
  fprintf(fp,"~%s[\n",fun_key);

  //================================================================================  
  //    II) User defined parameters 

  iuset = 1;itype = 1;
  fprintf(fp,"------------------------------------------------\n");
  fprintf(fp,"  User defined parameters with software overides\n");
  fprintf(fp,"------------------------------------------------\n");
  dict_print(fp,num_dict,dict,itype,iuset);

  //================================================================================  
  //   III) Default parameters 

  iuset = 0;itype = 1;
  fprintf(fp,"------------------------------------------------\n");
  fprintf(fp,"  Default parameters with software overides\n");
  fprintf(fp,"------------------------------------------------\n");
  dict_print(fp,num_dict,dict,itype,iuset);

  //================================================================================
  //   IV) End Meta key word  output

  fprintf(fp,"------------------------------------------------\n]\n");
  fprintf(fp,"================================================\n\n\n");

  //========================================================================
}//end routine 
//========================================================================



//==========================================================================
// Additional routine for GW
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void Config::read_lattice(GWBSEOPTS *gwbseopts){

  int i;
  double h1, h2, h3;
  FILE *fp;
  
  // open lattice.dat file. The name is hard-coded!
  fp = cfopen("lattice.dat","r");

  // first line: lattice constant
  //  fscanf(fp,"%lg",&(gwbseopts->latt));
  
  double lc;
  fscanf(fp,"%lg",&lc);
  gwbseopts->latt = lc;
  
  // read the rest
  for(i=0; i<3; i++){
    fscanf(fp,"%lg %lg %lg",&h1,&h2,&h3);
    gwbseopts->h[i*3]=h1*lc;
    gwbseopts->h[i*3+1]=h2*lc;
    gwbseopts->h[i*3+2]=h3*lc;
  }
  fclose(fp);
  
}//end routine 
//========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void Config::read_klist(GWBSEOPTS *gwbseopts){

  int nk;
  int i;
  double kx, ky, kz, wt;
  FILE *fp;
  
  // open klist.dat file. The name is hard-coded!
  fp = cfopen("klist.dat","r");
  fscanf(fp,"%d",&nk);

  // this is necessary to malloc 
  double ** kvec   = new double* [nk];
  for(int i=0; i<nk ; i++){
    kvec[i] = new double [3];
  }
  double *kwt = new double [nk];

  gwbseopts->kvec = kvec;
  gwbseopts->kwt = kwt;


  // good to check if nk is the same as gwbseopts.nkpt
  if(nk!=gwbseopts->nkpt){
    printf("The number of k points in the file klist.dat doesn't match with the number specified %d vs. %d.\n",nk,gwbseopts->nkpt);
    EXIT(1);
  }
  
  for(i=0; i<nk; i++){
    fscanf(fp,"%lg %lg %lg %lg",&kx,&ky,&kz,&wt);
    gwbseopts->kvec[i][0] = kx;
    gwbseopts->kvec[i][1] = ky;
    gwbseopts->kvec[i][2] = kz;
    gwbseopts->kwt[i] = wt;
  }
  fclose(fp);
  
}//end routine 
//========================================================================
