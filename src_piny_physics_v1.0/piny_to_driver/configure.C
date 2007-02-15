//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                         PI_MD:                                           
//             The future of simulation technology                          
//             ------------------------------------                         
//                     Module: control_sim_parms.c                          
//                                                                          
//                                                                          
// This subprogram reads in user simulation params and echoes them          
// and the default parameters to a file                                     
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#include "standard_include.h"
#include "../../include/configure.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_handle_entry.h"
#if CMK_PROJECTIONS_USE_ZLIB
#include "zlib.h"
#endif

//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::readConfig(char* input_name,int nstates_in, int nkf1, int nkf2, int nkf3, 
                        int maxIter_in,int ibinary_opt,int natm_nl_in, int fftopt_in,
                        int numPes_in, int natm_typ_in,int ees_eext_opt_in)
//===================================================================================
   {//begin routine
//===================================================================================

  int num_dict_fun;
  int num_dict_rho, num_dict_state, num_dict_pc;
  int num_dict_nl,  num_dict_gen,   num_dict_atm;
  DICT_WORD *dict_fun;
  DICT_WORD *dict_rho, *dict_state, *dict_pc;
  DICT_WORD *dict_nl, *dict_gen, *dict_atm;
  DICT_WORD word;            

  int nline;
  int nkey,nfun_key;
  char *fun_key,*fname;
  int ind_key;
  FILE *fp;
  int iii;

  fun_key = (char *)cmalloc(PINY_MAXWORD*sizeof(char),"Config::readConfig");
  fname   = (char *)cmalloc(1024*sizeof(char),"Config::readConfig");

//===================================================================================
// Tell everyone you are busy 

  PRINTF("  Reading charm parallel input from file : %s\n",input_name);

  if(PINY_MAXWORD!=MAXWORD || PINY_MAXLINE != MAXLINE){
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("   Incorrect word and line sizes\n");
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

//===================================================================================
// Set some parameters directly from PINY input

  nstates      = nstates_in;
  numFFTPoints = nkf1 * nkf2 * nkf3;
  maxIter      = maxIter_in;
  fftopt       = fftopt_in;
  sGrainSize   = nstates_in;
  natm_nl      = natm_nl_in;
  numPes       = numPes_in;
  natm_typ     = natm_typ_in;
  ees_eext_opt = ees_eext_opt_in;

//===================================================================================
// Set up the dictionaries

  useCommlib = 1;  //default value
  useTimeKeeper = 0;  //default value
  set_config_dict_fun  (&num_dict_fun  ,&dict_fun);
  set_config_dict_rho  (&num_dict_rho  ,&dict_rho);
  set_config_dict_state(&num_dict_state,&dict_state);
  set_config_dict_pc   (&num_dict_pc,   &dict_pc);
  set_config_dict_nl   (&num_dict_nl,   &dict_nl);
  set_config_dict_gen  (&num_dict_gen,  &dict_gen);
  set_config_dict_atm  (&num_dict_atm,  &dict_atm);

//===================================================================================
// Read the input file and fill the dictionaries with user input

  fp = cfopen(input_name,"r");

  nline = 1;  nfun_key=0;
  while(get_fun_key(fp,fun_key,&nline,&nfun_key,input_name)){
    get_fun_key_index(fun_key,num_dict_fun,dict_fun,nline,nfun_key,
                      input_name,&ind_key);
    nkey  = 0;
    while(get_word(fp,&word,&nline,&nkey,nfun_key,input_name)){
      switch(ind_key){
        case 1 : put_word_dict(&word,dict_rho,num_dict_rho,fun_key,nline,
                               nkey,nfun_key,input_name);break;
        case 2 : put_word_dict(&word,dict_state,num_dict_state,fun_key,nline,
                               nkey,nfun_key,input_name);break;
        case 3 : put_word_dict(&word,dict_pc,num_dict_pc,fun_key,nline,
                               nkey,nfun_key,input_name);break;
        case 4 : put_word_dict(&word,dict_nl,num_dict_nl,fun_key,nline,
                               nkey,nfun_key,input_name);break;
        case 5 : put_word_dict(&word,dict_gen,num_dict_gen,fun_key,nline,
                               nkey,nfun_key,input_name);break;
        case 6 : put_word_dict(&word,dict_atm,num_dict_atm,fun_key,nline,
                               nkey,nfun_key,input_name);break;
      }//end switch
    }// end while 
  }//end while

  fclose(fp);  

//===================================================================================
// Take the information out of the dictionary and put it in the class

  set_config_params_gen  (dict_gen,  dict_fun[5].keyword,input_name);
  int iflag = 1; if(useCommlib!=1){iflag=0;}  // change default commlib option

  set_config_params_rho  (dict_rho,  dict_fun[1].keyword,input_name,iflag);
  set_config_params_state(dict_state,dict_fun[2].keyword,input_name,iflag);
  set_config_params_pc   (dict_pc,   dict_fun[3].keyword,input_name);
  set_config_params_nl   (dict_nl,   dict_fun[4].keyword,input_name,iflag);
  set_config_params_atm  (dict_atm,  dict_fun[6].keyword,input_name);

  simpleRangeCheck(); // redundant checking

//===================================================================================
// Read some info from the state file and check for consistency with PINY

  int sizex,sizey,sizez,nPacked,minx,maxx;

  sprintf (fname, "%s/state1.out", dataPath);
  PRINTF("   Opening state file : %s\n",fname);
    readStateInfo(nPacked,minx,maxx,sizex,sizey,sizez,fname,ibinary_opt);
  PRINTF("   Closing state file : %s\n\n",fname);

  if(sizex!=nkf1 || sizey!=nkf2 || sizez !=nkf3){
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("   Incorrect FFT size in state files.\n");
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

  numFFTPoints = nkf1 * nkf2 * nkf3;
  low_x_size   = minx+1;
  high_x_size  = maxx-1;
  numData      = nPacked;

//===================================================================================
// Set rhoG and stateG chare array sizes

  int nplane_x     = minx+1;
  int nplane_x_rho = 2*minx+1;

  double temp      = (gExpandFact)*((double)nplane_x);
  nchareG          = ((int)temp);

  double temp_rho  = (gExpandFactRho)*((double)nplane_x_rho);
  nchareRhoG       = ((int)temp_rho);

  scalc_per_plane  = (nstates/sGrainSize)*(nstates/sGrainSize);

//===================================================================================
// Improve user parameters and/or try to optimize unset parameters

  guesstimateParmsConfig(sizez,dict_gen,dict_rho,dict_state,dict_pc,dict_nl,dict_atm);

//===================================================================================
// Final consistency checks

  Finale(nkf1,nkf2,nkf3,nplane_x,nplane_x_rho);

//===================================================================================
// Output your parameter choices to the screen

  sprintf(fname,"%s.out",input_name);
  fp = cfopen(fname,"w");
   write_cpaimd_config(fp,dict_rho,  num_dict_rho,  dict_fun[1].keyword);
   write_cpaimd_config(fp,dict_state,num_dict_state,dict_fun[2].keyword);
   write_cpaimd_config(fp,dict_pc,   num_dict_pc,   dict_fun[3].keyword);
   write_cpaimd_config(fp,dict_nl,   num_dict_nl,   dict_fun[4].keyword);
   write_cpaimd_config(fp,dict_gen,  num_dict_gen,  dict_fun[5].keyword);
   write_cpaimd_config(fp,dict_atm,  num_dict_atm,  dict_fun[6].keyword);
  fclose(fp);

//===================================================================================
// Free memory : 

  cfree(fun_key,"Config::readCconfig");  
  cfree(fname  ,"Config::readCconfig");  

  cfree(&dict_fun[1]  ,"Config::readCconfig");
  cfree(&dict_rho[1]  ,"Config::readCconfig");
  cfree(&dict_state[1],"Config::readCconfig");
  cfree(&dict_pc[1]   ,"Config::readCconfig");
  cfree(&dict_nl[1]   ,"Config::readCconfig");
  cfree(&dict_gen[1]  ,"Config::readCconfig");

//===================================================================================
// Tell Everyone you are done

  PRINTF("  Completed reading charm parallel input from file : %s\n",input_name);

//----------------------------------------------------------------------------------
   }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_fun  (int *num_dict  ,DICT_WORD **dict){
//==================================================================================
//  I) Malloc the dictionary                                              

  num_dict[0] = 6;
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
  //  1)~charm_conf_rho_def[ ]
    ind = 1;
    strcpy((*dict)[ind].error_mes," ");
    strcpy((*dict)[ind].keyword,"charm_conf_rho_def");
    strcpy((*dict)[ind].keyarg," ");
  //------------------------------------------------------------------------------
  //  2)~charm_conf_state_def[ ]
    ind = 2;
    strcpy((*dict)[ind].error_mes," ");
    strcpy((*dict)[ind].keyword,"charm_conf_state_def");
    strcpy((*dict)[ind].keyarg," ");
  //------------------------------------------------------------------------------
  //  3)~charm_conf_PC_def[ ]
    ind = 3;
    strcpy((*dict)[ind].error_mes," ");
    strcpy((*dict)[ind].keyword,"charm_conf_PC_def");
    strcpy((*dict)[ind].keyarg," ");
  //------------------------------------------------------------------------------
  //  4)~charm_conf_NL_def[ ]
    ind = 4;
    strcpy((*dict)[ind].error_mes," ");
    strcpy((*dict)[ind].keyword,"charm_conf_NL_def");
    strcpy((*dict)[ind].keyarg," ");
  //------------------------------------------------------------------------------
  //  5)~charm_conf_gen_def[ ]
    ind = 5;
    strcpy((*dict)[ind].error_mes," ");
    strcpy((*dict)[ind].keyword,"charm_conf_gen_def");
    strcpy((*dict)[ind].keyarg," ");
  //------------------------------------------------------------------------------
  //  6)~charm_conf_atm_def[ ]
    ind = 6;
    strcpy((*dict)[ind].error_mes," ");
    strcpy((*dict)[ind].keyword,"charm_conf_atm_def");
    strcpy((*dict)[ind].keyarg," ");
//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_rho  (int *num_dict ,DICT_WORD **dict){
//==================================================================================
//  I) Malloc the dictionary                                              

  num_dict[0] = 25;
  *dict = (DICT_WORD *)cmalloc(num_dict[0]*sizeof(DICT_WORD),"set_config_dict_rho")-1;

//=================================================================================
//  II) Initialize the user set option(did the user set the key word      

  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iuset    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iflag    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].key_type = 1;} //spec only once

//=================================================================================
// III) Set up the dictionary
  int ind;
  //-----------------------------------------------------------------------------
  //  1)\rhorHartpriority{}
    ind =   1;
    strcpy((*dict)[ind].keyword,"rhorHartpriority");
    strcpy((*dict)[ind].keyarg,"2000000"); 
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  //  2)\rhogHartpriority{}
    ind =   2;
    strcpy((*dict)[ind].keyword,"rhogHartpriority");
    strcpy((*dict)[ind].keyarg,"2000000");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  //  3)\useGHartInsRhoRP{}
    ind =   3;
    strcpy((*dict)[ind].keyword,"useGHartInsRhoRP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  //  4)\useGHartInsRHart{}
    ind =   4;
    strcpy((*dict)[ind].keyword,"useGHartInsRHart");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  //  5)\useRHartInsGHart{}
    ind =   5;
    strcpy((*dict)[ind].keyword,"useRHartInsGHart");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  //  6)\useCentroidMapRho{}
    ind =   6;
    strcpy((*dict)[ind].keyword,"useCentroidMapRho");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  //  7)\rhorpriority{}
    ind =   7;
    strcpy((*dict)[ind].keyword,"rhorpriority");
    strcpy((*dict)[ind].keyarg,"2000000");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  //  8)\rhogpriority{}
    ind =   8;
    strcpy((*dict)[ind].keyword,"rhogpriority");
    strcpy((*dict)[ind].keyarg,"2000000");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  //  9)\gExpandFactRho{}
    ind =   9;
    strcpy((*dict)[ind].keyword,"gExpandFactRho");
    strcpy((*dict)[ind].keyarg,"1.25");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 10)\lbdensity{}
    ind =  10;
    strcpy((*dict)[ind].keyword,"lbdensity");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 11)\rhoGHelpers{}
    ind =  11;
    strcpy((*dict)[ind].keyword,"rhoGHelpers");
    strcpy((*dict)[ind].keyarg,"1");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 12)\rhoRsubplanes{}
    ind =  12;
    strcpy((*dict)[ind].keyword,"rhoRsubplanes");
    strcpy((*dict)[ind].keyarg,"1");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 13)\useGIns0RhoRP{}
    ind =  13;
    strcpy((*dict)[ind].keyword,"useGIns0RhoRP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 14)\useGIns1RhoRP{}
    ind =  14;
    strcpy((*dict)[ind].keyword,"useGIns1RhoRP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 15)\useGIns2RhoRP{}
    ind =  15;
    strcpy((*dict)[ind].keyword,"useGIns2RhoRP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 16)\useGIns3RhoRP{}
    ind =  16;
    strcpy((*dict)[ind].keyword,"useGIns3RhoRP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 17)\useGByrdInsRhoRBP{}
    ind =  17;
    strcpy((*dict)[ind].keyword,"useGByrdInsRhoRBP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 18)\useRInsRhoGP{}
    ind =  18;
    strcpy((*dict)[ind].keyword,"useRInsRhoGP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 19)\useRInsIGXRhoGP{}
    ind =  19;
    strcpy((*dict)[ind].keyword,"useRInsIGXRhoGP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 20)\useRInsIGYRhoGP{}
    ind =  20;
    strcpy((*dict)[ind].keyword,"useRInsIGYRhoGP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 21)\useRInsIGZRhoGP{}
    ind =  21;
    strcpy((*dict)[ind].keyword,"useRInsIGZRhoGP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 22)\prioEextFFTMsg{}
    ind =  22;
    strcpy((*dict)[ind].keyword,"prioEextFFTMsg");
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 23)\rhoLineOrder{}
    ind =  23;
    strcpy((*dict)[ind].keyword,"rhoLineOrder");
    strcpy((*dict)[ind].keyarg,"skip");    
    strcpy((*dict)[ind].error_mes,"skip,none,random");
//----------------------------------------------------------------------------------
  // 24)\nchareHartAtmT{}
    ind =  24;
    strcpy((*dict)[ind].keyword,"nchareHartAtmT");
    strcpy((*dict)[ind].keyarg,"1");    
    strcpy((*dict)[ind].error_mes," number >=1 and <= natmtype");
  //-----------------------------------------------------------------------------
  // 25)\rhoSubPlaneBalance{}
    ind =  25;
    strcpy((*dict)[ind].keyword,"rhoSubPlaneBalance");
    strcpy((*dict)[ind].keyarg,"off");
    strcpy((*dict)[ind].error_mes,"on/off");
//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_params_rho (DICT_WORD *dict, char *fun_key, char *input_name,int iflag){
//===================================================================================

  int ind;
  int ierr;

//===================================================================================
// Set some values

  //-----------------------------------------------------------------------------
  //  1)\rhorHartpriority{}
    ind =   1;
    sscanf(dict[ind].keyarg,"%d",&rhorHartpriority);
    if(rhorHartpriority<=0){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  2)\rhogHartpriority{}
    ind =   2;
    sscanf(dict[ind].keyarg,"%d",&rhogHartpriority);
    if(rhogHartpriority<=0){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  3)\useGHartInsRhoRP{}
    ind =   3;
    parse_on_off(dict[ind].keyarg,&useGHartInsRhoRP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useGHartInsRhoRP=0;}
  //-----------------------------------------------------------------------------
  //  4)\useGHartInsRHart{}
    ind =   4;
    parse_on_off(dict[ind].keyarg,&useGHartInsRHart,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useGHartInsRHart=0;}
  //-----------------------------------------------------------------------------
  //  5)\useRHartInsGHart{}
    ind =   5;
    parse_on_off(dict[ind].keyarg,&useRHartInsGHart,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useRHartInsGHart=0;}
  //-----------------------------------------------------------------------------
  //  6)\useCentroidMapRho{}
    ind =   6;
    parse_on_off(dict[ind].keyarg,&useCentroidMapRho,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  7)\rhorpriority{}
    ind =   7;
    sscanf(dict[ind].keyarg,"%d",&rhorpriority);
    if(rhorpriority<=0){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  8)\rhogpriority{}
    ind =   8;
    sscanf(dict[ind].keyarg,"%d",&rhogpriority);
    if(rhogpriority<=0){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  9)\gExpandFactRho{}
    ind =   9;
    sscanf(dict[ind].keyarg,"%lg",&gExpandFactRho);
    if(gExpandFactRho<=0){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 10)\lbdensity{}
    ind =  10;
    parse_on_off(dict[ind].keyarg,&lbdensity,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 11)\rhoGHelpers{}
    ind =  11;
    sscanf(dict[ind].keyarg,"%d",&rhoGHelpers);
    if(rhoGHelpers<=0){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 12)\rhoRsubplanes{}
    ind =  12;
    sscanf(dict[ind].keyarg,"%d",&rhoRsubplanes);
    if(rhoRsubplanes<=0){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 13)\useGIns0RhoRP{}
    ind =  13;
    parse_on_off(dict[ind].keyarg,&useGIns0RhoRP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useGIns0RhoRP=0;}
  //-----------------------------------------------------------------------------
  // 14)\useGIns1RhoRP{}
    ind =  14;
    parse_on_off(dict[ind].keyarg,&useGIns1RhoRP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useGIns1RhoRP=0;}
  //-----------------------------------------------------------------------------
  // 15)\useGIns2RhoRP{}
    ind =  15;
    parse_on_off(dict[ind].keyarg,&useGIns2RhoRP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useGIns2RhoRP=0;}
  //-----------------------------------------------------------------------------
  // 16)\useGIns3RhoRP{}
    ind =  16;
    parse_on_off(dict[ind].keyarg,&useGIns3RhoRP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useGIns3RhoRP=0;}
  //-----------------------------------------------------------------------------
  // 17)\useGByrdInsRhoRBP{}
    ind =  17;
    parse_on_off(dict[ind].keyarg,&useGByrdInsRhoRBP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useGByrdInsRhoRBP=0;}
  //-----------------------------------------------------------------------------
  // 18)\useRInsRhoGP{}
    ind =  18;
    parse_on_off(dict[ind].keyarg,&useRInsRhoGP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useRInsRhoGP=0;}
  //-----------------------------------------------------------------------------
  // 19)\useRInsIGXRhoGP{}
    ind =  19;
    parse_on_off(dict[ind].keyarg,&useRInsIGXRhoGP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useRInsIGXRhoGP=0;}
  //-----------------------------------------------------------------------------
  // 20)\useRInsIGYRhoGP{}
    ind =  20;
    parse_on_off(dict[ind].keyarg,&useRInsIGYRhoGP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useRInsIGYRhoGP=0;}
  //-----------------------------------------------------------------------------
  // 21)\useRInsIGZRhoGP{}
    ind =  21;
    parse_on_off(dict[ind].keyarg,&useRInsIGZRhoGP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useRInsIGZRhoGP=0;}
  //-----------------------------------------------------------------------------
  // 22)\prioEextFFTMsg{}
    ind =  22;
    parse_on_off(dict[ind].keyarg,&prioEextFFTMsg,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 23)\rhoLineOrder{}
    ind =  23;
    ierr= 0;
    if(strcasecmp(dict[ind].keyarg,"none")==0)  {rhoLineOrder =-1; ierr++;}
    if(strcasecmp(dict[ind].keyarg,"skip")==0)  {rhoLineOrder = 0; ierr++;}
    if(strcasecmp(dict[ind].keyarg,"random")==0){rhoLineOrder = 1; ierr++;}
    if(ierr!=1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 24)\nchareHartAtmT{}
    ind =   24;
    sscanf(dict[ind].keyarg,"%d",&nchareHartAtmT);
    if(nchareHartAtmT<1 || nchareHartAtmT>natm_typ){
      keyarg_barf(dict,input_name,fun_key,ind);
    }//endif
  //-----------------------------------------------------------------------------
  // 25)\rhoSubPlaneBalance{}
    ind =  25;
    parse_on_off(dict[ind].keyarg,&rhoSubPlaneBalance,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}

//===================================================================================
// Clean up 

  if(rhoRsubplanes>1){
     PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
     PRINTF("   rhoRsubplanes=%d. Disabling rho rs commlib use\n",rhoRsubplanes);
     PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
     useGHartInsRHart  = 0; strcpy(dict[3].keyarg,"off");
     useRHartInsGHart  = 0; strcpy(dict[4].keyarg,"off");
     useGHartInsRhoRP  = 0; strcpy(dict[5].keyarg,"off");
     useGIns0RhoRP     = 0; strcpy(dict[13].keyarg,"off");
     useGIns1RhoRP     = 0; strcpy(dict[14].keyarg,"off");
     useGIns2RhoRP     = 0; strcpy(dict[15].keyarg,"off");
     useGIns3RhoRP     = 0; strcpy(dict[16].keyarg,"off");
     useGByrdInsRhoRBP = 0; strcpy(dict[17].keyarg,"off");
     useRInsRhoGP      = 0; strcpy(dict[18].keyarg,"off");
     useRInsIGXRhoGP   = 0; strcpy(dict[19].keyarg,"off");
     useRInsIGYRhoGP   = 0; strcpy(dict[20].keyarg,"off");
     useRInsIGZRhoGP   = 0; strcpy(dict[21].keyarg,"off");
  }//endif

//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_state(int *num_dict ,DICT_WORD **dict){
//==================================================================================
//  I) Malloc the dictionary                                              

  num_dict[0] = 25;
  *dict = (DICT_WORD *)cmalloc(num_dict[0]*sizeof(DICT_WORD),"set_config_dict_state")-1;

//=================================================================================
//  II) Initialize the user set option(did the user set the key word      

  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iuset    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iflag    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].key_type = 1;} //spec only once

//=================================================================================
// III) Set up the dictionary
  int ind;
  //-----------------------------------------------------------------------------
  //  1)\dataPath\{}
    ind=1;
    strcpy((*dict)[ind].keyword,"dataPath");
    strcpy((*dict)[ind].keyarg,"./");    
    strcpy((*dict)[ind].error_mes,"a directory tree");
  //-----------------------------------------------------------------------------
  //  2)\gBucketSize{}
    ind=2;
    strcpy((*dict)[ind].keyword,"gBucketSize");
    strcpy((*dict)[ind].keyarg,"5");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  //  3)\rBucketSize{}
    ind=3;
    strcpy((*dict)[ind].keyword,"rBucketSize");
    strcpy((*dict)[ind].keyarg,"5");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  //  4)\gStreamPeriod{}
    ind=4;
    strcpy((*dict)[ind].keyword,"gStreamPeriod");
    strcpy((*dict)[ind].keyarg,"2");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  //  5)\rStreamPeriod{}
    ind=5;
    strcpy((*dict)[ind].keyword,"rStreamPeriod");
    strcpy((*dict)[ind].keyarg,"2");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  //  6)\gExpandFact{}
    ind=6;
    strcpy((*dict)[ind].keyword,"gExpandFact");
    strcpy((*dict)[ind].keyarg,"1.25");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  //  7)\Gstates_per_pe{}
    ind=7;
    strcpy((*dict)[ind].keyword,"Gstates_per_pe");
    sprintf((*dict)[ind].keyarg,"%d",nstates);  
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  //  8)\Rstates_per_pe{}
    ind=8;
    strcpy((*dict)[ind].keyword,"Rstates_per_pe");
    sprintf((*dict)[ind].keyarg,"%d",nstates);
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  //  9)\stateOutputOn{}
    ind=9;
    strcpy((*dict)[ind].keyword,"stateOutputOn");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"off/on");
  //-----------------------------------------------------------------------------
  // 10)\psipriority{}
    ind=10;
    strcpy((*dict)[ind].keyword,"psipriority");
    strcpy((*dict)[ind].keyarg,"400000000");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 11)\prioFFTMsg{}
    ind=11;
    strcpy((*dict)[ind].keyword,"prioFFTMsg");
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 12)\rsfftpriority{}
    ind=12;
    strcpy((*dict)[ind].keyword,"rsfftpriority");
    strcpy((*dict)[ind].keyarg,"2000000");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 13)\gsfftpriority{}
    ind=13;
    strcpy((*dict)[ind].keyword,"gsfftpriority");
    strcpy((*dict)[ind].keyarg,"1000000");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 14)\rsifftpriority{}
    ind=14;
    strcpy((*dict)[ind].keyword,"rsifftpriority");
    strcpy((*dict)[ind].keyarg,"100000000");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 15)\gsifftpriority{}
    ind=15;
    strcpy((*dict)[ind].keyword,"gsifftpriority");
    strcpy((*dict)[ind].keyarg,"200000000");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 16)\conserveMemory{}
    ind=16;
    strcpy((*dict)[ind].keyword,"conserveMemory");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 17)\lbgspace{}
    ind=17;
    strcpy((*dict)[ind].keyword,"lbgspace");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 18)\doublePack{}
    ind=18;
    strcpy((*dict)[ind].keyword,"doublePack");
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 19)\useCuboidMap{}
    ind=19;
    strcpy((*dict)[ind].keyword,"useCuboidMap");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 20)\useCuboidMapRS{}
    ind=20;
    strcpy((*dict)[ind].keyword,"useCuboidMapRS");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 21)\useCentroidMap{}
    ind=21;
    strcpy((*dict)[ind].keyword,"useCentroidMap");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 22)\useGssInsRealP{}
    ind=22;
    strcpy((*dict)[ind].keyword,"useGssInsRealP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 23)\useMssInsGP{}
    ind=23;
    strcpy((*dict)[ind].keyword,"useMssInsGP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 24)\loadMapFiles{}
    ind=24;
    strcpy((*dict)[ind].keyword,"loadMapFiles");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 25)\dumpMapFiles{}
    ind=25;
    strcpy((*dict)[ind].keyword,"dumpMapFiles");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");

//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_params_state(DICT_WORD *dict, char *fun_key, char *input_name,int iflag){
//===================================================================================

  int ind;
  int ierr;

//===================================================================================
// Fill the class with data

  //-----------------------------------------------------------------------------
  //  1)\dataPath\{}
    ind=1;
    strcpy(dataPath, dict[ind].keyarg);
  //-----------------------------------------------------------------------------
  //  2)\gBucketSize{}
    ind=2;
    sscanf(dict[ind].keyarg,"%d",&gBucketSize);
    if(gBucketSize<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  3)\rBucketSize{}
    ind=3;
    sscanf(dict[ind].keyarg,"%d",&rBucketSize);
    if(rBucketSize<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  4)\gStreamPeriod{}
    ind=4;
    sscanf(dict[ind].keyarg,"%lg",&gStreamPeriod);
    if(gStreamPeriod<=0.0){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  5)\rStreamPeriod{}
    ind=5;
    sscanf(dict[ind].keyarg,"%lg",&rStreamPeriod);
    if(rStreamPeriod<=0.0){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  6)\gExpandFact{}
    ind=6;
    sscanf(dict[ind].keyarg,"%lg",&gExpandFact);
    if(gExpandFact<=0.0){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  7)\Gstates_per_pe{}
    ind=7;
    sscanf(dict[ind].keyarg,"%d",&Gstates_per_pe);
    if(Gstates_per_pe<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  8)\Rstates_per_pe{}
    ind=8;
    sscanf(dict[ind].keyarg,"%d",&Rstates_per_pe);
    if(Rstates_per_pe<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  9)\stateOutputOn{}
    ind=9;
    parse_on_off(dict[ind].keyarg,&stateOutputOn,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 10)\psipriority{}
    ind=10;
    sscanf(dict[ind].keyarg,"%d",&psipriority);
    if(psipriority<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 11)\prioFFTMsg{}
    ind=11;
    parse_on_off(dict[ind].keyarg,&prioFFTMsg,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 12)\rsfftpriority{}
    ind=12;
    sscanf(dict[ind].keyarg,"%d",&rsfftpriority);
    if(rsfftpriority<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 13)\gsfftpriority{}
    ind=13;
    sscanf(dict[ind].keyarg,"%d",&gsfftpriority);
    if(gsfftpriority<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 14)\rsifftpriority{}
    ind=14;
    sscanf(dict[ind].keyarg,"%d",&rsifftpriority);
    if(rsifftpriority<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 15)\gsifftpriority{}
    ind=15;
    sscanf(dict[ind].keyarg,"%d",&gsifftpriority);
    if(gsifftpriority<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 16)\conserveMemory{}
    ind=16;
    parse_on_off(dict[ind].keyarg,&conserveMemory,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 17)\lbgspace{}
    ind=17;
    parse_on_off(dict[ind].keyarg,&lbgspace,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 18)\doublePack{}
    ind=18;
    parse_on_off(dict[ind].keyarg,&doublePack,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 19)\useCuboidMap{}
    ind=19;
    parse_on_off(dict[ind].keyarg,&useCuboidMap,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 20)\useCuboidMapRS{}
    ind=20;
    parse_on_off(dict[ind].keyarg,&useCuboidMapRS,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 21)\useCentroidMap{}
    ind=21;
    parse_on_off(dict[ind].keyarg,&useCentroidMap,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 22)\useGssInsRealP{}
    ind=22;
    parse_on_off(dict[ind].keyarg,&useGssInsRealP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useGssInsRealP=0;}
  //-----------------------------------------------------------------------------
  // 23)\useMssInsGP{}
    ind=23;
    parse_on_off(dict[ind].keyarg,&useMssInsGP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useMssInsGP=0;}
  //-----------------------------------------------------------------------------
  // 24)\loadMapFiles{}
    ind=24;
    parse_on_off(dict[ind].keyarg,&loadMapFiles,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 25)\dumpMapFiles{}
    ind=25;
    parse_on_off(dict[ind].keyarg,&dumpMapFiles,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}

//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_pc (int *num_dict ,DICT_WORD **dict){
//==================================================================================
//  I) Malloc the dictionary                                              

  num_dict[0] = 31;
  *dict = (DICT_WORD *)cmalloc(num_dict[0]*sizeof(DICT_WORD),"set_config_dict_pc")-1;

//=================================================================================
//  II) Initialize the user set option(did the user set the key word      

  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iuset    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iflag    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].key_type = 1;} //spec only once

//=================================================================================
// III) Set up the dictionary
  int ind;
  //-----------------------------------------------------------------------------
  // 1)\usePairDirectSend{}
    ind=1;
    strcpy((*dict)[ind].keyword,"usePairDirectSend");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 2)\PCCollectTiles{}
    ind=2;
    strcpy((*dict)[ind].keyword,"PCCollectTiles");
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 3)\PCdelayBWSend{}
    ind=3;
    strcpy((*dict)[ind].keyword,"PCdelayBWSend");
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 4)\PCstreamBWout{}
    ind=4;
    strcpy((*dict)[ind].keyword,"PCstreamBWout");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 5)\PCstreamFWblock{}
    ind=5;
    strcpy((*dict)[ind].keyword,"PCstreamFWblock");
    sprintf((*dict)[ind].keyarg,"%d",0);
    strcpy((*dict)[ind].error_mes,"a number >= 0");
  //-----------------------------------------------------------------------------
  // 6)\useOrthoDirect{}
    ind=6;
    strcpy((*dict)[ind].keyword,"useOrthoDirect");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 7)\useOrthoHelpers{}
    ind=7;
    strcpy((*dict)[ind].keyword,"useOrthoHelpers");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 8)\useOrthoSection{}
    ind=8;
    strcpy((*dict)[ind].keyword,"useOrthoSection");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 9)\useOrthoSectionRed{}
    ind=9;
    strcpy((*dict)[ind].keyword,"useOrthoSectionRed");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 10)\lambdaGrainSize{}
    ind=10;
    strcpy((*dict)[ind].keyword,"lambdaGrainSize");
    sprintf((*dict)[ind].keyarg,"%d",nstates);
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 11)\PCSpanFactor{}
    ind=11;
    strcpy((*dict)[ind].keyword,"PCSpanFactor");
    strcpy((*dict)[ind].keyarg,"2");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 12)\OrthoRedSpanFactor{}
    ind=12;
    strcpy((*dict)[ind].keyword,"OrthoRedSpanFactor");
    strcpy((*dict)[ind].keyarg,"2");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 13)\OrthoMcastSpanFactor{}
    ind=13;
    strcpy((*dict)[ind].keyword,"OrthoMcastSpanFactor");
    strcpy((*dict)[ind].keyarg,"16");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 14)\sGrainSize{}
    ind=14;
    strcpy((*dict)[ind].keyword,"sGrainSize");
    sprintf((*dict)[ind].keyarg,"%d",nstates);
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 15)\gemmSplitFWk{}
    ind=15;
    strcpy((*dict)[ind].keyword,"gemmSplitFWk");
    strcpy((*dict)[ind].keyarg,"16");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 16)\gemmSplitFWm{}
    ind=16;
    strcpy((*dict)[ind].keyword,"gemmSplitFWm");
    strcpy((*dict)[ind].keyarg,"16");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 17)\gemmSplitBW{}
    ind=17;
    strcpy((*dict)[ind].keyword,"gemmSplitBW");
    strcpy((*dict)[ind].keyarg,"16");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 18)\gemmSplitOrtho{}
    ind=18;
    strcpy((*dict)[ind].keyword,"gemmSplitOrtho");
    strcpy((*dict)[ind].keyarg,"32");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 19)\orthoGrainSize{}
    ind=19;
    strcpy((*dict)[ind].keyword,"orthoGrainSize");
    sprintf((*dict)[ind].keyarg,"%d",nstates);
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 20)\orthoStride{}
    ind=20;
    strcpy((*dict)[ind].keyword,"orthoStride");
    strcpy((*dict)[ind].keyarg,"0");    
    strcpy((*dict)[ind].error_mes,"a number >= 0");
  //-----------------------------------------------------------------------------
  // 21)\useBWBarrier{}
    ind=21;
    strcpy((*dict)[ind].keyword,"useBWBarrier");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 22)\phantomSym{}
    ind=22;
    strcpy((*dict)[ind].keyword,"phantomSym");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 23)\lbpaircalc{}
    ind=23;
    strcpy((*dict)[ind].keyword,"lbpaircalc");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 24)\lambdapriority{}
    ind=24;
    strcpy((*dict)[ind].keyword,"lambdapriority");
    strcpy((*dict)[ind].keyarg,"300000000");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 25)\toleranceInterval{}
    ind=25;
    strcpy((*dict)[ind].keyword,"toleranceInterval");
    strcpy((*dict)[ind].keyarg,"1");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 26)\gSpaceSum{}
    ind=26;
    strcpy((*dict)[ind].keyword,"gSpaceSum");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 27)\numChunks{}
    ind=27;
    strcpy((*dict)[ind].keyword,"numChunks");
    strcpy((*dict)[ind].keyarg,"1");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 28)\numChunksSym{}
    ind=28;
    strcpy((*dict)[ind].keyword,"numChunksSym");
    strcpy((*dict)[ind].keyarg,"1");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 29)\numChunksAsym{}
    ind=29;
    strcpy((*dict)[ind].keyword,"numChunksAsym");
    strcpy((*dict)[ind].keyarg,"1");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 30)\prioBW{}
    ind=30;
    strcpy((*dict)[ind].keyword,"prioBW");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 31)\usePairEtoM{}
    ind=31;
    strcpy((*dict)[ind].keyword,"usePairEtoM");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");

//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_params_pc  (DICT_WORD *dict, char *fun_key, char *input_name){
//===================================================================================

  int ind,ierr;

//===================================================================================
// Fill me with joy.

  //-----------------------------------------------------------------------------
  // 1)\usePairDirectSend{}
    ind=1;
    parse_on_off(dict[ind].keyarg,&usePairDirectSend,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 2)\PCCollectTiles{}
    ind=2;
    parse_on_off(dict[ind].keyarg,&PCCollectTiles,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 3)\PCdelayBWSend{}
    ind=3;
    parse_on_off(dict[ind].keyarg,&PCdelayBWSend,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 4)\PCstreamBWout{}
    ind=4;
    parse_on_off(dict[ind].keyarg,&PCstreamBWout,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 5)\PCstreamFWblock{}
    ind=5;
    sscanf(dict[ind].keyarg,"%d",&PCstreamFWblock);
    if(PCstreamFWblock<0){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 6)\useOrthoDirect{}
    ind=6;
    parse_on_off(dict[ind].keyarg,&useOrthoDirect,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 7)\useOrthoHelpers{}
    ind=7;
    parse_on_off(dict[ind].keyarg,&useOrthoHelpers,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 8)\useOrthoSection{}
    ind=8;
    parse_on_off(dict[ind].keyarg,&useOrthoSection,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 9)\useOrthoSectionRed{}
    ind=9;
    parse_on_off(dict[ind].keyarg,&useOrthoSectionRed,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 10)\lambdaGrainSize{}
    ind=10;
    sscanf(dict[ind].keyarg,"%d",&lambdaGrainSize);
    if(lambdaGrainSize<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 11)\PCSpanFactor{}
    ind=11;
    sscanf(dict[ind].keyarg,"%d",&PCSpanFactor);
    if(PCSpanFactor<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 12)\OrthoRedSpanFactor{}
    ind=12;
    sscanf(dict[ind].keyarg,"%d",&OrthoRedSpanFactor);
    if(OrthoRedSpanFactor<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 13)\OrthoMcastSpanFactor{}
    ind=13;
    sscanf(dict[ind].keyarg,"%d",&OrthoMcastSpanFactor);
    if(OrthoMcastSpanFactor<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 14)\sGrainSize{}
    ind=14;
    sscanf(dict[ind].keyarg,"%d",&sGrainSize);
    if(sGrainSize<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 15)\gemmSplitFWk{}
    ind=15;
    sscanf(dict[ind].keyarg,"%d",&gemmSplitFWk);
    if(gemmSplitFWk<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 16)\gemmSplitFWm{}
    ind=16;
    sscanf(dict[ind].keyarg,"%d",&gemmSplitFWm);
    if(gemmSplitFWm<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 17)\gemmSplitBW{}
    ind=17;
    sscanf(dict[ind].keyarg,"%d",&gemmSplitBW);
    if(gemmSplitBW<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 18)\gemmSplitOrtho{}
    ind=18;
    sscanf(dict[ind].keyarg,"%d",&gemmSplitOrtho);
    if(gemmSplitOrtho<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 19)\orthoGrainSize{}
    ind=19;
    sscanf(dict[ind].keyarg,"%d",&orthoGrainSize);
    if(orthoGrainSize<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 20)\orthoStride{}
    ind=20;
    sscanf(dict[ind].keyarg,"%d",&orthoStride);
    if(orthoStride<0){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 21)\useBWBarrier{}
    ind=21;
    parse_on_off(dict[ind].keyarg,&useBWBarrier,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 22)\phantomSym{}
    ind=22;
    parse_on_off(dict[ind].keyarg,&phantomSym,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 23)\lbpaircalc{}
    ind=23;
    parse_on_off(dict[ind].keyarg,&lbpaircalc,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 24)\lambdapriority{}
    ind=24;
    sscanf(dict[ind].keyarg,"%d",&lambdapriority);
    if(lambdapriority<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 25)\toleranceInterval{}
    ind=25;
    sscanf(dict[ind].keyarg,"%d",&toleranceInterval);
    if(toleranceInterval<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 26)\gSpaceSum{}
    ind=26;
    parse_on_off(dict[ind].keyarg,&gSpaceSum,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 27)\numChunks{}
    ind=27;
    sscanf(dict[ind].keyarg,"%d",&numChunks);
    if(numChunks<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 28)\numChunksSym{}
    ind=28;
    sscanf(dict[ind].keyarg,"%d",&numChunksSym);
    if(numChunksSym<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 29)\numChunksAsym{}
    ind=29;
    sscanf(dict[ind].keyarg,"%d",&numChunksAsym);
    if(numChunksAsym<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 30)\prioBW{}
    ind=30;
    parse_on_off(dict[ind].keyarg,&prioBW,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 31)\usePairEtoM{}
    ind=31;
    parse_on_off(dict[ind].keyarg,&usePairEtoM,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}

//===================================================================================
// Clean up the user input

  if(numChunks>numChunksSym) {
     PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
     PRINTF("   Cleaning up numchunks > numchunksSym : %d %d\n",numChunks,numChunksSym);
     PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
     numChunksSym=numChunks;
     sprintf(dict[28].keyarg,"%d",numChunksSym);
  }//endif

  if(numChunks>numChunksAsym){
     PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
     PRINTF("   Cleaning up numchunks > numchunksAsym : %d %d\n",numChunks,numChunksAsym);
     PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
     numChunksAsym=numChunks;
     sprintf(dict[29].keyarg,"%d",numChunksAsym);
  }//endif

  if(lambdaGrainSize==nstates && orthoGrainSize!=nstates){
     PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
     PRINTF("   Cleaning up lambdaGrainSize==nstates && orthoGrainSize!=nstates\n");
     PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
     lambdaGrainSize=orthoGrainSize;
     sprintf(dict[10].keyarg,"%d",lambdaGrainSize);
  }//endif

//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_nl (int *num_dict ,DICT_WORD **dict){
//==================================================================================
//  I) Malloc the dictionary                                              

  num_dict[0] = 9;
  *dict = (DICT_WORD *)cmalloc(num_dict[0]*sizeof(DICT_WORD),"set_config_dict_nl")-1;

//=================================================================================
//  II) Initialize the user set option(did the user set the key word      

  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iuset    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iflag    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].key_type = 1;} //spec only once

//=================================================================================
// III) Set up the dictionary
  int ind;
  //-----------------------------------------------------------------------------
  // 1)\sfpriority{}
    ind=1;
    strcpy((*dict)[ind].keyword,"sfpriority");
    strcpy((*dict)[ind].keyarg,"10000000");
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 2)\prioNLFFTMsg{}
    ind=2;
    strcpy((*dict)[ind].keyword,"prioNLFFTMsg");
    strcpy((*dict)[ind].keyarg,"on");
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 3)\rsNLfftpriority{}
    ind=3;
    strcpy((*dict)[ind].keyword,"rsNLfftpriority");
    strcpy((*dict)[ind].keyarg,"2300000");
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 4)\gsNLfftpriority{}
    ind=4;
    strcpy((*dict)[ind].keyword,"gsNLfftpriority");
    strcpy((*dict)[ind].keyarg,"2500000");
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 5)\numSfGrps{}
    ind=5;
    strcpy((*dict)[ind].keyword,"numSfGrps");
    strcpy((*dict)[ind].keyarg,"1");
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 6)\numSfDups{}
    ind=6;
    strcpy((*dict)[ind].keyword,"numSfDups");
    strcpy((*dict)[ind].keyarg,"1");
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 7)\launchNLeesFromRho{}
    ind=7;
    int ierr= 0;
    strcpy((*dict)[ind].keyword,"launchNLeesFromRho");
    strcpy((*dict)[ind].keyarg,"rs");    
    strcpy((*dict)[ind].error_mes,"rs rhor rhog");
  //-----------------------------------------------------------------------------
  // 8)\useGssInsRealPP{}
    ind=8;
    strcpy((*dict)[ind].keyword,"useGssInsRealPP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 9)\useMssInsGPP{}
    ind=9;
    strcpy((*dict)[ind].keyword,"useMssInsGPP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");

//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_params_nl (DICT_WORD *dict, char *fun_key, char *input_name,int iflag){
//===================================================================================

  int ind, ierr;

//===================================================================================
// Fill 'er up
    
  //-----------------------------------------------------------------------------
  // 1)\sfpriority{}
    ind=1;
    sscanf(dict[ind].keyarg,"%d",&sfpriority);
    if(sfpriority<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 2)\prioNLFFTMsg{}
    ind=2;
    parse_on_off(dict[ind].keyarg,&prioNLFFTMsg,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 3)\rsNLfftpriority{}
    ind=3;
    sscanf(dict[ind].keyarg,"%d",&rsNLfftpriority);
    if(rsNLfftpriority<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 4)\gsNLfftpriority{}
    ind=4;
    sscanf(dict[ind].keyarg,"%d",&gsNLfftpriority);
    if(gsNLfftpriority<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 5)\numSfGrps{}
    ind=5;
    sscanf(dict[ind].keyarg,"%d",&numSfGrps);
    if(numSfGrps<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 6)\numSfDups{}
    ind=6;
    sscanf(dict[ind].keyarg,"%d",&numSfDups);
    if(numSfDups<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 7)\launchNLeesFromRho{}
    ind=7;
    if(strcasecmp(dict[ind].keyarg,"rs")==0)  {launchNLeesFromRho =0; ierr++;}
    if(strcasecmp(dict[ind].keyarg,"rhor")==0)  {launchNLeesFromRho = 1; ierr++;}
    if(strcasecmp(dict[ind].keyarg,"rhog")==0){launchNLeesFromRho = 2; ierr++;}
    if(ierr!=1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(launchNLeesFromRho<0){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 8)\useGssInsRealPP{}
    ind=8;
    parse_on_off(dict[ind].keyarg,&useGssInsRealPP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useGssInsRealPP=0;}
  //-----------------------------------------------------------------------------
  // 9)\useMssInsGPP{}
    ind=9;
    parse_on_off(dict[ind].keyarg,&useMssInsGPP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useMssInsGPP=0;}

//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_gen (int *num_dict ,DICT_WORD **dict){
//==================================================================================
//  I) Malloc the dictionary                                              

  num_dict[0] = 7;
  *dict = (DICT_WORD *)cmalloc(num_dict[0]*sizeof(DICT_WORD),"set_config_dict_gen")-1;

//=================================================================================
//  II) Initialize the user set option(did the user set the key word      

  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iuset    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iflag    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].key_type = 1;} //spec only once

//=================================================================================
// III) Set up the dictionary
  int ind;
  //-----------------------------------------------------------------------------
  // 1)\fftprogresssplit{}
    ind=1;
    strcpy((*dict)[ind].keyword,"fftprogresssplit");
#ifdef CMK_VERSION_BLUEGENE
    strcpy((*dict)[ind].keyarg,"20");    
#else
    strcpy((*dict)[ind].keyarg,"1000");    
#endif  
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 2)\fftprogresssplitReal{}
    ind=2;
    strcpy((*dict)[ind].keyword,"fftprogresssplitReal");
#ifdef CMK_VERSION_BLUEGENE
    strcpy((*dict)[ind].keyarg,"10");    
#else
    strcpy((*dict)[ind].keyarg,"1000");    
#endif  
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 3)\useCommlib{}
    ind=3;
    strcpy((*dict)[ind].keyword,"useCommlib");
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 4)\useGMulticast{}
    ind=4;
    strcpy((*dict)[ind].keyword,"useGMulticast");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 5)\numMulticastMsgs{}
    ind=5;
    strcpy((*dict)[ind].keyword,"numMulticastMsgs");
    strcpy((*dict)[ind].keyarg,"10");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 6)\useCommlibMulticast{}
    ind=6;
    strcpy((*dict)[ind].keyword,"useCommlibMulticast");
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"on/off");
//----------------------------------------------------------------------------------
  // 7)\useTimeKeeper{}
    ind=7;
    strcpy((*dict)[ind].keyword,"useTimeKeeper");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------

  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_params_gen (DICT_WORD *dict, char *fun_key, char *input_name){
//===================================================================================

  int ind,ierr;

//===================================================================================
// Fill 

  //-----------------------------------------------------------------------------
  // 1)\fftprogresssplit{}
    ind=1;
    sscanf(dict[ind].keyarg,"%d",&fftprogresssplit);
    if(fftprogresssplit<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 2)\fftprogresssplitReal{}
    ind=2;
    sscanf(dict[ind].keyarg,"%d",&fftprogresssplitReal);
    if(fftprogresssplitReal<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 3)\useCommlib{}
    ind=3;
    parse_on_off(dict[ind].keyarg,&useCommlib,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 4)\useGMulticast{}
    ind=4;
    parse_on_off(dict[ind].keyarg,&useGMulticast,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 5)\numMulticastMsgs{}
    ind=5;
    sscanf(dict[ind].keyarg,"%d",&numMulticastMsgs);
    if(numMulticastMsgs<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 6)\useCommlibMulticast{}
    ind=6;
    parse_on_off(dict[ind].keyarg,&useCommlibMulticast,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}

//----------------------------------------------------------------------------------
  // 7)\useTimeKeeper{}
    ind=7;
    parse_on_off(dict[ind].keyarg,&useTimeKeeper,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}

//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_atm (int *num_dict ,DICT_WORD **dict){
//==================================================================================
//  I) Malloc the dictionary                                              

  num_dict[0] = 3;
  *dict = (DICT_WORD *)cmalloc(num_dict[0]*sizeof(DICT_WORD),"set_config_dict_atm")-1;

//=================================================================================
//  II) Initialize the user set option(did the user set the key word      

  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iuset    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iflag    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].key_type = 1;} //spec only once

//=================================================================================
// III) Set up the dictionary
  int ind;
  //-----------------------------------------------------------------------------
  // 1)\localAtomBarrier{}
    ind=1;
    strcpy((*dict)[ind].keyword,"localAtomBarrier");
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 2)\localEnergyBarrier{}
    ind=2;
    strcpy((*dict)[ind].keyword,"localEnergyBarrier"); 
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 3)\atmOutputOn{}
    ind=3;
    strcpy((*dict)[ind].keyword,"atmOutputOn");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_params_atm (DICT_WORD *dict, char *fun_key, char *input_name){
//===================================================================================

  int ind, ierr; 

//===================================================================================
// Fill the class with data 

  //-----------------------------------------------------------------------------
  // 1)\localAtomBarrier{}
    ind=1;
    parse_on_off(dict[ind].keyarg,&localAtomBarrier,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 2)\localEnergyBarrier{}
    ind=2;
    parse_on_off(dict[ind].keyarg,&localEnergyBarrier,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 3)\atmOutputOn{}
    ind=3;
    parse_on_off(dict[ind].keyarg,&atmOutputOn,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}

//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================



//=============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=============================================================================
/**
 * Guesses decent values for configuration parameters based on user
 * set values, system size, and the number of processes available.
 *
 * For large pe runs, scheme is thus: 
 *    1. find power2 nchareG. (typically 32 or 64)
 *    2. determine numchunks and grainsize
 *      a. find numPes()/nchareG.
 *      b. set sGrainsize to nstates, set numchunks to 1
 *      c. set numgrains to (nstates/sGrainSize)^2
 *      
 *    3. if phantoms are in use, numchunks is same for asym and sym
 *    4. set rstates_per_pe and gstates_per_pe as high as they will
 *       go for numPes and this system (constrained by nstates,
 *       nplanes, nchareg
 *    5. if numPes is large enable gSpaceSum, cuboid mapping,
 *       centroid mapping, useDirectSend and probably a few other things
 *
 * Supported parameters: sGrainSize, gExpandFact, gExpandFactRho,
 * fftprogresssplit, fftprogresssplitReal, numSfGrps, numSfDups,
 * rhoGHelpers, numMulticastMsgs

 */
//=============================================================================
void Config::guesstimateParmsConfig(int sizez,DICT_WORD *dict_gen,DICT_WORD *dict_rho,
                            DICT_WORD *dict_state,DICT_WORD *dict_pc,
                            DICT_WORD *dict_nl,DICT_WORD *dict_atm){
//=============================================================================

    int sqrtpes    = (int) sqrt((double)numPes);
    int sqrtstates = (int) sqrt((double)nstates);
    int igo;

//=============================================================================
// gExpandFact not set 

    igo = dict_state[6].iuset;

    if(igo==0){ 
      if(numPes>low_x_size){
         int i=1;
         double mypow=1;
         while((mypow=pow(2.0, (double)i)) <= low_x_size){i++;}
         gExpandFact = mypow / (double) low_x_size;
         nchareG     = (int)( gExpandFact * (double) low_x_size);
         sprintf(dict_state[6].keyarg,"%g",gExpandFact);
      }//endif
    }//endif : gexpandfact not set

    nchareG  = (int)(gExpandFact*(double)low_x_size);

//=============================================================================
// numChunks, numChunksSym, numChunksAsym, sgrainsize are not set

    igo = (dict_pc[27].iuset+dict_pc[28].iuset+dict_pc[29].iuset+dict_pc[14].iuset);

    if(igo==0){
       int numGrains = nstates/sGrainSize;
       numGrains    *=numGrains;
       numChunks     = numPes/(nchareG*numGrains);
       if(numChunks<1){numChunks=1;}
       while(numChunks>16){
          sGrainSize = sGrainSize/2;
          numGrains  = nstates/sGrainSize;
          numGrains *= numGrains;
          numChunks  = numPes/(nchareG*numGrains);
       }//end while
       numChunksSym  = numChunks;
       numChunksAsym = numChunks;
       sprintf(dict_pc[14].keyarg,"%d",sGrainSize);
       sprintf(dict_pc[27].keyarg,"%d",numChunks);
       sprintf(dict_pc[28].keyarg,"%d",numChunksSym);
       sprintf(dict_pc[29].keyarg,"%d",numChunksAsym);
    }//endif : chunks not set

//=============================================================================
// gstates_per_pe not set
    
    igo = dict_state[7].iuset;

    if(numPes!=1 && igo==0){
      if(numPes<=128){
        Gstates_per_pe=nstates/4;
      }else{ 
	if(numPes>128 && numPes<=512){
          Gstates_per_pe=nstates/16;
	}else{
          Gstates_per_pe= nchareG*nstates/numPes;
	}//endif
      }//endif
      if(Gstates_per_pe==0){Gstates_per_pe=1;}
      sprintf(dict_state[7].keyarg,"%d",Gstates_per_pe);
    }//endif : Gstates_per_pe not set

//=============================================================================
// rstates_per_pe not set

    igo = dict_state[8].iuset;

    if(numPes!=1 && igo==0){
      if(numPes<=128){
        Rstates_per_pe=nstates/4;
      }else{
        Rstates_per_pe= sizez*nstates/numPes;
      }//endif
      if (Rstates_per_pe==0){
        Rstates_per_pe=1;
      }//endif
      sprintf(dict_state[8].keyarg,"%d",Rstates_per_pe);
    }//endif : Rstates_per_pe not set

//=============================================================================
// Rejiggering orthograinsize if dumb values employed :
//         This is lame and should be replace with something which finds
//         an even mod of any sGrainSize

    igo = dict_pc[19].iuset;

    if( (sGrainSize%orthoGrainSize !=0) || (sGrainSize==orthoGrainSize)){
      int iii = orthoGrainSize;
      if(orthoGrainSize>32){orthoGrainSize=32;}
      if(orthoGrainSize>sGrainSize){orthoGrainSize=sGrainSize;}
      if(igo==1){
         PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
         PRINTF("   Changing your choice of orthograinsize from %d\n",iii);
         PRINTF("   to %d\n",orthoGrainSize);
         PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      }//endif
      sprintf(dict_pc[19].keyarg,"%d",orthoGrainSize);
    }//endif

//=============================================================================
// Rejiggering lambdagrainsize if dumb values employed :
//            This is lame and should be replace with something which finds
//            an even mod of any (non prime) sGrainSize

    igo = dict_pc[10].iuset;

    if((sGrainSize%lambdaGrainSize !=0)|| (sGrainSize==lambdaGrainSize)){
      if(igo==1){
         PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
         PRINTF("   (sGrainSize%lambdaGrainSize !=0)|| (sGrainSize==lambdaGrainSize)");
         PRINTF("   lambdaGrainSize=orthoGrainSize = %d\n",orthoGrainSize);
         PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      }//endif
      lambdaGrainSize=orthoGrainSize;
      sprintf(dict_pc[10].keyarg,"%d",lambdaGrainSize);
    }//endif


//=============================================================================
// Gexpandfact rho

    igo = dict_rho[9].iuset;

    if(igo==0){ 
       if(numPes>low_x_size*4){
         gExpandFactRho   = 1.0+fabs((double) (numPes/2-low_x_size*4)/ (double)( numPes));
       }else{
         if(numPes>low_x_size){
           gExpandFactRho = 1.0+(double) sqrtpes/ (double)( numPes);
         }//endif
       }//endif
       sprintf(dict_rho[9].keyarg,"%g",gExpandFactRho);
    }//endif

//=============================================================================
// numsfgrps
    // number of groups to chop atom calc into
    // needs to grow with size and number of PEs
    // range 1->natm_nl 

    igo = dict_nl[5].iuset;

    if(igo==0){
       int atmstates=natm_nl*nstates;
       if(numPes<atmstates){
           double ratio=(double)numPes/(double)atmstates;
           numSfGrps=(int) (ratio*(double)natm_nl);
       }else{ 
           numSfGrps=natm_nl-1;
       }//endif
       if(numSfGrps==0){numSfGrps=1;}
       sprintf(dict_nl[5].keyarg,"%d",numSfGrps);
    }//endif

//=============================================================================
// numsfdups
       // numbers of duplicate caches to create needs to grow with
       // number of atom chunks mapped to multiple PEs 
       // which is map dependant... yikes. can't set this here.
       // Real range is 1 -> nstates

    igo = dict_nl[6].iuset;

    if(numSfDups==1){
       int atmstates=natm_nl*nstates;
       if(numPes<atmstates){
           double ratio=(double)numPes/(double)atmstates;
           numSfDups=(int) (ratio*(double)nstates);
       }else{ 
           numSfDups=nstates-1;
       }//endif
       if(numSfDups==0){numSfDups=1;}
       sprintf(dict_nl[6].keyarg,"%d",numSfDups);
    }//endif

//=============================================================================
// Fix rhoGHelpers : ranges from 1 to numPes/numChareRhoG

    igo = dict_rho[11].iuset;

    if(igo==0){
      int temp_rho  = (int) (gExpandFactRho*2.0*((double) low_x_size + 1.0));
      if(numPes>temp_rho){rhoGHelpers=numPes/temp_rho;}
      sprintf(dict_rho[11].keyarg,"%d",rhoGHelpers);
    }//endif

//============================================================================
   }//end routine
//============================================================================



//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::readStateInfo(int &nPacked,int &minx, int &maxx, int &nx, int &ny, int &nz,
                            const char *fromFile, int ibinary_opt) {
//===================================================================================
// Check for errors

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
    strcpy(localFile,fromFile);
    strcat(localFile,".gz");
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
    strcpy(localFile,fromFile);
    strcat(localFile,".gz");
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
      fread(&(nPacked),sizeof(int),n,fp);
      fread(&(nx),sizeof(int),n,fp);
      fread(&(ny),sizeof(int),n,fp);
      fread(&(nz),sizeof(int),n,fp);
      nktot=0;
      nplane0=0;
      for(int pNo=0;pNo<nPacked;pNo++) {
       double re,im; int x,y,z;
       fread(&(re),sizeof(double),n,fp);
       fread(&(im),sizeof(double),n,fp);
       fread(&(x),sizeof(int),n,fp);
       fread(&(y),sizeof(int),n,fp);
       fread(&(z),sizeof(int),n,fp);
       if(pNo==0){minx=x; maxx=x;}
       if(x<minx){minx=x;}
       if(x>maxx){maxx=x;}
       if(x==0){nplane0++;}
       nktot++;
       if(x==0 && y==0 && z==0 && doublePack)break;
      }//endfor
      fclose(fp);

  }//endif::binary

//===================================================================================
// Set a few parameters

  if(minx<0){minx+=nx;}
  n    = minx; 
  minx = maxx; 
  maxx = n;

  if(doublePack){nPacked=nktot+nplane0-1;}

//----------------------------------------------------------------------------------
   }//end routine
//===================================================================================


//===================================================================================
// Some simple range checking : needs some love
//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::simpleRangeCheck(){ 
//===================================================================================

//  rangeExit(launchNLeesFromRho,"launchNLeesFromRho",0);
  rangeExit(prioFFTMsg,"prioFFTMsg",1);
  rangeExit(stateOutputOn,"stateOutputOn",1);
  rangeExit(atmOutputOn,"atmOutputOn",1);
  rangeExit(localAtomBarrier,"localAtomBarrier",1);
  rangeExit(localEnergyBarrier,"localEnergyBarrier",1);
  rangeExit(lbpaircalc,"lbpaircalc",1);
  rangeExit(lbgspace,"lbgspace",1);
  rangeExit(lbdensity,"lbgspace",1);
  rangeExit(useGMulticast,"useGMulticast",1);
  rangeExit(useCommlibMulticast,"useCommlibMulticast",1);
  rangeExit(useCommlib,"useCommlib",1);
  rangeExit(usePairEtoM,"usePairEtoM",1);
  rangeExit(usePairDirectSend,"usePairDirectSend",1);
  rangeExit(PCCollectTiles,"PCCollectTiles",1);
  rangeExit(PCdelayBWSend,"PCdelayBWSend",1);    
  rangeExit(PCstreamBWout,"PCstreamBWout",1);
  rangeExit(doublePack,"doublePack",1);
  rangeExit(conserveMemory,"conserveMemory",1);
  rangeExit(fftprogresssplit,"fftprogresssplit",0);
  rangeExit(fftprogresssplitReal,"fftprogresssplitReal",0);
  rangeExit(rhoGHelpers,"rhoGHelpers",0);
  rangeExit(rhoRsubplanes,"rhoRsubplanes",0);
  rangeExit(numMulticastMsgs,"numMulticastMsgs",0);
  rangeExit(PCSpanFactor,"numMulticastMsgs",0);
  rangeExit(toleranceInterval,"toleranceInterval;",0);
  rangeExit(numChunks,"numChunks;",0);
  rangeExit(phantomSym,"phantomSym;",1);
  rangeExit(gSpaceSum,"gSpaceSum;",1);
  rangeExit(useCuboidMap,"useCuboidMap;",1);
  rangeExit(useCuboidMapRS,"useCuboidMapRS;",1);
  rangeExit(useCentroidMap,"useCentroidMap;",1);
  rangeExit(loadMapFiles,"loadMapFiles:",1);
  rangeExit(dumpMapFiles,"dumpMapFiles:",1);
  rangeExit(useCentroidMapRho,"useCentroidMapRho;",1);
  rangeExit(numChunksAsym,"numChunksAsym;",0);
  rangeExit(numChunksSym,"numChunksSym;",0);

//---------------------------------------------------------------------------------
  }//end routine
//===================================================================================

//============================================================================
// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Config::rangeExit(int param, char *name, int iopt){
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


//===================================================================================
// Consistency Checks on the input
//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
  void Config::Finale(int nkf1,int nkf2,int nkf3,int nplane_x,int nplane_x_rho){
//===================================================================================
// Code deficiency checks

    if(nkf1!=nkf2 || nkf1!=nkf3){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   Only Cubic boxes for now\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

    if(doublePack!= 1){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   Non-double Pack code is broken\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

//===================================================================================
// Mapping checks

#ifndef CMK_VERSION_BLUEGENE
    if(useCuboidMap || useCuboidMapRS){
       PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       PRINTF("   useCuboidMap requires CMK_VERSION_BLUEGENE\n");
       PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       EXIT(1);
    }//endif
#endif

    if(Gstates_per_pe<1 || Gstates_per_pe>nstates){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   The number of Gstates per pe must be >=1 < num states\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

    if(Rstates_per_pe<1 || Rstates_per_pe>nstates){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   The number of Rstates per pe must be >=1 < num states\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

//===================================================================================
// Chare size checks

    if(gExpandFact<1.0){
      PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      PRINTF("   Chare array expansion factor out of range %g\n",gExpandFact);
      PRINTF("   This probably could work but I'd check first.\n");
      PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    }//endif

    if(gExpandFactRho<1.0){
      PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      PRINTF("   RhoChare array expansion factor out of range %g\n",gExpandFactRho);
      PRINTF("   This probably could work but I'd check first.\n");
      PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    }//endif

    if(nchareG<nplane_x){
      PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      PRINTF("   Too few g-space chares %d %d\n",nplane_x,nchareG);
      PRINTF("   This probably could work but I'd check first.\n");
      PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    }//endif

    if(nchareRhoG<nplane_x_rho){
      PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      PRINTF("   Too few rhog-space chares %d %d\n",nplane_x_rho,nchareRhoG);
      PRINTF("   This probably could work but I'd check first.\n");
      PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    }//endif

//===================================================================================
// Ortho/PC checks

    if(gSpaceSum && !usePairDirectSend){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   gSpaceSum requires usePairDirectSend\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

    if(scalc_per_plane<1){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   Scalc per plane cannot be less then 1\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

    if (nstates % sGrainSize != 0){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   number of states should be divisible by S matrix grain-size\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

    if (sGrainSize %orthoGrainSize != 0){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   S matrix grain-size should be divisible by orthoGrainSize\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

    if (sGrainSize %lambdaGrainSize != 0){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   S matrix grain-size should be divisible by lambdaGrainSize\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

    if(usePairEtoM==1 && useCommlib!=1){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   EachToMany pairCalc requires Commlib!\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

    if(phantomSym && !gSpaceSum){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   Current implementation of phantomSym requires gSpaceSum\n");
      PRINTF("   The price of midnight hacking sessions.\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

    if(gemmSplitFWk > sGrainSize){
       PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       PRINTF("   gemmSplitFWk %d must not be greater than sGrainSize %d !\n",
               gemmSplitFWk, sGrainSize);
       PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       EXIT(1);
    }//endif

    if(gemmSplitFWm > sGrainSize){
       PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       PRINTF("   gemmSplitFWm %d must not be greater than sGrainSize %d !\n",
               gemmSplitFWm, sGrainSize);
       PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       EXIT(1);
    }//endif

    if(gemmSplitBW > sGrainSize){
       PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       PRINTF("   gemmSplitBW %d must not be greater than sGrainSize %d !\n",
               gemmSplitBW, sGrainSize);
       PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       EXIT(1);
    }//endif

//===================================================================================
// Density Controls

    if(useCommlibMulticast+useGMulticast!=1){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   You can't use both the g and commlib multicast\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

    if (nstates / numMulticastMsgs <= 0){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   Problem in the configuration of number of mcast msgs");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

    if(ees_eext_opt==0 && nchareHartAtmT>1){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   Multiple external/hartree atm-type chares only");
      PRINTF("   enabled for the ees external method");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }

//===================================================================================
// Nonlocal Controls 

    if(numSfGrps<1 || numSfGrps> natm_nl){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   The number of sf atm groups must be >=1 < natm_nl\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

    if(numSfDups<1||numSfDups>nstates){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   The number of sf dup groups must be >=1 < num states\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

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
