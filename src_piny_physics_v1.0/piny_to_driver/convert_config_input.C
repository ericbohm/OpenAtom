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

#define _COOL_CONVERSION_ON_

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
using namespace std;
#include "../../include/configure.h"

int main();

//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
int main(){
//===================================================================================

  Config config;
  int num_dict_fun;
  int num_dict_rho, num_dict_state, num_dict_pc;
  int num_dict_nl,  num_dict_gen,   num_dict_atm;
  DICT_WORD *dict_fun;
  DICT_WORD *dict_rho, *dict_state, *dict_pc;
  DICT_WORD *dict_nl, *dict_gen, *dict_atm;
  DICT_WORD word;            

  FILE *fp;
  char input_name[1024];
  char output_name[1024];


//===================================================================================
// read in the old keywords

  printf("Enter your old input file: "); scanf("%s",input_name);
  fp = fopen(input_name,"r");
    int eol   = (int)'\n';
    int nline = 0;
    int ch    = fgetc(fp);
    if(ch != EOF){ 
      while( (ch = fgetc(fp)) != EOF){
        if(ch == eol){nline++;}
      }//endif
    }//endif
  fclose(fp);   
  
  PINY_NAME *keyarg  = (PINY_NAME *)malloc(sizeof(PINY_NAME)*nline)-1;
  PINY_NAME *keyword = (PINY_NAME *)malloc(sizeof(PINY_NAME)*nline)-1;
  int *ifound        = (int *)malloc(sizeof(int)*nline)-1;


  fp = fopen(input_name,"r");
    for(int i=1;i<=nline;i++){
      ifound[i]=0;
      fscanf(fp,"%s %s\n",keyword[i],keyarg[i]); 
    }//endfor
  fclose(fp);     

//===================================================================================
// Set up the dictionaries

  config.useCommlib = 1;   //default value
  config.nstates    = 128; //default value

  config.set_config_dict_fun  (&num_dict_fun  ,&dict_fun);
  config.set_config_dict_rho  (&num_dict_rho  ,&dict_rho);
  config.set_config_dict_state(&num_dict_state,&dict_state);
  config.set_config_dict_pc   (&num_dict_pc,   &dict_pc);
  config.set_config_dict_nl   (&num_dict_nl,   &dict_nl);
  config.set_config_dict_gen  (&num_dict_gen,  &dict_gen);
  config.set_config_dict_atm  (&num_dict_atm,  &dict_atm);

//===================================================================================
// Load the keywords into the dictionary

  config.load_cpaimd_config(dict_rho,  num_dict_rho,  keyarg, keyword, nline, ifound);
  config.load_cpaimd_config(dict_state,num_dict_state,keyarg, keyword, nline, ifound);
  config.load_cpaimd_config(dict_pc,   num_dict_pc,   keyarg, keyword, nline, ifound);
  config.load_cpaimd_config(dict_nl,   num_dict_nl,   keyarg, keyword, nline, ifound);
  config.load_cpaimd_config(dict_atm,  num_dict_atm,  keyarg, keyword, nline, ifound);
  config.load_cpaimd_config(dict_gen,  num_dict_gen,  keyarg, keyword, nline, ifound); 

  for(int i=1;i<=nline;i++){
    if(ifound[i]!=1){
      printf("keyword %s found %d times\n",keyword[i],ifound[i]);
    }//endif
  }//endfor

//===================================================================================

  sprintf(output_name,"%s.out",input_name);
  fp = fopen(output_name,"w");
   config.write_cpaimd_config(fp,dict_rho,  num_dict_rho,  dict_fun[1].keyword);
   config.write_cpaimd_config(fp,dict_state,num_dict_state,dict_fun[2].keyword);
   config.write_cpaimd_config(fp,dict_pc,   num_dict_pc,   dict_fun[3].keyword);
   config.write_cpaimd_config(fp,dict_nl,   num_dict_nl,   dict_fun[4].keyword);
   config.write_cpaimd_config(fp,dict_gen,  num_dict_gen,  dict_fun[5].keyword);
   config.write_cpaimd_config(fp,dict_atm,  num_dict_atm,  dict_fun[6].keyword);
  fclose(fp);

  return 1;
//===================================================================================
  }//end routine
//===================================================================================



//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_fun  (int *num_dict  ,DICT_WORD **dict){
//==================================================================================
//  I) Malloc the dictionary                                              

  num_dict[0] = 6;
  *dict = (DICT_WORD *)malloc(num_dict[0]*sizeof(DICT_WORD))-1;

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

  num_dict[0] = 22;
  *dict = (DICT_WORD *)malloc(num_dict[0]*sizeof(DICT_WORD))-1;

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
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  //  4)\useGHartInsRHart{}
    ind =   4;
    strcpy((*dict)[ind].keyword,"useGHartInsRHart");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  //  5)\useRHartInsGHart{}
    ind =   5;
    strcpy((*dict)[ind].keyword,"useRHartInsGHart");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  //  6)\useCentroidMapRho{}
    ind =   6;
    strcpy((*dict)[ind].keyword,"useCentroidMapRho");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  //  7)\rhorpriority{}
    ind =   7;
    strcpy((*dict)[ind].keyword,"rhorpriority");
    strcpy((*dict)[ind].keyarg,"2000000");    

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
    (*dict)[ind].iflag = 1;
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
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 14)\useGIns1RhoRP{}
    ind =  14;
    strcpy((*dict)[ind].keyword,"useGIns1RhoRP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 15)\useGIns2RhoRP{}
    ind =  15;
    strcpy((*dict)[ind].keyword,"useGIns2RhoRP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 16)\useGIns3RhoRP{}
    ind =  16;
    strcpy((*dict)[ind].keyword,"useGIns3RhoRP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 17)\useGByrdInsRhoRBP{}
    ind =  17;
    strcpy((*dict)[ind].keyword,"useGByrdInsRhoRBP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 18)\useRInsRhoGP{}
    ind =  18;
    strcpy((*dict)[ind].keyword,"useRInsRhoGP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 19)\useRInsIGXRhoGP{}
    ind =  19;
    strcpy((*dict)[ind].keyword,"useRInsIGXRhoGP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 20)\useRInsIGYRhoGP{}
    ind =  20;
    strcpy((*dict)[ind].keyword,"useRInsIGYRhoGP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 21)\useRInsIGZRhoGP{}
    ind =  21;
    strcpy((*dict)[ind].keyword,"useRInsIGZRhoGP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 22)\prioEextFFTMsg{}
    ind =  22;
    strcpy((*dict)[ind].keyword,"prioEextFFTMsg");
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_state(int *num_dict ,DICT_WORD **dict){
//==================================================================================
//  I) Malloc the dictionary                                              

  num_dict[0] = 23;
  *dict = (DICT_WORD *)malloc(num_dict[0]*sizeof(DICT_WORD))-1;

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
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
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
    (*dict)[ind].iflag = 1;
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
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 17)\lbgspace{}
    ind=17;
    strcpy((*dict)[ind].keyword,"lbgspace");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 18)\doublePack{}
    ind=18;
    strcpy((*dict)[ind].keyword,"doublePack");
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 19)\useCuboidMap{}
    ind=19;
    strcpy((*dict)[ind].keyword,"useCuboidMap");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 20)\useCuboidMapRS{}
    ind=20;
    strcpy((*dict)[ind].keyword,"useCuboidMapRS");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 21)\useCentroidMap{}
    ind=21;
    strcpy((*dict)[ind].keyword,"useCentroidMap");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 22)\useGssInsRealP{}
    ind=22;
    strcpy((*dict)[ind].keyword,"useGssInsRealP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 23)\useMssInsGP{}
    ind=23;
    strcpy((*dict)[ind].keyword,"useMssInsGP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
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
  *dict = (DICT_WORD *)malloc(num_dict[0]*sizeof(DICT_WORD))-1;

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
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 2)\PCCollectTiles{}
    ind=2;
    strcpy((*dict)[ind].keyword,"PCCollectTiles");
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 3)\PCdelayBWSend{}
    ind=3;
    strcpy((*dict)[ind].keyword,"PCdelayBWSend");
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 4)\PCstreamBWout{}
    ind=4;
    strcpy((*dict)[ind].keyword,"PCstreamBWout");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
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
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 7)\useOrthoHelpers{}
    ind=7;
    strcpy((*dict)[ind].keyword,"useOrthoHelpers");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 8)\useOrthoSection{}
    ind=8;
    strcpy((*dict)[ind].keyword,"useOrthoSection");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 9)\useOrthoSectionRed{}
    ind=9;
    strcpy((*dict)[ind].keyword,"useOrthoSectionRed");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
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
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 22)\phantomSym{}
    ind=22;
    strcpy((*dict)[ind].keyword,"phantomSym");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 23)\lbpaircalc{}
    ind=23;
    strcpy((*dict)[ind].keyword,"lbpaircalc");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
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
    (*dict)[ind].iflag = 1;
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
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 31)\usePairEtoM{}
    ind=31;
    strcpy((*dict)[ind].keyword,"usePairEtoM");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;

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
  *dict = (DICT_WORD *)malloc(num_dict[0]*sizeof(DICT_WORD))-1;

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
    (*dict)[ind].iflag = 1;
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
    strcpy((*dict)[ind].keyword,"launchNLeesFromRho");
    strcpy((*dict)[ind].keyarg,"off");
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 8)\useGssInsRealPP{}
    ind=8;
    strcpy((*dict)[ind].keyword,"useGssInsRealPP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 9)\useMssInsGPP{}
    ind=9;
    strcpy((*dict)[ind].keyword,"useMssInsGPP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;

//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_gen (int *num_dict ,DICT_WORD **dict){
//==================================================================================
//  I) Malloc the dictionary                                              

  num_dict[0] = 6;
  *dict = (DICT_WORD *)malloc(num_dict[0]*sizeof(DICT_WORD))-1;

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
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 4)\useGMulticast{}
    ind=4;
    strcpy((*dict)[ind].keyword,"useGMulticast");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
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
    (*dict)[ind].iflag = 1;
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
  *dict = (DICT_WORD *)malloc(num_dict[0]*sizeof(DICT_WORD))-1;

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
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 2)\localEnergyBarrier{}
    ind=2;
    strcpy((*dict)[ind].keyword,"localEnergyBarrier"); 
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
  //-----------------------------------------------------------------------------
  // 3)\atmOutputOn{}
    ind=3;
    strcpy((*dict)[ind].keyword,"atmOutputOn");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
    (*dict)[ind].iflag = 1;
//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void Config::write_cpaimd_config(FILE *fp,DICT_WORD *dict, int num_dict, char *fun_key){
//========================================================================
//     I) Write out meta key word 

   fprintf(fp,"================================================\n");
   fprintf(fp,"cccccccccccccccccccccccccccccccccccccccccccccccc\n");
   fprintf(fp,"================================================\n");
   fprintf(fp,"~%s[\n",fun_key);
   for(int i = 1;i<=num_dict;i++){
     if(dict[i].iuset==1){
       fprintf(fp,"   \\%s{%s}\n",dict[i].keyword,dict[i].keyarg);
     }//endif
   }//endfo
   fprintf(fp,"]\n================================================\n\n\n");

//========================================================================
   }//end routine 
//========================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::load_cpaimd_config (DICT_WORD *dict,  int num_dict,  
                                 PINY_NAME *keyarg, PINY_NAME *keyword, 
                                 int nline, int *ifound){
//==================================================================================

  for(int i=1;i<=nline;i++){
    for(int j=1;j<=num_dict;j++){
      if(strcasecmp(keyword[i],dict[j].keyword)==0){
        ifound[i]++;
        dict[j].iuset=1;
        if(dict[j].iflag==0){
          strcpy(dict[j].keyarg,keyarg[i]);
	}else{
          strcpy(dict[j].keyarg,"on");
	  if(strcasecmp(keyword[i],"0")==0){strcpy(dict[j].keyarg,"off");}
	}//endif
      }//endif
    }//endfor
  }//endfor

//==================================================================================
  }//end routine
//==================================================================================
