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
#include "../../include/CPcharmParaInfo.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../class_defs/PINY_INIT/PhysicsParamTrans.h"
#if CMK_PROJECTIONS_USE_ZLIB
#include "zlib.h"
#endif

//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::readConfig(char* input_name,int nstates_in, int nkf1, int nkf2, int nkf3, 
                        int maxIter_in,int ibinary_opt,int natm_nl_in, int fftopt_in,
                        int numPes_in, int natm_typ_in,int ees_eext_opt_in,
                        int gen_wave_in,int ncoef, int cp_min_opt, int nchareRhoRHart,
                        int doublePack_in,int pi_beads, int nkpoint, int ntemper, int nspin)
//===================================================================================
   {//begin routine
//===================================================================================

  int num_dict_fun;
  int num_dict_rho, num_dict_state, num_dict_pc;
  int num_dict_nl,  num_dict_gen,   num_dict_map;
  int num_dict_nfreq;
  DICT_WORD *dict_fun;
  DICT_WORD *dict_rho, *dict_state, *dict_pc;
  DICT_WORD *dict_nl, *dict_gen, *dict_map;
  DICT_WORD *dict_nfreq;
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
// set the Uber parameters

  UberImax = pi_beads; //pi_beads : fixing > 1 implementation now 
  UberJmax = nkpoint;  //nkpoint  : fixing > 1 implementation now 
  UberKmax = ntemper;  //ntemper must be 1 for now
  UberMmax = nspin;    //nspin   spin not yet in here

  // Warn the folks when dicey things are going down
  if(nkpoint>1){
    PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("    Danger, Danger, nkpoint > 1 = %d\n",nkpoint);
    PRINTF("    Put on your debugging shoes and get ready to boogy\n");
    PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n");
  }//endif

  if(nspin>1){
    PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("    Danger, Danger, nspin > 1 = %d\n",nspin);
    PRINTF("    This is not yet supported in any way, shape or form\n");
    PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n");
    EXIT(1);
  }//endif

  if(ntemper>1){
    PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("    Danger, Danger, ntemper > 1 = %d\n",ntemper);
    PRINTF("    This is not yet supported in any way, shape or form\n");
    PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n");
  }//endif

//===================================================================================
// Tell everyone you are busy 

  PRINTF("  =============================================================\n");
  PRINTF("  Reading charm parallel input from file : %s\n",input_name);
  PRINTF("  -------------------------------------------------------------\n\n");

  if(PINY_MAXWORD!=MAXWORD || PINY_MAXLINE != MAXLINE){
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("   Incorrect word and line sizes\n");
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

//===================================================================================
// set the Uber parameters
  temperCycle=0; // this needs to be set from piny somehow
  UberImax = pi_beads; //pi_beads : fixing > 1 implementation now 
  UberJmax = nkpoint;  //nkpoint  : fixing > 1 implementation now 
  UberKmax = ntemper;  //ntemper 
  UberMmax = nspin;    //nspin   

  // Warn the folks when dicey things are going down
  if(pi_beads>1){
    PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("    Danger, Danger, pi_beads > 1 = %d\n",pi_beads);
    PRINTF("    Put on your debugging shoes and get ready to boogy\n");
    PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n");
  }//endif

  if(nkpoint>1){
    PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("    Danger, Danger, nkpoint > 1 = %d\n",nkpoint);
    PRINTF("    This is not yet supported in any way, shape or form in this version\n");
    PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n");
    EXIT(1);
  }//endif

  if(nspin>1){
    PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("    Danger, Danger, nspin > 1 = %d\n",nspin);
    PRINTF("    Put on your debugging shoes and get ready to boogy\n");
    PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n");
  }//endif

  if(ntemper>1){
    PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    PRINTF("    Danger, Danger, ntemper > 1 = %d\n",ntemper);
    PRINTF("    Put on your debugging shoes and get ready to boogy\n");
    PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n");
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
  gen_wave     = gen_wave_in;
  doublePack   = doublePack_in;

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
  set_config_dict_map  (&num_dict_map,  &dict_map);
  set_config_dict_nfreq(&num_dict_nfreq,&dict_nfreq);

//===================================================================================
// Read the input file and fill the dictionaries with user input

  fp = cfopen((const char *) input_name,"r");

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
        case 6 : put_word_dict(&word,dict_map, num_dict_map, fun_key, nline,
                               nkey, nfun_key, input_name); break;
        case 7 : put_word_dict(&word,dict_nfreq, num_dict_nfreq, fun_key, nline,
                               nkey, nfun_key, input_name); break;
      }//end switch
    }// end while 
  }//end while

  fclose(fp);  

//===================================================================================
// Take the information out of the dictionary and put it in the class

  set_config_params_gen  (dict_gen,  dict_fun[5].keyword,input_name);
  int iflag = 1; if(useCommlib!=1){iflag=0;}  // change default commlib option

  set_config_params_rho  (dict_rho,  dict_fun[1].keyword, input_name, iflag);
  set_config_params_state(dict_state,dict_fun[2].keyword, input_name, iflag);
  set_config_params_pc   (dict_pc,   dict_fun[3].keyword, input_name);
  set_config_params_nl   (dict_nl,   dict_fun[4].keyword, input_name, iflag);
  set_config_params_map  (dict_map,  dict_fun[6].keyword, input_name);
  set_config_params_nfreq(dict_nfreq,dict_fun[7].keyword, input_name);

  simpleRangeCheck(); // redundant checking

//===================================================================================
// Read some info from the state file and check for consistency with PINY

  int sizex,sizey,sizez,nPacked,minx,maxx;

  // check input directories if gen_wave is off
  if(gen_wave==0){
    for(int ispin=0;ispin<nspin;ispin++){
    for(int ikpt=0;ikpt<nkpoint;ikpt++){
    for(int ibead=0;ibead<pi_beads;ibead++){
    for(int itemper=0;itemper<ntemper;itemper++){
      sprintf (fname, "%s/Spin.%d_Kpt.%d_Bead.%d_Temper.%d/ChkDirExistOAtom",dataPath,ispin,ikpt,ibead,itemper);
      FILE *fpck = fopen(fname,"w");
      if(fpck==NULL){
        sprintf (fname, "%s/Spin.%d_Kpt.%d_Bead.%d_Temper.%d/ChkDirExistOAtom",dataPath,ispin,ikpt,ibead,itemper);
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        PRINTF("   Input directory, %s , is not present\n",fname);
        PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        EXIT(1);
      }//endif
      fclose(fpck);
    }}}}// end fors : input directory check
  }//endif : gen_wave off

  if(gen_wave==0){
    sprintf (fname, "%s/Spin.0_Kpt.0_Bead.0_Temper.0/state1.out", dataPath);
    PRINTF("  Opening state file : %s\n",fname);
  }//endif

  readStateInfo(nPacked,minx,maxx,sizex,sizey,sizez,fname,ibinary_opt,
                nkf1,nkf2,nkf3,ncoef);

  if(gen_wave==0){
    sprintf (fname, "%s/Spin.0_Kpt.0_Bead.0_Temper.0/state1.out", dataPath);
    PRINTF("  Closing state file : %s\n\n",fname);
  }//endif

  //Check for the output directories
  for(int ispin   = 0; ispin   < nspin   ; ispin++   ){
  for(int ikpt    = 0; ikpt    < nkpoint ; ikpt++    ){
  for(int ibead   = 0; ibead   < pi_beads; ibead++   ){
  for(int itemper = 0; itemper < ntemper ; itemper++ ){
    sprintf (fname, "%s/Spin.%d_Kpt.%d_Bead.%d_Temper.%d/ChkDirExistOAtom",dataPathOut,ispin,ikpt,ibead,itemper);
    FILE *fpck = fopen(fname,"w");
    if(fpck==NULL){
      sprintf (fname, "%s/Spin.%d_Kpt.%d_Bead.%d_Temper.%d",dataPathOut,ispin,ikpt,ibead,itemper);
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   Output directory, %s , is not present\n",fname);
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif
    fclose(fpck);
  }}}}// end fors : directory check
  
//===================================================================================
// FFT sizes 

  if(sizex!=nkf1 || sizey!=nkf2 || sizez !=nkf3){
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("   Incorrect FFT size in state files.\n");
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
  }//endif

  numFFTPoints = nkf1 * nkf2 * nkf3;
  nRplane_x    = nkf1;
  nGplane_x    = maxx-minx+1;
  numData      = nPacked;

//===================================================================================
// Set rhoG and stateG chare array sizes based on FFT stuff

  int nplane_x     = nGplane_x;
  int nplane_x_rho = 2*nGplane_x-1;
  PRINTF("  nplane = %d and nplane_rho = %d for the current system\n",
         nplane_x,nplane_x_rho);

  double temp      = (gExpandFact)*((double)nplane_x);
  nchareG          = ((int)temp);

  double temp_rho  = (gExpandFactRho)*((double)nplane_x_rho);
  nchareRhoG       = ((int)temp_rho);

  scalc_per_plane  = (nstates/sGrainSize)*(nstates/sGrainSize);

//===================================================================================
// Improve user parameters and/or try to optimize unset parameters

  guesstimateParmsConfig(sizez,dict_gen,dict_rho,dict_state,dict_pc,dict_nl,dict_map, 
                         nchareRhoRHart, nplane_x_rho, natm_typ);

//===================================================================================
// Final consistency checks

  Finale(nkf1,nkf2,nkf3,nplane_x,nplane_x_rho,cp_min_opt);

//===================================================================================
// Output your parameter choices to the screen

  // Tokenize the input cpaimd config file name to remove all directory paths
  char *cfgInputName  = new char[strlen(input_name)+1];
  char *cfgOutputName = NULL;
  char dirSeparator[] = "/"; ///< @warning: Will we ever run on Windows and get screwed?
  char *tokenized     = strtok( strcpy(cfgInputName,input_name), dirSeparator);
  while (tokenized != NULL)
  {
      cfgOutputName = tokenized;
      tokenized = strtok(NULL,dirSeparator);
  }
  // The cpaimd config output file is written in the current directory (and not where the input file is located)
  sprintf(fname,"%s.out",cfgOutputName);
  fp = cfopen((const char*) fname,"w");
   write_cpaimd_config(fp,dict_rho,  num_dict_rho,  dict_fun[1].keyword);
   write_cpaimd_config(fp,dict_state,num_dict_state,dict_fun[2].keyword);
   write_cpaimd_config(fp,dict_pc,   num_dict_pc,   dict_fun[3].keyword);
   write_cpaimd_config(fp,dict_nl,   num_dict_nl,   dict_fun[4].keyword);
   write_cpaimd_config(fp,dict_gen,  num_dict_gen,  dict_fun[5].keyword);
   write_cpaimd_config(fp,dict_map,  num_dict_map,  dict_fun[6].keyword);
   write_cpaimd_config(fp,dict_nfreq,num_dict_nfreq,dict_fun[7].keyword);
  fclose(fp);
  delete [] cfgInputName;
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
  cfree(&dict_map[1]  ,"Config::readCconfig");
  cfree(&dict_nfreq[1],"Config::readCconfig");

//===================================================================================
// Tell Everyone you are done

  PRINTF("  -------------------------------------------------------------\n");
  PRINTF("  Completed reading charm parallel input from file : %s\n",input_name);
  PRINTF("  =============================================================\n\n");


//----------------------------------------------------------------------------------
   }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_fun  (int *num_dict  ,DICT_WORD **dict){
//==================================================================================
//  I) Malloc the dictionary                                              

  num_dict[0] = 7;
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
  //  6)~charm_conf_map_def[ ]
    ind = 6;
    strcpy((*dict)[ind].error_mes," ");
    strcpy((*dict)[ind].keyword,"charm_conf_map_def");
    strcpy((*dict)[ind].keyarg," ");
  //------------------------------------------------------------------------------
  //  6)~charm_conf_nfreq_def[ ]
    ind = 7;
    strcpy((*dict)[ind].error_mes," ");
    strcpy((*dict)[ind].keyword,"charm_conf_nfreq_def");
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

  num_dict[0] = 26;
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
  //  6)\rhorpriority{}
    ind =   6;
    strcpy((*dict)[ind].keyword,"rhorpriority");
    strcpy((*dict)[ind].keyarg,"2000000");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  //  7)\rhogpriority{}
    ind =   7;
    strcpy((*dict)[ind].keyword,"rhogpriority");
    strcpy((*dict)[ind].keyarg,"2000000");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  //  8)\gExpandFactRho{}
    ind =   8;
    strcpy((*dict)[ind].keyword,"gExpandFactRho");
    strcpy((*dict)[ind].keyarg,"1.25");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 9)\lbdensity{}
    ind =  9;
    strcpy((*dict)[ind].keyword,"lbdensity");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 10)\rhoGHelpers{}
    ind =  10;
    strcpy((*dict)[ind].keyword,"rhoGHelpers");
    strcpy((*dict)[ind].keyarg,"1");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 11)\rhoRsubplanes{}
    ind =  11;
    strcpy((*dict)[ind].keyword,"rhoRsubplanes");
    strcpy((*dict)[ind].keyarg,"1");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 12)\useGIns0RhoRP{}
    ind =  12;
    strcpy((*dict)[ind].keyword,"useGIns0RhoRP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 13)\useGIns1RhoRP{}
    ind =  13;
    strcpy((*dict)[ind].keyword,"useGIns1RhoRP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 14)\useGIns2RhoRP{}
    ind =  14;
    strcpy((*dict)[ind].keyword,"useGIns2RhoRP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 15)\useGIns3RhoRP{}
    ind =  15;
    strcpy((*dict)[ind].keyword,"useGIns3RhoRP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 16)\useGByrdInsRhoRBP{}
    ind =  16;
    strcpy((*dict)[ind].keyword,"useGByrdInsRhoRBP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 17)\useRInsRhoGP{}
    ind =  17;
    strcpy((*dict)[ind].keyword,"useRInsRhoGP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 18)\useRInsIGXRhoGP{}
    ind =  18;
    strcpy((*dict)[ind].keyword,"useRInsIGXRhoGP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 19)\useRInsIGYRhoGP{}
    ind =  19;
    strcpy((*dict)[ind].keyword,"useRInsIGYRhoGP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 20)\useRInsIGZRhoGP{}
    ind =  20;
    strcpy((*dict)[ind].keyword,"useRInsIGZRhoGP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 21)\prioEextFFTMsg{}
    ind =  21;
    strcpy((*dict)[ind].keyword,"prioEextFFTMsg");
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 22)\rhoLineOrder{}
    ind =  22;
    strcpy((*dict)[ind].keyword,"rhoLineOrder");
    strcpy((*dict)[ind].keyarg,"skip");    
    strcpy((*dict)[ind].error_mes,"skip,none,random");
//----------------------------------------------------------------------------------
  // 23)\nchareHartAtmT{}
    ind =  23;
    strcpy((*dict)[ind].keyword,"nchareHartAtmT");
    strcpy((*dict)[ind].keyarg,"1");    
    strcpy((*dict)[ind].error_mes," number >=1 and <= natmtype");
  //-----------------------------------------------------------------------------
  // 24)\rhoSubPlaneBalance{}
    ind =  24;
    strcpy((*dict)[ind].keyword,"rhoSubPlaneBalance");
    strcpy((*dict)[ind].keyarg,"off");
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 25)\rhoGToRhoRMsgCombine{}
    ind =  25;
    strcpy((*dict)[ind].keyword,"rhoGToRhoRMsgCombine");
    strcpy((*dict)[ind].keyarg,"off");
    strcpy((*dict)[ind].error_mes,"on/off");
//----------------------------------------------------------------------------------

  // 23)\nchareVdW{}
    ind =  26;
    strcpy((*dict)[ind].keyword,"nchareVdW");
    strcpy((*dict)[ind].keyarg,"0");    
    strcpy((*dict)[ind].error_mes," number >=0 and <=  N_v^2");
  //-----------------------------------------------------------------------------

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
  //  6)\rhorpriority{}
    ind =   6;
    sscanf(dict[ind].keyarg,"%d",&rhorpriority);
    if(rhorpriority<=0){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  7)\rhogpriority{}
    ind =   7;
    sscanf(dict[ind].keyarg,"%d",&rhogpriority);
    if(rhogpriority<=0){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  8)\gExpandFactRho{}
    ind =   8;
    sscanf(dict[ind].keyarg,"%lg",&gExpandFactRho);
    if(gExpandFactRho<=0){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 9)\lbdensity{}
    ind =  9;
    parse_on_off(dict[ind].keyarg,&lbdensity,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 10)\rhoGHelpers{}
    ind =  10;
    sscanf(dict[ind].keyarg,"%d",&rhoGHelpers);
    if(rhoGHelpers<=0){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 11)\rhoRsubplanes{}
    ind =  11;
    sscanf(dict[ind].keyarg,"%d",&rhoRsubplanes);
    if(rhoRsubplanes<=0){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 12)\useGIns0RhoRP{}
    ind =  12;
    parse_on_off(dict[ind].keyarg,&useGIns0RhoRP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useGIns0RhoRP=0;}
  //-----------------------------------------------------------------------------
  // 13)\useGIns1RhoRP{}
    ind =  13;
    parse_on_off(dict[ind].keyarg,&useGIns1RhoRP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useGIns1RhoRP=0;}
  //-----------------------------------------------------------------------------
  // 14)\useGIns2RhoRP{}
    ind =  14;
    parse_on_off(dict[ind].keyarg,&useGIns2RhoRP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useGIns2RhoRP=0;}
  //-----------------------------------------------------------------------------
  // 15)\useGIns3RhoRP{}
    ind =  15;
    parse_on_off(dict[ind].keyarg,&useGIns3RhoRP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useGIns3RhoRP=0;}
  //-----------------------------------------------------------------------------
  // 16)\useGByrdInsRhoRBP{}
    ind =  16;
    parse_on_off(dict[ind].keyarg,&useGByrdInsRhoRBP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useGByrdInsRhoRBP=0;}
  //-----------------------------------------------------------------------------
  // 17)\useRInsRhoGP{}
    ind =  17;
    parse_on_off(dict[ind].keyarg,&useRInsRhoGP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useRInsRhoGP=0;}
  //-----------------------------------------------------------------------------
  // 18)\useRInsIGXRhoGP{}
    ind =  18;
    parse_on_off(dict[ind].keyarg,&useRInsIGXRhoGP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useRInsIGXRhoGP=0;}
  //-----------------------------------------------------------------------------
  // 19)\useRInsIGYRhoGP{}
    ind =  19;
    parse_on_off(dict[ind].keyarg,&useRInsIGYRhoGP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useRInsIGYRhoGP=0;}
  //-----------------------------------------------------------------------------
  // 20)\useRInsIGZRhoGP{}
    ind =  20;
    parse_on_off(dict[ind].keyarg,&useRInsIGZRhoGP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useRInsIGZRhoGP=0;}
  //-----------------------------------------------------------------------------
  // 21)\prioEextFFTMsg{}
    ind =  21;
    parse_on_off(dict[ind].keyarg,&prioEextFFTMsg,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 22)\rhoLineOrder{}
    ind =  22;
    ierr= 0;
    if(strcasecmp(dict[ind].keyarg,"none")==0)  {rhoLineOrder =-1; ierr++;}
    if(strcasecmp(dict[ind].keyarg,"skip")==0)  {rhoLineOrder = 0; ierr++;}
    if(strcasecmp(dict[ind].keyarg,"random")==0){rhoLineOrder = 1; ierr++;}
    if(ierr!=1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 23)\nchareHartAtmT{}
    ind =   23;
    sscanf(dict[ind].keyarg,"%d",&nchareHartAtmT);
    if(nchareHartAtmT<1 || nchareHartAtmT>natm_typ){
      keyarg_barf(dict,input_name,fun_key,ind);
    }//endif
  //-----------------------------------------------------------------------------
  // 24)\rhoSubPlaneBalance{}
    ind =  24;
    parse_on_off(dict[ind].keyarg,&rhoSubPlaneBalance,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 25)\rhoGToRhoRMsgCombine{}
    ind =  25;
    parse_on_off(dict[ind].keyarg,&rhoGToRhoRMsgComb,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 26)\nchareVdW{}
    ind =   26;
    sscanf(dict[ind].keyarg,"%d",&nchareVdW);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){nchareVdW=0;}
  //-----------------------------------------------------------------------------

  }//end routine

//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_state(int *num_dict ,DICT_WORD **dict){
//==================================================================================
//  I) Malloc the dictionary                                              

  num_dict[0] = 19;
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
    strcpy((*dict)[ind].keyarg,"./STATES");    
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
  //  7)\stateOutput{}
    ind=7;
    strcpy((*dict)[ind].keyword,"stateOutput");
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"off/on");
  //-----------------------------------------------------------------------------
  //  8)\psipriority{}
    ind=8;
    strcpy((*dict)[ind].keyword,"psipriority");
    strcpy((*dict)[ind].keyarg,"400000000");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  //  9)\prioFFTMsg{}
    ind=9;
    strcpy((*dict)[ind].keyword,"prioFFTMsg");
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 10)\rsfftpriority{}
    ind=10;
    strcpy((*dict)[ind].keyword,"rsfftpriority");
    strcpy((*dict)[ind].keyarg,"2000000");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 11)\gsfftpriority{}
    ind=11;
    strcpy((*dict)[ind].keyword,"gsfftpriority");
    strcpy((*dict)[ind].keyarg,"1000000");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 12)\rsifftpriority{}
    ind=12;
    strcpy((*dict)[ind].keyword,"rsifftpriority");
    strcpy((*dict)[ind].keyarg,"100000000");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 13)\gsifftpriority{}
    ind=13;
    strcpy((*dict)[ind].keyword,"gsifftpriority");
    strcpy((*dict)[ind].keyarg,"200000000");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 14)\conserveMemory{}
    ind=14;
    strcpy((*dict)[ind].keyword,"conserveMemory");
    strcpy((*dict)[ind].keyarg,"0");    
    strcpy((*dict)[ind].error_mes,"-1 to optimize for speed 0 for normal or 1 to minimize footprint");
  //-----------------------------------------------------------------------------
  // 15)\lbgspace{}
    ind=15;
    strcpy((*dict)[ind].keyword,"lbgspace");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 16)\doublePack{}
    ind=16;
    strcpy((*dict)[ind].keyword,"doublePack");
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"Deprecated keyword");
  //-----------------------------------------------------------------------------
  // 17)\useGssInsRealP{}
    ind=17;
    strcpy((*dict)[ind].keyword,"useGssInsRealP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 18)\useMssInsGP{}
    ind=18;
    strcpy((*dict)[ind].keyword,"useMssInsGP");
    if(useCommlib==0){strcpy((*dict)[ind].keyarg,"off");}
    if(useCommlib==1){strcpy((*dict)[ind].keyarg,"on");}
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 19)\dataPathOut\{}
    ind=19;
    strcpy((*dict)[ind].keyword,"dataPathOut");
    strcpy((*dict)[ind].keyarg,"./STATES_OUT");    
    strcpy((*dict)[ind].error_mes,"a directory tree");

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
  //  7)\stateOutput{}
    ind=7;
    parse_on_off(dict[ind].keyarg,&stateOutput,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(stateOutput==0){
        PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
        PRINTF("    You are in a state of danger! Your stateOutput is off! \n");
        PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n");
    }//
  //-----------------------------------------------------------------------------
  //  8)\psipriority{}
    ind=8;
    sscanf(dict[ind].keyarg,"%d",&psipriority);
    if(psipriority<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  9)\prioFFTMsg{}
    ind=9;
    parse_on_off(dict[ind].keyarg,&prioFFTMsg,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 10)\rsfftpriority{}
    ind=10;
    sscanf(dict[ind].keyarg,"%d",&rsfftpriority);
    if(rsfftpriority<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 11)\gsfftpriority{}
    ind=11;
    sscanf(dict[ind].keyarg,"%d",&gsfftpriority);
    if(gsfftpriority<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 12)\rsifftpriority{}
    ind=12;
    sscanf(dict[ind].keyarg,"%d",&rsifftpriority);
    if(rsifftpriority<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 13)\gsifftpriority{}
    ind=13;
    sscanf(dict[ind].keyarg,"%d",&gsifftpriority);
    if(gsifftpriority<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 14)\conserveMemory{}
    ind=14;
    sscanf(dict[ind].keyarg,"%d",&conserveMemory);
    if(conserveMemory<-1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 15)\lbgspace{}
    ind=15;
    parse_on_off(dict[ind].keyarg,&lbgspace,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 16)\doublePack{}
    ind=16;
    if(dict[ind].iuset==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 17)\useGssInsRealP{}
    ind=17;
    parse_on_off(dict[ind].keyarg,&useGssInsRealP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useGssInsRealP=0;}
  //-----------------------------------------------------------------------------
  // 18)\useMssInsGP{}
    ind=18;
    parse_on_off(dict[ind].keyarg,&useMssInsGP,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(iflag==0 && dict[ind].iuset==0){useMssInsGP=0;}
  //-----------------------------------------------------------------------------
  // 19)\dataPathOut\{}
    ind=19;
    strcpy(dataPathOut, dict[ind].keyarg);
    if(strcasecmp(dataPathOut,dataPath)==0){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   input %s and output %s datapath must differ \n",dataPath,dataPathOut);
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_pc (int *num_dict ,DICT_WORD **dict){
//==================================================================================
//  I) Malloc the dictionary                                              

  num_dict[0] = 32;
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
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 2)\PCCollectTiles{}
    ind=2;
    strcpy((*dict)[ind].keyword,"PCCollectTiles");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 3)\PCdelayBWSend{}
    ind=3;
    strcpy((*dict)[ind].keyword,"PCdelayBWSend");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 4)\PCstreamBWout{}
    ind=4;
    strcpy((*dict)[ind].keyword,"PCstreamBWout");
    strcpy((*dict)[ind].keyarg,"on");    
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
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 9)\useOrthoSectionRed{}
    ind=9;
    strcpy((*dict)[ind].keyword,"useOrthoSectionRed");
    strcpy((*dict)[ind].keyarg,"on");    
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
  // 20)\useBWBarrier{}
    ind=20;
    strcpy((*dict)[ind].keyword,"useBWBarrier");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 21)\phantomSym{}
    ind=21;
    strcpy((*dict)[ind].keyword,"phantomSym");
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 22)\lbpaircalc{}
    ind=22;
    strcpy((*dict)[ind].keyword,"lbpaircalc");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 23)\lambdapriority{}
    ind=23;
    strcpy((*dict)[ind].keyword,"lambdapriority");
    strcpy((*dict)[ind].keyarg,"300000000");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 24)\toleranceInterval{}
    ind=24;
    strcpy((*dict)[ind].keyword,"toleranceInterval");
    strcpy((*dict)[ind].keyarg,"1");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 25)\gSpaceSum{}
    ind=25;
    strcpy((*dict)[ind].keyword,"gSpaceSum");
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 26)\numChunks{}
    ind=26;
    strcpy((*dict)[ind].keyword,"numChunks");
    strcpy((*dict)[ind].keyarg,"1");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 27)\numChunksSym{}
    ind=27;
    strcpy((*dict)[ind].keyword,"numChunksSym");
    strcpy((*dict)[ind].keyarg,"1");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 28)\numChunksAsym{}
    ind=28;
    strcpy((*dict)[ind].keyword,"numChunksAsym");
    strcpy((*dict)[ind].keyarg,"1");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 29)\prioBW{}
    ind=29;
    strcpy((*dict)[ind].keyword,"prioBW");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 30)\usePairEtoM{}
    ind=30;
    strcpy((*dict)[ind].keyword,"usePairEtoM");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 31)\invsqr_tolerance{}
    ind=31;
    strcpy((*dict)[ind].keyword,"invsqr_tolerance");
    sprintf((*dict)[ind].keyarg,"%0.16g",1e-15);
    strcpy((*dict)[ind].error_mes,"a number > 0.0");
  //-----------------------------------------------------------------------------
  // 32)\invsqr_max_iter{}
    ind=32;
    strcpy((*dict)[ind].keyword,"invsqr_max_iter");
    sprintf((*dict)[ind].keyarg,"%d",1000);
    strcpy((*dict)[ind].error_mes,"a number >= 0");
  //-----------------------------------------------------------------------------


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
  // 20)\useBWBarrier{}
    ind=20;
    parse_on_off(dict[ind].keyarg,&useBWBarrier,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 21)\phantomSym{}
    ind=21;
    parse_on_off(dict[ind].keyarg,&phantomSym,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 22)\lbpaircalc{}
    ind=22;
    parse_on_off(dict[ind].keyarg,&lbpaircalc,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 23)\lambdapriority{}
    ind=23;
    sscanf(dict[ind].keyarg,"%d",&lambdapriority);
    if(lambdapriority<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 24)\toleranceInterval{}
    ind=24;
    sscanf(dict[ind].keyarg,"%d",&toleranceInterval);
    if(toleranceInterval<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 25)\gSpaceSum{}
    ind=25;
    parse_on_off(dict[ind].keyarg,&gSpaceSum,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 26)\numChunks{}
    ind=26;
    sscanf(dict[ind].keyarg,"%d",&numChunks);
    if(numChunks<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 27)\numChunksSym{}
    ind=27;
    sscanf(dict[ind].keyarg,"%d",&numChunksSym);
    if(numChunksSym<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 28)\numChunksAsym{}
    ind=28;
    sscanf(dict[ind].keyarg,"%d",&numChunksAsym);
    if(numChunksAsym<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 29)\prioBW{}
    ind=29;
    parse_on_off(dict[ind].keyarg,&prioBW,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 30)\usePairEtoM{}
    ind=30;
    parse_on_off(dict[ind].keyarg,&usePairEtoM,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 31)\invsqr_tolerance{}
    ind=31;
    sscanf(dict[ind].keyarg,"%lg",&invsqr_tolerance);
    if(invsqr_tolerance<=1e-100){
      //      CkPrintf("what the hell is wrong with %.15g \n",invsqr_tolerance);
      keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 32)\invsqr_max_iter{}
    ind=32;
    sscanf(dict[ind].keyarg,"%d",&invsqr_max_iter);
    if(invsqr_max_iter<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------


//===================================================================================
// Clean up the user input

  if(numChunks>numChunksSym) {
     PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
     PRINTF("   Cleaning up numchunks > numchunksSym : %d %d\n",numChunks,numChunksSym);
     PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n\n");
     numChunksSym=numChunks;
     sprintf(dict[28].keyarg,"%d",numChunksSym);
  }//endif

  if(numChunks>numChunksAsym){
     PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
     PRINTF("   Cleaning up numchunks > numchunksAsym : %d %d\n",numChunks,numChunksAsym);
     PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n\n");
     numChunksAsym=numChunks;
     sprintf(dict[29].keyarg,"%d",numChunksAsym);
  }//endif

  if(lambdaGrainSize==nstates && orthoGrainSize!=nstates){
     PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
     PRINTF("   Cleaning up lambdaGrainSize==nstates && orthoGrainSize!=nstates\n");
     PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n\n");
     lambdaGrainSize=orthoGrainSize;
     sprintf(dict[10].keyarg,"%d",lambdaGrainSize);
  }//endif

  if(PCstreamFWblock>0){
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    PRINTF("   PC forward streaming is broken\n");
    PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    EXIT(1);
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

  num_dict[0] = 8;
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
    strcpy((*dict)[ind].keyarg,"1000");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 2)\fftprogresssplitReal{}
    ind=2;
    strcpy((*dict)[ind].keyword,"fftprogresssplitReal");
    strcpy((*dict)[ind].keyarg,"1000");    
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 3)\useCommlib{}
    ind=3;
    strcpy((*dict)[ind].keyword,"useCommlib");
    strcpy((*dict)[ind].keyarg,"off");
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 4)\useGMulticast{}
    ind=4;
    strcpy((*dict)[ind].keyword,"useGMulticast");
    strcpy((*dict)[ind].keyarg,"on");
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 5)\numMulticastMsgs{}
    ind=5;
    strcpy((*dict)[ind].keyword,"numMulticastMsgs");
    strcpy((*dict)[ind].keyarg,"4");
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 6)\useCommlibMulticast{}
    ind=6;
    strcpy((*dict)[ind].keyword,"useCommlibMulticast");
    strcpy((*dict)[ind].keyarg,"off");
    strcpy((*dict)[ind].error_mes,"on/off");
//----------------------------------------------------------------------------------
  // 7)\useTimeKeeper{}
    ind=7;
    strcpy((*dict)[ind].keyword,"useTimeKeeper");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 8)\atmOutput{}
    ind=8;
    strcpy((*dict)[ind].keyword,"atmOutput");
    strcpy((*dict)[ind].keyarg,"on");    
    strcpy((*dict)[ind].error_mes,"on/off");

//----------------------------------------------------------------------------------
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
  //-----------------------------------------------------------------------------
  // 8)\atmOutput{}
    ind=8;
    parse_on_off(dict[ind].keyarg,&atmOutput,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(atmOutput==0){
        PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
        PRINTF("   You are in danger! Your atmOutput is off! \n");
        PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n\n");
    }//
 
//----------------------------------------------------------------------------------
  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_map (int *num_dict ,DICT_WORD **dict){
//==================================================================================
//  I) Malloc the dictionary                                              

  num_dict[0] = 20;
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
  // 1)\torusMap{}
    ind = 1;
    strcpy((*dict)[ind].keyword,"torusMap");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 2)\useCuboidMap{}
    ind = 2;
    strcpy((*dict)[ind].keyword,"useCuboidMap");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 3)\useCuboidMapRS{}
    ind = 3;
    strcpy((*dict)[ind].keyword,"useCuboidMapRS");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 4)\useCentroidMap{}
    ind = 4;
    strcpy((*dict)[ind].keyword,"useCentroidMap");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  // 5)\useCentroidMapRho{}
    ind = 5;
    strcpy((*dict)[ind].keyword,"useCentroidMapRho");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 6)\Gstates_per_pe{}
    ind = 6;
    strcpy((*dict)[ind].keyword,"Gstates_per_pe");
    sprintf((*dict)[ind].keyarg,"%d",nstates);  
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 7)\Rstates_per_pe{}
    ind = 7;
    strcpy((*dict)[ind].keyword,"Rstates_per_pe");
    sprintf((*dict)[ind].keyarg,"%d",nstates);
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 8)\loadMapFiles{}
    ind=8;
    strcpy((*dict)[ind].keyword,"loadMapFiles");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 9)\dumpMapFiles{}
    ind=9;
    strcpy((*dict)[ind].keyword,"dumpMapFiles");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 10)\useRhoExclusionMap{}
    ind=10;
    strcpy((*dict)[ind].keyword,"useRhoExclusionMap");
    strcpy((*dict)[ind].keyarg,"on");
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 11)\fakeTorus{}
    ind = 11;
    strcpy((*dict)[ind].keyword,"fakeTorus");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 12)\torusDimNX{}
    ind = 12;
    strcpy((*dict)[ind].keyword,"torusDimNX");
    sprintf((*dict)[ind].keyarg,"%d",1);  
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 13)\torusDimNY{}
    ind = 13;
    strcpy((*dict)[ind].keyword,"torusDimNY");
    sprintf((*dict)[ind].keyarg,"%d",1);  
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 14)\torusDimNZ{}
    ind = 14;
    strcpy((*dict)[ind].keyword,"torusDimNZ");
    sprintf((*dict)[ind].keyarg,"%d",1);  
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 15)\torusDimNT{}
    ind = 15;
    strcpy((*dict)[ind].keyword,"torusDimNT");
    sprintf((*dict)[ind].keyarg,"%d",1);  
    strcpy((*dict)[ind].error_mes,"a number > 0");
  //-----------------------------------------------------------------------------
  // 16)\useStrictCuboid{}
    ind = 16;
    strcpy((*dict)[ind].keyword,"useStrictCuboid");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // 17)\useReductionExclusionMap{}
    ind = 17;
    strcpy((*dict)[ind].keyword,"useReductionExclusionMap");
    strcpy((*dict)[ind].keyarg,"on");
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 18)\excludePE0{}
    ind = 18;
    strcpy((*dict)[ind].keyword,"excludePE0");
    strcpy((*dict)[ind].keyarg,"off");
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 19)\dumpMapCoordFiles{}
    ind=19;
    strcpy((*dict)[ind].keyword,"dumpMapCoordFiles");
    strcpy((*dict)[ind].keyarg,"off");    
    strcpy((*dict)[ind].error_mes,"on/off");
  //-----------------------------------------------------------------------------
  // 20)\forceMappingAxis{}
    ind = 20;
    strcpy((*dict)[ind].keyword,"forceMappingAxis");
    sprintf((*dict)[ind].keyarg,"%d",-1);  
    strcpy((*dict)[ind].error_mes,"-1  for no force, 0=X 1=Y 2=Z");
  //-----------------------------------------------------------------------------

  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_params_map (DICT_WORD *dict, char *fun_key, char *input_name){
//===================================================================================

  int ind, ierr; 

//===================================================================================
// Fill the class with data 

  //-----------------------------------------------------------------------------
  //  1)\torusMap{}
    ind = 1;
    parse_on_off(dict[ind].keyarg, &torusMap, &ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  2)\useCuboidMap{}
    ind = 2;
    parse_on_off(dict[ind].keyarg,&useCuboidMap,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  3)\useCuboidMapRS{}
    ind = 3;
    parse_on_off(dict[ind].keyarg,&useCuboidMapRS,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  4)\useCentroidMap{}
    ind = 4;
    parse_on_off(dict[ind].keyarg,&useCentroidMap,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  5)\useCentroidMapRho{}
    ind = 5;
    parse_on_off(dict[ind].keyarg,&useCentroidMapRho,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  6)\Gstates_per_pe{}
    ind = 6;
    sscanf(dict[ind].keyarg,"%d",&Gstates_per_pe);
    if(Gstates_per_pe<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  7)\Rstates_per_pe{}
    ind = 7;
    sscanf(dict[ind].keyarg,"%d",&Rstates_per_pe);
    if(Rstates_per_pe<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  8)\loadMapFiles{}
    ind = 8;
    parse_on_off(dict[ind].keyarg,&loadMapFiles,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  9)\dumpMapFiles{}
    ind = 9;
    parse_on_off(dict[ind].keyarg,&dumpMapFiles,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 10)\useRhoExclusionMap{}
    ind=10;
    parse_on_off(dict[ind].keyarg,&useRhoExclusionMap,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(dict[ind].iuset==0){useRhoExclusionMap=1;}
  //-----------------------------------------------------------------------------
  //  11)\fakeTorus{}
    ind = 11;
    parse_on_off(dict[ind].keyarg, &fakeTorus, &ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  12)\torusDimNX{}
    ind = 12;
    sscanf(dict[ind].keyarg,"%d",&torusDimNX);
    if(torusDimNX<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  13)\torusDimNY{}
    ind = 13;
    sscanf(dict[ind].keyarg,"%d",&torusDimNY);
    if(torusDimNY<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  14)\torusDimNZ{}
    ind = 14;
    sscanf(dict[ind].keyarg,"%d",&torusDimNZ);
    if(torusDimNZ<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  15)\torusDimNT{}
    ind = 15;
    sscanf(dict[ind].keyarg,"%d",&torusDimNT);
    if(torusDimNT<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  16)\useStrictCuboid{}
    ind = 16;
    parse_on_off(dict[ind].keyarg,&useStrictCuboid,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------

  //  17)\useReductionExclusionMap{}
    ind = 17;
    parse_on_off(dict[ind].keyarg,&useReductionExclusionMap,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(dict[ind].iuset==0){useReductionExclusionMap=1;}
  //-----------------------------------------------------------------------------
  //  18)\excludePE0{}
    ind = 18;
    parse_on_off(dict[ind].keyarg,&excludePE0,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
    if(dict[ind].iuset==0){excludePE0=0;}
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  // 19)\dumpMapCoordFiles{}
    ind = 19;
    parse_on_off(dict[ind].keyarg,&dumpMapCoordFiles,&ierr);
    if(ierr==1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  //  20)\forceMappingAxis{}
    ind = 20;
    sscanf(dict[ind].keyarg,"%d",&forceMappingAxis);
    if(forceMappingAxis>2){keyarg_barf(dict,input_name,fun_key,ind);}

  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_dict_nfreq (int *num_dict ,DICT_WORD **dict){
//==================================================================================
//  I) Malloc the dictionary                                              

  num_dict[0] = 8;
  *dict = (DICT_WORD *)cmalloc(num_dict[0]*sizeof(DICT_WORD),"set_config_dict_nfreq")-1;

//=================================================================================
//  II) Initialize the user set option(did the user set the key word      

  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iuset    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].iflag    = 0;}
  for(int i=1;i<=num_dict[0];i++){(*dict)[i].key_type = 1;} //spec only once

  // Setup the default values
  int nfreq_cpintegrate       = 400;   ///< CPINTEGRATE::CP_integrate_min_STD, CPINTEGRATE::CP_integrate_min_CG
  int nfreq_cplocal_hartext   = 4;     ///< CPLOCAL::CP_hart_eext_calc
  int nfreq_cplocal_eeshart   = 100;   ///< CPLOCAL::eesHartEextGchare
  int nfreq_cplocal_eesewald  = 100;   ///< CPLOCAL::eesEwaldGchare
  int nfreq_cpnonlocal_eke    = 400;   ///< CPNONLOCAL::CP_eke_calc
  int nfreq_cpnonlocal_eesfwd = 100;   ///< CPNONLOCAL::eesProjGchare, CPNONLOCAL::eesYlmOnD
  int nfreq_cpnonlocal_eesbk  = 100;   ///< CPNONLOCAL::eesPsiForcGspace
  int nfreq_xcfnctl           = 8;     ///< CPXCFNCTS::CP_exc_calc, CPXCFNCTS::CP_getGGAFunctional

//=================================================================================
// III) Set up the dictionary
  int ind;
  //-----------------------------------------------------------------------------
  // 1)\integrate{}
    ind = 1;
    strcpy((*dict)[ind].keyword,"integrate");
    sprintf((*dict)[ind].keyarg,"%d",nfreq_cpintegrate);
    strcpy((*dict)[ind].error_mes,"freq at which network progress should be called on BGL");
  //-----------------------------------------------------------------------------
  // 2)\localHartExt{}
    ind = 2;
    strcpy((*dict)[ind].keyword,"localHartExt");
    sprintf((*dict)[ind].keyarg,"%d",nfreq_cplocal_hartext);
    strcpy((*dict)[ind].error_mes,"freq at which network progress should be called on BGL");
  //-----------------------------------------------------------------------------
  // 3)\localEesHart{}
    ind = 3;
    strcpy((*dict)[ind].keyword,"localEesHart");
    sprintf((*dict)[ind].keyarg,"%d",nfreq_cplocal_eeshart);
    strcpy((*dict)[ind].error_mes,"freq at which network progress should be called on BGL");
  //-----------------------------------------------------------------------------
  // 4)\localEesEwald{}
    ind = 4;
    strcpy((*dict)[ind].keyword,"localEesEwald");
    sprintf((*dict)[ind].keyarg,"%d",nfreq_cplocal_eesewald);
    strcpy((*dict)[ind].error_mes,"freq at which network progress should be called on BGL");
  //-----------------------------------------------------------------------------
  // 5)\nonlocalEke{}
    ind = 5;
    strcpy((*dict)[ind].keyword,"nonlocalEke");
    sprintf((*dict)[ind].keyarg,"%d",nfreq_cpnonlocal_eke);
    strcpy((*dict)[ind].error_mes,"freq at which network progress should be called on BGL");
  //-----------------------------------------------------------------------------
  // 6)\nonlocalEesFwd{}
    ind = 6;
    strcpy((*dict)[ind].keyword,"nonlocalEesFwd");
    sprintf((*dict)[ind].keyarg,"%d",nfreq_cpnonlocal_eesfwd);
    strcpy((*dict)[ind].error_mes,"freq at which network progress should be called on BGL");
  //-----------------------------------------------------------------------------
  // 7)\nonlocalEesBk{}
    ind = 7;
    strcpy((*dict)[ind].keyword,"nonlocalEesBk");
    sprintf((*dict)[ind].keyarg,"%d",nfreq_cpnonlocal_eesbk);
    strcpy((*dict)[ind].error_mes,"freq at which network progress should be called on BGL");
  //-----------------------------------------------------------------------------
  // 8)\xcfnctl{}
    ind = 8;
    strcpy((*dict)[ind].keyword,"xcfnctl");
    sprintf((*dict)[ind].keyarg,"%d",nfreq_xcfnctl);
    strcpy((*dict)[ind].error_mes,"freq at which network progress should be called on BGL");
  }//end routine
//===================================================================================


//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::set_config_params_nfreq  (DICT_WORD *dict, char *fun_key, char *input_name){
//===================================================================================

  int ind,ierr;

//===================================================================================
  //-----------------------------------------------------------------------------
  // 1)\integrate{}
    ind = 1;
    sscanf(dict[ind].keyarg,"%d",&nfreq_cpintegrate);
    if(nfreq_cpintegrate<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 2)\localHartExt{}
    ind = 2;
    sscanf(dict[ind].keyarg,"%d",&nfreq_cplocal_hartext);
    if(nfreq_cplocal_hartext<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 3)\localEesHart{}
    ind = 3;
    sscanf(dict[ind].keyarg,"%d",&nfreq_cplocal_eeshart);
    if(nfreq_cplocal_eeshart<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 4)\localEesEwald{}
    ind = 4;
    sscanf(dict[ind].keyarg,"%d",&nfreq_cplocal_eesewald);
    if(nfreq_cplocal_eesewald<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 5)\nonlocalEke{}
    ind = 5;
    sscanf(dict[ind].keyarg,"%d",&nfreq_cpnonlocal_eke);
    if(nfreq_cpnonlocal_eke<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 6)\nonlocalEesFwd{}
    ind = 6;
    sscanf(dict[ind].keyarg,"%d",&nfreq_cpnonlocal_eesfwd);
    if(nfreq_cpnonlocal_eesfwd<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 7)\nonlocalEesBk{}
    ind = 7;
    sscanf(dict[ind].keyarg,"%d",&nfreq_cpnonlocal_eesbk);
    if(nfreq_cpnonlocal_eesbk<1){keyarg_barf(dict,input_name,fun_key,ind);}
  //-----------------------------------------------------------------------------
  // 8)\xcfnctl{}
    ind = 8;
    sscanf(dict[ind].keyarg,"%d",&nfreq_xcfnctl);
    if(nfreq_xcfnctl<1){keyarg_barf(dict,input_name,fun_key,ind);}
  }//end routine
//===================================================================================




//=============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=============================================================================
/**
 * Guesses decent values for configuration parameters based on user
 * set values, system size, and the number of cores available.
 *
 * For large pe runs, scheme is thus: 
 *    1. find nchareG which is a multiple of at least one torus dimension
 *      a. on BG you can just choose a convenient power of 2 and it
 *         will work well for most torus cases.
 *      b. BGW full machine runs have the 20 or 40 (VN) rack dimension
 *         water256 scales better using multiple of that in
 *         nchareG. YMMV.  This result is not reproduced by the
 *         guesstimate code, if you have access to >=20480 processors a
 *         little time tweaking the config file yourself is a good idea.
 *
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
 *       centroid mapping, useDirectSend and probably a few other
 *       things
 *    6. if nstates%2=0 and useCuboid then useStrictCuboid=1 unless
 *       user set it off in config.
 *     
 *      
 *
 *  for nonpower of 2 systems the scheme is revised:
 * 
 *    1. find a grainsize which meets the following constraints where:
 *      nstates%grainsize=remainder
 *
 *      a. the remainder is less than 8 (the smallest remotely
 *         justifiable choice for orthograinsize).  
 *      b. the grainsize has a factor small enough for either
 *         orthograinsize^2==numprocs, or orthograinsize^2==numprocs/2
 *         and enable orthohelpers.  with orthograinsize > remainder
 *      c. odd numbers should work, but will slightly degrade your
 *         ability to keep the double hummer pipeline loaded in which
 *          will reduce your flops so you should avoid them.
 *       
 * Supported parameters: sGrainSize, orthoGrainSize, gExpandFact,
 * gExpandFactRho, fftprogresssplit, fftprogresssplitReal, numSfGrps,
 * numSfDups, rhoGHelpers, numMulticastMsgs

 */
//=============================================================================
void Config::guesstimateParmsConfig(int sizez,DICT_WORD *dict_gen,DICT_WORD *dict_rho,
                            DICT_WORD *dict_state,DICT_WORD *dict_pc,
                            DICT_WORD *dict_nl,DICT_WORD *dict_map, 
			    int nchareRhoRHart, int nplane_x_rho, int natm_typ){
//=============================================================================
  if(fakeTorus)
    { //
      numPes=torusDimNX * torusDimNY * torusDimNZ * torusDimNT;
      CkPrintf("  Using fake torus node %d X %d X %d X %d numPes %d\n", torusDimNX, torusDimNY, torusDimNZ, torusDimNT, numPes);
    }

    // TODO get the number of instances.  This is a Glenn item as most
    // of the instance stuff will be on the piny side and determining
    // how many instances there are of each type, and then the total
    // sum, will presumably come directly from CPcharmParaInfoGrp.

    // in default case we have only 1 instance, so we just use the 0,0,0th
    // because it is simple.

    numInstances = UberImax * UberJmax * UberKmax * UberMmax; // numIntegrals * numKpoints * numTempers * numSpin;
    numPesPerInstance = numPes / numInstances;
    numPesPerInstance = (numPesPerInstance>0) ? numPesPerInstance : 1;

    int sqrtpes    = (int) sqrt((double)numPes);
    int sqrtstates = (int) sqrt((double)nstates);
    int igo;

//=============================================================================
// gExpandFact not set 

    igo = dict_state[6].iuset;

    if(igo==0){ 
      if(numPes>nGplane_x){
         int i=1;
         double mypow=1;
	 int targetNchare=nGplane_x;
	 // need a way to increase the targetNchare so that
	 // nChareG * (nstates/sGrainSize)^2 * numChunks = numPes
	 // without going overboard on the planes, grains, or chunks.
	 // in the general case we haven't determined any of these
	 // variables yet.  

	 // roughly speaking the desired nchare goes as sqrtpes)
	 if(nGplane_x<sqrtpes) targetNchare=sqrtpes / nGplane_x * nGplane_x;
         while((mypow=pow(2.0, (double)i)) <= targetNchare){i++;}

         gExpandFact = mypow / (double) nGplane_x;
         nchareG     = (int)( gExpandFact * (double) nGplane_x);
         sprintf(dict_state[6].keyarg,"%g",gExpandFact);
      }//endif
    }//endif : gexpandfact not set

    
    nchareG  = (int)(gExpandFact*(double)nGplane_x);
    CkPrintf("  nchareG now %d based on gExpandFact %.5g\n\n", nchareG, gExpandFact);
//=============================================================================
// numChunks, numChunksSym, numChunksAsym, sgrainsize are not set

    igo = (dict_pc[27].iuset+dict_pc[28].iuset+dict_pc[29].iuset+dict_pc[14].iuset);
    if(igo==0){
       int numGrains = nstates/sGrainSize;
       numGrains    *=numGrains;
       numChunks     = numPes/(nchareG*numGrains);
       if(numChunks<1){numChunks=1;}
       int remainder;
       while(numChunks>8){
	 if(!isPow2(nstates)) // not a power of two, find a good grainsize
	   // allowing for remainder
	   {
	     sGrainSize = sGrainSize/2;
	     remainder=approxFactor(nstates,sGrainSize, orthoGrainSize,numPes);
	   }
	 else
	   {  //
	     sGrainSize = sGrainSize/2;
	     remainder=approxFactor(nstates,sGrainSize, orthoGrainSize,numPes);
	   }
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
       sprintf(dict_pc[19].keyarg,"%d",orthoGrainSize);
       if(gemmSplitFWk>sGrainSize) gemmSplitFWk=sGrainSize;
       if(gemmSplitFWm>sGrainSize) gemmSplitFWm=sGrainSize;
       if(gemmSplitOrtho>orthoGrainSize) gemmSplitOrtho=orthoGrainSize;
       if(gemmSplitOrtho%2!=0) gemmSplitOrtho-=1;
       if(gemmSplitFWk%2!=0) gemmSplitFWk-=1;
       if(gemmSplitFWm%2!=0) gemmSplitFWm-=1;
       sprintf(dict_pc[15].keyarg,"%d",gemmSplitFWk);
       sprintf(dict_pc[16].keyarg,"%d",gemmSplitFWm);
       sprintf(dict_pc[17].keyarg,"%d",gemmSplitOrtho);
       dict_pc[19].iuset=1;
    }//endif : chunks not set
    else{
      // CkPrintf("guesstimate using user defined PC decomp\n");
    }

//=============================================================================
// gstates_per_pe not set
    
    igo = dict_map[6].iuset;

    if(numPesPerInstance!=1 && igo==0){
      if(numPesPerInstance<=128){
        Gstates_per_pe=nstates/numPesPerInstance;
      }else{ 
	if(numPesPerInstance>128 && numPesPerInstance<=512){
          Gstates_per_pe=nstates/16;
	}else{
          Gstates_per_pe= nchareG*nstates/numPesPerInstance;
	}//endif
      }//endif
      if(Gstates_per_pe==0){Gstates_per_pe=1;}
      sprintf(dict_map[6].keyarg,"%d",Gstates_per_pe);
    }//endif : Gstates_per_pe not set

//=============================================================================
// rstates_per_pe not set

    igo = dict_map[7].iuset;

    if(numPesPerInstance!=1 && igo==0){
      if(numPesPerInstance<=128){
        Rstates_per_pe=nstates/numPesPerInstance;
      }else{
        Rstates_per_pe= sizez*nstates/numPesPerInstance;
      }//endif
      if (Rstates_per_pe==0){
        Rstates_per_pe=1;
      }//endif
      sprintf(dict_map[7].keyarg,"%d",Rstates_per_pe);
    }//endif : Rstates_per_pe not set

//=============================================================================
// Rejiggering orthograinsize if dumb values employed :
//         This is lame and should be replace with something which finds
//         an even mod of any sGrainSize

    igo = dict_pc[19].iuset;

    if(igo==0){
      if( (sGrainSize%orthoGrainSize !=0) || (sGrainSize==orthoGrainSize)){
        int iii = orthoGrainSize;
        if(orthoGrainSize>32){orthoGrainSize=32;}
        if(orthoGrainSize>sGrainSize){orthoGrainSize=sGrainSize;}
        PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
        PRINTF("    Changing the default choice of orthograinsize\n");
        PRINTF("    from %d to %d\n", iii, orthoGrainSize);
        PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
        sprintf(dict_pc[19].keyarg,"%d",orthoGrainSize);
      }//endif
    }//endif


//=============================================================================
// Rejiggering lambdagrainsize if dumb values employed :
//            This is lame and should be replace with something which finds
//            an even mod of any (non prime) sGrainSize

    igo = dict_pc[10].iuset;

    if((sGrainSize%lambdaGrainSize !=0)|| (sGrainSize==lambdaGrainSize)){
      if(igo==1){
        PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
        PRINTF("    (sGrainSize%%lambdaGrainSize !=0)|| (sGrainSize==lambdaGrainSize)");
        //PRINTF("   lambdaGrainSize=orthoGrainSize = %d\n",orthoGrainSize);
        PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
      }//endif
      lambdaGrainSize=orthoGrainSize;
      sprintf(dict_pc[10].keyarg,"%d",lambdaGrainSize);
    }//endif


//=============================================================================
// Gexpandfact rho

    // In an ideal universe rhoRS, RhoGS, rhoRhart, rhoGhart are all
    // mapped to their own processor.  However, at the upper limit you
    // can only make rho pieces so small before communication overhead
    // swamps any possibly parallel overlap.  So we must cap the
    // decomposition. At the lower limit you want to choose parameters
    // so that the heaviest hitters, rhors and rhog, can be exclusion
    // mapped.  As numproc increases we need to gradually increase
    // gExpandFactRho, the number of subplanes, and the hart atm
    // type.  

    // Missing here is the logic for when ghelpers should be
    // increased.  Difficulty in that primarily lies in the fact that
    // it is unclear under what conditions ghelpers is useful when
    // using the modern EES method.  We're not even going to try to
    // optimize for non-EES, anyone using that creaky old stuff is on
    // their own. -EJB
    
    igo = dict_rho[9].iuset+dict_rho[11].iuset;
    int nchareRhoR      = sizez;
    if(igo==0){ 
      int numRS=nchareRhoR*rhoRsubplanes;
      int numGH=nchareRhoG*rhoGHelpers;
      if(numPesPerInstance < numRS*3){
	useReductionExclusionMap=0;
	
	strcpy(dict_map[17].keyarg,"off");
	CkPrintf("Disabling reduction exclusion map, too few processors to matter\n");
      }
      else
	{
	  useReductionExclusionMap=1;
	  strcpy(dict_map[17].keyarg,"on");
	  CkPrintf("Enabling reduction exclusion map\n");
	}

      bool notdone=true;
      int usedPes=numRS;
      int maxSubplanes=20; // arbitary choice, but probably good
      if(useReductionExclusionMap) usedPes+=nchareG;
      while(numPesPerInstance > usedPes && notdone)
	{

	  // trick here is to keep bumping up both subplanes and
	  // expandfact
	  int target=numPesPerInstance - usedPes;
	  if(numRS<target) target=numRS;
	  if(numPesPerInstance>numRS && nchareRhoR>nchareRhoG)
	    { // bring up expandfact to fill in 
	      if(target>nplane_x_rho && target <= numRS && gExpandFactRho< (double) nplane_x_rho)
		{

		  gExpandFactRho= (double) target / double (nplane_x_rho);
		  //		  CkPrintf("i1 rhoRsubplanes now %d gExpandFactRho now %g nchareHartAtmT now %d making numRS %d nchareRhoG %d\n",rhoRsubplanes, gExpandFactRho, nchareHartAtmT, numRS, nchareRhoG);
		  double temp_rho  = (gExpandFactRho)*((double)nplane_x_rho);
		  nchareRhoG       = ((int)temp_rho);
		  //notdone=false;
		  usedPes=nchareRhoG+numRS;
		  if(useReductionExclusionMap) usedPes+=nchareG;
		}
	    }
	  if (numPesPerInstance>usedPes && numPesPerInstance >nchareRhoR*(rhoRsubplanes+1) && (numPesPerInstance >= nchareRhoRHart* (rhoRsubplanes+1)* nchareHartAtmT) &&(rhoRsubplanes<maxSubplanes))
	    {
	      numRS= (++rhoRsubplanes)*nchareRhoR;
	      usedPes=nchareRhoG+numRS;
	      if(useReductionExclusionMap) usedPes+=nchareG;
	    }
	  else
	    {
	      notdone=false;
	    }
		   
	  if (numPesPerInstance>numRS && nchareRhoR<=nchareRhoG && gExpandFactRho< (double) nplane_x_rho)
	    { // keep ncharerhog close to numRS
	      int target=numPesPerInstance - numRS;
	      if(target>numRS) target=numRS;

	      if ( gExpandFactRho > double (nplane_x_rho))
		{
		  gExpandFactRho=(double)nplane_x_rho;
		}
	      if(gExpandFactRho>rhoRsubplanes) gExpandFactRho=(double) rhoRsubplanes;
	      double temp_rho  = (gExpandFactRho)*((double)nplane_x_rho);
	      //	      CkPrintf("i2 rhoRsubplanes now %d gExpandFactRho now %g nchareHartAtmT now %d making numRS %d nchareRhoG %d\n",rhoRsubplanes, gExpandFactRho, nchareHartAtmT, numRS, nchareRhoG);
	      nchareRhoG       = ((int)temp_rho);
	      usedPes=nchareRhoG+numRS;
	      if(useReductionExclusionMap) usedPes+=nchareG;
	    }	   
	  else
	    {
	      notdone=false;
	    }
	  if(ees_eext_opt==1 && (nchareHartAtmT< natm_typ) && 
	     (nchareRhoRHart* (rhoRsubplanes)* (nchareHartAtmT+1)<=numPesPerInstance/2))
	    {
	      nchareHartAtmT++;
	    }
       }//endwhile
      if(rhoRsubplanes>1)
	{
	  rhoLineOrder=-1;
	  strcpy(dict_rho[22].keyarg,"none");
	  rhoSubPlaneBalance=1;
	  strcpy(dict_rho[24].keyarg,"on");
	  rhoGToRhoRMsgComb=1;
	  strcpy(dict_rho[25].keyarg,"on");
	}

       sprintf(dict_rho[11].keyarg,"%d",rhoRsubplanes);
       sprintf(dict_rho[23].keyarg,"%d",nchareHartAtmT);
       sprintf(dict_rho[8].keyarg,"%g",gExpandFactRho);
       CkPrintf("rhoRsubplanes now %d gExpandFactRho now %g nchareHartAtmT now %d making numRS %d nchareRhoG %d\n",rhoRsubplanes, gExpandFactRho, nchareHartAtmT, numRS, nchareRhoG);
    }//endif
    double temp_rho  = (gExpandFactRho)*((double)nplane_x_rho);
    nchareRhoG       = ((int)temp_rho);



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
       if(numSfGrps==0 && natm_nl>0){numSfGrps=1;}
       sprintf(dict_nl[5].keyarg,"%d",numSfGrps);
    }//endif

    if(natm_nl==0){numSfGrps=0;}

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

    igo = dict_rho[10].iuset;

    if(igo==0){
      if(ees_eext_opt==1)
	{
	  //rhoGHelpers=rhoRsubplanes/2;
	  //leave rhoghelpers alone
	}
      else
	{
	  int temp_rho  = (int) (gExpandFactRho*2.0*((double) nGplane_x + 1.0));
	  if(numPes>temp_rho){rhoGHelpers=numPes/temp_rho;}
	}
      sprintf(dict_rho[10].keyarg,"%d",rhoGHelpers);
	  
    }//endif

//============================================================================
// Turn topo options off if torusMap is off

    igo = dict_map[1].iuset;

    if(igo==1 && torusMap==0) {
      sprintf(dict_map[1].keyarg, "%d", torusMap);
      useCuboidMap = 0;
      sprintf(dict_map[2].keyarg, "%d", useCuboidMap);
      useCuboidMapRS = 0;
      sprintf(dict_map[3].keyarg, "%d", useCuboidMapRS);
      useCentroidMap = 0;
      sprintf(dict_map[4].keyarg, "%d", useCentroidMap);
      //      useCentroidMapRho = 0; works ok on nontorus
      //      sprintf(dict_map[5].keyarg, "%d", useCentroidMapRho);
      useStrictCuboid = 0;
      sprintf(dict_map[16].keyarg, "%d", useStrictCuboid);
    }

//===================================================================================
// Clean up 

  if(rhoRsubplanes>1){
     PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
     PRINTF("   rhoRsubplanes=%d. Disabling rho rs commlib use\n",rhoRsubplanes);
     PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n\n");
     useGHartInsRHart  = 0; strcpy(dict_rho[3].keyarg,"off");
     useRHartInsGHart  = 0; strcpy(dict_rho[4].keyarg,"off");
     useGHartInsRhoRP  = 0; strcpy(dict_rho[5].keyarg,"off");
     useGIns0RhoRP     = 0; strcpy(dict_rho[12].keyarg,"off");
     useGIns1RhoRP     = 0; strcpy(dict_rho[13].keyarg,"off");
     useGIns2RhoRP     = 0; strcpy(dict_rho[14].keyarg,"off");
     useGIns3RhoRP     = 0; strcpy(dict_rho[15].keyarg,"off");
     useGByrdInsRhoRBP = 0; strcpy(dict_rho[16].keyarg,"off");
     useRInsRhoGP      = 0; strcpy(dict_rho[17].keyarg,"off");
     useRInsIGXRhoGP   = 0; strcpy(dict_rho[18].keyarg,"off");
     useRInsIGYRhoGP   = 0; strcpy(dict_rho[19].keyarg,"off");
     useRInsIGZRhoGP   = 0; strcpy(dict_rho[20].keyarg,"off");
  }//endif

  if(nchareHartAtmT>1  && rhoRsubplanes==1){
     PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
     PRINTF("   nchareHartAtmT=%d. Disabling rho hart commlib use\n",nchareHartAtmT);
     PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n\n");
     useGHartInsRHart  = 0; strcpy(dict_rho[3].keyarg,"off");
     useRHartInsGHart  = 0; strcpy(dict_rho[4].keyarg,"off");
     useGHartInsRhoRP  = 0; strcpy(dict_rho[5].keyarg,"off");
  }//endif

//----------------------------------------------------------------------------------


//============================================================================
   }//end routine
//============================================================================



//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
void Config::readStateInfo(int &nPacked,int &minx, int &maxx, int &nx, int &ny, int &nz,
                           const char *fromFile, int ibinary_opt,
                           int nkf1, int nkf2, int nkf3,int ncoef) {
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
  rangeExit(stateOutput,"stateOutput",1);
  rangeExit(atmOutput,"atmOutput",1);
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
  rangeExit(fftprogresssplit,"fftprogresssplit",0);
  rangeExit(fftprogresssplitReal,"fftprogresssplitReal",0);
  rangeExit(rhoGHelpers,"rhoGHelpers",0);
  rangeExit(rhoRsubplanes,"rhoRsubplanes",0);
  rangeExit(numMulticastMsgs,"numMulticastMsgs",0);
  rangeExit(PCSpanFactor,"numMulticastMsgs",0);
  rangeExit(toleranceInterval,"toleranceInterval;",0);
  //  rangeExit(invsqr_tolerance,"invsqr_tolerance;",0);
  rangeExit(invsqr_max_iter,"invsqr_max_iter;",0);
  rangeExit(numChunks,"numChunks;",0);
  rangeExit(phantomSym,"phantomSym;",1);
  rangeExit(gSpaceSum,"gSpaceSum;",1);
  rangeExit(useCuboidMap,"useCuboidMap;",1);
  rangeExit(useStrictCuboid,"useStrictCuboid;",1);
  rangeExit(useCuboidMapRS,"useCuboidMapRS;",1);
  rangeExit(useCentroidMap,"useCentroidMap;",1);
  rangeExit(loadMapFiles,"loadMapFiles:",1);
  rangeExit(dumpMapFiles,"dumpMapFiles:",1);
  rangeExit(dumpMapCoordFiles,"dumpMapCoordFiles:",1);
  rangeExit(useCentroidMapRho,"useCentroidMapRho;",1);
  rangeExit(numChunksAsym,"numChunksAsym;",0);
  rangeExit(numChunksSym,"numChunksSym;",0);

//---------------------------------------------------------------------------------
  }//end routine
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


//===================================================================================
// Consistency Checks on the input
//===================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//===================================================================================
  void Config::Finale(int nkf1,int nkf2,int nkf3,int nplane_x,int nplane_x_rho, int cp_min_opt){
//===================================================================================
// Code deficiency checks

#ifdef _CUBIC_BOXES_ONLY_
    if(nkf1!=nkf2 || nkf1!=nkf3){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   Only Cubic boxes for now\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif
#endif

    if(doublePack!= 1){
      PRINTF("\n  $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      PRINTF("     Non-double Pack code is under development!!\n");
      PRINTF("     Put on your debug shoes and boogy\n");
      PRINTF("  $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n\n");
    }//endif

//===================================================================================
// Mapping checks

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
      PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n\n");
    }//endif

    if(gExpandFactRho<1.0){
      PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      PRINTF("   RhoChare array expansion factor out of range %g\n",gExpandFactRho);
      PRINTF("   This probably could work but I'd check first.\n");
      PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n\n");
    }//endif

    if(nchareG<nplane_x){
      PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      PRINTF("   Too few g-space chares %d %d\n",nplane_x,nchareG);
      PRINTF("   This probably could work but I'd check first.\n");
      PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n\n");
    }//endif

    if(nchareRhoG<nplane_x_rho){
      PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      PRINTF("   Too few rhog-space chares %d %d\n",nplane_x_rho,nchareRhoG);
      PRINTF("   This probably could work but I'd check first.\n");
      PRINTF("   $$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n\n");
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


    if(PCCollectTiles && cp_min_opt!=1)
      {
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   Must NOT PCCollectTiles in dynamics\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
      }

    if((PCCollectTiles && PCstreamBWout)||(!PCCollectTiles && ! PCstreamBWout)) {
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   Must PCCollectTiles and PCstreamBWout are mutually exclusive.   Choose one (PCstreamBWout is usually faster).\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }

    if (nstates % sGrainSize > sGrainSize/2 ){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   Your remainder should be less than 1/2 of your grainsize\n   Or Load Imbalance will substantially degrade performance\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }
    if (nstates % sGrainSize != 0 && !gSpaceSum){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   Must enable gSpaceSum to support number of states not divisible\n    by S matrix grain-size\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif


    if (sGrainSize %orthoGrainSize != 0){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   S matrix grain-size %d must be divisible by orthoGrainSize %d\n",sGrainSize, orthoGrainSize);
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

    if (sGrainSize %lambdaGrainSize != 0){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   S matrix grain-size must be divisible by lambdaGrainSize\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

    if (nstates % sGrainSize >= orthoGrainSize ){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   Your remainder should be less than orthograinsize\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }
    if(usePairEtoM==1 && useCommlib!=1){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   EachToMany pairCalc requires Commlib!\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

    if(phantomSym && !gSpaceSum){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   Current implementation of phantomSym requires gSpaceSum\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

    if(gemmSplitFWk %2 !=0){
       PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       PRINTF("   gemmSplitFWk %d must be an even number !\n",
               gemmSplitFWk);
       PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       EXIT(1);
    }//endif

    if(gemmSplitOrtho %2 !=0){
       PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       PRINTF("   gemmSplitOrtho %d must be an even number !\n",
               gemmSplitOrtho);
       PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       EXIT(1);
    }//endif


    if((gemmSplitFWm %2 !=0) || gemmSplitFWm > sGrainSize){
       PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       PRINTF("   gemmSplitFWm %d must be an even number less than sGrainSize %d !\n",
               gemmSplitFWm, sGrainSize);
       PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       EXIT(1);
    }//endif

    if((gemmSplitBW %2 !=0) || gemmSplitBW > sGrainSize){
       PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       PRINTF("   gemmSplitBW %d must be an even number greater than sGrainSize %d !\n",
               gemmSplitBW, sGrainSize);
       PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       EXIT(1);
    }//endif

    int size1 = (nstates/sGrainSize)*(nstates/sGrainSize)*numChunksSym*nchareG;
    int size2 = (nstates/sGrainSize)*(nstates/sGrainSize)*numChunksAsym*nchareG;
    int size3 = (nstates/orthoGrainSize)*(nstates/orthoGrainSize);
   
    if(size1>numPes || size2 > numPes || size3 > numPes){
      PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
      PRINTF("    You only have %d processors\n",numPes);
      PRINTF("    You have asked for %d Symm  PC\n",size1);
      PRINTF("    You have asked for %d Asymm PC\n",size2);
      PRINTF("    You have asked for %d orthos  \n",size3);
      PRINTF("    nstates=%d sGrainSize=%d numChunks=%d %d orhoGrain=%d\n",
                 nstates,sGrainSize,numChunksSym,numChunksAsym,orthoGrainSize);
      PRINTF("  $$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    }//endif

//===================================================================================
// Density Controls

    if(useCommlibMulticast+useGMulticast!=1){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   Must turn on either useGMulticast or useCommlibMulticast\n");
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
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

    if((numSfGrps<1 && natm_nl>0)|| numSfGrps> natm_nl){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   The number of sf atm groups must be >=1 < natm_nl %d\n",numSfGrps);
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      EXIT(1);
    }//endif

    if((numSfDups<1 && natm_nl>0)|| numSfDups>nstates){
      PRINTF("   @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      PRINTF("   The number of sf dup groups must be >=1 < num states %d\n",numSfDups);
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

// PRE: input sGrainSize is upper bound, nstates & numpes >0
// POST: 
//     oGrainSize is a factor of sGrainSize
//     (nstates/ograinsize)^<=numPes
//     remainder= nstates%sGrainSize
//     remainder < sGrainSize/2
//     remainder < orthoGrainSize
//     remainder returned from function
int Config::approxFactor(int nstates,int &sGrainSize, int &oGrainSize,int numPes)
{


  // we refine grainSize to something with enough factors
  // for orthograinsize, orthograinsize must be factor of
  // sgrain.  So this breaks down to a factorization
  // problem on possible sGrainSizes.  sGrainSize can
  // generally be assumed to be 5 digits max, but we will
  // have to check a bunch of sGrainSizes, so doing it fast
  // will be rewarded.  Our factorization problem is
  // simplified by the fact that we don't care if
  // our factors are primes, but increased by the fact that
  // we need one with a factor that meets the
  // numograin~=numproc constraint.

  // We can actually leap more directly to hunting for the right
  // factor by operating in reverse.  Take the sqrt of the number of
  // processors and jigger that to find a factor of nstates
  int numograin=nstates/oGrainSize;
  int sqrtpes    = (int) sqrt((double)numPes);
  oGrainSize=nstates/sqrtpes;

  // our candidate is a lower bound on ograinsize, we can go larger
  // and leave some processors empty.
  // Empirically we know that finer granularity than ograinsize
  // 8 is extremely unlikely to improve performance.
  if (oGrainSize<8) oGrainSize=8;

  // Our actual constraint here is that ograinsize must be a
  // factor of sgrainsize, while minimizing the remainder


  // given our ograin lower bound and sgrain upper bound
  // we have an n^2 approach to step them towards each other until we
  // find a match.  
  int sGrainSizeLB=sGrainSize/2;
  int remainder=nstates%sGrainSize;
  int sGrainSizeUB=sGrainSize;
  CkPrintf("Approx given sGrainSize %d, nstates %d, numProcs %d \n",sGrainSize,nstates, numPes);
  for(  ; oGrainSize<=sGrainSize; oGrainSize++)
    {

      //TODO: this is accelerated by only testing multiples of oGrainSize
      int cand=sGrainSizeUB/oGrainSize*oGrainSize;
      remainder=nstates%cand;
      for(sGrainSize=cand; ((sGrainSize>sGrainSizeLB) &&

			    (remainder<sGrainSize/4) && (remainder<oGrainSize))
				      ;sGrainSize-=oGrainSize)
	{
	  CkPrintf("Approx trying oGrainSize %d sGrainSize %d \n",oGrainSize, sGrainSize);
	  remainder=nstates%sGrainSize;
	  if((remainder<sGrainSize/4) && (remainder<oGrainSize)) 
	    {
	      CkPrintf("Approx remainder %d \n",remainder);
      
	      return(remainder);
	    }
	}
    }
  // you are screwed
  return(-1);
}



