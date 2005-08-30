/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_sim_parms.c                          */
/*                                                                          */
/*                                                                          */
/* This subprogram reads in user simulation params and echoes them          */
/* and the default parameters to a file                                     */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../class_defs/allclass_gen.h"
#include "../class_defs/allclass_cp.h"
#include "../class_defs/allclass_mdintegrate.h"
#include "../class_defs/allclass_mdatoms.h"
#include "../class_defs/allclass_mdinter.h"
#include "../class_defs/allclass_mdintra.h"
#include "../class_defs/typedefs_par.h"
#include "../proto_defs/proto_sim_params_entry.h"
#include "../proto_defs/proto_sim_params_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_handle_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


void control_sim_params(MDINTEGRATE *mdintegrate, MDATOMS *mdatoms,
                        MDINTER *mdinter,GENERAL_DATA *general_data,
                        MDINTRA *bonded,CP *cp,CLASS_PARSE *class_parse,
                        CP_PARSE *cp_parse,FILENAME_PARSE *filename_parse)

/*=======================================================================*/

{/*begin routine*/ 

/*=======================================================================*/
/*          Local variable declarations                                  */

  int iii;
  int num_dict_fun,num_dict_list,num_dict_cp,num_dict_gen,num_dict_vpot; 
  int num_dict_run,num_dict_nhc,num_dict_vol,num_dict_write,num_dict_pimd;  
  int num_dict_velo,num_dict_msqd,num_dict_iikt_iso,num_dict_ickt_iso;
  int num_dict_rdf; 
                             /* Num: Number of words in the 
                                     dictionary of simulation
                                     key words                   */    
  DICT_WORD *dict_fun,*dict_list,*dict_cp,*dict_gen,*dict_vpot,*dict_run; 
  DICT_WORD *dict_nhc,*dict_vol,*dict_write,*dict_pimd; 
  DICT_WORD *dict_velo,*dict_msqd,*dict_iikt_iso,*dict_ickt_iso,*dict_rdf;
  DICT_WORD word;            
                             /* Str: Dictionary of key words
                                     key arguments, etc;
                                     Lth:num_dict                */ 
  FILE *fp;                  /* Fle: Simulation file pointer     */
  int nline;                 /* Num: Current line number in 
                                     simulation input file       */
  int nkey,nfun_key;         /* Num: Current key, func key       */
  char *fun_key;             /* Chr: Simulation Input File       */
  int ind_key;               /* Num: Index of functional keyword */
  double now_memory;         /* Num: Memory now                  */

/*          Local pointer declarations */
  char *input_name = filename_parse->input_name; 

/*========================================================================*/
/*    I) Write to the screen                                              */

  PRINTF("\n");
  PRINT_LINE_STAR;
  PRINTF("Reading simulation input file %s\n",input_name);
  PRINT_LINE_DASH; PRINTF("\n");

/*=======================================================================*/
/*   II) Set up dictionary and default parameters                        */
/*            (set_sim_dict.c)                                           */

  set_sim_dict_fun(&num_dict_fun,&dict_fun);
  set_sim_dict_list(&num_dict_list,&dict_list);
  set_sim_dict_cp(&num_dict_cp,&dict_cp);
  set_sim_dict_gen(&num_dict_gen,&dict_gen);
  set_sim_dict_vpot(&num_dict_vpot,&dict_vpot);
  set_sim_dict_run(&num_dict_run,&dict_run);
  set_sim_dict_nhc(&num_dict_nhc,&dict_nhc);
  set_sim_dict_vol(&num_dict_vol,&dict_vol);
  set_sim_dict_write(&num_dict_write,&dict_write);
  set_sim_dict_pimd(&num_dict_pimd,&dict_pimd);
  set_sim_dict_velo(&num_dict_velo,&dict_velo);
  set_sim_dict_msqd(&num_dict_msqd,&dict_msqd);
  set_sim_dict_iikt_iso(&num_dict_iikt_iso,&dict_iikt_iso);
  set_sim_dict_ickt_iso(&num_dict_ickt_iso,&dict_ickt_iso);
  set_sim_dict_rdf(&num_dict_rdf,&dict_rdf);

  fun_key     = (char *)cmalloc(MAXWORD*sizeof(char),"control_sim_params");

/*=======================================================================*/
/*  III) Malloc some memory */

  filename_parse->simname            = (char *)cmalloc(MAXWORD*sizeof(char),
                                       "control_sim_params");
  filename_parse->molsetname         = (char *)cmalloc(MAXWORD*sizeof(char),
                                       "control_sim_params");
  filename_parse->dnamei             = (char *)cmalloc(MAXWORD*sizeof(char),
                                       "control_sim_params");
  filename_parse->dnameci            = (char *)cmalloc(MAXWORD*sizeof(char),
                                       "control_sim_params");
  general_data->genfilenames.dname   = (char *)cmalloc(MAXWORD*sizeof(char),
                                       "control_sim_params");
  general_data->genfilenames.dnamei  = (char *)cmalloc(MAXWORD*sizeof(char),
                                       "control_sim_params");
  general_data->genfilenames.iname   = (char *)cmalloc(MAXWORD*sizeof(char),
                                       "control_sim_params");
  general_data->genfilenames.cpname  = (char *)cmalloc(MAXWORD*sizeof(char),
                                       "control_sim_params");
  general_data->genfilenames.cvname  = (char *)cmalloc(MAXWORD*sizeof(char),
                                       "control_sim_params");
  general_data->genfilenames.dnamec  = (char *)cmalloc(MAXWORD*sizeof(char),
                                       "control_sim_params");
  general_data->genfilenames.ccname  = (char *)cmalloc(MAXWORD*sizeof(char),
                                       "control_sim_params");
  general_data->genfilenames.cpparname = (char *)cmalloc(MAXWORD*sizeof(char),
                                       "control_sim_params");
  general_data->genfilenames.centname  = (char *)cmalloc(MAXWORD*sizeof(char),
                                       "control_sim_params");
  general_data->genfilenames.forcename = (char *)cmalloc(MAXWORD*sizeof(char),
                                       "control_sim_params");
  general_data->genfilenames.ksname  = (char *)cmalloc(MAXWORD*sizeof(char),
                                       "control_sim_params");
  general_data->genfilenames.elfname = (char *)cmalloc(MAXWORD*sizeof(char),
                                       "control_sim_params");
  cp->cppseudo.vxc_typ               = (char *)cmalloc(MAXWORD*sizeof(char),
                                       "control_sim_params");
  cp->cppseudo.ggax_typ              = (char *)cmalloc(MAXWORD*sizeof(char),
                                       "control_sim_params");
  cp->cppseudo.ggac_typ              = (char *)cmalloc(MAXWORD*sizeof(char),
                                       "control_sim_params");

  now_memory         = (12*sizeof(char)*MAXWORD)*1.0e-06;
  now_memory        += (10*sizeof(char)*MAXWORD)*1.0e-06;
  now_memory        += (80*sizeof(double))*1.0e-06;
  general_data->tot_memory += now_memory;

  PRINTF("Simulation param. allocation: %g Mbytes; Total memory: %g Mbytes\n",
          now_memory,general_data->tot_memory);

/*=======================================================================*/
/* IV) Open the data file and read in the information                    */

  fp = cfopen(input_name,"r");

  nline = 1;
  nfun_key=0;
  while(get_fun_key(fp,fun_key,&nline,&nfun_key,input_name)){
    get_fun_key_index(fun_key,num_dict_fun,dict_fun,nline,nfun_key,
                      input_name,&ind_key);

    nkey  = 0;
    while(get_word(fp,&word,&nline,&nkey,nfun_key,input_name)){
      switch(ind_key){
        case 1 : put_word_dict(&word,dict_list,num_dict_list,fun_key,nline,
                               nkey,nfun_key,input_name);break;
        case 2 : put_word_dict(&word,dict_cp,num_dict_cp,fun_key,nline,
                               nkey,nfun_key,input_name);break;
        case 3 : put_word_dict(&word,dict_gen,num_dict_gen,fun_key,nline,
                               nkey,nfun_key,input_name);break;
        case 4 : put_word_dict(&word,dict_vpot,num_dict_vpot,fun_key,nline,
                               nkey,nfun_key,input_name);break;
        case 5 : put_word_dict(&word,dict_run,num_dict_run,fun_key,nline,
                               nkey,nfun_key,input_name);break;
        case 6 : put_word_dict(&word,dict_nhc,num_dict_nhc,fun_key,nline,
                               nkey,nfun_key,input_name);break;
        case 7 : put_word_dict(&word,dict_vol,num_dict_vol,fun_key,nline,
                               nkey,nfun_key,input_name);break;
        case 8 : put_word_dict(&word,dict_write,num_dict_write,fun_key,nline,
                               nkey,nfun_key,input_name);break;
        case 9 : put_word_dict(&word,dict_pimd,num_dict_pimd,fun_key,nline,
                               nkey,nfun_key,input_name);break;
	case 10 : put_word_dict(&word,dict_velo,num_dict_velo,fun_key,nline,
				nkey,nfun_key,input_name);break;
	case 11 : put_word_dict(&word,dict_msqd,num_dict_msqd,fun_key,nline,
				nkey,nfun_key,input_name);break;
	case 12 : put_word_dict(&word,dict_iikt_iso,num_dict_iikt_iso,
                                fun_key,nline,nkey,nfun_key,input_name);break;
	case 13 : put_word_dict(&word,dict_ickt_iso,num_dict_ickt_iso,
                                fun_key,nline,nkey,nfun_key,input_name);break;
        case 14 : put_word_dict(&word,dict_rdf,num_dict_rdf,
                                fun_key,nline,nkey,nfun_key,input_name);break;
      }/*end switch*/
    }/* end while */
  }/*end while*/

  fclose(fp);  

/*=======================================================================*/
/*   IV) Stuff the information in the structures                         */
/*        (non-commuting order of calls)                                 */

  set_sim_params_gen (mdintegrate,mdatoms,mdinter,general_data,
                      bonded,cp,class_parse,cp_parse,filename_parse,
                      dict_gen,dict_fun[3].keyword);
  set_sim_params_list(mdintegrate,mdatoms,mdinter,general_data,
                      bonded,cp,class_parse,cp_parse,filename_parse,
                      dict_list,dict_fun[1].keyword);
  set_sim_params_cp  (mdintegrate,mdatoms,mdinter,general_data,
                      bonded,cp,class_parse,cp_parse,filename_parse,
                      dict_cp,dict_fun[2].keyword);
  set_sim_params_vpot(mdintegrate,mdatoms,mdinter,general_data,
                      bonded,cp,class_parse,cp_parse,filename_parse,
                      dict_vpot,dict_fun[4].keyword);
  set_sim_params_run (mdintegrate,mdatoms,mdinter,general_data,
                      bonded,cp,class_parse,cp_parse,filename_parse,
                      dict_run,dict_fun[5].keyword);
  set_sim_params_nhc (mdintegrate,mdatoms,mdinter,general_data,
                      bonded,cp,class_parse,cp_parse,filename_parse,
                      dict_nhc,dict_fun[6].keyword);
  set_sim_params_vol (mdintegrate,mdatoms,mdinter,general_data,
                      bonded,cp,class_parse,cp_parse,filename_parse,
                      dict_vol,dict_fun[7].keyword);
  set_sim_params_write(mdintegrate,mdatoms,mdinter,general_data,
                       bonded,cp,class_parse,cp_parse,filename_parse,
                       dict_write,dict_fun[8].keyword);
  set_sim_params_pimd(mdintegrate,mdatoms,mdinter,general_data,
                      bonded,cp,class_parse,cp_parse,filename_parse,
                      dict_pimd,dict_fun[9].keyword);


  set_sim_params_finale(mdintegrate,mdatoms,mdinter,general_data,
                        bonded,cp,class_parse,cp_parse,filename_parse);

/*=========================================================================*/
/*   V) Write out the simulation parameters to the simulation file        */

  fp = cfopen(filename_parse->simname,"w");

   write_simfile(fp,dict_list,num_dict_list,dict_fun[1].keyword);
   write_simfile(fp,dict_cp,num_dict_cp,dict_fun[2].keyword);
   write_simfile(fp,dict_gen,num_dict_gen,dict_fun[3].keyword);
   write_simfile(fp,dict_vpot,num_dict_vpot,dict_fun[4].keyword);
   write_simfile(fp,dict_run,num_dict_run,dict_fun[5].keyword);
   write_simfile(fp,dict_nhc,num_dict_nhc,dict_fun[6].keyword);
   write_simfile(fp,dict_vol,num_dict_vol,dict_fun[7].keyword);
   write_simfile(fp,dict_write,num_dict_write,dict_fun[8].keyword);
   write_simfile(fp,dict_pimd,num_dict_pimd,dict_fun[9].keyword);

  fclose(fp);

/*=========================================================================*/
/*   VI) Print out to the screen         */

  PRINTF("\n");
  PRINT_LINE_DASH;
  PRINTF("Completed reading simulation input file %s\n",input_name);
  PRINT_LINE_STAR;
  PRINTF("\n");

/*=========================================================================*/
/*   VII) Free the memory                                                  */

  cfree(fun_key,"control_sim_params");  
  cfree(filename_parse->input_name,"control_sim_params");
  cfree(filename_parse->simname,"control_sim_params");
  cfree(&dict_fun[1],"control_sim_params");
  cfree(&dict_list[1],"control_sim_params");
  cfree(&dict_cp[1],"control_sim_params");
  cfree(&dict_gen[1],"control_sim_params");
  cfree(&dict_vpot[1],"control_sim_params");
  cfree(&dict_run[1],"control_sim_params");
  cfree(&dict_nhc[1],"control_sim_params");
  cfree(&dict_vol[1],"control_sim_params");
  cfree(&dict_write[1],"control_sim_params");
  cfree(&dict_pimd[1],"control_sim_params");

/*========================================================================*/
    }/*end routine*/ 
/*==========================================================================*/







/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_simfile(FILE *fp,DICT_WORD *dict, int num_dict, char *fun_key)

/*==========================================================================*/
/*               Begin subprogram:                                          */
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

   int iuset,itype;

/*========================================================================*/
/*     I) Write out meta key word */

   fprintf(fp,"~%s[\n",fun_key);

/*========================================================================*/  
/*    II)User defined parameters */

   iuset = 1;itype = 1;
      fprintf(fp,"----------------------------------------\n");
      fprintf(fp,"user defined parameters\n");
      fprintf(fp,"----------------------------------------\n");
      dict_print(fp,num_dict,dict,itype,iuset);

/*========================================================================*/  
/*   III)Default parameters */

   iuset = 0;itype = 1;
      fprintf(fp,"----------------------------------------\n");
      fprintf(fp,"Default parameters\n");
      fprintf(fp,"----------------------------------------\n");
      dict_print(fp,num_dict,dict,itype,iuset);

/*========================================================================*/
/*   IV) End Meta key word */

   fprintf(fp,"\n]\n\n");

/*========================================================================*/
/*                   End Subprogram:                                      */
   }/*end routine*/ 
/*========================================================================*/


