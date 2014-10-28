//==========================================================================
//                 Input-output                                             
//             {Variables needed for mem allocation:                        
//                       file_len                           }               
//                                                                          

#ifndef _GENFILENAMES_
#define _GENFILENAMES_

class GENFILENAMES {

  //----------------
  public:
    int iwrite_screen;          // Num: Freq of screen writes          
    int iwrite_inst;            // Num: Freq of atm-pos config writes  
    int iwrite_dump;            // Num: Freq of atm/PW dump file writes          
    int iwrite_kseigs;          // Num: Freq of KS eigs writes         
    int iwrite_confv;           // Num: Freq of atm-vel conf writes    
    int iwrite_confp;           // Num: Freq of atm-pos conf writes    
    int iwrite_par_confp;       // Num: Freq of atm-pos partial conf writes 
    int iwrite_confc;           // Num: Freq of PW coef conf writes    
    int iwrite_path_cent;       // Num: Freq of centroid conf writes   
    int iwrite_atm_for;         // Num: Freq of atm force              
    int iwrite_elf;             // Num: Elf write freq
    int low_lim_par;            // Num: lower limit of partial conf write 
    int high_lim_par;           // Num: upper limit of partial conf write 
    int ifile_open;             // Opt: file open flag                      
    int iwrite_conf_binary;     // Opt: Write conf files in binary     
    int iwrite_units;           // Opt: Write screen output units      

    char *iname;                // Chr: Instananeous data file name    

    char *dnamei;               // Chr: Atm pos-coord input file name             
    char *dname;                // Chr: Atm pos-coord dump file name             
    char *cpname;               // Chr: Atm pos-coord conf file name         
    char *cvname;               // Chr: Atm vel-coord conf file name         
    char *forcename;            // Chr: Atm force-coord conf file name           
    char *cpparname;            // Chr: Atm partial pos-coord conf file name 
    char *atm_crd_dir_out;      // Chr: Dir for atm output coord files
    char *atm_crd_dir_in;       // Chr: Dir for atm input coord files

    char *dnamec;               // Chr: PW coef dump file name         
    char *ksname;               // Chr: KS eigs file name              
    char *elfname;              // Chr: ELF file name                  
    char *ccname;               // Chr: PW coef conf file name         
    char *centname;             // Chr: Centroid conf file name        

    //----------------
    //con-destruct:
    GENFILENAMES(){
      iwrite_screen    = 0;         
      iwrite_inst      = 0;           
      iwrite_dump      = 0;           
      iwrite_kseigs    = 0;         
      iwrite_confv     = 0;          
      iwrite_confp     = 0;          
      iwrite_par_confp = 0;      
      iwrite_confc     = 0;          
      iwrite_path_cent = 0;      
      iwrite_atm_for   = 0;        
      iwrite_elf       = 0;            
      low_lim_par      = 0;           
      high_lim_par     = 0;          
      ifile_open       = 0;            
      iwrite_conf_binary = 0;    
      iwrite_units     = 0;          
    };
    ~GENFILENAMES(){};

#ifdef PUP_ON
    //----------------
    //pupping
    void pup(PUP::er &p){
      //pupping ints
      p | iwrite_screen;
      p | iwrite_inst;
      p | iwrite_dump;
      p | iwrite_kseigs;
      p | iwrite_confv;
      p | iwrite_confp;
      p | iwrite_par_confp;
      p | iwrite_confc;
      p | iwrite_path_cent;
      p | iwrite_atm_for;
      p | iwrite_elf;
      p | low_lim_par;
      p | high_lim_par;
      p | ifile_open;
      p | iwrite_conf_binary;
      p | iwrite_units;
      //pupping char arrays
      pup1d_char(p,&iname,MAXWORD);
      pup1d_char(p,&dnamei,MAXWORD);
      pup1d_char(p,&dname,MAXWORD);
      pup1d_char(p,&dnamec,MAXWORD);
      pup1d_char(p,&cpname,MAXWORD);
      pup1d_char(p,&ksname,MAXWORD);
      pup1d_char(p,&elfname,MAXWORD);
      pup1d_char(p,&cvname,MAXWORD);
      pup1d_char(p,&ccname,MAXWORD);
      pup1d_char(p,&cpparname,MAXWORD);
      pup1d_char(p,&centname,MAXWORD);
      pup1d_char(p,&forcename,MAXWORD);
      pup1d_char(p,&atm_crd_dir_in,MAXWORD);
      pup1d_char(p,&atm_crd_dir_out,MAXWORD);
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif       
    } // end pup
#endif

    void state_class_out(){

      char fileName [255];
      sprintf (fileName, "%d_genfilenames.state", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");

      fprintf(fp,"iwrite_screen %d\n",iwrite_screen);
      fprintf(fp,"iwrite_inst %d\n",iwrite_inst);
      fprintf(fp,"iwrite_dump %d\n",iwrite_dump);
      fprintf(fp,"iwrite_kseigs %d\n",iwrite_kseigs);
      fprintf(fp,"iwrite_confv %d\n",iwrite_confv);
      fprintf(fp,"iwrite_confp %d\n",iwrite_confp);
      fprintf(fp,"iwrite_par_confp %d\n",iwrite_par_confp);
      fprintf(fp,"iwrite_confc %d\n",iwrite_confc);
      fprintf(fp,"iwrite_path_cent %d\n",iwrite_path_cent);
      fprintf(fp,"iwrite_atm_for %d\n",iwrite_atm_for);
      fprintf(fp,"iwrite_elf %d\n",iwrite_elf);
      fprintf(fp,"low_lim_par %d\n",low_lim_par);
      fprintf(fp,"high_lim_par %d\n",high_lim_par);
      fprintf(fp,"ifile_open %d\n",ifile_open);
      fprintf(fp,"iwrite_conf_binary %d\n",iwrite_conf_binary);
      fprintf(fp,"iwrite_units %d\n",iwrite_units);
      //pupping char arrays
      fprintf(fp,"iname %s\n",iname);
      fprintf(fp,"dnamei %s\n",dnamei);
      fprintf(fp,"dname %s\n",dname);
      fprintf(fp,"dnamec %s\n",dnamec);
      fprintf(fp,"cpname %s\n",cpname);
      fprintf(fp,"ksname %s\n",ksname);
      fprintf(fp,"elfname %s\n",elfname);
      fprintf(fp,"cvname %s\n",cvname);
      fprintf(fp,"ccname %s\n",ccname);
      fprintf(fp,"cpparname %s\n",cpparname);
      fprintf(fp,"centname %s\n",centname);
      fprintf(fp,"forcename %s\n",forcename);

      fclose(fp);
    } // end routine

}; // GENFILENAMES;

#ifdef PUP_ON
PUPmarshall(GENFILENAMES);
#endif


#endif

//==========================================================================
