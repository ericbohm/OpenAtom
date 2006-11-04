//==========================================================================
//                  Simulation options                                      
//             {Variables needed for mem allocation:}                       

#ifndef _GENSIMOPTS_
#define _GENSIMOPTS_

class GENSIMOPTS {

 //----------------
 public:
  int istart;                 // Opt: startup type
  int istart_cp;              // Opt: startup type
  int md;                     // Opt: Classical MD                         
  int pimd;                   // Opt: Path integral MD                     
  int minimize;               // Opt: Classical minimization               
  int cp;                     // Opt: Full CP                              
  int cp_pimd;                // Opt: Full CP                              
  int cp_wave;                // Opt: CP, wave function only               
  int cp_wave_pimd;           // Opt: CP, wave function only               
  int cp_wave_min_pimd;       // Opt: CP min, wave function only           
  int cp_min;                 // Opt: Full CP minimization                 
  int cp_wave_min;            // Opt: CP min, wave function only           
  int cp_any_on;              // Opt: Is cp of any type ``on''
  int debug;                  // Opt: Internal use-backdoor chks           
  int debug_pimd;             // Opt: Internal use-backdoor chks for pimd  
  int debug_cp;               // Opt: Internal use-backdoor chks for CP    
  int debug_cp_pimd;          // Opt: Internal use-backdoor chks for CP    
  int pi_beads;               // Num: # of path integral descritizations   
  int pi_md_typ;              // Opt: Staging or normal modes              
  int initial_spread_opt;     // Opt: Spread coordinates for pimd          
  int anneal_opt;             // Opt: Do simulated annealing (on/off)
  int hess_calc;              // Opt: Calculate the atomic hessian 
  int fftopt;


  double ann_rate;            // Num: Annealing rate                       
  double ann_start_temp;      // Num: Annealing start temperature          
  double ann_target_temp;     // Num: Annealing final temperature          

 //----------------
 //con-destruct:
   GENSIMOPTS(){
     istart           = 0;
     istart_cp        = 0;
     md               = 0;
     pimd             = 0;
     minimize         = 0;
     cp               = 0;
     cp_pimd          = 0;
     cp_wave          = 0;
     cp_wave_pimd     = 0;
     cp_wave_min_pimd = 0;
     cp_min           = 0;
     cp_wave_min      = 0;
     cp_any_on        = 0;
     debug            = 0;
     debug_pimd       = 0;
     debug_cp         = 0;
     debug_cp_pimd    = 0;
     pi_beads         = 0;
     pi_md_typ        = 0;
     initial_spread_opt= 0;
     fftopt           = 0;
     anneal_opt       = 0;
     hess_calc        = 0;
     ann_rate         = 0;
     ann_start_temp   = 0;
     ann_target_temp  = 0;
   };
  ~GENSIMOPTS(){};

#ifdef PUP_ON
 //----------------
 //pupping
  void pup(PUP::er &p){
    //pupping ints
        p | istart;
        p | istart_cp;
        p | md;
        p | pimd;
        p | minimize;
        p | cp;
        p | cp_pimd;
        p | cp_wave;
        p | cp_wave_pimd;
        p | cp_wave_min_pimd;
        p | cp_min;
        p | cp_wave_min;
        p | cp_any_on;
        p | debug;
        p | debug_pimd;
        p | debug_cp;
        p | debug_cp_pimd;
        p | pi_beads;
        p | pi_md_typ;
        p | initial_spread_opt;
        p | anneal_opt;
        p | hess_calc;
	p | fftopt;
    //pupping dbles
        p | ann_rate;
        p | ann_start_temp;
        p | ann_target_temp;
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif        
  } // end pup
#endif

  void state_class_out(){

     char fileName [255];
     sprintf (fileName, "%d_gensimopts.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");

    // ints
     fprintf(fp,"istart %d\n",istart);
     fprintf(fp,"istart_cp %d\n",istart_cp);
     fprintf(fp,"md %d\n",md);
     fprintf(fp,"pimd %d\n",pimd);
     fprintf(fp,"minimize %d\n",minimize);
     fprintf(fp,"cp %d\n",cp);
     fprintf(fp,"cp_pimd %d\n",cp_pimd);
     fprintf(fp,"cp_wave %d\n",cp_wave);
     fprintf(fp,"cp_wave_pimd %d\n",cp_wave_pimd);
     fprintf(fp,"cp_wave_min_pimd %d\n",cp_wave_min_pimd);
     fprintf(fp,"cp_min %d\n",cp_min);
     fprintf(fp,"cp_wave_min %d\n",cp_wave_min);
     fprintf(fp,"cp_any_on %d\n",cp_any_on);
     fprintf(fp,"debug %d\n",debug);
     fprintf(fp,"debug_pimd %d\n",debug_pimd);
     fprintf(fp,"debug_cp %d\n",debug_cp);
     fprintf(fp,"debug_cp_pimd %d\n",debug_cp_pimd);
     fprintf(fp,"pi_beads %d\n",pi_beads);    
     fprintf(fp,"pi_md_typ %d\n",pi_md_typ);
     fprintf(fp,"initial_spread_opt %d\n",initial_spread_opt);
     fprintf(fp,"anneal_opt %d\n",anneal_opt);
     fprintf(fp,"hess_calc %d\n",hess_calc);
    // dbles
     fprintf(fp,"ann_rate %g\n",ann_rate);
     fprintf(fp,"ann_start_temp %g\n",ann_start_temp);
     fprintf(fp,"ann_target_temp %g\n",ann_target_temp);
   fclose(fp);
  }// end routine


}; // GENSIMOPTS;

#ifdef PUP_ON
PUPmarshall(GENSIMOPTS);
#endif

#endif

//==========================================================================
