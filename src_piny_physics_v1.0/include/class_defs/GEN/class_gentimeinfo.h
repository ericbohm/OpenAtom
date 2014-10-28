//==========================================================================
//                  Timestep info                                           
//             {Variables needed for mem allocation:}                       

#ifndef _GENTIMEINFO_
#define _GENTIMEINFO_

class GENTIMEINFO{

  //----------------
  public:
    int ntime,itime;            // Num: Time-minimization steps           
    int ix_respa;               // Opt: Extended system respa opt         
    int int_res_ter,int_res_tra;// Opt: Intra-inter respa opts            
    int nres_ter,nres_tra;      // Num: # intra-inter respa steps          
    int int_res_tor,nres_tor;   // Opt: Torsion respa opts                
    int int_res_pimd,nres_pimd; // Opt: Path integral respa options       
    int iget_pe_real_inter_freq;//Opt: Freq calc of realspace inter PE     
    int exit_flag;              // Num: Used to tell code to exit under   
    //      annealing and auto exit opt
    double dt;                  // Num: Time step                         

    //----------------
    //con-destruct:
    GENTIMEINFO(){
      ntime        = 0;
      itime        = 0;
      ix_respa     = 0;
      int_res_ter  = 0;
      int_res_tra  = 0;
      nres_ter     = 0;
      nres_tra     = 0;
      int_res_tor  = 0;
      nres_tor     = 0;
      int_res_pimd = 0;
      nres_pimd    = 0;
      exit_flag    = 0;
      dt           = 0;
      iget_pe_real_inter_freq = 0;
    };
    ~GENTIMEINFO(){};

#ifdef PUP_ON
    //----------------
    //pupping
    void pup(PUP::er &p){
      //pupping ints
      p | ntime;
      p | itime;
      p | ix_respa;
      p | int_res_ter;
      p | int_res_tra;
      p | nres_ter;
      p | nres_tra;
      p | int_res_tor;
      p | nres_tor;
      p | int_res_pimd;
      p | nres_pimd;
      p | iget_pe_real_inter_freq;
      p | exit_flag;
      //pupping dbles
      p |  dt;
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif      
    } // end pup
#endif

    void state_class_out(){

      char fileName [255];
      sprintf (fileName, "%d_gentimeinfo.state", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");

      fprintf(fp,"ntime %d\n",ntime);
      fprintf(fp,"itime %d\n",itime);
      fprintf(fp,"ix_respa %d\n",ix_respa);
      fprintf(fp,"int_res_ter %d\n",int_res_ter);
      fprintf(fp,"int_res_tra %d\n",int_res_tra);
      fprintf(fp,"nres_ter %d\n",nres_ter);
      fprintf(fp,"nres_tra %d\n",nres_tra);
      fprintf(fp,"int_res_tor %d\n",int_res_tor);
      fprintf(fp,"nres_tor %d\n",nres_tor);
      fprintf(fp,"int_res_pimd %d\n",int_res_pimd);
      fprintf(fp,"nres_pimd %d\n",nres_pimd);
      fprintf(fp,"iget_pe_real_inter_freq %d\n",iget_pe_real_inter_freq);
      fprintf(fp,"exit_flag %d\n",exit_flag);
      //dbles
      fprintf(fp," dt %g\n", dt);
      fclose(fp);
    }// end routine

}; // GENTIMEINFO;

#ifdef PUP_ON
PUPmarshall(GENTIMEINFO);
#endif

#endif

//==========================================================================
