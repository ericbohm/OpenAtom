//==========================================================================
//                  Tempering_ctrl Class
//                  
//==========================================================================

#ifndef _GENTEMPERING_CTRL_
#define _GENTEMPERING_CTRL_

class GENTEMPERING_CTRL {
  //----------------
  public:
    int npara_temps;    /* number of puppies              */
    int nvt,npt,nst;    /* type of tempering              */
    int rsmpl_opt;      /* resample velocities on switch(1/0) */
    int switch_steps;   /* steps between switches         */
    int history_frq;    /* how often to output history    */
    int ipt_restart;    /* restart PT run or from scratch */

    double *t_ext;      /* master temperature list */
    double *p_ext;      /* master pressure list    */
    double *s_ext;      /* master surface tension  */

    char *history_name; /* history file name to be appended        */
    char *wgt_name;     /* frenkel weight file name to be appended */
    char *troyer_name;  /* troyer f(T) file name to be appended    */
    char *output_directory; /* output directory name  */

    //----------------
    //con-destruct:
    GENTEMPERING_CTRL(){
      npara_temps  = 1;
      nvt=npt=nst  = 0;
      rsmpl_opt    = 1; 
      switch_steps = 100;
      history_frq  = 10;
      ipt_restart  = 0;
      t_ext        = NULL;
      p_ext        = NULL;
      s_ext        = NULL;
      history_name = NULL;
      wgt_name     = NULL;
      troyer_name  = NULL;
      output_directory = NULL;
    };
    ~GENTEMPERING_CTRL(){};

#ifdef PUP_ON
    //----------------
    //pupping
    void pup(PUP::er &p){
      //pupping ints
      p | npara_temps;
      p | nvt; p | npt; p | nst;
      p | rsmpl_opt;
      p | switch_steps;
      p | history_frq;
      p | ipt_restart;
      if(npara_temps>1){
        pup1d_dbl(p,&t_ext,npara_temps);
        pup1d_dbl(p,&p_ext,npara_temps);
        pup1d_dbl(p,&s_ext,npara_temps);
        pup1d_char(p,&history_name,MAXWORD);
        pup1d_char(p,&wgt_name,MAXWORD);
        pup1d_char(p,&troyer_name,MAXWORD);
      }/*endif*/
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking()){
        state_class_out ();
      }/*endif*/
#endif        
    } // end pup
#endif

    void state_class_out(){
      char fileName [255];
      sprintf (fileName, "%d_gentempering_ctrl.state", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");
      fprintf(fp,"Glenn is tired and doesn't want to write the output\n");
      fclose(fp);
    }// end routine

}; // GENTEMPERING_CTRL;

#ifdef PUP_ON
PUPmarshall(GENTEMPERING_CTRL);
#endif

#endif
//==========================================================================
