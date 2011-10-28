//==========================================================================
//    Class used to generate initial guess for wavefunction 
//             
//             


#ifndef _GEN_WAVE_
#define _GEN_WAVE_

class GEN_WAVE{
  public:
   int natm_typ_cp,nsplin,nab_initio;
   int cp_lda, cp_lsda;
   double dg,gmin,gmax;

   int *n_ang;
   int *iatm_state_str;
   int *iatm_state_end;
   int *iatm_state_excess_str;
   int *iatm_state_excess_end;
   int *iatm_state_excess_off;
   int *iatm_atm_typ_cp;

   double *gpsi00;
   double **gpsi_now;
   double ***gpsi0,***gpsi1,***gpsi2,***gpsi3;

   GEN_WAVE();
  ~GEN_WAVE();
   void fill_gw_gpsi(CPATOM_MAPS * ,CPCOEFFS_INFO *,
   		     int ,double ,double ,int *iatm_atm_typ, int natm_typ,
                     NAME *vps_name, int cp_lda_in, int cp_lsda_in, int occupation_file_set);

   void create_coefs(int *k_x,int *k_y,int *k_z,
                     int gSpaceSize,int nstate_in,complex *gspace_coefs,
                     double *xfull, double *yfull, double *zfull, int kpoint_ind);

   void splin_btrans(double *g,double **gpsi0, double **gpsi1,
                     double **gpsi2, double **gpsi3,
                     double *gpsi00, int *n_ang_out, char *fname_ps);

   void bess_trans(double *v_rphi,int nr,double dr,double *r,
                   double *fv_rphi,double *g,int iang,double *gzero);

   void get_gpsi(double g, double **gpsi0,double **gpsi1,
                 double **gpsi2,double **gpsi3,double *gpsi_now,int n_ang);

   void fit_spline(double *c0i,double *c1i,double *c2i,
                   double *c3i,double *xi);

   void get_ylm(double xk,double yk,double zk,double g,
                double *ylmr,double *ylmi,CPYLM_CONS *ylm_cons);

static void read_occupation_numbers(double *occ_up,double *occ_dn,
                                    int nstate_up, int nstate_dn,int cp_lda_tmp,
                                    char *occupation_file,int *uniform_flag);

static void read_kpoints(int , char *kpt_file_name, int istart);

//-------------------------------------------------------------------------
#ifdef PUP_ON
 //pupping
  void pup(PUP::er &p){
    //pupping ints
      p | natm_typ_cp;
      p | nsplin;
      p | nab_initio;
      p | cp_lda;
      p | cp_lsda;
    //pupping dbles
      p | dg;
      p | gmin;
      p | gmax;
    //pupping arrays
      if(natm_typ_cp>0){
        pup1d_int(p,&n_ang,natm_typ_cp);
      }//endif

      if(nab_initio>0){
        pup1d_int(p,&iatm_state_str,nab_initio);
        pup1d_int(p,&iatm_state_end,nab_initio);
        pup1d_int(p,&iatm_state_excess_str,nab_initio);
        pup1d_int(p,&iatm_state_excess_end,nab_initio);
        pup1d_int(p,&iatm_state_excess_off,nab_initio);
        pup1d_int(p,&iatm_atm_typ_cp,nab_initio);
      }//endif

      if(natm_typ_cp>0 && nsplin>0){

        pup1d_dbl(p,&gpsi00,natm_typ_cp);
        pup2d_dbl(p,&gpsi_now,natm_typ_cp,3,"classgen_wave");
        
        double *x0 = (double *)cmalloc(3*natm_typ_cp*nsplin*sizeof(double),"genwave_pup");
        double *x1 = (double *)cmalloc(3*natm_typ_cp*nsplin*sizeof(double),"genwave_pup");
        double *x2 = (double *)cmalloc(3*natm_typ_cp*nsplin*sizeof(double),"genwave_pup");
        double *x3 = (double *)cmalloc(3*natm_typ_cp*nsplin*sizeof(double),"genwave_pup");
        if(p.isPacking()){
          int iii=0;
          for(int i=1;i<=natm_typ_cp;i++){
          for(int j=1;j<=3;j++){
          for(int k=1;k<=nsplin;k++){
            x0[iii] = gpsi0[i][j][k];
            x1[iii] = gpsi1[i][j][k];
            x2[iii] = gpsi2[i][j][k];
            x3[iii] = gpsi3[i][j][k];
            iii++;
	  }}}
	}//endif

        PUParray(p,x0,3*nsplin*natm_typ_cp);
        PUParray(p,x1,3*nsplin*natm_typ_cp);
        PUParray(p,x2,3*nsplin*natm_typ_cp);
        PUParray(p,x3,3*nsplin*natm_typ_cp);

        if(p.isUnpacking()){
          gpsi0 = (double ***) cmalloc(natm_typ_cp*sizeof(double **),"gen_wave_pup")-1;
          gpsi1 = (double ***) cmalloc(natm_typ_cp*sizeof(double **),"gen_wave_pup")-1;
          gpsi2 = (double ***) cmalloc(natm_typ_cp*sizeof(double **),"gen_wave_pup")-1;
          gpsi3 = (double ***) cmalloc(natm_typ_cp*sizeof(double **),"gen_wave_pup")-1;
          for(int i=1; i<= natm_typ_cp; i++){
            gpsi0[i] = (double **) cmalloc(3*sizeof(double *),"gen_wave_pup")-1;
            gpsi1[i] = (double **) cmalloc(3*sizeof(double *),"gen_wave_pup")-1;
            gpsi2[i] = (double **) cmalloc(3*sizeof(double *),"gen_wave_pup")-1;
            gpsi3[i] = (double **) cmalloc(3*sizeof(double *),"gen_wave_pup")-1;
            for(int j=1; j<=3; j++){
              gpsi0[i][j] = (double *) cmalloc(nsplin*sizeof(double ),"gen_wave_pup")-1;
              gpsi1[i][j] = (double *) cmalloc(nsplin*sizeof(double ),"gen_wave_pup")-1;
              gpsi2[i][j] = (double *) cmalloc(nsplin*sizeof(double ),"gen_wave_pup")-1;
              gpsi3[i][j] = (double *) cmalloc(nsplin*sizeof(double ),"gen_wave_pup")-1;
            }//endfor
          }//endfor
          int iii=0;
          for(int i=1;i<=natm_typ_cp;i++){
          for(int j=1;j<=3;j++){
          for(int k=1;k<=nsplin;k++){
            gpsi0[i][j][k] = x0[iii];
            gpsi1[i][j][k] = x1[iii];
            gpsi2[i][j][k] = x2[iii];
            gpsi3[i][j][k] = x3[iii];
            iii++;
	  }}}
	}//endif

        cfree(x0,"genwave_pup");
        cfree(x1,"genwave_pup");
        cfree(x2,"genwave_pup");
        cfree(x3,"genwave_pup");

      }//endif
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif     
  } // end pup
#endif
//----------------------------------------------------------------------
  void state_class_out(){
     int i,j,k;
     char fileName [255];
     sprintf (fileName, "%d_gen_wave.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");

    // ints
     fprintf(fp,"natm_typ_cp %d\n",natm_typ_cp);
     fprintf(fp,"dg %lg\n",dg);
     if(natm_typ_cp>0){
      for(i=1;i<=natm_typ_cp;i++){fprintf(fp,"n_ang[%d] %d\n",i,n_ang[i]);}
      for(i=1;i<=natm_typ_cp;i++){fprintf(fp,"gpsi00[%d] %lg\n",i,gpsi00[i]);}

      for(i=1; i<= natm_typ_cp;i++){
      for(j=1; j<= 3;j++){
	fprintf(fp,"gpsi_now[%d][%d] %lg \n",i,j,gpsi_now[i][j]);
      }}

      for(i=1; i<= natm_typ_cp;i++){
      for(j=1; j<= 3;j++){
      for(k=1; k<= nsplin; k++){
	fprintf(fp,"gpsi0[%d][%d][%d] %lg\n",i,j,k,gpsi0[i][j][k]);
	fprintf(fp,"gpsi1[%d][%d][%d] %lg\n",i,j,k,gpsi1[i][j][k]);
	fprintf(fp,"gpsi2[%d][%d][%d] %lg\n",i,j,k,gpsi2[i][j][k]);
	fprintf(fp,"gpsi3[%d][%d][%d] %lg\n",i,j,k,gpsi3[i][j][k]);
      }}}//endfor

   }//endif
   fclose(fp);
  }// end routine 

//-------------------------------------------------------------------------
  }; //GEN_WAVE
//==========================================================================

#ifdef PUP_ON
PUPmarshall(GEN_WAVE);
#endif

#endif
//==========================================================================
