/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                      class_mdinteract.h                                  */
/*                                                                          */
/*    Class definition for intermolecular interaction options and details    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#ifndef _MDINTERACT_
#define _MDINTERACT_

class MDINTERACT {
 public:
  int natm_typ;                /* Num: Number of unique atom types */
  int nsplin;                  /* Num: # pts in spline of interact typ*/
  int iswit_vdw;                     /* Opt: Switch inter shift on/off */
  int dielectric_opt;          /* Flag to turn on distant dependent 
                                  dielectric constant                 */
  int ishave_opt;              /* Num: Option to shave skin in force_npol.c */

  int nter_typ;                /* Num: # inter-atm interaction 
                                       types=natm_typ*(natm_typ+1)/2  */
  int nsplin_tot;              /* Num: nsplin*nter_typ                */
  int ninter;
  int ninter_unique;           /* Num: Number of unique interactions */
  int nsplin_tot_unique;       /* Num: nsplin*ninter_unique           */
  int lj_coul_explicit_opt;
  int *inter_map_index;        /* Num: length ninter maps into spline
                                       arrays for unique cv0 cdv0 */
  int **inter_map_mat;         /* Map : matrix form of inter_map_index */
  int **inter_map_mat0;        /* Map : map to non-unique interaction types */
  double dielectric_rheal;     /* Num: Dielectric constant healing length */ 
  double dielectric_cut;       /* Num: Distant Dependent Dielectric cutoff*/ 
  double dielectric_eps;       /* Num: Distant Dependent Dielectric constant*/

  double brnch_root_skin;      /* Num: Extra Skin needed to correctly obtain
                                       brnch-brnch pairs > 
                                       cut_root-2*brnch_root_cut           */
  double skin;                 /* Num: Skin depth used to generate ver_list 
                                       skin=skin_comp+brnch_root_skin      */
  double brnch_root_cut;       /* Num: An atom bonded to only one other atom 
                                       with eq_bond < brnch_root_cut is a 
                                       brnch */
  double rheal_res;            /* Num: RESPA healing length           */ 
  double pten_inter_guess;     /* Num: Estimate of intermolecular pressure */
  double pten_kin_guess;       /* Num: Estimate of kinetic pressure */

  double clong;                /* Num: Long range correction to 
                                       C_6/r^6 potentials             */
  double clong_res;            /* Num: RESPA long range correction to 
                                       C_6/r^6 pots                   */
  double spread_now, spread;   /* Num: Spread of the path integral beads */
  double spread_lnk;
  double brnch_root_skin_lnk;
  double cutoff_max;


  double *cutoff;              /* Lst: Inter-atm interact cutoff dist;
                                 Lth: nter_typ                        */ 
  double *cutoff_res;          /* Lst: Inter-atm interact RESPA 
                                       cutoff distance;
                                  Lth: nter_typ                       */ 
  double *cutti;               /* Lst: inner cutoff         
                                  Lth: nter_typ                       */ 
  double *cutskin;             /* Lst: Inter-atm interaction 
                                       cutoff+skin  distance;  
                                  Lth: nter_typ                       */ 
  double *cutskin_res;         /* Lst: Inter-atm interaction 
                                       RESPA cutoff+skin distance;  
                                  Lth: nter_typ                       */ 
  double *cutskin_root;        /* Lst: Inter-atm interaction
                                       cutoff+skin+branch root skin distance;
                                  Lth: nter_typ                       */
  double *cutskin_root_res;    /* Lst: Inter-atm interaction
                                       RESPA cutoff+skin+branch root skin
                                       distance;
                                  Lth: nter_typ                       */
  double *vcut_coul;           /* Lst: Value of coulomb potential at 
                                  interaction cutoff distance;        
                                  Lth: nter_typ                       */ 
  double *rmin_spl;            /* Lst: Min interact distance splined 
                                       of each atm-interact typ;
                                  Lth: nter_typ                       */ 
  double *dr_spl,*dri_spl;     /* Lst: Spacing of data pts used in 
                                       spline of each atm-interact typ
                                  Lth: nter_typ                       */ 
  double *cv0;                 /* Lst: Spline coef of atm-interact pot
                                           Lth: nsplin_tot_unique     */ 
  double *cdv0;               /* Lst: Spline of  atm-interact force
                                             Lth: nsplin_tot_unique     */ 
  double *cv0_c;              /* Lst: Spline coef of couloumb
                                                       -interact potential
                                                   Lth: nsplin_c      */ 
  double *cdv0_c;             /* Lst: Spline coef of 
                                                     couloumb-interact force;
                                                     Lth: nsplin_c      */
  double *lj_sigma;           /* Lst: sigmas  Lth : natm_typ;       */
  double *lj_epsilon;         /* Lst: sigmas  Lth : natm_typ;       */
  double *cut_atm_typ;        /* Lst: cut/2 over each atom type     */
  double *cut_atm_typ_res;    /* Lst: cut/2 over each atom type     */

//----------------------------------------------------------------------------
// Default constructor/destructor

   MDINTERACT(){
    natm_typ   = 0;
    nsplin     = 0;            
    iswit_vdw  = 0;         
    dielectric_opt = 0;    
    ishave_opt = 0;   
    nter_typ   = 0;         
    nsplin_tot = 0;       
    ninter     = 0;
    ninter_unique = 0;
    nsplin_tot_unique = 0;
    lj_coul_explicit_opt = 0;
   }
  ~MDINTERACT(){}

//----------------------------------------------------------------------------
// Pack/Unpack

#ifdef PUP_ON
  void pup(PUP::er &p){
    // PUP ints
    p | natm_typ;
    p | nsplin;       
    p | iswit_vdw;        
    p | dielectric_opt;   
    p | ishave_opt;       
    p | nter_typ;         
    p | nsplin_tot;       
    p | ninter;
    p | ninter_unique;    
    p | nsplin_tot_unique;
    p | lj_coul_explicit_opt;


    // PUP doubles

    p | dielectric_rheal; 
    p | dielectric_cut;   
    p | dielectric_eps;   
    p | brnch_root_skin;  
    p | skin;             
    p | brnch_root_cut;   
    p | rheal_res;        
    p | pten_inter_guess; 
    p | pten_kin_guess;   
    p | clong;            
    p | clong_res;        
    p | spread_now; 
    p | spread;
    p | spread_lnk;
    p | brnch_root_skin_lnk;
    p | cutoff_max;

   // PUP arrays

  if(lj_coul_explicit_opt == 0){ 
   if(ninter > 0){
    pup1d_int(p,&inter_map_index,ninter); 
   }

   pup2d_int(p,&inter_map_mat,natm_typ,natm_typ); 
   pup2d_int(p,&inter_map_mat0,natm_typ,natm_typ); 

   if(nter_typ > 0){
    pup1d_dbl(p,&cutoff,nter_typ);
    pup1d_dbl(p,&cutoff_res,nter_typ);
    pup1d_dbl(p,&cutti,nter_typ);
    pup1d_dbl(p,&cutskin,nter_typ);
    pup1d_dbl(p,&cutskin_res,nter_typ); 
    pup1d_dbl(p,&cutskin_root,nter_typ);
    pup1d_dbl(p,&cutskin_root_res,nter_typ);
   }
   if(ninter_unique > 0){
    pup1d_dbl(p,&vcut_coul,ninter_unique);
    pup1d_dbl(p,&rmin_spl,ninter_unique);
    pup1d_dbl(p,&dr_spl,ninter_unique);
    pup1d_dbl(p,&dri_spl,ninter_unique); 
   }
   if(nsplin_tot_unique > 0){
    pup1d_dbl(p,&cv0,nsplin_tot_unique);
    pup1d_dbl(p,&cdv0,nsplin_tot_unique);
    pup1d_dbl(p,&cv0_c,nsplin_tot_unique);
    pup1d_dbl(p,&cdv0_c,nsplin_tot_unique);
   }
   }else{    // lj_coul_explicit
    pup1d_dbl(p,&lj_sigma,natm_typ);
    pup1d_dbl(p,&lj_epsilon,natm_typ);
    pup1d_dbl(p,&cut_atm_typ,natm_typ);
    pup1d_dbl(p,&cut_atm_typ_res,natm_typ);
   }//endif
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif
  } /* end pack unpack */
#endif

//----------------------------------------------------------------------------
// Print out state of class

  void state_class_out(){
    
     char fileName [255];
     sprintf (fileName, "%d_mdinteract.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");
    // Print ints
    fprintf(fp,"mdinteract:   natm_typ %d\n", natm_typ);           
    fprintf(fp,"mdinteract:   nsplin %d\n", nsplin);           
    fprintf(fp,"mdinteract:   iswit_vdw %d\n", iswit_vdw);        
    fprintf(fp,"mdinteract:   dielectric_opt %d\n", dielectric_opt);   
    fprintf(fp,"mdinteract:   ishave_opt %d\n", ishave_opt);       
    fprintf(fp,"mdinteract:   nter_typ %d\n", nter_typ);         
    fprintf(fp,"mdinteract:   nsplin_tot %d\n", nsplin_tot);       
    fprintf(fp,"mdinteract:   ninter %d\n",ninter);
    fprintf(fp,"mdinteract:   ninter_unique %d\n",ninter_unique);
    fprintf(fp,"mdinteract:   lj_coul_explicit_opt %d\n",lj_coul_explicit_opt);

    fprintf(fp,"mdinteract:   dielectric_rheal %.12g\n", dielectric_rheal);   
    fprintf(fp,"mdinteract:   dielectric_cut %.12g\n", dielectric_cut);   
    fprintf(fp,"mdinteract:   dielectric_eps %.12g\n", dielectric_eps);   
    fprintf(fp,"mdinteract:   brnch_root_skin %.12g\n", brnch_root_skin);  
    fprintf(fp,"mdinteract:   skin %.12g\n", skin);             
    fprintf(fp,"mdinteract:   brnch_root_cut %.12g\n", brnch_root_cut);   
    fprintf(fp,"mdinteract:   rheal_res %.12g\n", rheal_res);        
    fprintf(fp,"mdinteract:   pten_inter_guess %.12g\n", pten_inter_guess); 
    fprintf(fp,"mdinteract:   pten_kin_guess %.12g\n", pten_kin_guess);   
    fprintf(fp,"mdinteract:   clong %.12g\n", clong);            
    fprintf(fp,"mdinteract:   clong_res %.12g\n", clong_res);        
    fprintf(fp,"mdinteract:   spread_now %.12g\n", spread_now); 
    fprintf(fp,"mdinteract:   spread %.12g\n", spread);
    fprintf(fp,"mdinteract:   spread_lnk %.12g\n", spread_lnk);
    fprintf(fp,"mdinteract:   brnch_root_skin_lnk %.12g\n", brnch_root_skin_lnk);
    fprintf(fp,"mdinteract:   cutoff_max %.12g\n", cutoff_max);

    // Print Arrays

    int i;
  if(lj_coul_explicit_opt == 0){  // SPLINE the potential 
    for(i=1;i<=ninter;i++){
      fprintf(fp,"mdinteract:  inter_map_index[%d] %d\n",i,inter_map_index[i]);}

    for(i=1;i<=nter_typ;i++){
      fprintf(fp,"mdinteract:  cutoff[%d] %.12g\n",i,cutoff[i]);}
    for(i=1;i<=nter_typ;i++){
      fprintf(fp,"mdinteract:  cutoff_res[%d] %.12g\n",i,cutoff_res[i]);}
    for(i=1;i<=nter_typ;i++){
      fprintf(fp,"mdinteract:  cutti[%d] %.12g\n",i,cutti[i]);}
    for(i=1;i<=nter_typ;i++){
      fprintf(fp,"mdinteract:  cutskin[%d] %.12g\n",i,cutskin[i]);}
    for(i=1;i<=nter_typ;i++){
      fprintf(fp,"mdinteract:  cutskin_res[%d] %.12g\n",i,cutskin_res[i]);}
    for(i=1;i<=nter_typ;i++){
      fprintf(fp,"mdinteract:  cutskin_root[%d] %.12g\n",i,cutskin_root[i]);}
    for(i=1;i<=nter_typ;i++){
      fprintf(fp,"mdinteract:  cutskin_root_res[%d] %.12g\n",
                               i,cutskin_root_res[i]);}

    for(i=1;i<=nter_typ;i++){
      fprintf(fp,"mdinteract:  vcut_coul[%d] %.12g\n",i,vcut_coul[i]);}
    for(i=1;i<=nter_typ;i++){
      fprintf(fp,"mdinteract:  rmin_spl[%d] %.12g\n",i,rmin_spl[i]);}
    for(i=1;i<=nter_typ;i++){
      fprintf(fp,"mdinteract:  dr_spl[%d] %.12g\n",i,dr_spl[i]);}
    for(i=1;i<=nter_typ;i++){
      fprintf(fp,"mdinteract:  dri_spl[%d] %.12g\n",i,dri_spl[i]);}
    for(i=1; i<= nsplin_tot_unique; i++){
      fprintf(fp,"mdinteract: cv0[%d] %.12g \n",i,cv0[i]);}
    for(i=1; i<= nsplin_tot_unique; i++){
      fprintf(fp,"mdinteract: cdv0[%d] %.12g\n",i,cdv0[i]);}
    for(i=1; i<= nsplin_tot_unique; i++){
      fprintf(fp,"mdinteract: cv0_c[%d] %.12g\n",i,cv0_c[i]);}
    for(i=1; i<= nsplin_tot_unique; i++){
      fprintf(fp,"mdinteract: cdv0_c[%d] %.12g\n",i,cdv0_c[i]);}
  }else{   // lj_coul_explicit
    for(i=1; i<= natm_typ; i++){
      fprintf(fp,"mdinteract: lj_sigma[%d] %.12g\n",i,lj_sigma[i]);}
    for(i=1; i<= natm_typ; i++){
      fprintf(fp,"mdinteract: lj_epsilon[%d] %.12g\n",i,lj_epsilon[i]);}
    for(i=1; i<= natm_typ; i++){
      fprintf(fp,"mdinteract: cut_atm_typ[%d] %.12g\n",i,cut_atm_typ[i]);}
    for(i=1; i<= natm_typ; i++){
      fprintf(fp,"mdinteract: cut_atm_typ_res[%d] %.12g\n",i,cut_atm_typ_res[i]);}
   }//end if lj_coul_explicit_opt
    fclose(fp);

  }// end member function 

    
}; // end class definition

#ifdef PUP_ON
PUPmarshall(MDINTERACT);
#endif

#endif

