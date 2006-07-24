//==========================================================================
//                  CP coeffs and states                                    
//             {Variables needed for mem allocation:                        
//                   ncoef_l,ncoef,(0,ncoef),nstate_up,nstate_dn}         
//                                                                          
//==========================================================================

#ifndef _CPCOEFFS_INFO_
#define _CPCOEFFS_INFO_

class CPCOEFFS_INFO {

 //---------------------------------------------------------------------------
 public:
  int cp_any_on;              // Opt: Is cp of any type ``on''
  int pi_beads;               // Num: # of path integral beads        
  int nstate_up,nstate_dn;    // Num: # spin up/dn states             
  int nstate2;                // Num : nstate*nstate
  int icmass_unif;            // Opt: Uniform coeff PW masses         
  int ncoef;                  // Num: # of PW coef on small 
                              //         spherical cutoff grid dens_cp_box   
  int ncoef_l;                // Num: # of coef on large 
                              //         spherical cutoff grid sparse box    
  int ncoef_l_dens_cp_box;    // Num: # of coef on small dense rho
                              //      spherical cutoff cp_grid (hmat_cp)  
                              //      on dens_cp_box  
  double ecut_psi;            // Num: Energy cutoff small sparse grid    
  double ecut_mass;           // Num: mass preconditioning
  double tau_mass;            // Num: mass tau
  double ecut;                // Num: Energy cutoff large sparse grid    
  double ecut_dens_cp_box;    // Num: Energy cutoff cp_box               
  double *occ_up,*occ_dn;     // Lst: orbital occupation numbers     
                              // Lth: nstate_up                     
 //---------------------------------------------------------------------------
 //con-destruct:
   CPCOEFFS_INFO(){
    cp_any_on       = 0;
    pi_beads        = 0;               
    nstate_up       = 0;
    nstate_dn       = 0;    
    nstate2         = 0;    
    icmass_unif     = 0;            
    ncoef           = 0;
    ncoef_l         = 0;                
    ncoef_l_dens_cp_box = 0;    
   };
  ~CPCOEFFS_INFO(){};

 //---------------------------------------------------------------------------
#ifdef PUP_ON
 //pupping
  void pup(PUP::er &p){
    //pupping ints
      p | cp_any_on;
      p | pi_beads;
      p | nstate_up;
      p | nstate_dn;
      p | nstate2;
      p | icmass_unif;
      p | ncoef;
      p | ncoef_l;
      p | ncoef_l_dens_cp_box;
    //pupping dbles
      p | ecut_psi;
      p | ecut_mass;
      p | tau_mass;
      p | ecut;
      p | ecut_dens_cp_box;

  // PUP Arrays
      if(cp_any_on==1){
       pup1d_dbl(p,&occ_up,nstate_up);
       if(nstate_dn>0){pup1d_dbl(p,&occ_dn,nstate_dn);}
      }//endif
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif           
  } // end pup
#endif

 //---------------------------------------------------------------------------
  void state_class_out(){
    int i;

    char fileName [255];
    sprintf (fileName, "%d_cpcoeffs_info.state", CkMyPe());
    FILE *fp;  fp = fopen(fileName,"w");

    // ints
    fprintf(fp,"cp_any_on %d\n",cp_any_on);
    fprintf(fp,"pi_beads %d\n",pi_beads);
    fprintf(fp,"nstate_up %d\n",nstate_up);
    fprintf(fp,"nstate_dn %d\n",nstate_dn);
    fprintf(fp,"nstate2 %d\n",nstate2);
    fprintf(fp,"icmass_unif %d\n",icmass_unif);
    fprintf(fp,"ncoef %d\n",ncoef);
    fprintf(fp,"ncoef_l %d\n",ncoef_l);
    fprintf(fp,"ncoef_l_dens_cp_box %d\n",ncoef_l_dens_cp_box);
   // dbles
    fprintf(fp,"tau_mass %g\n",tau_mass);
    fprintf(fp,"ecut_mass %g\n",ecut_mass);
    fprintf(fp,"ecut_psi  %g\n",ecut_psi);
    fprintf(fp,"ecut %g\n",ecut);
    fprintf(fp,"ecut_dens_cp_box %g\n",ecut_dens_cp_box);
    if(cp_any_on==1){
      for(i=1;i<=nstate_up;i++) fprintf(fp,"occ_up[%d] %.12g\n",i,occ_up[i]);
      for(i=1;i<=nstate_dn;i++) fprintf(fp,"occ_dn[%d] %.12g\n",i,occ_dn[i]);
    }// endif
   fclose(fp);

  }// end routine
 //---------------------------------------------------------------------------

}; //CPCOEFFS_INFO
//==========================================================================

#ifdef PUP_ON
PUPmarshall(CPCOEFFS_INFO);
#endif

#endif
//==========================================================================
