//==========================================================================
//                  CP simulation options                                   
//             {Variables needed for mem allocation:}                       
//                                                                          

#ifndef _CPOPTS_
#define _CPOPTS_

class CPOPTS{

 //----------------
 public:
  int cp_any_on;              // Opt: Is cp of any type ``on''
  int nstate;                 // Num : number of states = MAX(up,dn)
  int nstate2;                // Num : nstate*nstate
  int cp_lda,cp_lsda;         // Opt: Lda-lsda opt                   
  int cp_sic;                 // Opt: Sic opt                           
  int cp_gga;                 // Opt: Gga opt                           
  int cp_nonint;              // Opt: Non-interacting elec opt       
  int cp_norb;                // Opt: Norb integration               
  int cp_ptens_calc;          // Opt: Calculate CP pressure tensor   
  int cp_init_orthog;         // Opt: Initially orthogonalize WFs    
  int cp_gs,cp_low,cp_normalize;
                              // Opt: Orthogonalization option for 
                              //         minimization                   
  int cp_ngrid_skip;          // Num: How many grid points to skip
                              //      when writing out grid-based quantities 
  int cp_hess_calc;           // Opt: Calculate the diagonal Hessian
                              //      for minimization               
  int zero_cp_vel;            // Opt: Zero all coefficient velocities 
  int iwrite_coef_binary;     // binary write option for coef file 
  int iread_coef_binary;      // binary read option for coef file 
  int icheck_perd_size;       // Opt: Check atm dists under CBCs 
  int icheck_dual_size;       // Opt: Check atm dists under CBCs 

  int cp_becke;               // cp GGA flags 
  int cp_pw91x;
  int cp_fila_1x; 
  int cp_fila_2x;
  int cp_pbe_x;
  int cp_brx89;
  int cp_brx2k;
  int cp_lyp;  
  int cp_lypm1;  
  int cp_pw91c;
  int cp_pbe_c;
  int cp_tau1_c;
  int cp_debug_xc;
  int cp_dual_grid_opt;        //Opt: on = 1 off = 0 
  int cp_isok_opt;             //Opt: Use Gaussian isokinetic for CP coeffs 

  int ks_rot_on;              // Opt: Perioidic rotation to 
                              //         KS states opt                   
  int n_ks_rot;               // Num: Freq of rotation to
                              //         KS states                        
  int cp_ke_dens_on;          // Opt: Calculate the electron kinetic
                              //         energy density                  
  int cp_tau_functional;      // Opt: A tau-dependent functional is 
                              //         being used in this calculation  
  int cp_elf_calc_frq;        // Num: Compute the electron localization
                              //         function evey ... steps         
  int cp_laplacian_on;        // Num: Compute the Laplacian of the density 


  double te_ext;              // Num: Temperature of PW coef          
  double tol_edge_dist;       //  Num: Tolerance on max atom dist under CBC
  double tol_coef;            // Num: cp and cp_wave fcoef tolerance
  double cp_hess_cut,cp_hess_tau;// Num: Hessian cutoff and Hessian tau 
  double pseud_hess_loc;      // Num: Local contribution to Hessian
  double cp_tau_nhc;

 //----------------
 //con-destruct:
   CPOPTS(){
      cp_any_on       = 0;
      nstate          = 0;               
      nstate2         = 0;              
      cp_lda          = 0;
      cp_lsda         = 0;       
      cp_sic          = 0;               
      cp_gga          = 0;               
      cp_nonint       = 0;            
      cp_norb         = 0;              
      cp_ptens_calc   = 0;        
      cp_init_orthog  = 0;       
      cp_gs           = 0;
      cp_low          = 0;
      cp_normalize    = 0;
      cp_ngrid_skip   = 0;        
      cp_hess_calc    = 0;         
      zero_cp_vel     = 0;          
      iwrite_coef_binary = 0;   
      iread_coef_binary  = 0;    
      icheck_perd_size   = 0;     
      icheck_dual_size   = 0;     
      cp_becke        = 0;     
      cp_pw91x        = 0;
      cp_fila_1x      = 0; 
      cp_fila_2x      = 0;
      cp_pbe_x        = 0;
      cp_brx89        = 0;
      cp_brx2k        = 0;
      cp_lyp          = 0;  
      cp_lypm1        = 0;  
      cp_pw91c        = 0;
      cp_pbe_c        = 0;
      cp_tau1_c       = 0;
      cp_debug_xc     = 0;
      cp_dual_grid_opt= 0; 
      cp_isok_opt     = 0;  
      ks_rot_on       = 0;              
      n_ks_rot        = 0;               
      cp_ke_dens_on   = 0;          
      cp_tau_functional = 0;      
      cp_elf_calc_frq = 0;        
      cp_laplacian_on = 0;        
   };
  ~CPOPTS(){};
  
#ifdef PUP_ON
 //----------------
 //pupping :
  void pup(PUP::er &p){
    //pupping ints
      p | cp_any_on;
      p | nstate;
      p | nstate2;
      p | cp_lda;
      p | cp_lsda;
      p | cp_sic;
      p | cp_gga;
      p | cp_nonint;
      p | cp_norb;
      p | cp_ptens_calc;
      p | cp_init_orthog;
      p | cp_gs;
      p | cp_low;
      p | cp_normalize;
      p | cp_ngrid_skip;
      p | cp_hess_calc;
      p | zero_cp_vel;
      p | iwrite_coef_binary;
      p | iread_coef_binary;
      p | icheck_perd_size;
      p | icheck_dual_size;
      p | cp_becke;
      p | cp_pw91x;
      p | cp_fila_1x; 
      p | cp_fila_2x;
      p | cp_pbe_x;
      p | cp_brx89;
      p | cp_brx2k;
      p | cp_lyp;  
      p | cp_lypm1;  
      p | cp_pw91c;
      p | cp_pbe_c;
      p | cp_tau1_c;
      p | cp_debug_xc;
      p | cp_dual_grid_opt;
      p | cp_isok_opt;
      p | ks_rot_on;
      p | n_ks_rot;
      p | cp_ke_dens_on;
      p | cp_tau_functional;
      p | cp_elf_calc_frq;
      p | cp_laplacian_on;
    //pupping dbles
      p | te_ext;
      p | tol_edge_dist;
      p | tol_coef;
      p | cp_hess_cut;
      p | cp_hess_tau;
      p | pseud_hess_loc;
      p | cp_tau_nhc;
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif           
  } // end pup
#endif

  void state_class_out(){
     char fileName [255];
     sprintf (fileName, "%d_cpopts.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");
   // ints
     fprintf(fp,"cp_any_on %d\n",cp_any_on);
     fprintf(fp,"nstate %d\n",nstate);
     fprintf(fp,"nstate2 %d\n",nstate2);
     fprintf(fp,"cp_lda %d\n",cp_lda);
     fprintf(fp,"cp_lsda %d\n",cp_lsda);
     fprintf(fp,"cp_sic %d\n",cp_sic);
     fprintf(fp,"cp_gga %d\n",cp_gga);
     fprintf(fp,"cp_nonint %d\n",cp_nonint);
     fprintf(fp,"cp_norb %d\n",cp_norb);
     fprintf(fp,"cp_ptens_calc %d\n",cp_ptens_calc);
     fprintf(fp,"cp_init_orthog %d\n",cp_init_orthog);
     fprintf(fp,"cp_gs %d\n",cp_gs);
     fprintf(fp,"cp_low %d\n",cp_low);
     fprintf(fp,"cp_normalize %d\n",cp_normalize);
     fprintf(fp,"cp_ngrid_skip %d\n",cp_ngrid_skip);
     fprintf(fp,"cp_hess_calc %d\n",cp_hess_calc);
     fprintf(fp,"zero_cp_vel %d\n",zero_cp_vel);
     fprintf(fp,"iwrite_coef_binary %d\n",iwrite_coef_binary);
     fprintf(fp,"iread_coef_binary %d\n",iread_coef_binary);
     fprintf(fp,"icheck_perd_size %d\n",icheck_perd_size);
     fprintf(fp,"icheck_dual_size %d\n",icheck_dual_size);
     fprintf(fp,"cp_becke %d\n",cp_becke);
     fprintf(fp,"cp_pw91x %d\n",cp_pw91x);
     fprintf(fp,"cp_fila_1x %d\n",cp_fila_1x); 
     fprintf(fp,"cp_fila_2x %d\n",cp_fila_2x);
     fprintf(fp,"cp_pbe_x %d\n",cp_pbe_x);
     fprintf(fp,"cp_brx89 %d\n",cp_brx89);
     fprintf(fp,"cp_brx2k %d\n",cp_brx2k);
     fprintf(fp,"cp_lyp %d\n",cp_lyp);  
     fprintf(fp,"cp_lypm1 %d\n",cp_lypm1);  
     fprintf(fp,"cp_pw91c %d\n",cp_pw91c);
     fprintf(fp,"cp_pbe_c %d\n",cp_pbe_c);
     fprintf(fp,"cp_tau1_c %d\n",cp_tau1_c);
     fprintf(fp,"cp_debug_xc %d\n",cp_debug_xc);
     fprintf(fp,"cp_dual_grid_opt %d\n",cp_dual_grid_opt);
     fprintf(fp,"cp_isok_opt %d\n",cp_isok_opt);
     fprintf(fp,"ks_rot_on %d\n",ks_rot_on);
     fprintf(fp,"n_ks_rot %d\n",n_ks_rot);
     fprintf(fp,"cp_ke_dens_on %d\n",cp_ke_dens_on);
     fprintf(fp,"cp_tau_functional %d\n",cp_tau_functional);
     fprintf(fp,"cp_elf_calc_frq %d\n",cp_elf_calc_frq);
     fprintf(fp,"cp_laplacian_on %d\n",cp_laplacian_on);
    //dbles
     fprintf(fp,"te_ext %g\n",te_ext);
     fprintf(fp,"tol_edge_dist %g\n",tol_edge_dist);
     fprintf(fp,"tol_coef %g\n",tol_coef);
     fprintf(fp,"cp_hess_cut %g\n",cp_hess_cut);
     fprintf(fp,"cp_hess_tau %g\n",cp_hess_tau);
     fprintf(fp,"pseud_hess_loc %g\n",pseud_hess_loc);
     fprintf(fp,"cp_tau_nhc %g\n",cp_tau_nhc);
    fclose(fp);
    
  }// end routine 


}; // CPOPTS

#ifdef PUP_ON
PUPmarshall(CPOPTS);
#endif


#endif

//==========================================================================
