//==========================================================================
//                Pseudopotential interaction info                          
//             {Variables needed for mem allocation:                        
//                 vxctyplen,nsplin_g,nsplin_g_tot,n_ang_max,               
//                 num_nl_lst       }                                       
//                                                                          
//==========================================================================

#ifndef _CPPSEUDO_
#define _CPPSEUDO_

#include "../class_defs/CP/class_psnonlocal.h"

//==========================================================================
class CPPSEUDO{

//-------------------------------------------------------------------------
 public:
  int cp_any_on;              // Opt: Is cp of any type ``on''
  int ees_nonloc_on;          // Opt : nonlocal ees method
  int ees_eext_on;            // Opt : eext ees method
  int cp_ptens_calc;          // Opt : pressure tensor on
  int natm_typ;               // Num : num atom types
  int natm_tot;               // Num : num atom types
  int n_ang_max;              // Num: Max # of angular momentum 
                              //      channels in an e-atm pseudopot 
  int n_ang_max1;             // Num: Max # of angular momentum +1
  int n_ang_max_gh;           // Num: Max # of angular momentum 
                              //      channels in an e-atm pseudopot 
                              //       only Gauss-Hermite type
  int n_ang_max_kb;           // Num: Max # of angular momentum 
                              //         channels in an e-atm pseudopot 
                              //         KB and Goedecker type
  int n_rad_max;              // Num: Max # radial channels in any
                              //                        e-atm pseudopot 
  int nsplin_g;               // Num: # of g-space spline pts in any 
                              //         e-atm interaction channel      
  int nsplin_g_tot;           // Num: Total # of g-space spline pts  
  int num_nl_lst;             // Num: Total # of atms involved in the
                              //         pseudopot angular momentum 
                              //         channels =sum np_nl(j)         
  int nl_cut_on;              // Opt: Non-local cutoff scheme opt    
  int n_interp_ps;            //Num: Order of interpolation for pseudo stuff
  int n_interp_pme_dual;      //Num: Order of interpolation
  int ngrid_eext_a;           //Num : External eext stuff
  int ngrid_eext_b;
  int ngrid_eext_c;
  int nka_eext;
  int nkb_eext;
  int nkc_eext;
  int natm_eext_max;
  int natm_typ_nl;            //Num: number of nonlocal atom types  
                              //    this is a subset of iatm_atm_typ
  int natm_typ_gh;            //Num: number of nonlocal atom types  
  int ngh;                    // Num: Number of gauss-hermite points for
                              //        each nonlocal atom type 
  int ngh_tot;
  int norm_size;
  int nlist;
  int np_loc_cp_box;          // Num: # of atms in small cp box        
  int np_nonloc_cp_box;       // Num: # of atms in small cp box        
  int np_nonloc_cp_box_kb;    // Num: # of kleinman-bylander nonlocal atoms
  int np_nonloc_cp_box_gh;    // Num: # of gauss-hermite nonlocal atoms 
  int natm_nonloc;              // Num: s-wave non-local atoms KB type 

  double fft_size_scale_ps;   // Num: fft size scale factor for pseudo ees method
  double gmin_true;           // Num: Mag of smallest g-vector       
  double gmin_spl;            // Num: Min mag of g-vec in spline     
  double gmax_spl;            // Num: Max mag of g-vec in spline     
  double dg_spl;              // Num: Spacing between g-vectors pts
                              //      in spline of e-atm pseudos     
  double gga_cut;             // Num: Gradient cutoff value          
  double nlvps_skin;          // Num: Nonlocal pseudopotential list
                              //      skin length                     
  double alpha_conv_dual;     // Num: convergence factor for long     
                              //  and short range break up for dual grid
                                

  int *n_ang;                 // Lst: # of angular momentum channels 
                              //         in each e-atm pseudopot;
                              //    Lth: natm_typ                        
  int *nrad_0,*nrad_1,*nrad_2,*nrad_3;// Lst: # rad channels in each
                              //         l channel in each e-atm pseudopot;
                              //    Lth: natm_typ                       
  int *nrad_max_l;            // Lst: max # rad channels in each
                              //      l channel;
                              //    Lth: natm_typ                       
  int *loc_opt;               // Lst: Angular momentum channel 
                              //      chosen as local; Lth: natm_typ  
  int *ivps_label;            // Lst: Type label of e-atm pseudopots;
                              //    Lth: natm_typ                        
  int *np_nl;                 // Lst: # of atms in each  angular 
                              //       momentum  channel except GAUSS-HERMITE
                              //    Lth:  (n_ang_max+1)                  
  int *np_nl_gh;              // Lst: # of gauss-hermite atms in each 
                              //         angular momentum  channel;
                              //   Lth:  (n_ang_max+1)
  int **np_nl_rad_str;        // Lst : where each rad channel strs    
  int **np_nl_rad_end;        // Lst : where each rad channel ends    
  int *ip_nl;                 // Lst: index of atms involved in each 
                              //         angular momentum channel;
                              //    Lth: num_nl_lst<natm_tot*(n_ang_max+1)
  int *ip_nl_rev;             // Lst: index of atms involved in each 
                              //         angular momentum channel;
                              //    Lth: num_nl_lst<natm_tot*(n_ang_max+1)
  int *ip_nl_gh;              // Lst: index of atms involved in each 
                              //         angular momentum channel
                              //         GAUSS-HERMITE
                              //    Lth: num_nl_lst<natm_tot*(n_ang_max+1)
  int *ip_nl_rev_gh;          // Lst: index of atms involved in each 
                              //         angular momentum channel
                              //         GAUSS-HERMITE
                              //    Lth: num_nl_lst<natm_tot*(n_ang_max+1)
  int *iatm_nonloc;           // Lst: index of atms with s-KB type nonlocality

  int *map_nl;                // Lst: order atoms in decreasing order
                              //         of # of radial channels          

  int *natm_eext;             // Lst : number of atoms of each type
  int **map_eext;             // Lst : map of atoms of each type

  int *ip_loc_cp_box;         // Lst: index of atms on small cp box   
                              // Lth: natm_tot                        


  double *vps0,*vps1,*vps2,*vps3;// Lst: Spline coef of pseudopot 
                                 //   Lth: nsplin_g_tot                  
  double *dvps0,*dvps1,*dvps2,*dvps3;// Lst: Spline coef of pseudopot
                                    //   Lth: nsplin_g_tot                 
  double *gzvps;              // Lst: g=0 term of local part of 
                              //         e-atm pseudopot; Lth: natm_typ  
  double *gzvps0;             // Lst: g=0 term of l=0 channel of 
                              //    e-atm pseudopot; Lth: natm_typ       
  double *vpsnorm;            // Lst: Norm of KB type non-local 
                              //    pseudopot; Lth: natm_typ*(n_ang_max+1)  
  double *rcut_nl;            // Lst: Non-local cutoff distance;
                              //    Lth: (natm_typ)                      
  double *rgh,*wgh;           // Lst: Gauss-Hermite nodes and weights 
                              //    Lth: ngh  
  double *q_pseud;            // charge associated with pseudo  
                              // Lth: natm_typ 
  double *q_typ;              // charge associated with the atoms
                              // Lth: natm_typ 
  char *vxc_typ;              // Chr: Exchange-correlation type 
                              //    Lth: MAXWORD                        
  char *ggax_typ;             // Chr: GGA-Exchange-correlation type 
                              //    Lth: MAXWORD                        
  char *ggac_typ;             // Chr: GGA-Exchange-correlation type 
                              //   Lth: MAXWORD                        
  PSNONLOCAL nonlocal;
//-------------------------------------------------------------------------
//con-destruct:
   CPPSEUDO(){
     cp_any_on           = 0;
     ees_nonloc_on       = 0;
     ees_eext_on         = 0;
     cp_ptens_calc       = 0;
     natm_typ            = 0;
     natm_tot            = 0;
     n_ang_max           = 0;         
     n_ang_max1          = 0;        
     n_ang_max_gh        = 0;      
     n_ang_max_kb        = 0;      
     n_rad_max           = 0;         
     nsplin_g            = 0;          
     nsplin_g_tot        = 0;      
     num_nl_lst          = 0;        
     nl_cut_on           = 0;         
     n_interp_ps         = 0; 
     ngrid_eext_a        = 0;
     ngrid_eext_b        = 0;
     ngrid_eext_c        = 0;
     nka_eext            = 0;         
     nkb_eext            = 0;
     nkc_eext            = 0;
     n_interp_pme_dual   = 0; 
     natm_typ_nl         = 0;       
     natm_typ_gh         = 0;       
     ngh                 = 0;               
     ngh_tot             = 0;
     norm_size           = 0;
     nlist               = 0;
     np_loc_cp_box       = 0;     
     np_nonloc_cp_box    = 0;     
     np_nonloc_cp_box_kb = 0;  
     np_nonloc_cp_box_gh = 0;  
     natm_nonloc         = 0;
     fft_size_scale_ps   = 1.0;
   };
  ~CPPSEUDO(){};

//-------------------------------------------------------------------------
#ifdef PUP_ON
 //pupping
  void pup(PUP::er &p){
    //pupping ints
      p | cp_any_on;      
      p | ees_nonloc_on;
      p | ees_eext_on;
      p | cp_ptens_calc;      
      p | natm_typ;
      p | natm_tot;
      p | n_ang_max;
      p | n_ang_max1;
      p | n_ang_max_gh;
      p | n_ang_max_kb;
      p | n_rad_max;
      p | nsplin_g;
      p | nsplin_g_tot;
      p | num_nl_lst;
      p | nl_cut_on;
      p | n_interp_ps;
      p | ngrid_eext_a;
      p | ngrid_eext_b;
      p | ngrid_eext_c;
      p | nka_eext;
      p | nkb_eext;
      p | nkc_eext;
      p | n_interp_pme_dual;
      p | natm_typ_nl;
      p | natm_typ_gh;
      p | ngh;
      p | np_loc_cp_box;
      p | np_nonloc_cp_box;
      p | np_nonloc_cp_box_kb;
      p | np_nonloc_cp_box_gh;
      p | ngh_tot;
      p | norm_size;
      p | nlist;
      p | natm_nonloc;
      p | natm_eext_max;
    //pupping dbles
      p | fft_size_scale_ps;
      p | gmin_true;
      p | gmin_spl;
      p | gmax_spl;
      p | dg_spl;
      p | gga_cut;
      p | nlvps_skin;
      p | alpha_conv_dual;
    
    //pupping arrays
      if(cp_any_on==1){

        if(natm_typ>0){
          pup1d_int(p,&n_ang,natm_typ);
          pup1d_int(p,&nrad_0,natm_typ);
          pup1d_int(p,&nrad_1,natm_typ);
          pup1d_int(p,&nrad_2,natm_typ);
          pup1d_int(p,&nrad_3,natm_typ);
          pup1d_int(p,&nrad_max_l,natm_typ);
          pup1d_int(p,&loc_opt,natm_typ);
          pup1d_int(p,&ivps_label,natm_typ);
          pup1d_int(p,&natm_eext,(natm_typ+1));
          pup2d_int(p,&map_eext,natm_typ,natm_eext_max);

          pup1d_dbl(p,&q_pseud,natm_typ);
          pup1d_dbl(p,&q_typ,natm_typ);
          pup1d_dbl(p,&gzvps,natm_typ);
          pup1d_dbl(p,&gzvps0,natm_typ);
          pup1d_dbl(p,&rcut_nl,natm_typ);

  int **map_eext;             // Lst : map of atoms of each type

        }//endif
        if(norm_size>0){
          pup1d_dbl(p,&vpsnorm,norm_size);
        }//endif
        if(n_ang_max1>0){
          pup1d_int(p,&np_nl,n_ang_max1);
          pup1d_int(p,&np_nl_gh,n_ang_max1);
	}// endif
        if(nlist>0){
          pup1d_int(p,&ip_nl,nlist);
          pup1d_int(p,&ip_nl_rev,nlist);
          pup1d_int(p,&ip_nl_gh,nlist);
          pup1d_int(p,&ip_nl_rev_gh,nlist);
	}//endif
        if(natm_nonloc>0){
          pup1d_int(p,&iatm_nonloc,natm_nonloc);
        }
        if(natm_tot>0){
          pup1d_int(p,&map_nl,natm_tot);
          pup1d_int(p,&ip_loc_cp_box,natm_tot);
	}//endif
        if((n_ang_max1>0)&&(n_rad_max>0)){
          pup2d_int(p,&np_nl_rad_str,n_ang_max1,n_rad_max);
          pup2d_int(p,&np_nl_rad_end,n_ang_max1,n_rad_max);
        }// endif
        if(nsplin_g_tot>0){
          pup1d_dbl(p,&vps0,nsplin_g_tot);
          pup1d_dbl(p,&vps1,nsplin_g_tot);
          pup1d_dbl(p,&vps2,nsplin_g_tot);
          pup1d_dbl(p,&vps3,nsplin_g_tot);
          if(cp_ptens_calc==1){
            pup1d_dbl(p,&dvps0,nsplin_g_tot);
            pup1d_dbl(p,&dvps1,nsplin_g_tot);
            pup1d_dbl(p,&dvps2,nsplin_g_tot);
            pup1d_dbl(p,&dvps3,nsplin_g_tot);
          }//endif
        }//endif

        if(ngh>0){
          pup1d_dbl(p,&rgh,ngh);
          pup1d_dbl(p,&wgh,ngh_tot);
        }//

        pup1d_char(p,&vxc_typ,MAXWORD);
        pup1d_char(p,&ggax_typ,MAXWORD);
        pup1d_char(p,&ggac_typ,MAXWORD);
        nonlocal.pup(p);
      }// endif : cp is on
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())state_class_out ();
#endif     
  } // end pup routine
#endif

//----------------------------------------------------------------------
  void state_class_out(){
     int i;
     char fileName [255];
     sprintf (fileName, "%d_cppseudo.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");

    // ints
     fprintf(fp,"cp_any_on %d\n",cp_any_on);
     fprintf(fp,"cp_ptens_calc %d\n",cp_ptens_calc);
     fprintf(fp,"natm_tot %d\n",natm_tot);
     fprintf(fp,"natm_typ %d\n",natm_typ);
     fprintf(fp,"n_ang_max %d\n",n_ang_max);
     fprintf(fp,"n_ang_max1 %d\n",n_ang_max1);
     fprintf(fp,"n_ang_max_gh %d\n",n_ang_max_gh);
     fprintf(fp,"n_ang_max_kb %d\n",n_ang_max_kb);
     fprintf(fp,"n_rad_max %d\n",n_rad_max);
     fprintf(fp,"nsplin_g %d\n",nsplin_g);
     fprintf(fp,"nsplin_g_tot %d\n",nsplin_g_tot);
     fprintf(fp,"num_nl_lst %d\n",num_nl_lst);
     fprintf(fp,"nl_cut_on %d\n",nl_cut_on);
     fprintf(fp,"n_interp_pme_dual %d\n",n_interp_pme_dual);
     fprintf(fp,"natm_typ_nl %d\n",natm_typ_nl);
     fprintf(fp,"natm_typ_gh %d\n",natm_typ_gh);
     fprintf(fp,"ngh %d\n",ngh);
     fprintf(fp,"ngh_tot %d\n",ngh_tot);
     fprintf(fp,"norm_size %d\n",norm_size);
     fprintf(fp,"nlist %d\n",nlist);
     fprintf(fp,"np_loc_cp_box %d\n",np_loc_cp_box);
     fprintf(fp,"np_nonloc_cp_box %d\n",np_nonloc_cp_box);
     fprintf(fp,"np_nonloc_cp_box_kb %d\n",np_nonloc_cp_box_kb);
     fprintf(fp,"np_nonloc_cp_box_gh %d\n",np_nonloc_cp_box_gh);
     if(cp_any_on==1){
      for(i=1;i<=natm_typ;i++){fprintf(fp,"n_ang[%d] %d\n",i,n_ang[i]);}
      for(i=1;i<=natm_typ;i++){fprintf(fp,"nrad_0[%d] %d\n",i,nrad_0[i]);}
      for(i=1;i<=natm_typ;i++){fprintf(fp,"nrad_1[%d] %d\n",i,nrad_1[i]);}
      for(i=1;i<=natm_typ;i++){fprintf(fp,"nrad_2[%d] %d\n",i,nrad_2[i]);}
      for(i=1;i<=natm_typ;i++){fprintf(fp,"nrad_3[%d] %d\n",i,nrad_3[i]);}
      for(i=1;i<=natm_typ;i++){fprintf(fp,"nrad_max_l[%d] %d\n",
                                   i,nrad_max_l[i]);}
      for(i=1;i<=natm_typ;i++){fprintf(fp,"loc_opt[%d] %d\n",i,loc_opt[i]);}
      for(i=1;i<=natm_typ;i++){fprintf(fp,"ivps_label[%d] %d\n",
                                   i,ivps_label[i]);}
      for(i=1;i<=n_ang_max1;i++){fprintf(fp,"np_nl[%d] %d\n",i,np_nl[i]);}
      for(i=1;i<=n_ang_max1;i++){fprintf(fp,"np_nl_gh[%d] %d\n",
                                   i,np_nl_gh[i]);}
      for(i=1;i<=natm_typ;i++){fprintf(fp,"q_pseud[%d] %g\n",i,q_pseud[i]);}
      fprintf(fp,"vxc_typ %s\n",vxc_typ);
      fprintf(fp,"ggax_typ %s\n",ggax_typ);
      fprintf(fp,"ggac_typ %s\n",ggac_typ);
     }//endif
   fclose(fp);

  }// end routine 

//-------------------------------------------------------------------------
  }; //CPPSEUDO
//==========================================================================

#ifdef PUP_ON
PUPmarshall(CPPSEUDO);
#endif

#endif
//==========================================================================
