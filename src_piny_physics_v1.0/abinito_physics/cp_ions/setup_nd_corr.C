//=============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=============================================================================
//                 Ewald correction sums
//=============================================================================
/** \file name cp_rspace_ion.C
 ** \brief The physics routines that setup correction to 
    the g-space atom ewald energy for wires, surfaces and clusters
 */
//==========================================================================



//=============================================================================
//=============================================================================
/** \brief Setup correction to g-space part of atomic Ewald sum for surface:  
      Add more gc for small ga and gb when periodicity = 1
*/
//=============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=============================================================================
void setup_2d_ewd_corner_corr(GENEWALD *genewald, GENCELL *gencell){
//=============================================================================
// Unpack 

  int iperd      = genewald->iperd;
  int nfix_para  = genewald->nfix_para;
  int nfix_perp  = genewald->nfix_perp;
  int nka_max    = genewald->nka_max;
  int nkb_max    = genewald->nkb_max;
  int nkc_max    = genewald->nkc_max;

  double alp_ewd = gencell->alp_ewd;
  double vol     = gencell->vol;
  double *hmati  = gencell->hmati;
  double *hmat   = gencell->hmat;

//=============================================================================
// Output

  PRINTF("\n");
  PRINT_LINE_STAR;
  PRINTF("Setting up the surface ewald sum corrections\n");
  PRINT_LINE_DASH; PRINTF("\n");
 
//=============================================================================
// Check the parameters

  if(nfix_para>nka_max || nfix_para>nkb_max){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("Surface periodicity fix value %d too large > %d %d %d\n",
	    nfix_para,nka_max,nkb_max,nkc_max);
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     FFLUSH(stdout); EXIT(1);
  }//endif

  if(iperd!=2){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("Surface periodicity fix invoked for perd=%d\n",iperd);
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     FFLUSH(stdout); EXIT(1);
  }//endif

//=============================================================================
// Count the k-vectors : 

   int kcmax   = nkc_max;
   int nkc_fix = nfix_perp*kcmax;

   int ic = 0;
   for(int ka=0;ka<=nfix_para;ka++){
     int kbmin = (ka==0 ? 0 : -nfix_para);
     for(int kb=kbmin;kb<=nfix_para;kb++){
       if(ka!=0 && kb!=0){
         for(int kc=-nkc_fix;kc<=-(kcmax+1);kc++){ic++;}
       }//endif
       for(int kc=(kcmax+1);kc<=nkc_fix;kc++){ic++;}
     }//endfor
   }//endfor

//=============================================================================
// Malloc them 

   int nktot_fix       = ic;
   int *ka_fix         = (int *)   cmalloc(nktot_fix*sizeof(int))-1;
   int *kb_fix         = (int *)   cmalloc(nktot_fix*sizeof(int))-1;
   int *kc_fix         = (int *)   cmalloc(nktot_fix*sizeof(int))-1;
   double *kernel_fix  = (double *)cmalloc(nktot_fix*sizeof(double))-1;

   PRINTF("You are using nfix_perp %d nfix_para %d and nktot_fix %d\n",
           nfix_perp, nfix_para, nktot_fix);
   PRINTF("to correct the reciprocal space ewald sum for surfaces\n");

//=============================================================================
// Store the k-vectors

   double tpivol = 2.0*M_PI/vol;
   double falp2  = 4.0*alp_ewd;

   ic = 0;
   for(int ka=0;ka<=nfix_para;ka++){
     int kbmin = (ka==0 ? 0 : -nfix_para);
     for(int kb=kbmin;kb<=nfix_para;kb++){
       if(ka!=0 && kb!=0){
         for(int kc=-nkc_fix;kc<=-(kcmax+1);kc++){
           ic++;
           double aka    = (double)ka;
           double akb    = (double)kb;
           double akc    = (double)kc;
           double xk     = (aka*hmati[1]+akb*hmati[2])*tpi;
           double yk     = (aka*hmati[4]+akb*hmati[5])*tpi;
           double zk     =  akc*hmati[9]*tpi;
           double g2     = xk*xk+yk*yk+zk*zk;
           double gs     = sqrt(xk*xk+yk*yk);
           double pre_g2 = exp(-g2/falp2)
           double phase  = cos(M_PI*akc)
           double pre_gs = exp(-gs*hmat[9]);

           ka_fix[ic]     = ka;
           kb_fix[ic]     = kb;
           kc_fix[ic]     = kc;
           kernel_fix[ic] = tpivol*(pre_g2-phase*pre_gs);
         }//endfor
       }//endif
       for(int kc=(kcmax+1);kc<=nkc_fix;kc++){
          ic++;
          double aka    = (double)ka;
          double akb    = (double)kb;
          double akc    = (double)kc;
          double xk     = (aka*hmati[1]+akb*hmati[2])*tpi;
          double yk     = (aka*hmati[4]+akb*hmati[5])*tpi;
          double zk     =  akc*hmati[9]*tpi;
          double g2     = xk*xk+yk*yk+zk*zk;
          double gs     = sqrt(xk*xk+yk*yk);
          double pre_g2 = exp(-g2/falp2)
          double phase  = cos(M_PI*akc)
          double pre_gs = exp(-gs*hmat[9]);

          ka_fix[ic]     = ka;
          kb_fix[ic]     = kb;
          kc_fix[ic]     = kc;
          kernel_fix[ic] = tpivol*(pre_g2-phase*pre_gs);
       }//endfor
     }//endfor
   }//endfor

//=============================================================================
// Pack up

   genewald->nktot_fix     = nktot_fix;
   genewald->nktot_mid_fix = nktot_fix;
   genewald->ka_fix        = ka_fix;
   genewald->kb_fix        = kb_fix;
   genewald->kc_fix        = kc_fix;
   genewald->kernel_fix    = kernel_fix;

//=============================================================================
// Output 

  PRINTF("\n");
  PRINT_LINE_DASH;
  PRINTF("Completed the surface ewald correction setup\n");
  PRINT_LINE_STAR;
  PRINTF("\n");

//-----------------------------------------------------------------------------
  }//end routine
//=============================================================================


//=============================================================================
/** \brief Setup correction to the g-space part of the atomic Ewald for wires:  
      Add more kc for small ka and kb 
      Add more kb for small ka and kc 
*/
//=============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=============================================================================
void setup_1d_ewd_corner_corr(GENEWALD *genewald, GENCELL *gencell){
//=============================================================================
// Unpack

  int iperd      = genewald->iperd;
  int nfix_para  = genewald->nfix_para;
  int nfix_perp  = genewald->nfix_perp;
  int nka_max    = genewald->nka_max;
  int nkb_max    = genewald->nkb_max;
  int nkc_max    = genewald->nkc_max;

  double alp_ewd = gencell->alp_ewd;
  double vol     = gencell->vol;
  double *hmati  = gencell->hmati;
  double *hmat   = gencell->hmat;

//=============================================================================
// Output

  PRINTF("\n");
  PRINT_LINE_STAR;
  PRINTF("Setting up the wire ewald sum corrections\n");
  PRINT_LINE_DASH; PRINTF("\n");
 
//=============================================================================
// Check the parameters

   if(nfix_para>nka_max || nfix_para>nkb_max || nfix_para>nkc_max){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("Wire periodicity fix value %d too large > %d %d %d\n",
	    nfix_para,nka_max,nkb_max,nkc_max);
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     FFLUSH(stdout); EXIT(1);
   }//endif

   if(iperd!=1){
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     PRINTF("Wire periodicity fix invoked for perd=%d \n",iperd);
     PRINTF("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     FFLUSH(stdout); EXIT(1);
   }//endif

//=============================================================================
// Count the k-vectors

   kbmax   = nkb_max;
   nkb_fix = nfix_perp*kbmax;

   kcmax   = nkc_max;
   nkc_fix = nfix_perp*kcmax;

   int ic = 0;
   for(int ka=0;ka<=nfix_para;ka++){
     int kbmin = (ka==0 ? 0 : -nfix_para);
     for(int kb=kbmin;kb<=nfix_para;kb++){
       if(ka!=0 && kb!=0){
         for(int kc=-nkc_fix;kc<=-(kcmax+1);kc++){ic++;}
       }//endif
       for(int kc=(kcmax+1);kc<=nkc_fix;kc++){ic++;}
     }//endfor
   }//endfor

   nktot_mid_fix = ic;

   for(int ka=0;ka<=nka_fix;ka++){
     for(int kc=-nkc_fix;kc<=nkc_fix;kc++){
       if(ka!=0){
         for(int kb=-nkb_fix;kb<=-(kbmax+1);kb++){ic++;}
       }//endif
       for(int kb=(kbmax+1);kb<=nkb_fix;kb++){ic++;}
     }//endfor
   }//endfor

//=============================================================================
// Malloc them 

   nktot_fix   = ic;
   ka_fix      = (int *)   cmalloc(nktot_fix*sizeof(int))-1;
   kb_fix      = (int *)   cmalloc(nktot_fix*sizeof(int))-1;
   kc_fix      = (int *)   cmalloc(nktot_fix*sizeof(int))-1;
   kernel_fix  = (double *)cmalloc(nktot_fix*sizeof(double))-1;

   PRINTF("You are using nfix_perp %d nfix_para %d and nktot_fix %d\n",
           nfix_perp, nfix_para, nktot_fix);
   PRINTF("to correct the reciprocal space ewald sum for wires\n");

//=============================================================================
// Store the k-vectors

   ic = 0;
   for(int ka=0;ka<=nfix_para;ka++){
     int kbmin = (ka==0 ? 0 : -nfix_para);
     for(int kb=kbmin;kb<=nfix_para;kb++){
       if(ka!=0 && kb!=0){
         for(int kc=-nkc_fix;kc<=-(kcmax+1);kc++){
           ic++;
           ka_fix[ic]   = ka;
           kb_fix[ic]   = kb;
           kc_fix[ic]   = kc;
           kernel_fix[ic] = ;
         }//endfor
       }//endif
       for(int kc=(kcmax+1);kc<=nkc_fix;kc++){ic++;}
     }//endfor
   }//endfor

   for(int ka=0;ka<=nka_fix;ka++){
     for(int kc=-nkc_fix;kc<=nkc_fix;kc++){
       if(ka!=0){
         for(int kb=-nkb_fix;kb<=-(kbmax+1);kb++){
           ic++;
           ka_fix[ic]   = ka;
           kb_fix[ic]   = kb;
           kc_fix[ic]   = kc;
           kernel_fix[ic] = ;
         }//endif
       }//endif
       for(int kb=(kbmax+1);kb<=nkb_fix;kb++){
         ic++;
         ka_fix[ic]   = ka;
         kb_fix[ic]   = kb;
         kc_fix[ic]   = kc;
         kernel_fix[ic] = ;
       }//endfor
     }//endfor
   }//endfor

//=============================================================================
// Pack up

   genewald->nktot_fix     = nktot_fix;
   genewald->nktot_mid_fix = nktot_mid_fix;
   genewald->ka_fix        = ka_fix;
   genewald->kb_fix        = kb_fix;
   genewald->kc_fix        = kc_fix;
   genewald->kernel_fix    = kernel_fix;

//=============================================================================
// Output

  PRINTF("\n");
  PRINT_LINE_DASH;
  PRINTF("Completed wire ewald correction setup\n");
  PRINT_LINE_STAR;
  PRINTF("\n");

//-----------------------------------------------------------------------------
  }//end routine
//=============================================================================
