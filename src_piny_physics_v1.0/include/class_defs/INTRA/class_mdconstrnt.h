//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                            class_mdconstrnt.h                            
//                                                                          
//                Class definition for bonds constraints                    
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#ifndef _MDCONSTRNT_
#define _MDCONSTRNT_

class MDCONSTRNT{
  public:
    int natm_tot;
    int nmol_typ;
    int nfreeze;
    int pimd_freez_typ;         // Opt: type of freezing 
    int iconstrnt;              /* Opt: Atm constraints operative      */
    int max_iter;               /* Num: maximum iterations             */
    int nres_tot;               /* Num : Total number of residues in system */

    double tolshake;            /* Num: Atm shake tolerence            */
    double tolratl;             /* Num: Atm rattle tolerence           */

    int *freeze_flag;           /* Lst: 1/0 if atoms is/is not frozen
Lth: natm_tot                     */
    int *freeze_map;            /* Map: indices of frozen atoms 
Lth: nfreeze  */
    int *atom_label;            /* Num: label of atom for freezing of length
                                   natm_tot.
Options: Standard(0), Backbone(1), Sidechain(2) */ 
    int *icons_jmol_typ;     /* Lst: Flag to indicate if there are constraints
                                in  the jmol_typth molecule type;         
Lth: nmol_typ           */
    int *icons_jres_jmol_typ; /* Lst :  nres_tot*/

    //============================================================================
    // Default constructor/destructor

    MDCONSTRNT(){
      natm_tot       = 0;
      nmol_typ       = 0;
      pimd_freez_typ = 0;
      nfreeze        = 0;
      iconstrnt      = 0;             
      max_iter       = 0;              
      nres_tot       = 0;              
      tolshake       = 1.e-05;
      tolratl        = 1.e-05;
    }
    ~MDCONSTRNT(){}

    //============================================================================
    // Pack/Unpack

#ifdef PUP_ON
    void pup(PUP::er &p){
      p | natm_tot;
      p | nmol_typ;
      p | pimd_freez_typ;
      p | nfreeze;
      p | iconstrnt;
      p | max_iter;
      p | nres_tot;
      p | tolshake;
      p | tolratl;
      if(natm_tot>0){
        pup1d_int(p,&freeze_flag,natm_tot);
        pup1d_int(p,&atom_label,natm_tot);
      }// endif
      if(nfreeze>0){
        pup1d_int(p,&freeze_map,nfreeze);
      }// endif
      if(nmol_typ>0){
        pup1d_int(p,&icons_jmol_typ,nmol_typ);
      }//endif
      if(nres_tot>0){
        pup1d_int(p,&icons_jres_jmol_typ,nres_tot);
      }
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif 
    }// end Pack/Unpack 
#endif

    //============================================================================
    // Print out state of class 

    void state_class_out(){

      int i;
      char fileName [255];
      sprintf (fileName, "%d_mdconstrnt.state", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");

      fprintf(fp,"mdcons: natm_tot %d\n",natm_tot);
      fprintf(fp,"mdcons: nmol_typ %d\n",nmol_typ);
      fprintf(fp,"mdcons: nfreeze %d\n",nfreeze);
      fprintf(fp,"mdcons: iconstrnt %d\n",iconstrnt);
      fprintf(fp,"mdcons: max_iter %d\n",max_iter);
      fprintf(fp,"mdcons: tolshake %g\n",tolshake);
      fprintf(fp,"mdcons: tolratl %g\n",tolratl);

      for(i=1;i<=natm_tot;i++){fprintf(fp,"mdcons: freeze_flag[%d] %d\n",
          i,freeze_flag[i]);}

      for(i=1;i<=natm_tot;i++){fprintf(fp,"mdcons: atom_label[%d] %d\n",
          i,atom_label[i]);}
      for(i=1;i<=nfreeze;i++){fprintf(fp,"mdcons: freeze_map[%d] %d\n",
          i,freeze_map[i]);}
      for(i=1;i<=nmol_typ;i++){fprintf(fp,"mdcons: icons_jmol_typ[%d] %d\n",
          i,icons_jmol_typ[i]);}
      if(iconstrnt == 1){
        for(i=1;i<=nres_tot;i++){fprintf(fp,"mdcons:icons_jres_jmol_typ[%d] %d\n",
            i,icons_jres_jmol_typ[i]);}
      }
      fclose(fp);
    }// end member function 

    //-------------------------------------------------------------------------
}; // end class definition
//==========================================================================

#ifdef PUP_ON
PUPmarshall(MDCONSTRNT);
#endif

#endif
//==========================================================================
