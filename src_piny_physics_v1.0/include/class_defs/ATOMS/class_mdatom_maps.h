//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                      class_mdatom_maps.h                                 
//                                                                          
//         Class definition for atom maps relating to connectivity          
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#ifndef _MDATOM_MAPS_
#define _MDATOM_MAPS_

class MDATOM_MAPS {
  public:
    int natm_tot;                /* Num: Total number of atoms */
    int nmol_typ;                /* Num: # of mol typs                     */
    int nres_typ;                /* Num: # of res typs                     */
    int natm_typ;                /* Num: # of atm typs                     */  

    int nres_typ_max;
    int nres_max;    
    int nres_tot;

    NAME *atm_typ;               /* Lst: Atm types as strgs; Lth: natm_typ */
    NAME *res_typ;               /* Lst: Res types as strgs; Lth: nres_typ */
    NAME *mol_typ;               /* Lst: Mol types as strgs; Lth: nmol_typ */

    int *nmol_jmol_typ;          /* Lst: # of mol of each Mol type  
Lth: nmol_typ                          */
    int *nres_1mol_jmol_typ;     /* Lst: number of residues in 1 molecule
                                    of the jmol_typth  molecule type;         
Lth: nmol_typ  */ 
    int *jatm_jmol_typ_strt;     /* Lst: start index in atm list where jmol_typth
                                    atoms are listed; Lth: nmol_typ */
    int *jres_jmol_typ_strt;     /* Lst: start index in residue lists 
                                    where jmol_typth's residues are listed
Lth: nmol_typ */
    int *ires_typ_jres_jmol_typ; /* Map: jth residue of the jth molecule type
                                    to index of the residue type
Lth: sum(nres_1mol_jmol_typ) see below */
    int *jatm_jres_1mol_jmol_typ_strt;
    /* Lst: start index in atm list where 
       the atoms of the jth res
       of the 1st molecule of the 
       jth molecule type is 
       listed; Lth: sum(nres_1mol_jmol_typ) */
    int *natm_1mol_jmol_typ;     /* Lst: number of atms in 1 molecule
                                    of the jmol_typth  molecule type;         
Lth: nmol_typ           */ 
    int *natm_jres_jmol_typ;     /* Lst: number of atms in the jth residue
                                    of the jmol_typth  molecule type;         
Lth: sum(nres_1mol_jmol_typ) see below   */ 
    int *nfree_1mol_jmol_typ;     /* Lst: number of degrees of freedom
                                     in 1 molecule of the jmol_typth 
                                     molecule type; Lth: nmol_typ           */
    int *nfree_jres_jmol_typ;     /* Lst: number of degrees of freedom
                                     in the jth residue of the jmol_typth 
                                     molecule type; 
Lth: sum(nres_1mol_jmol_typ)           */
    int *iatm_mol_typ;           /* Map: Atm ind -> mol_typ ind;     
Ex: 25th atm is of 27th mol_typ; 
Lth: natm_tot           */
    int *iatm_res_typ;           /* Map: Atm ind -> res_typ ind;     
Ex: 25th atm is in the 27th res_typ; 
Lth: natm_tot           */
    int *iatm_atm_typ;           /* Map: Atm ind -> atm_typ ind;
Ex: 25th atm is of 27th atm_typ; 
Lth: natm_tot           */
    int *iatm_mol_num;           /* Map: Atm ind -> mol ind 
Ex: 25th atm is the 27th molecule of 
mol typ specified by iatm_mol_typ;  
Lth: natm_tot           */
    int *iatm_res_num;           /* Map: Atm ind -> res num
Ex: 25th atm is in the 27th residue of 
mol typ specified by iatm_mol_typ;  
Lth: natm_tot           */
    //============================================================================
    // Default constructor/destructor

    MDATOM_MAPS(){

      natm_tot = 0;                
      nmol_typ = 0;               
      nres_typ = 0;               
      natm_typ = 0;               

      nres_typ_max = 0;
      nres_max = 0;    
      nres_tot = 0;

    }
    ~MDATOM_MAPS(){}

    //=============================================================================
    // Calculate some constants

    void calc_nres_tot_max(){
      nres_tot = 0;  nres_max = 1;
      for(int i=2;i<=nmol_typ;i++){
        nres_tot += MAX(nres_1mol_jmol_typ[i],1);
        nres_max = MAX(nres_1mol_jmol_typ[i],nres_max);
      }/* endfor */
    }/* end member function */

    //=============================================================================
    // Pack/Unpack

#ifdef PUP_ON
    void pup(PUP::er &p){

      // PUP integers

      p | natm_tot;
      p | nmol_typ;   
      p | nres_typ;   
      p | natm_typ;   

      p | nres_typ_max;
      p | nres_tot;
      p | nres_max;    

      // PUP the arrays

      if(natm_typ>0){
        pup1d_name(p,&atm_typ,natm_typ);   
      }

      if(nres_typ>0){
        pup1d_name(p,&res_typ,nres_typ);   
      }

      if(nmol_typ>0){
        pup1d_name(p,&mol_typ,nmol_typ);   
        pup1d_int(p,&nmol_jmol_typ,nmol_typ);
        pup1d_int(p,&nres_1mol_jmol_typ,nmol_typ);
        pup1d_int(p,&jatm_jmol_typ_strt,nmol_typ); 
        pup1d_int(p,&jres_jmol_typ_strt,nmol_typ);  
        pup1d_int(p,&natm_1mol_jmol_typ,nmol_typ);
        pup1d_int(p,&nfree_1mol_jmol_typ,nmol_typ); 
      }

      if(nres_tot>0){
        pup1d_int(p,&ires_typ_jres_jmol_typ,nres_tot); 
        pup1d_int(p,&jatm_jres_1mol_jmol_typ_strt,nres_tot); 
        pup1d_int(p,&natm_jres_jmol_typ,nres_tot);
        pup1d_int(p,&nfree_jres_jmol_typ,nres_tot);
      }

      if(natm_tot>0){
        pup1d_int(p,&iatm_mol_typ,natm_tot);  
        pup1d_int(p,&iatm_res_typ,natm_tot); 
        pup1d_int(p,&iatm_atm_typ,natm_tot);
        pup1d_int(p,&iatm_mol_num,natm_tot);  
        pup1d_int(p,&iatm_res_num,natm_tot); 
      }
#ifdef _PARALLEL_DEBUG_        
      if (p.isUnpacking())
        state_class_out ();
#endif     
    }/* end pack/unpack */
#endif

    //============================================================================
    // Print out the state of the class

    void state_class_out(){

      char fileName [255];
      sprintf (fileName, "%d_mdclatom_maps.state", CkMyPe());
      FILE *fp;  fp = fopen(fileName,"w");

      fprintf(fp,"mdclatom_maps:  natm_tot %d\n",natm_tot);
      fprintf(fp,"mdclatom_maps:  nmol_typ %d\n",nmol_typ);
      fprintf(fp,"mdclatom_maps:  nres_typ %d\n",nres_typ);
      fprintf(fp,"mdclatom_maps:  natm_typ %d\n",natm_typ);
      fprintf(fp,"mdclatom_maps:  nres_typ_max %d\n",nres_typ_max);
      fprintf(fp,"mdclatom_maps:  nres_tot %d\n",nres_tot);
      fprintf(fp,"mdclatom_maps:  nres_max %d\n",nres_max);

      int i;

      for(i = 1; i <= natm_typ;i++){fprintf(fp,"mdclatom_maps: atm_typ[%d] %s \n",i,atm_typ[i]);}
      for(i = 1; i <= nres_typ;i++){fprintf(fp,"mdclatom_maps: res_typ[%d] %s\n",i,res_typ[i]);}
      for(i = 1; i <= nmol_typ;i++){fprintf(fp,"mdclatom_maps: mol_typ[%d] %s\n",i,mol_typ[i]);}

      for(i = 1; i <= nmol_typ;i++){fprintf(fp,"mdclatom_maps: nmol_jmol_typ[%d] %d\n",i,nmol_jmol_typ[i]);}
      for(i = 1; i <= nmol_typ;i++){fprintf(fp,"mdclatom_maps: nres_1mol_jmol_typ[%d] %d\n",i,nres_1mol_jmol_typ[i]);}
      for(i = 1; i <= nres_tot;i++){fprintf(fp,"mdclatom_maps: jatm_jmol_typ_strt[%d] %d\n",i,jatm_jmol_typ_strt[i]);}
      for(i = 1; i <= nres_tot;i++){fprintf(fp,"mdclatom_maps: jres_jmol_typ_strt[%d] %d\n",i,jres_jmol_typ_strt[i]);}

      for(i = 1; i <= nres_tot;i++)
      {fprintf(fp,"mdclatom_maps: ires_typ_jres_jmol_typ[%d] %d\n",i,ires_typ_jres_jmol_typ[i]);}
      for(i = 1; i <=  nres_tot;i++) 
      {fprintf(fp,"mdclatom_maps: jatm_jres_1mol_jmol_typ_strt[%d] %d\n",i,jatm_jres_1mol_jmol_typ_strt[i]);}
      for(i = 1; i <= nmol_typ;i++){fprintf(fp,"mdclatom_maps: natm_1mol_jmol_typ[%d] %d\n",i,natm_1mol_jmol_typ[i]);}
      for(i = 1; i <= nres_tot;i++){fprintf(fp,"mdclatom_maps: natm_jres_jmol_typ[%d] %d\n",i,natm_jres_jmol_typ[i]);}

      for(i = 1; i <= nmol_typ;i++){fprintf(fp,"mdclatom_maps: nfree_1mol_jmol_typ[%d] %d\n",i,nfree_1mol_jmol_typ[i]);}
      for(i = 1; i <= nres_tot;i++){fprintf(fp,"mdclatom_maps: nfree_jres_jmol_typ[%d] %d\n",i,nfree_jres_jmol_typ[i]);}
      for(i = 1; i <= natm_tot;i++){fprintf(fp,"mdclatom_maps: iatm_mol_typ[%d] %d\n",i,iatm_mol_typ[i]);}

      for(i = 1; i <= natm_tot;i++){fprintf(fp,"mdclatom_maps: iatm_res_typ[%d] %d\n",i,iatm_res_typ[i]);}
      for(i = 1; i <= natm_tot;i++){fprintf(fp,"mdclatom_maps: iatm_atm_typ[%d] %d\n",i,iatm_atm_typ[i]);}
      for(i = 1; i <= natm_tot;i++){fprintf(fp,"mdclatom_maps: iatm_mol_num[%d] %d\n",i,iatm_mol_num[i]);}
      for(i = 1; i <= natm_tot;i++){fprintf(fp,"mdclatom_maps: iatm_res_num[%d] %d\n",i,iatm_res_num[i]);}

      fclose(fp);
    }/* end member function */

};/* end class definition */;

#ifdef PUP_ON
PUPmarshall(MDATOM_MAPS);
#endif

#endif
