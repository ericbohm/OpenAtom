//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                      class_cpatom_maps.h                                 
//                                                                          
//         Class definition for atom maps relating to connectivity          
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#ifndef _CPATOM_MAPS_
#define _CPATOM_MAPS_

class CPATOM_MAPS {
 public:

  int natm_tot;
  int natm_nl;                 /* Num: Number of nonlocal atoms          */
  int natm_gh_nl;              /* Num: Number of nonlocal atoms Gauss-Hermite*/
  int natm_typ_nl;             /* Num: Number of nonlocal atom types     */
  int nab_initio;              /* Num: Number of QM atoms                */

  int *cp_vlnc_true_up;        /* lst: cp atoms cp_valence_true_up > 0 */
  int *cp_vlnc_true_dn;        /* lst: cp atoms cp_valence_true_dn > 0 */
  int *cp_vlnc_up,*cp_vlnc_dn; /* lst: cp atoms cp_valence > 0           */
  int *cp_atm_flag;            /* lst: list of cp atoms :lth natm_tot    */
  int *iatm_atm_typ_nl;        /* Map: Atm ind -> cp nonlocal atm_typ ind;
                                  Ex: 5th cp atm is of 3rd nonlocal atm_typ; 
                                  Lth: number of atoms with nonlocal pseuds*/
  int *iatm_atm_typ_nl_rev;    /* Map: Atm ind -> cp nonlocal atm_typ ind;
                                  Ex: 5th cp atm is of 7th global atm_typ; */
  int *imap_atm_typ_nl;        /* Map: Atm ind -> cp nonlocal atm_typ ind;
                                  Ex: nonlocal type 1 is global type 3     */
  int *imap_atm_typ_nl_gh;     /* Map: Atm ind -> cp gauss-hermite (gh)
                                      nonlocal atm_typ ind;
                                  Ex: gh nonlocal type 3 is gh type 2      */

  int *cp_atm_lst;             /* Map: indices of cp atoms used in gen_wave*/
                               /* also used with dual option               */ 


//----------------------------------------------------------------------------
// Default constructor/destructor
 
  CPATOM_MAPS(){
    natm_tot   =0;
    natm_nl    =0;
    natm_gh_nl =0;
    natm_typ_nl=0;
    nab_initio =0;
  }
  ~CPATOM_MAPS(){}

//----------------------------------------------------------------------------
#ifdef PUP_ON
 //----------------
 //pupping
  void pup(PUP::er &p){
    //pupping ints

   p | natm_tot;       
   p | natm_nl;       
   p | natm_gh_nl;    
   p | natm_typ_nl;   
   p | nab_initio;     

   // Pupping the integer arrays

   if(nab_initio>0){
    pup1d_int(p,&cp_vlnc_up,nab_initio);
    pup1d_int(p,&cp_vlnc_dn,nab_initio);
    pup1d_int(p,&cp_vlnc_true_up,nab_initio);
    pup1d_int(p,&cp_vlnc_true_dn,nab_initio); 
    pup1d_int(p,&cp_atm_lst,nab_initio);
   }//endif
   if(natm_tot>0){
    pup1d_int(p,&cp_atm_flag,natm_tot);          
   }//endif
   if(natm_nl>0){
    pup1d_int(p,&iatm_atm_typ_nl,natm_nl);
    pup1d_int(p,&iatm_atm_typ_nl_rev,natm_nl);
   }//endif
   if(natm_typ_nl>0){
    pup1d_int(p,&imap_atm_typ_nl,natm_typ_nl);
   }//endif
   if(natm_gh_nl>0){
    pup1d_int(p,&imap_atm_typ_nl_gh,natm_gh_nl);
   }//endif
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif        
  }// end PUP
#endif


//============================================================================
// Print out the state of the class

  void state_class_out(){

    int i;
    char fileName [255];
    sprintf (fileName, "%d_cpatom_maps.state", CkMyPe());
    FILE *fp;  fp = fopen(fileName,"w");

   fprintf(fp,"natm_nl %d\n",natm_nl);
   fprintf(fp,"natm_gh_nl %d\n",natm_gh_nl);
   fprintf(fp,"natm_typ_nl %d\n",natm_typ_nl);
   fprintf(fp,"nab_initio %d\n",nab_initio);

    for(i=1;i<=natm_nl;i++){
         fprintf(fp,"mdclatom_maps: iatm_atm_typ_nl[%d] %d\n",
                    i,iatm_atm_typ_nl[i]);}
    for(i=1;i<=natm_nl;i++){
         fprintf(fp,"mdclatom_maps: iatm_atm_typ_nl_rev[%d] %d\n",
                    i,iatm_atm_typ_nl_rev[i]);}
    for(i=1;i<=natm_typ_nl;i++){
         fprintf(fp,"mdclatom_maps: imap_atm_typ_nl[%d] %d\n",
                   i,imap_atm_typ_nl[i]);}
    for(i=1;i<=natm_gh_nl;i++){
         fprintf(fp,"mdclatom_maps: imap_atm_typ_nl_gh[%d] %d\n",
                   i,imap_atm_typ_nl_gh[i]);}
    for(i=1;i<=nab_initio;i++){
         fprintf(fp,"mdclatom_maps: cp_atm_lst[%d] %d\n",
                   i,cp_atm_lst[i]);}

 }// End print out state of class 


 }; // end class

#ifdef PUP_ON
PUPmarshall(CPATOM_MAPS);
#endif

#endif
//==========================================================================
