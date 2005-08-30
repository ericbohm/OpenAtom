//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                      class_mdsurf.h                                      
//                                                                          
//    Class definition for surface potentials                               
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

#ifndef _MDSURFACE_
#define _MDSURFACE_

class MDSURFACE {
 public:
  int    isurf_on;                  // Opt : Surface potential on/off  
  int    natm_typ;                  // Num : Number of atom types      
  int    nsplin_surf;               // Num : spline points             
  int    nsplin_tot;                // Num : natm_typ*nsplin           
  
  double zheal;                     // Num : Healing length            
  double surface_height;            // Num : Height of the surface     

  double *surface_pot;              // Lst : surface pot               
  double *surface_forc;             // Lst : surface force             
                                    //     Lth: nsplin_tot             
  double *zcut_off;                 // Lst : surface pot cutoff        
                                    //     Lth: natm_typ               
  double *dz_spl;                   // Lst : surface pot dz            
                                    //     Lth: natm_typ               
  double *dzi_spl;                  // Lst : surface pot dz            
                                    //     Lth: natm_typ               
  double *zmin_spl;                 // Lst : surface pot zmin          
                                    //     Lth: natm_typ               
  char   *surface_type;             // Type of the surface             

//---------------------------------------------------------------------------
// Default constructor/destructor

   MDSURFACE(){
     isurf_on    = 0;         
     natm_typ    = 0;         
     nsplin_surf = 0;      
     nsplin_tot  = 0;       
     surface_type = (char *)cmalloc(MAXWORD*sizeof(char),"constr:surface");
     strcpy(surface_type,"null");
   }
  ~MDSURFACE(){}

//---------------------------------------------------------------------------
// Pack/Unpack

#ifdef PUP_ON
  void pup(PUP::er &p){

   // PUP ints

   p | isurf_on; 

   if(isurf_on==1){
   
     p | natm_typ; 
     p | nsplin_surf; 
     p | nsplin_tot;  
     // PUP doubles
     p | zheal;     
     p | surface_height;
     // PUP arrays
     pup1d_dbl(p,&zcut_off,natm_typ); 
     pup1d_dbl(p,&dz_spl,natm_typ);  
     pup1d_dbl(p,&dzi_spl,natm_typ);  
     pup1d_dbl(p,&zmin_spl,natm_typ);  
     pup1d_dbl(p,&surface_pot,nsplin_tot);  
     pup1d_dbl(p,&surface_forc,nsplin_tot); 
     pup1d_char(p,&surface_type,MAXWORD);
   }// endif
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif
  } // end pack unpack
#endif

//---------------------------------------------------------------------------
// Print out state of class

  void state_class_out(){
    
     int i;
     char fileName [255];
     sprintf (fileName, "%d_mdsurface.state", CkMyPe());
     FILE *fp;  fp = fopen(fileName,"w");

     fprintf(fp,"mdsurface: isurf_on %d\n",isurf_on);
     if(isurf_on==1){
       fprintf(fp,"mdsurface: natm_typ %d\n",natm_typ);
       fprintf(fp,"mdsurface: nsplin_surf %d\n",nsplin_surf);
       fprintf(fp,"mdsurface: nsplin_tot %d\n",nsplin_tot);

       fprintf(fp,"mdsurface: zheal %.12g\n",zheal);
       fprintf(fp,"mdsurface: surface_height %.12g\n",surface_height);

       for(i=1;i<=natm_typ;i++){fprintf(fp,"mdsurf: zcut_off[%d] %g\n",
                                 i,zcut_off[i]);}
       for(i=1;i<=natm_typ;i++){fprintf(fp,"mdsurf: dz_spl[%d] %g\n",
                                 i,dz_spl[i]);}
       for(i=1;i<=natm_typ;i++){fprintf(fp,"mdsurf: dzi_spl[%d] %g\n",
                                 i,dzi_spl[i]);}
       for(i=1;i<=natm_typ;i++){fprintf(fp,"mdsurf: zmin_spl[%d] %g\n",
                                 i,zmin_spl[i]);}

       fprintf(fp,"mdsurface: surface_type %s\n",surface_type);
     }//endif

     fclose(fp);
  }// end member function

//---------------------------------------------------------------------------
  }; // end class definition
//==========================================================================


#ifdef PUP_ON
PUPmarshall(MDSURFACE);
#endif

#endif
