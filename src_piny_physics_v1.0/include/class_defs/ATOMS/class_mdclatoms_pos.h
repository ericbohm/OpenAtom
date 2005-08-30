//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                                                                          
//                      class_mdclatoms_pos.h                               
//                                                                          
//         Class definition for classical atom positions, velocities,       
//         forces, mode forces, hessian elements.                           
//                                                                          
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================



//==========================================================================
#ifndef _MDCLATOMS_POS_
#define _MDCLATOMS_POS_

class MDCLATOMS_POS {

 public:
  int natm_tot;                /* Num: Total number of atoms */
  double *x,*y,*z;             /* Lst: Atm positions;      Lth: natm_tot */
  double *vx,*vy,*vz;          /* Lst: Atm velocity;       Lth: natm_tot */

//--------------------------------------------------------------------------
// Default constructor/destructor

  MDCLATOMS_POS(){
   natm_tot = 0;
  }
  ~MDCLATOMS_POS(){}

//--------------------------------------------------------------------------
// Pack/Unpack

#ifdef PUP_ON
  void pup(PUP::er &p){

    p | natm_tot;

    // PUP the arrays

    if(natm_tot>0){
     pup1d_dbl(p,&x,natm_tot);
     pup1d_dbl(p,&y,natm_tot);
     pup1d_dbl(p,&z,natm_tot);

     pup1d_dbl(p,&vx,natm_tot);
     pup1d_dbl(p,&vy,natm_tot);
     pup1d_dbl(p,&vz,natm_tot);
    }// endif
#ifdef _PARALLEL_DEBUG_        
    if (p.isUnpacking())
     state_class_out ();
#endif     
  }// end pack/unpack
#endif


//----------------------------------------------------------------------------
// Write out the state of the class

  void state_class_out(){

   int i;
   char fileName [255];
   sprintf (fileName, "%d_mdclatoms_pos.state", CkMyPe());
   FILE *fp;  fp = fopen(fileName,"w");

    //Print the ints
   fprintf(fp,"mdclatoms_pos:  natm_tot %d\n",natm_tot);

   for(i=1;i<=natm_tot;i++){fprintf(fp,"mdclatoms_pos:x[%d] %.9g\n",i,x[i]);}
   for(i=1;i<=natm_tot;i++){fprintf(fp,"mdclatoms_pos:y[%d] %.9g\n",i,y[i]);}
   for(i=1;i<=natm_tot;i++){fprintf(fp,"mdclatoms_pos:z[%d] %.9g\n",i,z[i]);}
   for(i=1;i<=natm_tot;i++){fprintf(fp,"mdclatoms_pos:vx[%d] %.9g\n",i,vx[i]);}
   for(i=1;i<=natm_tot;i++){fprintf(fp,"mdclatoms_pos:vy[%d] %.9g\n",i,vy[i]);}
   for(i=1;i<=natm_tot;i++){fprintf(fp,"mdclatoms_pos:vz[%d] %.9g\n",i,vz[i]);}

   fclose(fp);
  } // end member function 
//-----------------------------------------------------------------------------

 }; // end class definition
//--------------------------------------------------------------------------

#ifdef PUP_ON
PUPmarshall(MDCLATOMS_POS);
#endif

#endif
//==========================================================================
