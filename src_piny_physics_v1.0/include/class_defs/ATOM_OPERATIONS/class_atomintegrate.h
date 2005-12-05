//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
//                          Atom Integration                                    
//==========================================================================

#ifndef _ATOMINTEGRATE_
#define _ATOMINTEGRATE_

class ATOMINTEGRATE{

 //---------------------------------------------------------------------------
 public:

 //---------------------------------------------------------------------------
 //con-destruct:

   ATOMINTEGRATE(){};
  ~ATOMINTEGRATE(){};

 //---------------------------------------------------------------------------
 // Functions

   static void ctrl_atom_integrate(int ,int ,int ,int ,int,int,Atom *,AtomNHC *,
                                   int,double *,double *,double*,int *,int);

   static void integrate_1st_half_step(int ,int ,int ,Atom *,AtomNHC *);
   static void integrate_2nd_half_step(int,int ,int ,int ,Atom *,AtomNHC *,
                                       double *,double *,double*);

   static void integrate_nve_1st_half(int ,Atom *);
   static void integrate_nve_2nd_half(int ,int ,Atom *,double *);
   static void integrate_nvt_1st_half(int ,int ,Atom *,AtomNHC *);
   static void integrate_nvt_2nd_half(int ,int ,int ,Atom *,AtomNHC *,
                                      double *,double *,double *);

   static void computeENHC(int ,int,AtomNHC *,double *,double *);
   static void computeEkin(int , Atom *, double *);

   static void applyNHC      (int ,int ,Atom *,AtomNHC *,double ,int,int);
   static void get_forc_NHC0 (int ,Atom *,AtomNHC *);
   static void get_forc_NHC  (int ,int ,AtomNHC *);
   static void evolve_vNHCM  (int ,int ,AtomNHC *,double );
   static void evolve_vNHC   (int ,int ,AtomNHC *,double ,double );
   static void evolve_vAtmNHC(int ,Atom *,AtomNHC *,double );
   static void evolve_pNHC   (int ,int , AtomNHC *,double );
   static void set_yosh      (int ,double ,double *,double *,double *,double *);

//---------------------------------------------------------------------------
   }; //ATOMINTEGRATE
//==========================================================================

#endif
