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

    static void ctrl_atom_integrate(int ,int ,int ,int, int, int ,int,int,Atom *,AtomNHC *,
        int,double *,double *,double*,int *,int,int,int,int,int);

    static void integrate_1st_half_step(int ,int ,int ,Atom *,AtomNHC *,int,int,int);
    static void integrate_2nd_half_step(int,int ,int ,int ,Atom *,AtomNHC *,
        double *,double *,double*,int,int,int);

    static void integrate_nve_1st_half(int ,Atom *,int,int,int);
    static void integrate_nve_2nd_half(int ,int ,Atom *,double *,int,int,int);
    static void integrate_nvt_1st_half(int ,int ,Atom *,AtomNHC *,int,int,int);
    static void integrate_nvt_2nd_half(int ,int ,int ,Atom *,AtomNHC *,
        double *,double *,double *,int,int,int);
    static void integrate_isonvt_1st_half(int ,int ,Atom *,AtomNHC *,int,int,int);
    static void integrate_isonvt_2nd_half(int ,int ,int ,Atom *,AtomNHC *,
        double *,double *,double *,int,int,int);

    static void computeENHC(int ,int,AtomNHC *,double *,double *,double *,int,int,int,int);
    static void computeEkin(int , Atom *, double *,int,int,int);

    static void applyNHC      (int ,int ,Atom *,AtomNHC *,double ,int,int,int,int,int);
    static void get_forc_NHC0 (int ,Atom *,AtomNHC *,int,int,int);
    static void get_forc_NHC  (int ,int ,AtomNHC *,int,int,int);
    static void evolve_vNHCM  (int ,int ,AtomNHC *,double ,int,int,int);
    static void evolve_vNHC   (int ,int ,AtomNHC *,double ,double ,int,int,int);
    static void evolve_vAtmNHC(int ,Atom *,AtomNHC *,double ,int,int,int);
    static void evolve_pNHC   (int ,int , AtomNHC *,double ,int,int,int);
    static void set_yosh      (int ,double ,double *,double *,double *,double *);

    static void scaleIso(int ,Atom *,AtomNHC *,int,int,int);
    static void applyIsoVel(int ,Atom *,AtomNHC *,double ,int,int,int);
    static void applyIsoNHC(int ,int ,Atom *,AtomNHC *,double ,int ,int ,int,int,int);
    static void evolve_vAtmIsoNHC(int, Atom *,AtomNHC *,double ,int,int,int);
    static void evolve_pIsoNHC(int, int ,AtomNHC *,double,int,int,int);

    //---------------------------------------------------------------------------
}; //ATOMINTEGRATE
//==========================================================================

#endif
