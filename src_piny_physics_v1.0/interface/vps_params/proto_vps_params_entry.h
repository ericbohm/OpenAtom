/*---------------------------------------------------------------*/
/* Control_vps_params.c */

void control_vps_params(CPPSEUDO *,GENCELL *,FILENAME_PARSE *,
    SPLINE_PARSE *,int *,int ,NAME [],
    double *,int, int,int ,int , 
    double);

void make_cp_atom_list(CPATOM_MAPS *,CPPSEUDO *,int, int, double *, int *);

void create_non_local_list(CPPSEUDO *,int,int*,int);

void set_ylm_cons(CPYLM_CONS *);

void nlEesSetIter(CPPSEUDO *);

void control_grimme_params(FILENAME_PARSE *,int *, int ,NAME *,
  	                   double *,int ,double *, double *);

/*---------------------------------------------------------------*/





