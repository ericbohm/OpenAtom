/*--------------------------------------------------------------------*/
/* A copy of this is in proto_search_local.h */

void surf_coef(DICT_WORD *, char *, char *, DATA_BASE_ENTRIES *, 
    CATM_LAB *, int );

/*--------------------------------------------------------------------*/

void assign_base_surf(DATA_BASE_ENTRIES *, int , int *, int , 
    double *, double *, int *, double *, double *,
    int *, int ); 

void set_surf_splin(double *, double *, int *, MDSURFACE *, int ); 

void spline_surf(double , double , double *, double *, int , 
    double *, double *, double , double , int , double ); 

void vlj_12_3_surf_bin(int , double *, double *, double *, double , 
    double , double *, double *);  

void vlj_9_3_surf_bin(int , double *, double *, double *, double , 
    double , double *, double *); 

void set_potsurf_dict(DICT_WORD *[], int *, int ); 

void vnull_surf_bin(int ,double *, double *); 






