// No piny structure or class definitions allowed in argument lists here

void setput_nd_eext_corrs(int , int *, int *, int *, double *);

void create_clus_corr(int , int *, int *, int *, double *, double *, 
                      int, double *, double *, double **,double, double, double *);

void create_wire_corr(int , int *, int *, int *, double *, double *, 
                      int, double *, double *,double **, double, double, double *);

void create_surf_corr_dummy(int , int *, int *, int *, double *, double *, 
                            int, double *, double *,double **,double,double,double *);

void create_surf_corr(int , int * ,int *, int *, double *, double *, double *);

void create_perd_corr_data(int ,int , int ,int *, int *, double **, double **,double ***, 
                           double *, double *);

void set_large_kvectors(double , double *, int **, int **,int **, int *,int , int , int);
