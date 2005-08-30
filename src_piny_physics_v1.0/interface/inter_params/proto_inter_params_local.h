/*---------------------------------------------------------------*/
/* control_inter_params.c : inter_coef needed in proto_search */

void set_inter_splin(double [],double [],double [],double [],double [],
                     double [],double [],double [],double [],double [],
                     double ,int ,int [],
                     SPLINE_PARSE *,MDINTERACT *,
                     double [],double [],
                     int ,int ,int);


void set_inter_params(DICT_WORD [],NAME ,char *,
                      int *,double *,double *,
            	      double *,double *,double *,double *,
                      double *,double *,double *,double *,double *,double *,
                      double *,double *,double *,double *,double);

void get_clong(int ,int ,int ,
		int [],double [],int [],
                double *,double *,double [],
		double [],int , double );

void get_clong_lj_explicit(int ,int ,int [],
                           double [],double [],double [],double [],
	                   double *,double *,int , double );


void inter_coef(DICT_WORD *,char [],char [], DATA_BASE_ENTRIES *,
                     CATM_LAB *,int );

void assign_base_inter(DATA_BASE_ENTRIES *,int,int *,
                       int, double *, double *,double *,double *,
                       double *,double *,double *, double *,double *,double *,
                       int *,double *,double *,double *,int *,int,CATM_LAB *,
                       CATM_LAB *);


void sort_inter_params(double *,double *,double *,double *,double *,
                       double *,double *,double *,double *,double *,
                       double *,double *,double *,
                       int *,int  *,int *,int );

void switchij(double *,double *);
void iswitchij(int *,int *);

/*---------------------------------------------------------------*/
/* Set_inter_dict.c : copy needed in proto_serach*/

void set_potinter_dict(DICT_WORD *[],int *, int );

void set_pot_explicit_lj_dict(DICT_WORD *[],int *, int );

void set_explicit_lj_params(DICT_WORD [],
                            char *,char *,double *,double *,double *,double *);

/*---------------------------------------------------------------*/
void spline_vdv(double , double , 
		double *, double *, 
		int , double *, double *,double , double ,
		double , double , double , double , 
                double , double , double , double , 
		double , double , int ,int , int ,
                double, int, double, double,double);

double glj(double , double , double , double *);

double gwill(double , double , double , double , double , double , 
	     double *);

double gcoul(double , double , double *, int , double );

double gaziz(double , double , double , double , double , double , 
              double , double , double , double *);


void vnull_bin(int , double [], double [],double []);

void vlj_bin(int , double [], double [],double [], double [],
             double [],double , double );


void vwill_bin(int ,double [],double [],double [],
               double [],double [],double ,double ,
               double ,double ,double );


void vwill_lj_bin(int ,double [],double [],double [],
               double [],double [],double ,double ,
               double ,double ,double ,double ,double );


void vaziz_bin(int ,double [],double [],double [],
               double [],double [],double ,double ,double ,
               double ,double ,double ,double ,double );


void vcoul_bin(int ,double [],double [],double [], double [],double [],
               double,int, double, int, double, double,double);














