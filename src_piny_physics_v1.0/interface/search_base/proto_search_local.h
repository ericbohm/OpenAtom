/*---------------------------------------------------------------------*/
/*     search_base.c                                                   */

void sort_atm_list(CATM_LAB *,int ,int *, WILD *,int );

void match_data_base(CATM_LAB *, CATM_LAB *, int , int ,int , int *,
    int *,int *,int *,int *,int ,WILD *,int );

void match_wild_base(CATM_LAB *, CATM_LAB *, int , int , int , int , int *,
    int *,int *,int *,int ,WILD *,int );

void comp_atmlst(int *,CATM_LAB *,WILD *,int ,int ,int );

void match_atmlst(int *,CATM_LAB *,int ,CATM_LAB *,int , WILD *,int );

void check_mult_base(CATM_LAB *,int *,int ,int , int , char *, WILD *);

/*---------------------------------------------------------------------*/
/* defined in various routines AND OTHER PROTOFILES  */

void set_potvps_dict(DICT_WORD *[],int *, int );

void set_grimme_dict(DICT_WORD *[],int *, int );

void set_pot_explicit_lj_dict(DICT_WORD *[],int *, int );

void set_potinter_dict(DICT_WORD *[],int *, int );

void set_potbond_dict(DICT_WORD *[],int *,int);

void set_potbend_dict(DICT_WORD *[],int *,int);

void set_potbend_bnd_dict(DICT_WORD *[],int *,int);

void set_pottors_dict(DICT_WORD *[],int *,int);

void set_potonfo_dict(DICT_WORD *[],int *,int);

void set_potsurf_dict(DICT_WORD *[], int *, int ); 

void inter_coef(DICT_WORD *,char [],char [], DATA_BASE_ENTRIES *,
    CATM_LAB *,int );

void bond_coef(DICT_WORD *,char *,char *,DATA_BASE_ENTRIES *,CATM_LAB *,
    int );

void bend_coef(DICT_WORD *,char *,char *,DATA_BASE_ENTRIES *,CATM_LAB *,
    int );

void tors_coef(DICT_WORD *,char *,char *,DATA_BASE_ENTRIES *,CATM_LAB *,
    int );

void onfo_coef(DICT_WORD *,char *,char *,DATA_BASE_ENTRIES *,CATM_LAB *,
    int );

void bend_bnd_coef(DICT_WORD *,char *,char *,DATA_BASE_ENTRIES *,
    CATM_LAB *, int );


void surf_coef(DICT_WORD *, char *, char *, DATA_BASE_ENTRIES *, 
    CATM_LAB *, int );

/*---------------------------------------------------------------------*/
