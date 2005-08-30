void control_set_cp_ewald(GENSIMOPTS *,GENCELL *,CPCOEFFS_INFO *, CPOPTS *,
                          GENEWALD *, CPEWALD *,CP_PARSE *,
                          double *,double *,double *, int ,int ,
                          double *,int ,MDPART_MESH *,MDECOR *, 
                          int ,int ,int);

void get_coul_clus_corr(GENEWALD *,GENCELL *,int ,int ,
                        double ,double ,int ,int ,
                        double *,int ,int);


void get_coul_2D_corr(GENEWALD *,GENCELL *,int ,int ,int,double *,int );

void get_coul_1D_corr(GENEWALD *,GENCELL *,int ,int ,int,double *,int );

void init_cp_dual_pme_bw(GENEWALD *,CPEWALD *,int );

void init_cp_atom_pme_bw(GENEWALD * ,CPEWALD *);

