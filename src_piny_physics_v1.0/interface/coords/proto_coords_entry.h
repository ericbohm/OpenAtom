void set_atm_NHC(GENENSOPTS *,GENSTATEPOINT *,
                 GENSIMOPTS *,MDCLATOMS_INFO *,
                 MDCLATOMS_PIMD *,MDGHOST_ATOMS *,
                 MDATOM_MAPS *,MDCONSTRNT *,
                 MDTHERM_INFO *,MDTHERM_INFO *,
                 MDBARO *,MDPAR_RAHMAN *,
                 CLASS_PARSE *,double *,int ,int );

void read_hmat(MDINTEGRATE *, MDATOMS *,MDINTER *,
               GENERAL_DATA *,FILENAME_PARSE *,int ,int ,
               double *,int *);

void check_box_center(GENCELL *,int );

void check_cell(GENCELL *,int ,double ,char *);
