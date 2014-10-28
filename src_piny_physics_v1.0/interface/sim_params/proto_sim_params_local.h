/*----------------------------------------------------------------------*/
/* control_sim_params.c */

void write_simfile(FILE *,DICT_WORD *, int, char *); 


/*----------------------------------------------------------------------*/
/* set_sim_params.c */

void set_sim_params_list(MDINTEGRATE *, MDATOMS *, MDINTER *,
    GENERAL_DATA *,MDINTRA *,CP *,CLASS_PARSE *,
    CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_cp(MDINTEGRATE *, MDATOMS *, MDINTER *,
    GENERAL_DATA *,MDINTRA *,CP *,CLASS_PARSE *,
    CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_gen(MDINTEGRATE *, MDATOMS *, MDINTER *,
    GENERAL_DATA *,MDINTRA *,CP *,CLASS_PARSE *,
    CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_vpot(MDINTEGRATE *, MDATOMS *, MDINTER *,
    GENERAL_DATA *,MDINTRA *,CP *,CLASS_PARSE *,
    CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_run(MDINTEGRATE *, MDATOMS *, MDINTER *,
    GENERAL_DATA *,MDINTRA *,CP *,CLASS_PARSE *,
    CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_nhc(MDINTEGRATE *, MDATOMS *, MDINTER *,
    GENERAL_DATA *,MDINTRA *,CP *,CLASS_PARSE *,
    CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_vol(MDINTEGRATE *, MDATOMS *, MDINTER *,
    GENERAL_DATA *,MDINTRA *,CP *,CLASS_PARSE *,
    CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_write(MDINTEGRATE *, MDATOMS *, MDINTER *,
    GENERAL_DATA *,MDINTRA *,CP *,CLASS_PARSE *,
    CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_pimd(MDINTEGRATE *, MDATOMS *, MDINTER *,
    GENERAL_DATA *,MDINTRA *,CP *,CLASS_PARSE *,
    CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_fun(MDINTEGRATE *, MDATOMS *, MDINTER *,
    GENERAL_DATA *,MDINTRA *,CP *,CLASS_PARSE *,
    CP_PARSE *,FILENAME_PARSE *,DICT_WORD *,char *);
void set_sim_params_finale(MDINTEGRATE *, MDATOMS *, MDINTER *,
    GENERAL_DATA *,MDINTRA *,CP *,CLASS_PARSE *,
    CP_PARSE *,FILENAME_PARSE *);
void set_sim_params_temper(GENTEMPERING_CTRL *,GENSIMOPTS *,DICT_WORD *dict,char *,char *);

/*----------------------------------------------------------------------*/
/* set_sim_dict.c */

void set_sim_dict_temper(int *,DICT_WORD *[]);
void set_sim_dict_list(int *,DICT_WORD *[]);
void set_sim_dict_cp(int *,DICT_WORD *[]);
void set_sim_dict_gen(int *,DICT_WORD *[]);
void set_sim_dict_vpot(int *,DICT_WORD *[]);
void set_sim_dict_run(int *,DICT_WORD *[]);
void set_sim_dict_nhc(int *,DICT_WORD *[]);
void set_sim_dict_vol(int *,DICT_WORD *[]);
void set_sim_dict_write(int *,DICT_WORD *[]);
void set_sim_dict_pimd(int *,DICT_WORD *[]);
void set_sim_dict_fun(int *,DICT_WORD *[]);
void set_sim_dict_velo(int *,DICT_WORD *[]);
void set_sim_dict_msqd(int *,DICT_WORD *[]);
void set_sim_dict_iikt_iso(int *,DICT_WORD *[]);
void set_sim_dict_ickt_iso(int *,DICT_WORD *[]);
void set_sim_dict_rdf(int *,DICT_WORD *[]);

/*----------------------------------------------------------------------*/



