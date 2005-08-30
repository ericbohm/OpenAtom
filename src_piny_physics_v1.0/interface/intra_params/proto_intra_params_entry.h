/*-------------------------------------------------------------------------*/
/* control_intra_params.c */

void control_intra_params(double *,MDCLATOMS_INFO *,MDCLATOMS_PIMD *,
                          CPATOM_MAPS *, 
                          MDATOM_MAPS *,MDINTRA *,
			  FILENAME_PARSE *,
			  FREE_PARSE *,CLASS_PARSE *,NULL_INTER_PARSE *,
                          GENSIMOPTS *, int );

/*-------------------------------------------------------------------------*/
/* close_intra_params.c */

void close_intra_params(MDCLATOMS_INFO *,MDCLATOMS_PIMD *, CPATOM_MAPS *,
                        MDATOM_MAPS *,  BUILD_INTRA *, MDINTRA *, 
                        NULL_INTER_PARSE *, double *, GENSIMOPTS *);

/*-------------------------------------------------------------------------*/
