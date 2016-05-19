/*************************************************************/
/*                                                           */
/*Procedure de classification en un nombre de classes donnees*/
/*                                                           */
/*************************************************************/

void mixturemodels(int nb_individus,int nb_dimensions,individu_t **individus,
		   int *nb_clusters,int nb_iterations_max,int nb_clusters_selection);

/*************************************************************/
/*                                                           */
/*Procedure de classification en un nombre de classes donnees*/
/*                                                           */
/*************************************************************/

int kmixturemodels(int nb_individus,int nb_dimensions,individu_t **individus,
		   int nb_clusters,groupe_t *clusters,int nb_iterations_max,
		   double *meilleure_vraisemblance,double densite_loi_uniforme);

/***************************************************************/
/*                                                             */
/*Procedure de calcul des log-vraisemblances de tous les points*/
/*par rapport a tous les clusters                              */
/*                                                             */
/***************************************************************/

void calcul_vraisemblances_points(int nb_individus,int nb_dimensions,individu_t **individus,
				  int nb_clusters,groupe_t *clusters);

/**************************************************/
/*                                                */
/*Procedure de calcul du mu moyen pour initialiser*/
/*la maximisation de la vraisemblance             */
/*                                                */
/**************************************************/

void calcul_mu_moyen(int nb_points,double **coordonnees,int nb_dimensions,
		     double *mu);

