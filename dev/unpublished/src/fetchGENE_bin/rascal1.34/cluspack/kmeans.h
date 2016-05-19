/****************************************************/
/*                                                  */
/*Procedure de realisation de clustering par k-means*/
/*en partant d'un debut de clustering               */ 
/*                                                  */
/****************************************************/

int kmeans(int nb_individus,int nb_dimensions,individu_t **individus,int nb_classes);

/***************************************************/
/*                                                 */
/*Procedure d'affectation des individus aux groupes*/
/*et de calcul des nouveaux centres de gravite     */
/*                                                 */
/***************************************************/

int affectation_individuskmeans(int nb_dimensions,int nb_individus,int nb_classes,
				int *nb_individus_par_classe,individu_t **individus,
				double **centres,double *inertie_intra_classe);
