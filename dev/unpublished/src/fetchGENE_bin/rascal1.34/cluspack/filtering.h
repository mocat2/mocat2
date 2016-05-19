/********************************************/
/*                                          */
/*Procedure de standardisation des individus*/
/*                                          */
/********************************************/

void standardisation_individus(int nb_individus,int nb_dimensions,individu_t **individus);

/**************************************************************/
/*                                                            */
/*Procedure de calcul des similarites a partir des coordonnees*/
/*en utilisant les correlations                               */
/*                                                            */
/**************************************************************/

void coordonnees_to_similarites(int nb_individus,int nb_dimensions,individu_t **individus);

/*********************************************/
/*                                           */
/*Procedure de filtrage des individus suivant*/ 
/*une valeur seuil et le type des donnees    */
/*                                           */
/*********************************************/

void filtre_individus(int *nb_individus,int nb_dimensions,individu_t **individus,
		      int type_donnees,double distance_min,int nb_clusters_selection);
