/********************************************************************************/
/*                                                                              */
/*Procedure de decoupage d'arbre avec le test anova sur les densites de presence*/
/*                                                                              */
/********************************************************************************/

void decoupage_arbre_DPC(int nb_individus,int nb_dimensions,individu_t **individus,
			 noeud_t *racine,int *nb_clusters,int type_donnees,
			 int nb_clusters_selection,int donnees_normalisees,int type_densite,
			 int nb_simulations);
