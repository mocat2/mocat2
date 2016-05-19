/**************************************************************/
/*                                                            */
/*Procedure de test pour savoir si l'on doit scinder un groupe*/
/*par la methode DPC (Density of Points Clustering)           */ 
/*                                                            */
/**************************************************************/

int test_DPC(int nb_individus,int nb_dimensions,groupe_t *groupe_pere,
	     groupe_t *groupe1,groupe_t *groupe2,double **coordonnees_individus,
	     double *cotes_hyperpave,double *valeurs_minimum,double risque,
	     int type_donnees,int type_densite,int donnees_normalisees,
	     int nb_simulations);
