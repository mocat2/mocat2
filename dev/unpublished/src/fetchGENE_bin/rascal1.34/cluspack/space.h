/*************************************************************/
/*                                                           */
/*Procedure de calcul des exposants d'Holder                 */
/*Methode : on forme un hyperpave avec pour longueur         */
/*la distances entre les deux coordonnees                    */
/*les plus eloignees sur chaque axe                          */
/*                                                           */
/*************************************************************/

void calcul_exposants_Holder_points(int nb_dimensions,int nb_individus_total,
				    double **coordonnees_tous_individus,int nb_individus,
				    double **coordonnees_individus,double *exposants_Holder,
				    double *cotes_hyperpave,double *valeurs_minimum,
				    int mode_verbose,int test);

