/********************************************************************/
/*                                                                  */
/*Procedure de recherche de coupe minimum normalisee de Xing et Karp*/
/*                                                                  */
/********************************************************************/

int normalized_cut(groupe_t *groupe,groupe_t *groupe1,groupe_t *groupe2);


/********************************************************/
/*                                                      */
/*Procedure renvoyant pour un vecteur donne la valeur   */
/*de la fonction objectif associee a la coupe normalisee*/
/*                                                      */
/********************************************************/

double evaluation_coupe_normalisee_minimum(int nb_dimensions,double *y,double *D,double **DminusW);
