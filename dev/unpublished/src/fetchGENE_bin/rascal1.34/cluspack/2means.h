#define GROUPE1 0
#define GROUPE2 1

/******************************************************/
/*                                                    */
/*Procedure de realisation d'un clustering d'un groupe*/ 
/*en deux par la methode des k-means avec k=2         */
/*                                                    */
/******************************************************/

int two_means(int nb_dimensions,groupe_t *groupe,groupe_t *groupe1,groupe_t *groupe2,
	      int type_donnees);

/***************************************************/
/*                                                 */
/*Procedure d'affectation des individus aux groupes*/
/*et de calcul des nouveaux centres de gravite     */
/*                                                 */
/***************************************************/

int affectation_individus2means(int nb_dimensions,int nb_individus,individu_t **individus,
				double *centre1,double *centre2,int *groupes_des_individus,
				double *inertie_intra_classe,int type_donnees,double *centre);
