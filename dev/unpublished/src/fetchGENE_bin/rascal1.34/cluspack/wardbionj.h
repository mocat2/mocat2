/******************************************/
/*                                        */
/*Procedure de classification hierarchique*/
/*                                        */
/******************************************/

void ward_sur_un_arbre_bionj(int nb_sequences,noeud_t **adresses_noeuds,
			     noeud_t *noeud_fictif,int *compteur_indices_dissimilarite,
			     double *historique_indices_dissimilarite,
			     noeud_t **historique_noeuds);

/*********************************************************/
/*                                                       */
/*Procedure d'initialisation des indices de dissimilarite*/
/*                                                       */
/*********************************************************/

void initialisation_indices_dissimilarite(int nb_individus,double **indices_dissimilarite_anciens,
					  noeud_t **adresses_noeuds);

/***************************************************/
/*                                                 */
/*Procedure de mise a jour des indices de Ward     */
/*lors de la classification hierarchique ascendante*/
/*                                                 */
/***************************************************/

void mise_a_jour_indices_Ward(int groupe1,int groupe2,noeud_t **adresses_noeuds,
			      int nb_groupes_courants,int *numeros_groupes_courants,
			      double **indices_Ward_anciens,double **indices_Ward_nouveaux);
