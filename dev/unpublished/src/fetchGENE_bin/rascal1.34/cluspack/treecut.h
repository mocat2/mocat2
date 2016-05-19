/*******************************************************************************/
/*                                                                             */
/*Procedure de decouverte du nombre de clusters par Secator a partir d'un arbre*/
/*                                                                             */
/*******************************************************************************/

void secator_tree(int nb_individus,noeud_t *racine,double *seuil_dissimilarite,int *nb_clusters);

/********************************************************/
/*                                                      */
/*Procedure de decoupage de l'arbre en coupant quand la */
/*dissimilarite est superieure au seuil de dissimilarite*/
/*                                                      */
/********************************************************/

void decoupage_arbre_seuil_dissimilarite(int nb_motifs,noeud_t *noeud,double seuil_dissimilarite,
					 int *nb_clusters);

/******************************************************************/
/*                                                                */
/*Procedure de recherche des premiers noeuds valides d'une branche*/
/*i.e les premiers a avoir une perte d'inertie superieur au seuil */
/*                                                                */
/******************************************************************/

void cherche_premiers_noeuds_valides(noeud_t *noeud,double seuil_dissimilarite,
				     int *nb_noeuds_trouves,noeud_t **noeuds_trouves);

/*************************************************************************/
/*                                                                       */
/*Procedure de recherche de toutes les feuilles a partir d'un noeud donne*/
/*                                                                       */
/*************************************************************************/

void cherche_feuilles(noeud_t *noeud,int *nb_feuilles,noeud_t **feuilles);

/**********************************************************************************/
/*                                                                                */
/*Procedure de decouverte du seuil de dissimilarite a partir du nombre de clusters*/
/*                                                                                */
/**********************************************************************************/

void nb_clusters_to_dissimilarity_threshold(int nb_individus,int nb_clusters,noeud_t *racine,
					    double *seuil_dissimilarite);

/*******************************************************************/
/*                                                                 */
/*Procedure recursive de recuperation des dissimilarites d'un arbre*/
/*                                                                 */
/*******************************************************************/

void recupere_dissimilarites(noeud_t *noeud,int *compteur,double *dissimilarites);
