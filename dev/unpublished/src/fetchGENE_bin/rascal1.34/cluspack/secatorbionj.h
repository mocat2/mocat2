/****************************************/
/*                                      */
/*Procedure de tri par ordre decroissant*/
/*                                      */
/****************************************/

void tri_rapide_local(double *valeurs,noeud_t **noeuds,int gauche,int droite);

/*******************************************/
/*                                         */
/*sous-procedure de tri_rapide             */ 
/*                                         */
/*******************************************/

void echanger_local(double *valeurs,noeud_t **noeuds,int element1,int element2);

/*****************************************************************/
/*                                                               */
/*Procedure de classification hierarchique base sur l'arbre bionj*/
/*avec decouverte du nombre de groupes par Secator               */
/*                                                               */
/*****************************************************************/

void secator_bionj(int nb_individus,noeud_t *racine,int *nb_clusters,int weighting,
		   int nb_clusters_selection,int nb_clusters_selected);

/**************************************************/
/*                                                */
/*Procedure de recherche des adresses des feuilles*/
/*                                                */
/**************************************************/

void cherche_adresses_noeuds(noeud_t *noeud,noeud_t **adresses_noeuds);

/**************************************************************/
/*                                                            */
/*Procedure de construction de la solution dans le cas general*/
/*                                                            */
/**************************************************************/

void clusters_building(int nb_individus,double seuil_pertes_inertie,
		       noeud_t *racine,noeud_t **adresses_noeuds);
