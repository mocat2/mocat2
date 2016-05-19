/**********************************************************************************/
/*                                                                                */
/*Procedure d'ecriture du fichier.dst des distances en vue du calcul de phylogenie*/
/*                                                                                */
/**********************************************************************************/

void ecriture_fichier_distances_sequences(char *nom_fichier,int nb_individus,
					  individu_t **individus);

/******************************************************************/
/*                                                                */
/*Procedure d'ecriture du fichier resultat presentant les clusters*/
/*                                                                */
/******************************************************************/

void ecriture_fichier_clusters(int nb_individus,int nb_dimensions,individu_t **individus,
			       int nb_clusters,int nb_clusters_selection,int individus_orphelins,
			       int write_coordinates,char *nom_fichier);

/**********************************************************/
/*                                                        */
/*Procedure d'ecriture du fichier de l'alignement multiple*/
/*avec les sequences reordonnees suivant leurs groupes    */
/*                                                        */
/**********************************************************/

void ecriture_fichier_clusters_alignement(char *nom_fichier,int longueur_alignement,
					  int nb_individus,individu_t **individus,
					  char **sequences,int individus_orphelins,
					  int nb_clusters);
