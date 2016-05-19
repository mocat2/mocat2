/***********************************************************/
/*                                                         */
/*Procedure de lecture du fichier d'un arbre phylogenetique*/
/*(au format bionj)                                        */
/*                                                         */
/***********************************************************/

void lecture_arbre_phylogenetique(char *nom_fichier,char **arbre);

/**********************************************************/
/*                                                        */
/*Procedure de lecture du fichier d'un alignement multiple*/
/*                                                        */
/**********************************************************/

void lecture_fichier_alignement(char *nom_alignement,int *longueur_alignement,int *nb_individus,
				individu_t ***individus,char ***sequences,int type_donnees);

/***************************************/
/*                                     */
/*Procedure de lecture d'un fichier.tfa*/
/*                                     */
/***************************************/

void lecture_fichier_tfa(char *nom_alignement,int *longueur_alignement,
			 int *nb_individus,char ***sequences,individu_t ***individus,int type_donnees);

/***************************************/
/*                                     */
/*Procedure de lecture d'un fichier.msf*/
/*                                     */
/***************************************/

void lecture_fichier_msf(char *nom_alignement,int *longueur_alignement,
			 int *nb_individus,char ***sequences,individu_t ***individus);

/********************************************/
/*                                          */
/*Procedure de lecture du fichier contenant */
/*les coordonnees des individus a classifier*/
/*                                          */
/********************************************/

void lecture_fichier_coordonnees(int *nb_individus,int *nb_dimensions,individu_t ***individus,
				 char *nom_fichier_entree);
