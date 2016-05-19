#define FEUILLE 0
#define NOEUD_ORDINAIRE 1

/*********************************************************************************/
/*                                                                               */
/*Procedure de classification hierarchique par la methode des voisins reciproques*/ 
/*                                                                               */
/*********************************************************************************/

void classification_ward(int nb_individus,int nb_dimensions,individu_t **individus,
			 int type_donnees,noeud_t **racine);
