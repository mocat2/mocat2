#define TROUVE 1
#define PAS_TROUVE 0

/******************************************************************************************/
/*                                                                                        */
/*Procedure construisant l'arbre bionj a partir d'un alignement de sequence dans sequences*/
/*                                                                                        */
/******************************************************************************************/

void creation_arbre_bionj(int nb_individus,individu_t **individus,
			  char *fichier_entree,noeud_t **racine);

/**************************************************/
/*                                                */
/*Procedure d'association des individus aux noeuds*/
/*                                                */
/**************************************************/

void association_individus_noeuds(int nb_individus,individu_t **individus,noeud_t *noeud);

/*************************************************************/
/*                                                           */
/*Procedure d'enracinement de l'arbre sur la sequence fictive*/
/*                                                           */
/*************************************************************/

void enracinement_arbre(int nb_individus,noeud_t **racine);

/**********************************************************************/
/*                                                                    */
/*Procedure de renommage des noeuds internes car apres restructuration*/
/*le nombre de noeuds n'est plus le meme                              */
/*                                                                    */
/**********************************************************************/

void renommage_noeuds_internes(int nb_individus,noeud_t *noeud,int *nb_noeuds_internes,
			       int *nb_feuilles);

/****************************************************************************/
/*                                                                          */
/*Procedure de restructuration de l'arbre et de calcul de la nouvelle racine*/
/*                                                                          */
/****************************************************************************/

void restructuration_arbre(noeud_t *noeud,noeud_t *nouveau_pere,double nouvelle_distance);

/************************************************************/
/*                                                          */
/*Procedure de recherche de la sequence fictive dans l'arbre*/
/*                                                          */
/************************************************************/

int cherche_noeud_sequence_fictive(int nb_individus,noeud_t *noeud,
				   noeud_t **noeud_sequence_fictive);

/*****************************************************/
/*                                                   */
/*Procedure de construction de l'arbre phylogenetique*/
/*                                                   */
/*****************************************************/

void construction_arbre(int nb_individus,char *arbre,noeud_t **racine);

/*******************************************************************************/
/*                                                                             */
/*Procedure de copie d'une portion de chaine de caracteres delimitee de maniere*/
/*stricte par deux pointeurs sur char dans une autre chaine de caracteres      */
/*                                                                             */
/*******************************************************************************/

void copie_partie_string(char *debut,char *fin,char *destination);

/******************************************************************************/
/*                                                                            */
/*Procedure de remplacement d'une partie d'une chaine de caracteres, delimitee*/
/*de maniere large par debut et fin, par une autre chaine de caracteres       */
/*                                                                            */
/******************************************************************************/

void remplace_partie_string(char **chaine,char *debut,char *fin,char *chaine_de_remplacement);

/*******************************************************************/
/*                                                                 */
/*Procedure d'affichage de l'arbre (essentiellement pour debuggage)*/
/*                                                                 */
/*******************************************************************/

void affichage_arbre(noeud_t *noeud);


