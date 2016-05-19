#define EXPOSANT_HOLDER_SEUIL 0.5

/*************************************************************/
/*                                                           */
/*classification des dissimilarites pour trouver le nombre de*/
/*dissimilarites elevees et donc le nombre de groupes        */
/*(avec l'exposant de Holder)                                */
/*                                                           */
/*************************************************************/

void clustering_dissimilarites_Holder(int nb_dissimilarites,double *historique_dissimilarites,
				      int *nb_dissimilarites_elevees,int resolution);

/******************************************************************************/
/*                                                                            */
/*Procedure de calcul des exposants d'Holder sur les valeurs de dissimilarites*/
/*                                                                            */
/******************************************************************************/

void calcul_exposants_Holder(int nb_dissimilarites,double *historique_dissimilarites,
			     int *nb_exposants_Holder,double *exposants_Holder);

/**********************************************************************/
/*                                                                    */
/*Procedure de separation en deux groupes des valeurs de dissimilarite*/
/*par calcul du minimum d'inertie intraclasse                         */
/*                                                                    */
/**********************************************************************/

void separation_dissimilarites_bete(int nb_dissimilarites,double *historique_dissimilarites,
				    int *nb_dissimilarites_premier_groupe);

/************************************************************************/
/*                                                                      */
/*Procedure de separation en deux groupes des valeurs de dissimilarite  */
/*par calcul du minimum d'inertie intraclasse mais uniquement aux points*/
/*d'exposants d'Holder faibles                                          */
/*                                                                      */
/************************************************************************/

void separation_dissimilarites_evolue(int nb_exposants_Holder,int nb_dissimilarites,
				      double *historique_dissimilarites,
				      int *nb_dissimilarites_premier_groupe,
				      double *exposants_Holder);
