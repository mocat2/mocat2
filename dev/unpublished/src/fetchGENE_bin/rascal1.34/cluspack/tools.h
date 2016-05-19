#define PI 3.1415927

/******************************************************/
/*                                                    */
/*Procedure de calcul de la distance entre deux points*/
/*                                                    */
/******************************************************/

double calcul_distance_carre(int nb_dimensions,double *point1,double *point2);

/******************************************************/
/*                                                    */
/*Procedure de calcul de la distance entre deux points*/
/*                                                    */
/******************************************************/

double calcul_distance(int nb_dimensions,double *point1,double *point2);

/************************************************************************/
/*                                                                      */
/*Procedure de calcul de la fonction de repartition de la loi de Student*/
/*                                                                      */
/************************************************************************/

double F_student(double n,double t);

/********************************************************************************/
/*                                                                              */
/*Procedure de calcul de la fonction de repartition de la loi de Fisher-Snedecor*/
/*                                                                              */
/********************************************************************************/

double F_fisher_snedecor(double v1,double v2,double f);

/*************************************************************************************/
/*                                                                                   */
/*Procedure de calcul de la fonction de repartition de la loi normale_centree_reduite*/
/*                                                                                   */
/*************************************************************************************/

double F_normale_centree_reduite(double u);

/************************************************************/
/*                                                          */
/*Procedure de calcul de la fonction de repartition du Chi2,*/
/*(sa probabilite de depasser la valeur x)                  */
/*                                                          */
/************************************************************/

double chi_square(double x,int k);

/***************************************************************************/
/*                                                                         */
/*Procedure de regression lineaire sur les mesures (renvoie le coefficient)*/
/*                                                                         */
/***************************************************************************/

double regression_lineaire(int n,double **mesures);

/*****************************************************/
/*                                                   */
/*renvoie la valeur absolue de a                     */
/*                                                   */
/*****************************************************/
double valeur_absolue(double a);
 
/*****************************************************/
/*                                                   */
/*renvoie le minimum de a et de b                    */
/*                                                   */
/*****************************************************/
int min(int a,int b);

/****************************************/
/*                                      */
/*Procedure de tri par ordre decroissant*/
/*                                      */
/****************************************/

void tri_rapide(double *valeurs,int gauche,int droite);

/******************************/
/*                            */
/*sous-procedure de tri_rapide*/ 
/*                            */
/******************************/

void echanger_valeurs_pour_tri(double *valeurs,int element1,int element2);

/****************************************/
/*                                      */
/*Procedure de tri par ordre decroissant*/
/*                                      */
/****************************************/

void tri_rapide_valeurs_positions(double *valeurs,int *positions,int gauche,int droite);

/******************************/
/*                            */
/*sous-procedure de tri_rapide*/ 
/*                            */
/******************************/

void echanger_valeurs_positions_pour_tri(double *valeurs,int *positions,int element1,int element2);

/****************************************************************/
/*                                                              */
/*Procedure de test d'ajustement avec la distribution gaussienne*/
/*par Smirnov-Kolmogorov                                        */
/*                                                              */
/****************************************************************/

int test_ajustement_distribution_gaussienne(int nb_valeurs,double *valeurs,double moyenne,
					    double ecart_type);

/*****************************************************/
/*                                                   */
/*Procedure de calcul de la factorielle d'un entier n*/
/*                                                   */
/*****************************************************/

int calcul_factorielle(int n);

/***************************************************************/
/*                                                             */
/*Procedure de calcul de la factorielle "impaire" d'un entier n*/
/*                                                             */
/***************************************************************/

int calcul_factorielle_impaire(int n);
