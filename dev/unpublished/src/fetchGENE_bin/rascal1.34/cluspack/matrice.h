/*******************************************************************************/
/*                                                                             */
/*Procedure de calcul de la norme au carre d'un vecteur pour une matrice donnee*/
/*                                                                             */
/*******************************************************************************/

double calcul_norme(int nb_dimensions,double *vecteur,double **matrice);

/****************************************************/
/*                                                  */
/*Procedure d'estimation de la matrice de covariance*/
/*                                                  */
/****************************************************/

void calcul_matrice_covariance(int nb_dimensions,int nb_points,double **points,
			       double **matrice_covariance);

/*******************************************************************/
/*                                                                 */
/*Procedure d'inversion de matrice qui renvoie le determinant      */
/*de la matrice a inverser ou 0 s'il y a un depassement de capacite*/
/*                                                                 */
/*******************************************************************/

double inversion_matrice(int nb_dimensions,double **matrice);

/*****************************************************/
/*                                                   */
/*Procedure d'echange de deux lignes dans une matrice*/
/*                                                   */
/*****************************************************/

void echange_lignes(int nb_dimensions,double **matrice,int ligne1,int ligne2);

