#include "principal.h"
#include "matrice.h"

/****************************************************/
/*                                                  */
/*Procedure d'estimation de la matrice de covariance*/
/*                                                  */
/****************************************************/

void calcul_matrice_covariance(int nb_dimensions,int nb_points,double **points,
			       double **matrice_covariance)
{
  /*declaration de variables*/
  int i,j,k;
  double *moyennes,diff;
  /*fin declaration de variables*/
  
  /*allocation memoire*/
  moyennes=(double *)malloc(sizeof(double)*nb_dimensions);
  /*fin allocation memoire*/

  /*calcul des moyennes*/
  for(j=0;j<nb_dimensions;j++)
    {
      moyennes[j]=0;
    }
  for(i=0;i<nb_points;i++)
    {
      for(j=0;j<nb_dimensions;j++)
	{
	  moyennes[j]+=points[i][j];
	}
    }
  for(j=0;j<nb_dimensions;j++)
    {
      moyennes[j]/=(double)nb_points;
    }

  for(i=0;i<nb_dimensions;i++)
    {
      for(j=0;j<nb_dimensions;j++)
	{
	  matrice_covariance[i][j]=0;
	}
    }	  
  for(k=0;k<nb_points;k++)
    {
      for(i=0;i<nb_dimensions;i++)
	{
	  diff=points[k][i]-moyennes[i];
	  for(j=i;j<nb_dimensions;j++)
	    {
	      matrice_covariance[i][j]+=diff*(points[k][j]-moyennes[j]);
	    }
	}	 
    }
   for(i=0;i<nb_dimensions;i++)
    {
      for(j=i;j<nb_dimensions;j++)
	{
	  matrice_covariance[i][j]/=(double)(nb_points-1);
	  
	  matrice_covariance[j][i]=matrice_covariance[i][j];
	}
    }	

  /*desallocation memoire*/
  free(moyennes);
  /*fin desallocation memoire*/
}


/*******************************************************************/
/*                                                                 */
/*Procedure d'inversion de matrice qui renvoie le determinant      */
/*de la matrice a inverser ou 0 s'il y a un depassement de capacite*/
/*                                                                 */
/*******************************************************************/

double inversion_matrice(int nb_dimensions,double **matrice)
{
  /*declaration de variables*/
  int i,j,k,pivot;
  double **matrice_inverse,valeur;
  double determinant;
  /*fin declaration de variables*/

  /*allocation memoire*/
  matrice_inverse=(double **)malloc(sizeof(double *)*nb_dimensions);
  for(i=0;i<nb_dimensions;i++)
    {
      matrice_inverse[i]=(double *)malloc(sizeof(double)*nb_dimensions);
    }
  /*fin allocation memoire*/

  for(i=0;i<nb_dimensions;i++)
    {
      for(j=0;j<nb_dimensions;j++)
	{
	  matrice_inverse[i][j]=0;
	}
      matrice_inverse[i][i]=1.0;
    }

  determinant=0.0;
  for(i=0;i<nb_dimensions;i++)
    { 
      pivot=-1;
      for(j=i;j<nb_dimensions;j++)
	{
	  if(fabs(matrice[j][i])>ZERO_LIMIT)
	    {
	      pivot=j;
	      break;
	    }
	}
      
      if(pivot==-1)
	{
	  return 0.0;
	}
      echange_lignes(nb_dimensions,matrice,i,pivot);
      echange_lignes(nb_dimensions,matrice_inverse,i,pivot);
 
      valeur=matrice[i][i];

      for(j=0;j<nb_dimensions;j++)
	{
	  matrice[i][j]/=valeur;
	  matrice_inverse[i][j]/=valeur;
	}
      
      determinant+=log(fabs(valeur));
      
      for(j=0;j<nb_dimensions;j++)
	{
	  if(j!=i)
	    { 
	      valeur=matrice[j][i];
	      for(k=0;k<nb_dimensions;k++)
		{
		  matrice[j][k]-=valeur*matrice[i][k];
		  matrice_inverse[j][k]-=valeur*matrice_inverse[i][k];
		}
	    }	
	}

    }
  
  for(i=0;i<nb_dimensions;i++)
    {
      for(j=0;j<nb_dimensions;j++)
	{
	  matrice[i][j]=matrice_inverse[i][j];
	}
    }

  /*desallocation memoire*/
  for(i=0;i<nb_dimensions;i++)
    {
      free(matrice_inverse[i]);
    }
  free(matrice_inverse);
  /*fin desallocation memoire*/
  
  return determinant;
}

/*****************************************************/
/*                                                   */
/*Procedure d'echange de deux lignes dans une matrice*/
/*                                                   */
/*****************************************************/

void echange_lignes(int nb_dimensions,double **matrice,int ligne1,int ligne2)
{
  /*declaration de variables*/
  int i;
  double val_temp;
  /*fin declaration de variables*/

  for(i=0;i<nb_dimensions;i++)
    {
      val_temp=matrice[ligne1][i];
      matrice[ligne1][i]=matrice[ligne2][i];
      matrice[ligne2][i]=val_temp;
    }
}

/*******************************************************************************/
/*                                                                             */
/*Procedure de calcul de la norme au carre d'un vecteur pour une matrice donnee*/
/*                                                                             */
/*******************************************************************************/

double calcul_norme(int nb_dimensions,double *vecteur,double **matrice)
{
  /*declaration de variables*/
  int i,j;
  double norme,*vecteuri;
  /*fin declaration de variables*/
  
  /*allocation memoire*/
  vecteuri=(double *)malloc(sizeof(double)*nb_dimensions);
  /*fin allocation memoire*/
  
  /*calcul du produit de gauche*/
  for(i=0;i<nb_dimensions;i++)
    {
      vecteuri[i]=0;
      for(j=0;j<nb_dimensions;j++)
	{
	  vecteuri[i]+=vecteur[j]*matrice[j][i];
	}
    }

  /*calcul du produit final*/
  norme=0;
  for(i=0;i<nb_dimensions;i++)
    {
      norme+=vecteuri[i]*vecteur[i];
    }

  /*desallocation memoire*/
  free(vecteuri);
  /*fin desallocation memoire*/

  return norme;
}

