#include "principal.h"
#include "symeig.h"
#include "tools.h"
#include "normalized_2cut.h"

/********************************************************************/
/*                                                                  */
/*Procedure de recherche de coupe minimum normalisee de Xing et Karp*/
/*                                                                  */
/********************************************************************/

int normalized_cut(groupe_t *groupe,groupe_t *groupe1,groupe_t *groupe2)
{
  /*declaration de variables*/
  int i,j,nb_individus,*positions,meilleure_solution;
  double **DminusW,**A,*D,*W,**V,*vecteur_propre,meilleur_score,score,*y;
  double *coordonnees_triees;
  /*fin declaration de variables*/

  nb_individus=groupe->nb_individus;

  /*allocation memoire*/
  DminusW=(double **)malloc(sizeof(double *)*nb_individus);
  A=(double **)malloc(sizeof(double *)*nb_individus);
  V=(double **)malloc(sizeof(double *)*nb_individus);
  for(i=0;i<nb_individus;i++)
    {
      DminusW[i]=(double *)malloc(sizeof(double)*nb_individus);
      A[i]=(double *)malloc(sizeof(double)*nb_individus);
      V[i]=(double *)malloc(sizeof(double)*nb_individus);
    }
  D=(double *)malloc(sizeof(double)*nb_individus);
  W=(double *)malloc(sizeof(double)*nb_individus);
  y=(double *)malloc(sizeof(double)*nb_individus);
  vecteur_propre=(double *)malloc(sizeof(double)*nb_individus);
  coordonnees_triees=(double *)malloc(sizeof(double)*nb_individus);
  positions=(int *)malloc(sizeof(int)*nb_individus);
  /*fin allocation memoire*/

  for(i=0;i<nb_individus;i++)
    {
      D[i]=0;
      for(j=0;j<nb_individus;j++)
	{
	  D[i]+=(double)groupe->individus[i]->valeurs_traitees[groupe->individus[j]->id];
	  DminusW[i][j]=0;
	}
    }

  for(i=0;i<nb_individus;i++)
    {
      DminusW[i][i]=D[i];
      for(j=i;j<nb_individus;j++)
	{
	  DminusW[i][j]-=(double)groupe->individus[i]->valeurs_traitees[groupe->individus[j]->id];
	  A[i][j]=DminusW[i][j];
	  DminusW[j][i]=DminusW[i][j];
	  A[j][i]=A[i][j];
	}
    }

  for(i=0;i<nb_individus;i++)
    {
      for(j=0;j<nb_individus;j++)
	{
	  A[i][j]/=sqrt(D[i]);
	}
    }

  for(i=0;i<nb_individus;i++)
    {
      for(j=0;j<nb_individus;j++)
	{
	  A[j][i]/=sqrt(D[i]);
	}
    }

  /*on cherche les vecteurs propres*/
  symeig(A,nb_individus,V,W);

  /*on cherche le vecteur propre a la deuxieme plus petite valeur propre*/
  for(i=0;i<nb_individus;i++)
    {
      positions[i]=i;
    }

  /*on trie dans l'ordre decroissant les valeurs propres*/
  tri_rapide_valeurs_positions(W,positions,0,nb_individus-1);
  
  for(i=0;i<nb_individus;i++)
    {
      vecteur_propre[i]=V[i][positions[nb_individus-2]]/sqrt(D[i]);
      coordonnees_triees[i]=vecteur_propre[i];
    }

  for(i=0;i<nb_individus;i++)
    {
      positions[i]=i;
    }

  /*on trie dans l'ordre decroissant les coordonnees du vecteur propre*/
  tri_rapide_valeurs_positions(coordonnees_triees,positions,0,nb_individus-1);

  for(i=0;i<nb_individus;i++)
    {
      y[i]=-1.0;
    }

  meilleur_score=-VALEUR_ENORME;
  for(i=0;i<nb_individus;i++)
    {
      y[positions[i]]=1.0;
      score=evaluation_coupe_normalisee_minimum(nb_individus,y,D,DminusW);
      printf("score : %f\n",score);
      if(score>meilleur_score)
	{
	  meilleur_score=score;
	  meilleure_solution=i;
	  printf("i : %d \n",i);
	}
    }

  for(i=0;i<nb_individus;i++)
    {
      y[i]=-1.0;
    }
  for(i=0;i<=meilleure_solution;i++)
    {
      y[positions[i]]=1.0;
    }
  
  groupe1->nb_individus=0;
  groupe2->nb_individus=0;
  for(i=0;i<nb_individus;i++)
    {
      if(y[i]<0)
	{
	  groupe1->individus[groupe1->nb_individus]=groupe->individus[i];
	  groupe1->nb_individus++;
	}
      else
	{
	  groupe2->individus[groupe2->nb_individus]=groupe->individus[i];
	  groupe2->nb_individus++;
	}
    }
  printf("groupe1 : %d ; groupe2 : %d\n",groupe1->nb_individus,groupe2->nb_individus);

  /*desallocation memoire*/
  for(i=0;i<nb_individus;i++)
    {
      free(DminusW[i]);
      free(A[i]);
      free(V[i]);
    }
  free(DminusW);
  free(A);
  free(D);
  free(V);
  free(W);
  free(vecteur_propre);
  free(y);
  free(coordonnees_triees);
  free(positions);
  /*fin desallocation memoire*/

  return OUI;
}

/********************************************************/
/*                                                      */
/*Procedure renvoyant pour un vecteur donne la valeur   */
/*de la fonction objectif associee a la coupe normalisee*/
/*                                                      */
/********************************************************/

double evaluation_coupe_normalisee_minimum(int nb_dimensions,double *y,double *D,double **DminusW)
{
  /*declaration de variables*/
  int i,j;
  double res,numerateur,denominateur,*vecteur;
  /*fin declaration de variables*/
  
  /*allocation memoire*/
  vecteur=(double *)malloc(sizeof(double)*nb_dimensions);
  /*fin allocation memoire*/

  denominateur=0;
  for(i=0;i<nb_dimensions;i++)
    {
      denominateur+=y[i]*y[i]*D[i];
    }

  for(i=0;i<nb_dimensions;i++)
    {
      vecteur[i]=0;
      for(j=0;j<nb_dimensions;j++)
	{
	  vecteur[i]+=y[j]*DminusW[j][i];
	}
    }

  numerateur=0;
  for(i=0;i<nb_dimensions;i++)
    {
      numerateur+=vecteur[i]*y[i];
    }

  res=numerateur/denominateur;

  /*desallocation memoire*/
  free(vecteur);
  /*fin desallocation memoire*/

  return res;
}

