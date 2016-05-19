#include "principal.h"
#include "weighting.h"

/*********************************************/
/*                                           */
/*Procedure de calcul des poids des sequences*/
/*                                           */
/*********************************************/

void calcul_poids_sequences(noeud_t *noeud,double poids)
{
  /*declaration des variables*/
  int nb_feuilles_gauche,nb_feuilles_droite;
  /*fin declaration des variables*/

  if(noeud->copains[1]!=NULL)
    {
      nb_feuilles_gauche=0;
      nb_feuilles_droite=0;
      calcul_nombre_feuilles(noeud->copains[1],&nb_feuilles_gauche);
      calcul_nombre_feuilles(noeud->copains[2],&nb_feuilles_droite);
      if(noeud->copains[1]==NULL)
	{
	  /*on est a la racine*/
	  calcul_poids_sequences(noeud->copains[1],0);
	  calcul_poids_sequences(noeud->copains[2],0);
	}
      else
	{
	  calcul_poids_sequences(noeud->copains[1],poids+noeud->distances[0]/
				 (double)nb_feuilles_gauche);
	  calcul_poids_sequences(noeud->copains[2],poids+noeud->distances[0]/
				 (double)nb_feuilles_droite);
	}
    }
  else
    {
      noeud->poids=poids+noeud->distances[0];
    }
}

/***************************************************************/
/*                                                             */
/*Procedure de calcul du nombre de feuilles a partir d'un noeud*/
/*                                                             */
/***************************************************************/

void calcul_nombre_feuilles(noeud_t *noeud,int *nb_feuilles)
{
  if(noeud->copains[1]==NULL)
    {
      (*nb_feuilles)++;
    }
  else
    {
      calcul_nombre_feuilles(noeud->copains[1],nb_feuilles);
      calcul_nombre_feuilles(noeud->copains[2],nb_feuilles);
    }
  
}



