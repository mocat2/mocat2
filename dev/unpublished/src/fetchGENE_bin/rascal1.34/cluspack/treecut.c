#include "principal.h"
#include "treecut.h"
#include "dissimilarities_clustering.h"

/*******************************************************************************/
/*                                                                             */
/*Procedure de decouverte du nombre de clusters par Secator a partir d'un arbre*/
/*                                                                             */
/*******************************************************************************/

void secator_tree(int nb_individus,noeud_t *racine,double *seuil_dissimilarite,int *nb_clusters)
{
  /*declaration de variables*/
  int compteur;  
  double *dissimilarites;
  /*fin declaration de variables*/

  /*allocation memoire*/
  dissimilarites=(double *)malloc(sizeof(double)*(nb_individus-1));
  /*fin allocation memoire*/

  compteur=0;
  recupere_dissimilarites(racine,&compteur,dissimilarites);

  /*tri des dissimilarites par ordre decroissant*/
  tri_rapide(dissimilarites,0,compteur-1); 

  /*decouverte du nombre de clusters*/
  clustering_dissimilarites_Holder(nb_individus-1,dissimilarites,nb_clusters,0);

  *seuil_dissimilarite=dissimilarites[*nb_clusters-2];

  /*desallocation memoire*/
  free(dissimilarites);
  /*fin desallocation memoire*/
}

/**********************************************************************************/
/*                                                                                */
/*Procedure de decouverte du seuil de dissimilarite a partir du nombre de clusters*/
/*                                                                                */
/**********************************************************************************/

void nb_clusters_to_dissimilarity_threshold(int nb_individus,int nb_clusters,noeud_t *racine,
					    double *seuil_dissimilarite)
{
  /*declaration de variables*/
  int i,compteur;
  double *dissimilarites;
  /*fin declaration de variables*/

  /*allocation memoire*/
  dissimilarites=(double *)malloc(sizeof(double)*(nb_individus-1));
  /*fin allocation memoire*/

  compteur=0;
  recupere_dissimilarites(racine,&compteur,dissimilarites);

  /*tri des dissimilarites par ordre decroissant*/
  tri_rapide(dissimilarites,0,compteur-1);

  *seuil_dissimilarite=dissimilarites[nb_clusters-2];

  /*desallocation memoire*/
  free(dissimilarites);
  /*fin desallocation memoire*/
}

/*******************************************************************/
/*                                                                 */
/*Procedure recursive de recuperation des dissimilarites d'un arbre*/
/*                                                                 */
/*******************************************************************/

void recupere_dissimilarites(noeud_t *noeud,int *compteur,double *dissimilarites)
{
  if(noeud->copain1!=NULL)
    {
      dissimilarites[*compteur]=noeud->dissimilarite;
      (*compteur)++;
      recupere_dissimilarites(noeud->copain1,compteur,dissimilarites);
      recupere_dissimilarites(noeud->copain2,compteur,dissimilarites);
    }
}

/********************************************************/
/*                                                      */
/*Procedure de decoupage de l'arbre en coupant quand la */
/*dissimilarite est superieure au seuil de dissimilarite*/
/*                                                      */
/********************************************************/

void decoupage_arbre_seuil_dissimilarite(int nb_motifs,noeud_t *noeud,
					 double seuil_dissimilarite,int *nb_clusters) 
{
  /*declaration des variables*/
  int i,j,nb_anciennes_branches,nb_nouvelles_branches,nb_noeuds_trouves,changement;
  int nb_feuilles;
  noeud_t **anciennes_branches,**nouvelles_branches,**noeuds_trouves,**feuilles;
  /*fin declaration des variables*/
  
  /*allocation memoire*/
  anciennes_branches=(noeud_t **)malloc(sizeof(noeud_t *)*nb_motifs);
  nouvelles_branches=(noeud_t **)malloc(sizeof(noeud_t *)*nb_motifs);
  noeuds_trouves=(noeud_t **)malloc(sizeof(noeud_t *)*nb_motifs);
  feuilles=(noeud_t **)malloc(sizeof(noeud_t *)*nb_motifs);
  /*fin allocation memoire*/

  nb_anciennes_branches=1;
  anciennes_branches[0]=noeud;
  
  while(1)
    {
      nb_nouvelles_branches=0;
      changement=NON;
      
      /*on cherche les premiers noeuds valides a partir de chaque branche*/
      for(i=0;i<nb_anciennes_branches;i++)
	{ 
	  nb_noeuds_trouves=0;

	  cherche_premiers_noeuds_valides(anciennes_branches[i],seuil_dissimilarite,
					  &nb_noeuds_trouves,noeuds_trouves);

	  if(nb_noeuds_trouves==0)
	    {
	      nouvelles_branches[nb_nouvelles_branches]=anciennes_branches[i];
	      nb_nouvelles_branches++;
	    }
	  else
	    {  
	      changement=OUI;
	      for(j=0;j<nb_noeuds_trouves;j++)
		{
		  nouvelles_branches[nb_nouvelles_branches]=noeuds_trouves[j]->copain1;
		  nb_nouvelles_branches++;
		  
		  nouvelles_branches[nb_nouvelles_branches]=noeuds_trouves[j]->copain2;
		  nb_nouvelles_branches++;
		}
	    }
	}
      
      if(changement==NON)
	{
	  break;
	}
      else
	{
	  nb_anciennes_branches=nb_nouvelles_branches;
	  for(i=0;i<nb_anciennes_branches;i++)
	    {
	      anciennes_branches[i]=nouvelles_branches[i];
	    }
	}      
    }

  *nb_clusters=0;
  /*on constitue les groupes a partir des feuilles de chaque branche*/
  for(i=0;i<nb_nouvelles_branches;i++)
    {
      nb_feuilles=0;
      cherche_feuilles(nouvelles_branches[i],&nb_feuilles,feuilles);
      
      for(j=0;j<nb_feuilles;j++)
	{
	  feuilles[j]->individu->cluster=*nb_clusters;
	}
      (*nb_clusters)++;
    }
  
  /*desallocation memoire*/
  free(anciennes_branches);
  free(nouvelles_branches);
  free(noeuds_trouves);
  free(feuilles);
  /*fin desallocation memoire*/

}

/******************************************************************/
/*                                                                */
/*Procedure de recherche des premiers noeuds valides d'une branche*/
/*i.e les premiers a avoir une perte d'inertie superieur au seuil */
/*                                                                */
/******************************************************************/

void cherche_premiers_noeuds_valides(noeud_t *noeud,double seuil_dissimilarite,
				     int *nb_noeuds_trouves,noeud_t **noeuds_trouves)
{
  if(noeud->dissimilarite>seuil_dissimilarite-ZERO_LIMIT)
    {
      noeuds_trouves[*nb_noeuds_trouves]=noeud;
      (*nb_noeuds_trouves)++;
    }
  else
    {
      if(noeud->copain1!=NULL)
	{
	  cherche_premiers_noeuds_valides(noeud->copain1,seuil_dissimilarite,
					  nb_noeuds_trouves,noeuds_trouves);
	  cherche_premiers_noeuds_valides(noeud->copain2,seuil_dissimilarite,
					  nb_noeuds_trouves,noeuds_trouves);
	}
    }
}

/*************************************************************************/
/*                                                                       */
/*Procedure de recherche de toutes les feuilles a partir d'un noeud donne*/
/*                                                                       */
/*************************************************************************/

void cherche_feuilles(noeud_t *noeud,int *nb_feuilles,noeud_t **feuilles)
{
  if(noeud->copain1==NULL)
    {
      feuilles[*nb_feuilles]=noeud;
      (*nb_feuilles)++;
    }
  else
    {
      cherche_feuilles(noeud->copain1,nb_feuilles,feuilles);
      cherche_feuilles(noeud->copain2,nb_feuilles,feuilles);
    }
}


