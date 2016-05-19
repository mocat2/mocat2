#include "principal.h"
#include "secatorbionj.h"
#include "dissimilarities_clustering.h"
#include "weighting.h"
#include "treecut.h"

/*****************************************************************/
/*                                                               */
/*Procedure de classification hierarchique base sur l'arbre bionj*/
/*avec decouverte du nombre de groupes par Secator               */
/*                                                               */
/*****************************************************************/

void secator_bionj(int nb_individus,noeud_t *racine,int *nb_clusters,int weighting,
		   int nb_clusters_selection,int nb_clusters_selected)
{
  /*declaration des variables*/
  int i,j,k,l,compteur_indices_dissimilarite,nb_exposants_Holder,resolution,nb_feuilles;
  double *historique_indices_dissimilarite,seuil_pertes_inertie,*RMSD_groupes;
  double *exposants_Holder,**pourcentages_identite;
  noeud_t **adresses_noeuds,**historique_noeuds,**feuilles;
  /*fin declaration des variables*/
  
  resolution=0;

  /*allocation memoire*/
  adresses_noeuds=(noeud_t **)malloc(sizeof(noeud_t)*(2*nb_individus-1));
  historique_indices_dissimilarite=(double *)malloc(sizeof(double)*nb_individus);
  historique_noeuds=(noeud_t **)malloc(sizeof(noeud_t *)*nb_individus);
  feuilles=(noeud_t **)malloc(sizeof(noeud_t *)*nb_individus);
  /*fin allocation memoire*/

  /*on recupere les adresses de tous les noeuds*/
  cherche_adresses_noeuds(racine,adresses_noeuds);

  if(weighting==OUI)
    {
      /*weighting des sequences s'il y a lieu*/
      calcul_poids_sequences(racine,0);
    }
  else
    {
      for(i=0;i<nb_individus;i++)
	{
	  adresses_noeuds[i]->poids=1.0;
	}
    }

  /*maintenant on peut deroote l'arbre comme on a les poids*/
  racine->copains[1]->copains[0]=racine->copains[2];
  racine->copains[2]->copains[0]=racine->copains[1];

  /*classification hierachique des sequences en respectant la topologie de l'arbre*/
  ward_sur_un_arbre_bionj(nb_individus,adresses_noeuds,racine,&compteur_indices_dissimilarite,
			  historique_indices_dissimilarite,historique_noeuds);


  /*affichage de l'arbre*/
  /*affichage_arbre(&noeud_fictif,OUI);*/

  /*tri par ordre decroissant des pertes d'inertie inter-classe*/
  tri_rapide_local(historique_indices_dissimilarite,historique_noeuds,0,
		   compteur_indices_dissimilarite-1);

  if(fabs(historique_indices_dissimilarite[0]-
	  historique_indices_dissimilarite[nb_individus-2])<ZERO_LIMIT)
    {
      *nb_clusters=1;
      nb_feuilles=0;
      cherche_feuilles(racine,&nb_feuilles,feuilles);
      for(i=0;i<nb_individus;i++)
	{
	  feuilles[i]->individu->cluster=0;
	}
    }
  else
    {
      if(nb_clusters_selection==SECATOR)
	{
	  /*dissimilarities clustering using Holder exponent*/
	  clustering_dissimilarites_Holder(nb_individus-1,historique_indices_dissimilarite,
					   nb_clusters,resolution);
	}
      else if(nb_clusters_selection==FIXED)
	{
	  *nb_clusters=nb_clusters_selected;
	}
      
      /*clusters's building*/
      clusters_building(nb_individus,historique_indices_dissimilarite[*nb_clusters-2],
			racine,adresses_noeuds);
    }  

  /*desallocation memoire*/
  for(i=0;i<2*nb_individus-1;i++)
    {
      free(adresses_noeuds[i]);
    }
  free(adresses_noeuds);
  free(historique_indices_dissimilarite);
  free(historique_noeuds);
  free(feuilles);
  /*fin desallocation memoire*/
  
}

/********************************************************/
/*                                                      */
/*Procedure de recherche des adresses de tous les noeuds*/
/*                                                      */
/********************************************************/

void cherche_adresses_noeuds(noeud_t *noeud,noeud_t **adresses_noeuds)
{
  adresses_noeuds[noeud->numero]=noeud;

  if(noeud->copains[1]!=NULL)
    {
      cherche_adresses_noeuds(noeud->copains[1],adresses_noeuds);
      cherche_adresses_noeuds(noeud->copains[2],adresses_noeuds);
    }
}

/**************************************************************/
/*                                                            */
/*Procedure de construction de la solution dans le cas general*/
/*                                                            */
/**************************************************************/

void clusters_building(int nb_individus,double seuil_pertes_inertie,
		       noeud_t *racine,noeud_t **adresses_noeuds)
{
  /*declaration des variables*/
  int i,j,nb_anciennes_branches,nb_nouvelles_branches,nb_noeuds_trouves,changement,nb_feuilles;
  noeud_t **anciennes_branches,**nouvelles_branches,**noeuds_trouves,**feuilles;
  /*fin declaration des variables*/
  
  /*allocation memoire*/
  anciennes_branches=(noeud_t **)malloc(sizeof(noeud_t *)*nb_individus);
  nouvelles_branches=(noeud_t **)malloc(sizeof(noeud_t *)*nb_individus);
  noeuds_trouves=(noeud_t **)malloc(sizeof(noeud_t *)*nb_individus);
  feuilles=(noeud_t **)malloc(sizeof(noeud_t *)*nb_individus);
  /*fin allocation memoire*/

  nb_anciennes_branches=1;
  anciennes_branches[0]=racine;
  
  while(1)
    {
      nb_nouvelles_branches=0;
      changement=NON;
      
      /*on cherche les premiers noeuds valides a partir de chaque branche*/
      for(i=0;i<nb_anciennes_branches;i++)
	{ 
	  nb_noeuds_trouves=0;

	  cherche_premiers_noeuds_valides(anciennes_branches[i],seuil_pertes_inertie,
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

  /*on each bough leaves receive their clusters numbers*/
  for(i=0;i<nb_nouvelles_branches;i++)
    {
      nb_feuilles=0;
      cherche_feuilles(nouvelles_branches[i],&nb_feuilles,feuilles);
      
      for(j=0;j<nb_feuilles;j++)
	{
	  feuilles[j]->individu->cluster=i;
	}
    }
  
  /*desallocation memoire*/
  free(anciennes_branches);
  free(nouvelles_branches);
  free(noeuds_trouves);
  free(feuilles);
  /*fin desallocation memoire*/
}

/****************************************/
/*                                      */
/*Procedure de tri par ordre decroissant*/
/*                                      */
/****************************************/

void tri_rapide_local(double *valeurs,noeud_t **noeuds,int gauche,int droite)
{ 
  /*declaration des variables*/
  int element_suivant;
  int indice_separateur;
  /*fin de declaration des variables*/
  
  if(gauche>=droite)
    {
      return;
    }
  
  indice_separateur=gauche;
  element_suivant=gauche+1;
  
  while(element_suivant<=droite)
    {
      if(valeurs[element_suivant]>valeurs[indice_separateur])
	{
	  echanger_local(valeurs,noeuds,indice_separateur+1,element_suivant);
	  echanger_local(valeurs,noeuds,indice_separateur,indice_separateur+1);
	  indice_separateur++;
	}
      element_suivant++;
    }
  
  tri_rapide_local(valeurs,noeuds,gauche,indice_separateur-1);
  tri_rapide_local(valeurs,noeuds,indice_separateur+1,droite);
}

/*******************************************/
/*                                         */
/*sous-procedure de tri_rapide             */ 
/*                                         */
/*******************************************/

void echanger_local(double *valeurs,noeud_t **noeuds,int element1,int element2)
{
  /*declaration des variables*/
  double temp_val;
  noeud_t *temp_noeud;
  /*fin de declaration des variables*/

  temp_val=valeurs[element1];
  valeurs[element1]=valeurs[element2];
  valeurs[element2]=temp_val;

  temp_noeud=noeuds[element1];
  noeuds[element1]=noeuds[element2];
  noeuds[element2]=temp_noeud;
}



