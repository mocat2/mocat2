#include "principal.h"
#include "dpc.h"
#include "treecut.h"
#include "graphpc.h"

/********************************************************************************/
/*                                                                              */
/*Procedure de decoupage d'arbre avec le test anova sur les densites de presence*/
/*                                                                              */
/********************************************************************************/

void decoupage_arbre_DPC(int nb_individus,int nb_dimensions,individu_t **individus,
			 noeud_t *racine,int *nb_clusters,int type_donnees,
			 int nb_clusters_selection,int donnees_normalisees,int type_densite,
			 int nb_simulations)
{
  /*declaration de variables*/
  int i,j,au_moins_un_noeud_secable,nb_noeuds_actifs,test,nb_fois_quon_separe;
  int noeud_courant,compteur;
  double **coordonnees_individus,*cotes_hyperpave,*valeurs_minimum,*valeurs_maximum,risque;
  double qualite;
  noeud_t **noeuds,**feuilles;
  groupe_t groupe_pere,groupe1,groupe2;
  /*fin declaration de variables*/

  /*allocation memoire*/
  noeuds=(noeud_t **)malloc(sizeof(noeud_t *)*nb_individus);
  feuilles=(noeud_t **)malloc(sizeof(noeud_t *)*nb_individus);
  groupe_pere.individus=(individu_t **)malloc(sizeof(individu_t *)*nb_individus);
  groupe1.individus=(individu_t **)malloc(sizeof(individu_t *)*nb_individus);
  groupe2.individus=(individu_t **)malloc(sizeof(individu_t *)*nb_individus);
  coordonnees_individus=(double **)malloc(sizeof(double *)*nb_individus);
  for(i=0;i<nb_individus;i++)
    {
      coordonnees_individus[i]=(double *)malloc(sizeof(double)*nb_dimensions);
    }
  cotes_hyperpave=(double *)malloc(sizeof(double)*nb_dimensions);
  valeurs_minimum=(double *)malloc(sizeof(double)*nb_dimensions);
  valeurs_maximum=(double *)malloc(sizeof(double)*nb_dimensions);
  /*fin allocation memoire*/

  /*on recupere toutes les feuilles*/
  compteur=0;
  cherche_feuilles(racine,&compteur,feuilles);
  for(i=0;i<nb_individus;i++)
    {
      for(j=0;j<nb_dimensions;j++)
	{
	  coordonnees_individus[i][j]=feuilles[i]->individu->valeurs_traitees[j];
	}
    }

  qualite=0;
  risque=0.0001;
  if(type_donnees==ALIGNEMENT)
    {
      risque=0.000001;
    }
  else if(type_densite==DENSITE2)
    {
      risque=0.001;
    }
  
    if(type_donnees==ALIGNEMENT)
    { 
      for(i=0;i<nb_dimensions;i++)
	{
	  valeurs_minimum[i]=0;
	  cotes_hyperpave[i]=1.0;
	}
    }
  else
    {
      /*on cherche les valeurs minimum et maximum pour chaque axe*/
      for(i=0;i<nb_dimensions;i++)
	{
	  valeurs_minimum[i]=coordonnees_individus[0][i];
	  valeurs_maximum[i]=coordonnees_individus[0][i];
	}
      
      for(i=0;i<nb_individus;i++)
	{
	  for(j=0;j<nb_dimensions;j++)
	    {
	      if(valeurs_minimum[j]>coordonnees_individus[i][j])
		{
		  valeurs_minimum[j]=coordonnees_individus[i][j];
		}
	      if(valeurs_maximum[j]<coordonnees_individus[i][j])
		{
		  valeurs_maximum[j]=coordonnees_individus[i][j];
		}
	    }
	}
      
      /*on calcule la longueur des cotes de l'hyperpave*/
      for(i=0;i<nb_dimensions;i++)
	{
	  cotes_hyperpave[i]=valeurs_maximum[i]-valeurs_minimum[i];
	}
    }
    
  nb_noeuds_actifs=1;
  noeuds[0]=racine;
  noeuds[0]->secable=OUI;
  *nb_clusters=0;
  while(1)
    {
      au_moins_un_noeud_secable=NON;
      for(i=0;i<nb_noeuds_actifs;i++)
	{
	  if(noeuds[i]->secable==OUI)
	    {
	      au_moins_un_noeud_secable=OUI;
	      noeud_courant=i;
	      break;
	    }
	}
      if(au_moins_un_noeud_secable==NON)
	{
	  break;
	}
      else
	{ 
	  /*on cherche les feuilles de la branche du noeud i*/
	  groupe_pere.nb_individus=0;
	  cherche_feuilles(noeuds[noeud_courant],&(groupe_pere.nb_individus),feuilles);
	  
	  for(i=0;i<groupe_pere.nb_individus;i++)
	    {
	      groupe_pere.individus[i]=feuilles[i]->individu;
	    }

	  if(groupe_pere.nb_individus<=2)
	    {
	      noeuds[noeud_courant]->secable=NON;
	      for(j=0;j<groupe_pere.nb_individus;j++)
		{
		  groupe_pere.individus[j]->cluster=*nb_clusters;
		}
	      (*nb_clusters)++;
	    }
	  else
	    {
	      /*on cherche les feuilles de la branche 1 du noeud i*/
	      groupe1.nb_individus=0;
	      cherche_feuilles(noeuds[noeud_courant]->copain1,&(groupe1.nb_individus),feuilles);
	      
	      for(i=0;i<groupe1.nb_individus;i++)
		{
		  groupe1.individus[i]=feuilles[i]->individu;
		}
	      
	      /*on cherche les feuilles de la branche 2 du noeud i*/
	      groupe2.nb_individus=0;
	      cherche_feuilles(noeuds[noeud_courant]->copain2,&(groupe2.nb_individus),feuilles);
	      
	      for(i=0;i<groupe2.nb_individus;i++)
		{
		  groupe2.individus[i]=feuilles[i]->individu;
		}
	      
	      if(nb_clusters_selection==GRAPHPC)
		{
		  /*test de division possible par graphpc*/
		  test=test_GRAPHPC(nb_individus,nb_dimensions,&groupe_pere,&groupe1,&groupe2);
		}
	      else if(nb_clusters_selection==DPC)
		{
		  /*test de division possible par DPC*/
		  nb_fois_quon_separe=0;
		  for(j=0;j<nb_simulations;j++)
		    {
		      /*on test si la separation de groupe_pere en groupe1 et
			groupe2 est pertinente*/
		      test=test_DPC(nb_individus,nb_dimensions,&groupe_pere,&groupe1,&groupe2,
				coordonnees_individus,cotes_hyperpave,valeurs_minimum,risque,
				    type_donnees,type_densite,donnees_normalisees,nb_simulations);
		      if(test==OUI)
			{
			  nb_fois_quon_separe++;
			}
		    }
		  if(nb_fois_quon_separe>nb_simulations/2)
		    {
		      test=OUI;
		      qualite+=(double)nb_fois_quon_separe/(double)nb_simulations;
		    }
		  else
		    {
		      test=NON;
		      qualite+=(double)(nb_simulations-nb_fois_quon_separe)/(double)nb_simulations;
		    }
		}
	      
	      if(test==OUI)
		{
		  noeuds[nb_noeuds_actifs]=noeuds[noeud_courant]->copain1;
		  noeuds[nb_noeuds_actifs+1]=noeuds[noeud_courant]->copain2;
		  nb_noeuds_actifs+=2;
		  noeuds[noeud_courant]->secable=NON;
		  noeuds[noeud_courant]->copain1->secable=OUI;
		  noeuds[noeud_courant]->copain2->secable=OUI;
		}
	      else
		{
		  for(j=0;j<groupe_pere.nb_individus;j++)
		    {
		      groupe_pere.individus[j]->cluster=*nb_clusters;
		    }
		  (*nb_clusters)++;
		  noeuds[noeud_courant]->secable=NON;
		}
	    }
	}
    }

  /*desallocation memoire*/
  free(noeuds);
  free(feuilles);
  free(groupe_pere.individus);
  free(groupe1.individus);
  free(groupe2.individus); 
  for(i=0;i<nb_individus;i++)
    {
      free(coordonnees_individus[i]);
    }
  free(coordonnees_individus);
  free(valeurs_minimum); 
  free(valeurs_maximum);
  free(cotes_hyperpave);
  /*fin desallocation memoire*/
}
