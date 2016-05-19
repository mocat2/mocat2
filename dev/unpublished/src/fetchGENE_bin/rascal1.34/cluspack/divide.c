#include "principal.h"
#include "divide.h"
#include "2means.h"
#include "dpc.h"
#include "space.h"
#include "tools.h"
#include "normalized_2cut.h"

/**************************************/
/*                                    */
/*Clustering par divisions successives*/
/*                                    */
/**************************************/

void divide(int nb_individus,int nb_dimensions,int *nb_clusters,individu_t **individus,
	    int type_donnees,int clustering_method,int nb_clusters_selection,
	    int donnees_normalisees,int type_densite)
{
  /*declaration des variables*/
  int i,j,k,resultat_test_groupe1,resultat_test_groupe2,compteur,nb_fois_quon_separe;
  int nb_simulations,succes_classification;
  groupe_t *premier_groupe,*groupe_courant,*groupe1,*groupe2,*nouveau_groupe;
  double **coordonnees_individus,*cotes_hyperpave,*valeurs_minimum,*valeurs_maximum;
  double qualite,risque;
  /*fin declaration des variables*/

  /********************/
  /*                  */
  /*Corps du programme*/
  /*                  */
  /********************/

  /*allocation memoire*/
  premier_groupe=(groupe_t *)malloc(sizeof(groupe_t));
  premier_groupe->individus=(individu_t **)malloc(sizeof(individu_t *)*nb_individus);
  premier_groupe->precedent=NULL;
  premier_groupe->suivant=NULL;

  groupe1=(groupe_t *)malloc(sizeof(groupe_t));
  groupe1->individus=(individu_t **)malloc(sizeof(individu_t *)*nb_individus);

  groupe2=(groupe_t *)malloc(sizeof(groupe_t));
  groupe2->individus=(individu_t **)malloc(sizeof(individu_t *)*nb_individus);
  coordonnees_individus=(double **)malloc(sizeof(double *)*nb_individus);
  for(i=0;i<nb_individus;i++)
    {
      coordonnees_individus[i]=(double *)malloc(sizeof(double)*nb_dimensions);
    }
  cotes_hyperpave=(double *)malloc(sizeof(double)*nb_dimensions);
  valeurs_minimum=(double *)malloc(sizeof(double)*nb_dimensions);
  valeurs_maximum=(double *)malloc(sizeof(double)*nb_dimensions);
  /*fin allocation memoire*/

  for(i=0;i<nb_individus;i++)
    {
      for(j=0;j<nb_dimensions;j++)
	{
	  coordonnees_individus[i][j]=individus[i]->valeurs_traitees[j];
	}
    }

  if(nb_clusters_selection==DPC)
    {
      qualite=0;
      risque=0.01;
    
      if(type_donnees==ALIGNEMENT)
	{
	  risque=0.0001;
	}
      else if(type_densite==DENSITE2)
	{
	  risque=0.01;
	}  
      
      nb_simulations=10;
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
    }

  /*on initialise le premier groupe*/
  premier_groupe->nb_individus=nb_individus;
  for(i=0;i<nb_individus;i++)
    {
      premier_groupe->individus[i]=individus[i];
    }
  premier_groupe->insecable=NON;

  *nb_clusters=1;
  while(1)
    {
      groupe_courant=premier_groupe;
      while(groupe_courant!=NULL)
	{
	  if(groupe_courant->insecable==NON)
	    {
	      break;
	    }
	  else
	    {
	      groupe_courant=groupe_courant->suivant;
	    }
	}
      if(groupe_courant==NULL)
	{
	  break;
	}

      /*on regarde si dans le groupe courant il y a au moins deux individus differents*/
      /*auquel cas on peut appliquer l'algo 2-means*/
      if(clustering_method==KMEANS)
	{
	  /*2-means*/
	  succes_classification=
	    two_means(nb_dimensions,groupe_courant,groupe1,groupe2,type_donnees);
	}
      else if(clustering_method==NORMALIZED_CUT)
	{
	  /*normalized 2cut*/
	  succes_classification=normalized_cut(groupe_courant,groupe1,groupe2);
	}
      
      if((groupe1->nb_individus<0)||(groupe2->nb_individus<0))
	{
	  groupe_courant->insecable=OUI;
	}
      else
	{
	  if(succes_classification==NON)
	    {
	      resultat_test_groupe1=NON;
	      resultat_test_groupe2=NON;
	    }
	  else if(nb_clusters_selection==GRAPHPC)
	    {
	      resultat_test_groupe1=test_GRAPHPC(nb_individus,nb_dimensions,
						 groupe_courant,groupe1,groupe2,
						 coordonnees_individus);
	    }
	  else if(nb_clusters_selection==DPC)
	    {
	      nb_fois_quon_separe=0;
	      for(i=0;i<nb_simulations;i++)
		{
		  resultat_test_groupe1=
		    test_DPC(nb_individus,nb_dimensions,groupe_courant,groupe1,
			     groupe2,coordonnees_individus,cotes_hyperpave,
			     valeurs_minimum,risque,type_donnees,type_densite,
			     donnees_normalisees,1);
		  if(resultat_test_groupe1==OUI)
		    {
		      nb_fois_quon_separe++;
		    }
		}
	      
	      if(nb_fois_quon_separe>(int)floor(0.6*(double)nb_simulations))
		{
		  resultat_test_groupe1=OUI;
		  qualite+=(double)nb_fois_quon_separe/(double)nb_simulations;
		}
	      else
		{
		  resultat_test_groupe1=NON;
		  qualite+=(double)(nb_simulations-nb_fois_quon_separe)/
		    (double)nb_simulations;
		}
	    }
	  
	  if(resultat_test_groupe1==OUI)
	    {
	      (*nb_clusters)++;
	    }
	  printf("current number of classes : %d\n",*nb_clusters);
	  
	  resultat_test_groupe2=resultat_test_groupe1;
	  
	  if((resultat_test_groupe1==OUI)||(resultat_test_groupe2==OUI))
	    { 
	      /*allocation memoire*/
	      nouveau_groupe=(groupe_t *)malloc(sizeof(groupe_t));
	      nouveau_groupe->individus=(individu_t **)malloc(sizeof(individu_t *)*
							      groupe1->nb_individus);
	      /*fin allocation memoire*/
	      
	      nouveau_groupe->precedent=premier_groupe;
	      nouveau_groupe->suivant=premier_groupe->suivant;
	      if(premier_groupe->suivant!=NULL)
		{
		  premier_groupe->suivant->precedent=nouveau_groupe;
		}
	      premier_groupe->suivant=nouveau_groupe;
	      
	      nouveau_groupe->nb_individus=groupe1->nb_individus;
	      for(i=0;i<nouveau_groupe->nb_individus;i++)
		{
		  nouveau_groupe->individus[i]=groupe1->individus[i];
		}
	      
	      groupe_courant->nb_individus=groupe2->nb_individus;
	      for(i=0;i<groupe_courant->nb_individus;i++)
		{
		  groupe_courant->individus[i]=groupe2->individus[i];
		}
	      
	      if(nouveau_groupe->nb_individus<3)
		{
		  nouveau_groupe->insecable=OUI;
		}
	      else
		{
		  nouveau_groupe->insecable=NON;
		}
	      
	      if(groupe_courant->nb_individus<3)
		{
		  groupe_courant->insecable=OUI;
		}
	      else
		{
		  groupe_courant->insecable=NON;
		}
	    }
	  else
	    {
	      groupe_courant->insecable=OUI;
	    }
	}
    }

  /*om compte le nombre de groupes*/
  groupe_courant=premier_groupe;
  *nb_clusters=0;
  
  while(groupe_courant!=NULL)
    {
      for(i=0;i<groupe_courant->nb_individus;i++)
	{
	  groupe_courant->individus[i]->cluster=*nb_clusters;
	}
      groupe_courant=groupe_courant->suivant;
      (*nb_clusters)++;
    }

  /*desallocation memoire*/
  free(groupe1->individus);
  free(groupe1);
  free(groupe2->individus);
  free(groupe2);

  groupe_courant=premier_groupe->suivant;
  if(groupe_courant==NULL)
    {
      free(premier_groupe->individus);
      free(premier_groupe);
    }
  else
    {
      while(1)
	{
	  free(groupe_courant->precedent->individus);
	  free(groupe_courant->precedent);

	  if(groupe_courant->suivant!=NULL)
	    {
	      groupe_courant=groupe_courant->suivant;
	    }
	  else
	    { 
	      free(groupe_courant->individus);
	      free(groupe_courant);
	      break;
	    }
	}
    }
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

