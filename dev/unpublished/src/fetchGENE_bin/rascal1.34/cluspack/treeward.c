#include "principal.h"
#include "treeward.h"

/********************************************************************/
/*                                                                  */
/*Procedure de classification par la methode des voisins reciproques*/ 
/*                                                                  */
/********************************************************************/

void classification_ward(int nb_individus,int nb_dimensions,individu_t **individus,
			 int type_donnees,noeud_t **racine)
{
  /*declaration des variables*/
  int i,j,k,l,nb_racines,graine,plus_proche_voisin,ancienne_graine,compteur_numeros;
  double **dissimilarites_courantes,dissimilarite_minimum;
  double *poids,*historique_dissimilarites;
  double dissimilarite_minimum_paria,seuil_dissimilarite;
  couple_t meilleur_regroupement;
  noeud_t *feuilles,**noeuds,*nouveau_noeud;
  /*fin de declaration des variables*/

  /*allocation de la memoire*/
  dissimilarites_courantes=(double **)malloc(sizeof(double *)*nb_individus);
  for(i=0;i<nb_individus;i++)
    {
      dissimilarites_courantes[i]=(double *)malloc(sizeof(double)*nb_individus);
    }
  poids=(double *)malloc(sizeof(double)*nb_individus);
  feuilles=(noeud_t *)malloc(sizeof(noeud_t)*nb_individus);
  noeuds=(noeud_t **)malloc(sizeof(noeud_t *)*nb_individus);

  historique_dissimilarites=(double *)malloc(sizeof(double)*nb_individus);
  /*fin allocation de la memoire*/

  /*initialisation des dissimilarites*/
  for(i=0;i<nb_individus;i++)
    {
      poids[i]=1.0;
      dissimilarites_courantes[i][i]=0;
    }

  compteur_numeros=0;
  /*initialisation des dissimilarites*/
  for(i=0;i<nb_individus-1;i++)
    {
      for(j=i+1;j<nb_individus;j++)
	{       
	  if(type_donnees==DISTANCES)
	    {
	      dissimilarites_courantes[i][j]=pow((double)(individus[i]->valeurs_traitees[j]),
						 2.0)*poids[i]*poids[j]/(poids[i]+poids[j]);
	      dissimilarites_courantes[j][i]=dissimilarites_courantes[i][j];
	    }
	  else
	    {
	      dissimilarites_courantes[i][j]=0.0;
	      for(k=0;k<nb_dimensions;k++)
		{
		  dissimilarites_courantes[i][j]+=
		    pow((double)(individus[i]->valeurs_traitees[k]-
				 individus[j]->valeurs_traitees[k]),2.0);
		}
	      
	      dissimilarites_courantes[i][j]*=poids[i]*poids[j]/(poids[i]+poids[j]);
	      dissimilarites_courantes[j][i]=dissimilarites_courantes[i][j];
	    }
	}
    }
    
   /*initialisation de la construction de l'arbre*/
   for(i=0;i<nb_individus;i++)
     {
       feuilles[i].copain1=NULL;
       feuilles[i].copain2=NULL;
       feuilles[i].numero=i;
       feuilles[i].qualite=FEUILLE;
       feuilles[i].distances[1]=0;
       feuilles[i].distances[2]=0;
       feuilles[i].dissimilarite=0;
       feuilles[i].individu=individus[i];
       noeuds[i]=&(feuilles[i]);
     }
  
   /*debut de la classification hierarchique proprement dite*/
   graine=0;
   for(nb_racines=nb_individus;nb_racines>1;nb_racines--)
    {
      dissimilarite_minimum=VALEUR_ENORME;
      dissimilarite_minimum_paria=VALEUR_ENORME;
      if(graine==0)
	{
	  plus_proche_voisin=1;
	}
      else
	{
	  plus_proche_voisin=0;
	}
      dissimilarite_minimum=dissimilarites_courantes[graine][plus_proche_voisin];
      if(nb_racines==2)
	{
	  meilleur_regroupement.i=0;
	  meilleur_regroupement.j=1;
	}
      else
	{
	  for(i=0;i<nb_racines;i++)
	    {
	      if((i!=graine)&&(dissimilarites_courantes[graine][i]<dissimilarite_minimum))
		{
		  plus_proche_voisin=i;
		  dissimilarite_minimum=dissimilarites_courantes[graine][i];
		}
	    }
	  ancienne_graine=graine;
	  graine=plus_proche_voisin;
	  if(graine==0)
	    {
	      plus_proche_voisin=1;
	    }
	  else
	    {
	      plus_proche_voisin=0;
	    }
	  dissimilarite_minimum=dissimilarites_courantes[plus_proche_voisin][graine];

	  while(1)
	    {
	      for(i=0;i<nb_racines;i++)
		{
		  if((i!=graine)&&(dissimilarites_courantes[graine][i]<dissimilarite_minimum))
		    {
		      plus_proche_voisin=i;
		      dissimilarite_minimum=dissimilarites_courantes[graine][i];
		    }
		}
	      if(plus_proche_voisin==ancienne_graine)
		{
		  break;
		}
	      else
		{
		  ancienne_graine=graine;
		  graine=plus_proche_voisin;
		  if(graine==0)
		    {
		      plus_proche_voisin=1;
		    }
		  else
		    {
		      plus_proche_voisin=0;
		    }
		  dissimilarite_minimum=dissimilarites_courantes[plus_proche_voisin][graine];
		}
	    }
	
	  if(graine<plus_proche_voisin)
	    {
	      meilleur_regroupement.i=graine;
	      meilleur_regroupement.j=plus_proche_voisin;
	    }
	  else
	    {
	      meilleur_regroupement.i=plus_proche_voisin;
	      meilleur_regroupement.j=graine;
	      graine=plus_proche_voisin;
	    }
	}

      historique_dissimilarites[nb_racines-2]=dissimilarite_minimum;
	
      /*construction de l'arbre*/
      nouveau_noeud=(noeud_t *)malloc(sizeof(noeud_t));
     
      nouveau_noeud->copain1=noeuds[meilleur_regroupement.i];
      nouveau_noeud->copain2=noeuds[meilleur_regroupement.j];
      nouveau_noeud->numero=compteur_numeros;
      compteur_numeros++;

      nouveau_noeud->qualite=NOEUD_ORDINAIRE;
      nouveau_noeud->dissimilarite=dissimilarite_minimum;
      nouveau_noeud->distances[1]=
	dissimilarite_minimum-noeuds[meilleur_regroupement.i]->dissimilarite;
      nouveau_noeud->distances[2]=
	dissimilarite_minimum-noeuds[meilleur_regroupement.j]->dissimilarite;
      noeuds[meilleur_regroupement.i]->copains[0]=nouveau_noeud;
      noeuds[meilleur_regroupement.j]->copains[0]=nouveau_noeud;
      noeuds[meilleur_regroupement.i]=nouveau_noeud;
      for(i=meilleur_regroupement.j;i<nb_racines-1;i++)
	{
	  noeuds[i]=noeuds[i+1];
	}
      
      /*mise a jour des dissimilarites*/
      for(i=0;i<meilleur_regroupement.i;i++)
	{
	  dissimilarites_courantes[i][meilleur_regroupement.i]=
	    ((poids[meilleur_regroupement.i]+poids[i])*
	     dissimilarites_courantes[i][meilleur_regroupement.i]+
	     (poids[meilleur_regroupement.j]+poids[i])*
	     dissimilarites_courantes[i][meilleur_regroupement.j]-
	     (poids[i])*
	     dissimilarites_courantes[meilleur_regroupement.i][meilleur_regroupement.j])/
	    (poids[i]+poids[meilleur_regroupement.i]+poids[meilleur_regroupement.j]);
	  dissimilarites_courantes[meilleur_regroupement.i][i]=
	    dissimilarites_courantes[i][meilleur_regroupement.i];
	}
      for(i=meilleur_regroupement.i+1;i<meilleur_regroupement.j;i++)
	{
	  dissimilarites_courantes[i][meilleur_regroupement.i]=
	    ((poids[meilleur_regroupement.i]+poids[i])*
	     dissimilarites_courantes[i][meilleur_regroupement.i]+
	     (poids[meilleur_regroupement.j]+poids[i])*
	     dissimilarites_courantes[i][meilleur_regroupement.j]-
	     (poids[i])*
	     dissimilarites_courantes[meilleur_regroupement.i][meilleur_regroupement.j])/
	    (poids[i]+poids[meilleur_regroupement.i]+
	     poids[meilleur_regroupement.j]);
	  dissimilarites_courantes[meilleur_regroupement.i][i]=
	    dissimilarites_courantes[i][meilleur_regroupement.i];
	}
      for(i=meilleur_regroupement.j+1;i<nb_racines;i++)
	{
	  dissimilarites_courantes[i][meilleur_regroupement.i]=
	    ((poids[meilleur_regroupement.i]+poids[i])*
	     dissimilarites_courantes[i][meilleur_regroupement.i]+
	     (poids[meilleur_regroupement.j]+poids[i])*
	     dissimilarites_courantes[i][meilleur_regroupement.j]-
	     (poids[i])*
	     dissimilarites_courantes[meilleur_regroupement.i][meilleur_regroupement.j])/
	    (poids[i]+poids[meilleur_regroupement.i]+
	     poids[meilleur_regroupement.j]);
	  dissimilarites_courantes[meilleur_regroupement.i][i]=
	    dissimilarites_courantes[i][meilleur_regroupement.i];
	}
      dissimilarites_courantes[meilleur_regroupement.i][meilleur_regroupement.i]=0;
      for(i=meilleur_regroupement.j;i<nb_racines-1;i++)
	{
	  for(j=0;j<nb_racines;j++)
	    {
	      dissimilarites_courantes[i][j]=dissimilarites_courantes[i+1][j];
	    }
	}
      for(i=meilleur_regroupement.j;i<nb_racines-1;i++)
	{
	  for(j=0;j<nb_racines-1;j++)
	    {
	      dissimilarites_courantes[j][i]=dissimilarites_courantes[j][i+1];
	    }
	}
      
      /*fin de la mise a jour des dissimilarites*/
      poids[meilleur_regroupement.i]+=poids[meilleur_regroupement.j];
      for(i=meilleur_regroupement.j;i<nb_racines-1;i++)
	{
	  poids[i]=poids[i+1];
	}
    }
   /*fin de la classification hierarchique*/

   *racine=noeuds[0];

   /*desallocation de la memoire*/
   for(i=0;i<nb_individus;i++)
     {
       free(dissimilarites_courantes[i]);
     }
   free(dissimilarites_courantes);
  
   free(poids);
   free(noeuds);
   free(historique_dissimilarites);
   /*fin desallocation de la memoire*/
}
