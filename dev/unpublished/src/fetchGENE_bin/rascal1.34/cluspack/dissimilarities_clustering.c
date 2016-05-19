#include "principal.h"
#include "tools.h"
#include "dissimilarities_clustering.h"

/*************************************************************/
/*                                                           */
/*classification des dissimilarites pour trouver le nombre de*/
/*dissimilarites elevees et donc le nombre de groupes        */
/*(avec l'exposant de Holder)                                */
/*                                                           */
/*************************************************************/

void clustering_dissimilarites_Holder(int nb_dissimilarites,double *historique_dissimilarites,
				      int *nb_clusters,int resolution)
{
  /*declaration de variables*/
  int i,j,nb_exposants_Holder,nb_divisions,**debuts_intervalles,nb_dimensions_fractales_faibles;
  int au_moins_un_exposant_Holder_unitaire,nb_exposants_Holder_faibles,nb_coupures_candidates;
  int *coupures_candidates,coupure_selectionnee,coupure_par_defaut;
  int nb_dissimilarites_elevees;
  double *exposants_Holder,*exposants_Holder_ordonnes;
  /*fin declaration de variables*/

  /*allocation memoire*/
  exposants_Holder=(double *)malloc(sizeof(double)*nb_dissimilarites);
  exposants_Holder_ordonnes=(double *)malloc(sizeof(double)*nb_dissimilarites);
  coupures_candidates=(int *)malloc(sizeof(int)*nb_dissimilarites);
  /*fin allocation memoire*/

  calcul_exposants_Holder(nb_dissimilarites,historique_dissimilarites,
			  &nb_exposants_Holder,exposants_Holder);

  au_moins_un_exposant_Holder_unitaire=NON;
  for(i=0;i<nb_exposants_Holder;i++)
    {
      if(exposants_Holder[i]>EXPOSANT_HOLDER_SEUIL)
	{
	  au_moins_un_exposant_Holder_unitaire=OUI;
	  break;
	}
    }
  
  if(au_moins_un_exposant_Holder_unitaire==OUI)
    {
      /*calcul du nombre de divisions pour la dimension fractale*/
      nb_divisions=(int)floor(log((double)nb_exposants_Holder)/log(2.0));
      
      /*allocation memoire*/
      debuts_intervalles=(int **)malloc(sizeof(int *)*nb_divisions);
      for(i=0;i<nb_divisions;i++)
	{
	  debuts_intervalles[i]=(int *)malloc(sizeof(int)*((int)pow(2.0,(double)i)+1));
	}
      /*fin allocation memoire*/
      
      /*clustering proprement dit des dissimilarites*/
      separation_dissimilarites_evolue(nb_exposants_Holder,nb_dissimilarites,
				       historique_dissimilarites,&nb_dissimilarites_elevees,
				       exposants_Holder);      
      
      nb_coupures_candidates=0;
      coupures_candidates[0]=0;

      if(resolution!=0)
	{
	  for(i=0;i<nb_exposants_Holder-1;i++)
	    {
	      if((exposants_Holder[i]<=EXPOSANT_HOLDER_SEUIL)&&
		 (exposants_Holder[i+1]>EXPOSANT_HOLDER_SEUIL))
		{
		  coupures_candidates[nb_coupures_candidates]=i;
		  if((i+1)==nb_dissimilarites_elevees)
		    {
		      coupure_par_defaut=nb_coupures_candidates;
		    }
		  nb_coupures_candidates++;
		}
	    }

	  printf("coupure_par_defaut : %d\n",coupure_par_defaut);

	  coupure_selectionnee=resolution+coupure_par_defaut;
	  if(coupure_selectionnee<0)
	    {
	      nb_dissimilarites_elevees=coupures_candidates[0]+1;
	    }
	  else if(coupure_selectionnee>=nb_coupures_candidates)
	    {
	      nb_dissimilarites_elevees=coupures_candidates[nb_coupures_candidates-1]+1;
	    }
	  else
	    {
	      nb_dissimilarites_elevees=coupures_candidates[coupure_selectionnee]+1;
	    }
	}
    }
  else
    {
      /*on separe les valeurs elevees des valeurs faibles*/
      separation_dissimilarites_bete(nb_dissimilarites,historique_dissimilarites,
				     &nb_dissimilarites_elevees);
    }
  
  *nb_clusters=nb_dissimilarites_elevees+1;

  /*desallocation memoire*/
  free(exposants_Holder);
  free(exposants_Holder_ordonnes);
  free(coupures_candidates);
  /*fin desallocation memoire*/
}

/******************************************************************************/
/*                                                                            */
/*Procedure de calcul des exposants d'Holder sur les valeurs de dissimilarites*/
/*                                                                            */
/******************************************************************************/

void calcul_exposants_Holder(int nb_dissimilarites,double *historique_dissimilarites,
			     int *nb_exposants_Holder,double *exposants_Holder)
{
  /*declaration des variables*/
  int i,j,nb_mesures,compteur_exposants_Holder,nb_exposants_Holder_faibles;
  double **valeurs_pour_regression,mesure_max;
  /*fin declaration des variables*/

  /*allocation memoire*/
  valeurs_pour_regression=(double **)malloc(sizeof(double *)*(nb_dissimilarites-1));
  for(i=0;i<nb_dissimilarites-1;i++)
    {
      valeurs_pour_regression[i]=(double *)malloc(sizeof(double)*2);
    }
  /*fin allocation memoire*/

  mesure_max=historique_dissimilarites[0]-historique_dissimilarites[nb_dissimilarites-1];
  
  compteur_exposants_Holder=0;
  for(i=0;i<nb_dissimilarites-1;i++)
    {
      nb_mesures=(int)floor(log((double)(nb_dissimilarites-i))/log(2.0));
      for(j=nb_mesures;j>0;j--)
	{
	  if(((historique_dissimilarites[i]-
	      historique_dissimilarites[i+(int)pow(2.0,(double)j)-1])/mesure_max)<ZERO_LIMIT)
	    {
	      nb_mesures--;
	    }
	  else
	    {
	      valeurs_pour_regression[nb_mesures-j][1]=
		log((double)(historique_dissimilarites[i]-
		     historique_dissimilarites[i+(int)pow(2.0,(double)j)-1])/mesure_max);
	      valeurs_pour_regression[nb_mesures-j][0]=-log(2.0)*(double)(nb_mesures-j);
	    }
	}

      if(nb_mesures>1)
	{
	  exposants_Holder[compteur_exposants_Holder]=
	    regression_lineaire(nb_mesures,valeurs_pour_regression);
	  if(exposants_Holder[compteur_exposants_Holder]>1.0)
	    {
	      exposants_Holder[compteur_exposants_Holder]=1.0;
	    }
	  compteur_exposants_Holder++;
	}
    }

  *nb_exposants_Holder=compteur_exposants_Holder;

  /*desallocation memoire*/
  for(i=0;i<nb_dissimilarites-1;i++)
    {
      free(valeurs_pour_regression[i]);
    }
  free(valeurs_pour_regression);
  /*fin desallocation memoire*/
}

/**********************************************************************/
/*                                                                    */
/*Procedure de separation en deux groupes des valeurs de dissimilarite*/
/*par calcul du minimum d'inertie intraclasse                         */
/*                                                                    */
/**********************************************************************/

void separation_dissimilarites_bete(int nb_dissimilarites,double *historique_dissimilarites,
				    int *nb_dissimilarites_premier_groupe)
{
  /*declaration des variables*/
  int i,j,nb_dissimilarites_valeur_minimum,nb_dissimilarites1,nb_dissimilarites2;
  double valeur_courante,valeur_minimum,centre_gravite_tot,centre_gravite1,centre_gravite2;
  /*fin declaration des variables*/

  valeur_minimum=-1;
  nb_dissimilarites_valeur_minimum=1;
  for(i=1;i<nb_dissimilarites;i++)
    {
      nb_dissimilarites1=i;
      nb_dissimilarites2=nb_dissimilarites-nb_dissimilarites1;
      
      centre_gravite1=0;
      for(j=0;j<nb_dissimilarites1;j++)
	{
	  centre_gravite1+=historique_dissimilarites[j];
	}
      centre_gravite1/=(double)nb_dissimilarites1;
    
      centre_gravite2=0;
      for(j=nb_dissimilarites1;j<nb_dissimilarites;j++)
	{
	  centre_gravite2+=historique_dissimilarites[j];
	}
      centre_gravite2/=(double)nb_dissimilarites2;

      valeur_courante=0;

      for(j=0;j<nb_dissimilarites1;j++)
	{
	  valeur_courante+=((double)(historique_dissimilarites[j]-centre_gravite1))*
	    ((double)(historique_dissimilarites[j]-centre_gravite1));
	}
      for(j=nb_dissimilarites1;j<nb_dissimilarites;j++)
	{
	  valeur_courante+=((double)(historique_dissimilarites[j]-centre_gravite2))*
	    ((double)(historique_dissimilarites[j]-centre_gravite2));
	}

      if((valeur_minimum<0)||(valeur_courante<valeur_minimum))
	{
	  valeur_minimum=valeur_courante;
	  nb_dissimilarites_valeur_minimum=nb_dissimilarites1;
	}
    }

  *nb_dissimilarites_premier_groupe=nb_dissimilarites_valeur_minimum;
}

/************************************************************************/
/*                                                                      */
/*Procedure de separation en deux groupes des valeurs de dissimilarite  */
/*par calcul du minimum d'inertie intraclasse mais uniquement aux points*/
/*d'exposants d'Holder faibles                                          */
/*                                                                      */
/************************************************************************/

void separation_dissimilarites_evolue(int nb_exposants_Holder,int nb_dissimilarites,
				      double *historique_dissimilarites,
				      int *nb_dissimilarites_premier_groupe,
				      double *exposants_Holder)
{
  /*declaration des variables*/
  int i,j,nb_dissimilarites_valeur_minimum,nb_dissimilarites1,nb_dissimilarites2,compteur;
  double valeur_courante,valeur_minimum,centre_gravite_tot,centre_gravite1,centre_gravite2;
  double *inerties_intra_classe;
  /*fin declaration des variables*/

  /*allocation memoire*/
  inerties_intra_classe=(double *)malloc(sizeof(double)*nb_dissimilarites);
  /*fin allocation memoire*/

  valeur_minimum=-1;
  nb_dissimilarites_valeur_minimum=1;
  compteur=0;
  for(i=1;i<nb_exposants_Holder;i++)
    {
      if(i>nb_dissimilarites/2)
	{
	  break;
	}
      if((exposants_Holder[i-1]<=EXPOSANT_HOLDER_SEUIL)&&
	 (exposants_Holder[i]>EXPOSANT_HOLDER_SEUIL))
	{
	  nb_dissimilarites1=i;
	  nb_dissimilarites2=nb_dissimilarites-nb_dissimilarites1;
	  
	  centre_gravite1=0;
	  for(j=0;j<nb_dissimilarites1;j++)
	    {
	      centre_gravite1+=historique_dissimilarites[j];
	    }
	  centre_gravite1/=(double)nb_dissimilarites1;
	  
	  centre_gravite2=0;
	  for(j=nb_dissimilarites1;j<nb_dissimilarites;j++)
	    {
	      centre_gravite2+=historique_dissimilarites[j];
	    }
	  centre_gravite2/=(double)nb_dissimilarites2;
	  
	  valeur_courante=0;
	  
	  for(j=0;j<nb_dissimilarites1;j++)
	    {
	      valeur_courante+=((double)(historique_dissimilarites[j]-centre_gravite1))*
		((double)(historique_dissimilarites[j]-centre_gravite1));
	    }
	  for(j=nb_dissimilarites1;j<nb_dissimilarites;j++)
	    {
	      valeur_courante+=((double)(historique_dissimilarites[j]-centre_gravite2))*
		((double)(historique_dissimilarites[j]-centre_gravite2));
	    }
	  
	  if((valeur_minimum<0)||(valeur_courante<valeur_minimum))
	    {
	      valeur_minimum=valeur_courante;
	      nb_dissimilarites_valeur_minimum=nb_dissimilarites1;
	    }
	  /*printf("inertie intra classe en %d : %.2f\n",i,valeur_courante);*/
	  inerties_intra_classe[i]=valeur_courante;
	  compteur++;
	}
    }

  *nb_dissimilarites_premier_groupe=nb_dissimilarites_valeur_minimum;

  /*desallocation memoire*/
  free(inerties_intra_classe);
  /*fin desallocation memoire*/
}


