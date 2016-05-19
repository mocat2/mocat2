#include "principal.h"

/********************************************/
/*                                          */
/*Procedure de standardisation des individus*/
/*                                          */
/********************************************/

void standardisation_individus(int nb_individus,int nb_dimensions,individu_t **individus)
{
  /*declaration de variables*/
  int i,j;
  double moyenne,ecart_type;
  /*fin declaration de variables*/
  
  for(i=0;i<nb_individus;i++)
    {
      moyenne=0;
      for(j=0;j<nb_dimensions;j++)
	{
	  moyenne+=individus[i]->valeurs_traitees[j];
	}
      moyenne/=(double)nb_dimensions;

      ecart_type=0;
      for(j=0;j<nb_dimensions;j++)
	{
	  ecart_type+=(individus[i]->valeurs_traitees[j]-moyenne)*(individus[i]->valeurs_traitees[j]-moyenne);
	}
      ecart_type=sqrt(ecart_type/(double)nb_dimensions);

      if(ecart_type<ZERO_LIMIT)
	{
	  printf("Standard deviation of %s is too low : %f\n",individus[i]->nom,ecart_type);
	  exit(1);
	}

      for(j=0;j<nb_dimensions;j++)
	{
	  individus[i]->valeurs_traitees[j]=(individus[i]->valeurs_traitees[j]-moyenne)/ecart_type;
	}
    }
}


/**************************************************************/
/*                                                            */
/*Procedure de calcul des similarites a partir des coordonnees*/
/*en utilisant les correlations                               */
/*                                                            */
/**************************************************************/

void coordonnees_to_similarites(int nb_individus,int nb_dimensions,individu_t **individus)
{
  /*declaration de variables*/
  int i,j,k;
  double **correlations,*ecarts_type,*moyennes;
  /*fin declaration de variables*/
  
  /*allocation memoire*/
  correlations=(double **)malloc(sizeof(double *)*nb_individus);
  for(i=0;i<nb_individus;i++)
    {
      correlations[i]=(double *)malloc(sizeof(double)*nb_individus);
    }
  moyennes=(double *)malloc(sizeof(double)*nb_individus);
  ecarts_type=(double *)malloc(sizeof(double)*nb_individus);
  /*fin allocation memoire*/
  
  for(i=0;i<nb_individus;i++)
    {
      moyennes[i]=0;
      for(j=0;j<nb_dimensions;j++)
	{
	  moyennes[i]+=individus[i]->valeurs_traitees[j];
	}
      moyennes[i]/=(double)nb_dimensions;
    }

  for(i=0;i<nb_individus;i++)
    {
      ecarts_type[i]=0;
      for(j=0;j<nb_dimensions;j++)
	{
	  ecarts_type[i]+=(individus[i]->valeurs_traitees[j]-moyennes[i])*
	    (individus[i]->valeurs_traitees[j]-moyennes[i]);
	}
      ecarts_type[i]=sqrt(ecarts_type[i]/(double)(nb_dimensions-1));
    }

  for(i=0;i<nb_individus;i++)
    {
      correlations[i][i]=1.0;
      for(j=i+1;j<nb_individus;j++)
	{
	  correlations[i][j]=0;
	  for(k=0;k<nb_dimensions;k++)
	    {
	      correlations[i][j]+=(individus[i]->valeurs_traitees[k]-moyennes[i])*
		(individus[j]->valeurs_traitees[k]-moyennes[j]);
	    }
	  correlations[i][j]/=(double)nb_dimensions;
	  correlations[j][i]=correlations[i][j];
	}
    }

  for(i=0;i<nb_individus;i++)
    {
      individus[i]->valeurs_traitees=
	(double *)realloc(individus[i]->valeurs_traitees,sizeof(double)*nb_individus);
      for(j=0;j<nb_individus;j++)
	{
	  individus[i]->valeurs_traitees[j]=correlations[i][j];
	  if(individus[i]->valeurs_traitees[j]<0)
	    {
	      individus[i]->valeurs_traitees[j]=0.0;
	    } 	  
	}
    }

  /*allocation memoire*/
  for(i=0;i<nb_individus;i++)
    {
      free(correlations[i]);
    }
  free(correlations);
  free(ecarts_type);
  free(moyennes);
  /*fin allocation memoire*/
  
}

/*********************************************/
/*                                           */
/*Procedure de filtrage des individus suivant*/ 
/*une valeur seuil et le type des donnees    */
/*                                           */
/*********************************************/

void filtre_individus(int *nb_individus,int nb_dimensions,individu_t **individus,
		      int type_donnees,double distance_min,int nb_clusters_selection)
{
  /*declaration de variables*/
  int i,j,k,*individus_elimines,nb_individus_elimines,compteur;
  double distance;
  /*fin declaration de variables*/

  /*allocation memoire*/
  individus_elimines=(int *)malloc(sizeof(int)*(*nb_individus));
  /*fin allocation memoire*/

  nb_individus_elimines=0;
 
  for(i=0;i<*nb_individus;i++)
    {
      individus_elimines[i]=NON;
      individus[i]->nb_individus_similaires=0;
      individus[i]->individus_similaires=NULL;
    }
  for(i=0;i<*nb_individus;i++)
    {
      for(j=i+1;j<*nb_individus;j++)
	{
	  /*COORDONNEES*/
	  if(type_donnees==COORDONNEES)
	    {
	      distance=0;
	      for(k=0;k<nb_dimensions;k++)
		{
		  distance+=(individus[i]->valeurs_traitees[k]-
			     individus[j]->valeurs_traitees[k])*
		    (individus[i]->valeurs_traitees[k]-
		     individus[j]->valeurs_traitees[k]);
		}
	      distance=sqrt(distance);
	    }
	  /*ALIGNEMENT*/
	  else if(type_donnees==ALIGNEMENT)
	    {
	      if(nb_clusters_selection==SECATOR)
		{
		  distance=individus[i]->valeurs_traitees[j];
		}
	      else
		{
		  distance=1.0-individus[i]->valeurs_traitees[j];
		}
	    }
	  /*DISTANCES*/
	  else if(type_donnees==DISTANCES)
	    {
	      distance=individus[i]->valeurs_traitees[j];
	    }
	  if((distance<distance_min)&&(individus_elimines[i]==NON)&&
	     (individus_elimines[j]==NON))
	    {
	      nb_individus_elimines++;
	      individus_elimines[j]=OUI;
	      individus[i]->nb_individus_similaires++;
	      individus[i]->individus_similaires=
		(individu_t **)realloc(individus[i]->individus_similaires,sizeof(individu_t *)
				       *individus[i]->nb_individus_similaires);
	      
	      individus[i]->individus_similaires[individus[i]->nb_individus_similaires-1]=
		individus[j];
	    }
	}
    }

  compteur=0;
  for(i=0;i<*nb_individus;i++)
    {
      if(individus_elimines[i]==NON)
	{
	  individus[compteur]=individus[i];
	  compteur++;
	}
    }
  *nb_individus-=nb_individus_elimines;

  /*desallocation memoire*/
  free(individus_elimines);
  /*fin desallocation memoire*/
}


