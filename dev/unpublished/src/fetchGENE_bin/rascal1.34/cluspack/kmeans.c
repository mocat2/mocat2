#include "principal.h"
#include "kmeans.h"
#include "tools.h"

/****************************************************/
/*                                                  */
/*Procedure de realisation de clustering par k-means*/
/*en partant d'un debut de clustering               */ 
/*                                                  */
/****************************************************/

int kmeans(int nb_individus,int nb_dimensions,individu_t **individus,int nb_classes)
{
  /*declaration des variables*/
  int i,j,*nb_individus_par_classe,*solution,compteur_iterations1,compteur_iterations2;
  int nb_kmeans,nb_iterations_max,*points_selectionnes,resultat;
  double **centres,inertie_intra_classe_precedente,inertie_intra_classe_courante;
  double meilleure_inertie_intra_classe;
  /*fin eclaration des variables*/

  /*allocation memoire*/
  centres=(double **)malloc(sizeof(double *)*nb_classes);
  for(i=0;i<nb_classes;i++)
    {
      centres[i]=(double *)malloc(sizeof(double)*nb_dimensions);
    }
  nb_individus_par_classe=(int *)malloc(sizeof(int)*nb_classes);
  solution=(int *)malloc(sizeof(int)*nb_individus);
  points_selectionnes=(int *)malloc(sizeof(int)*nb_individus);
  /*fin allocation memoire*/

  nb_kmeans=3;
  nb_iterations_max=20;
  meilleure_inertie_intra_classe=NOT_DEFINED;
  
  for(compteur_iterations1=0;compteur_iterations1<nb_kmeans;compteur_iterations1++)
    {
      for(i=0;i<nb_classes;i++)
	{
	  points_selectionnes[i]=-1;
	  while(points_selectionnes[i]==-1)
	    {
	      points_selectionnes[i]=
		(int)floor((double)nb_individus-1.0)*((double)rand()/(double)RAND_MAX);
	      for(j=0;j<i;j++)
		{
		  if(points_selectionnes[i]==points_selectionnes[j])
		    {
		      points_selectionnes[i]=-1;
		      break;
		    }
		}
	    }
	  
	  for(j=0;j<nb_dimensions;j++)
	    {
	      centres[i][j]=individus[points_selectionnes[i]]->valeurs_traitees[j];
	    }
	}      
      /*on affecte les individus a leurs groupes les plus proches*/
      resultat=affectation_individuskmeans(nb_dimensions,nb_individus,nb_classes,nb_individus_par_classe,
					   individus,centres,&inertie_intra_classe_precedente);

      if(resultat==ECHEC)
	{
	  break;
	}

      compteur_iterations2=0;
      while(1)
	{
	  
	  /*on affecte les individus a leurs groupes les plus proches*/
	  resultat=affectation_individuskmeans(nb_dimensions,nb_individus,nb_classes,
					       nb_individus_par_classe,individus,centres,
					       &inertie_intra_classe_courante);
	  
	  if(resultat==ECHEC)
	    {
	      break;
	    }
	  
	  if((inertie_intra_classe_courante-inertie_intra_classe_precedente<ZERO_LIMIT)&&
	     (-inertie_intra_classe_courante+inertie_intra_classe_precedente<ZERO_LIMIT))
	    {
	      break;
	    }
	  else
	    {
	      inertie_intra_classe_precedente=inertie_intra_classe_courante;
	    }

	  compteur_iterations2++;
	  
	  if(compteur_iterations2>nb_iterations_max)
	    {
	      break;
	    }
	}

      if((resultat!=ECHEC)&&((meilleure_inertie_intra_classe=NOT_DEFINED)||
	 (inertie_intra_classe_courante<meilleure_inertie_intra_classe)))
	{
	  for(i=0;i<nb_individus;i++)
	    {
	      solution[i]=individus[i]->cluster;
	    }
	  meilleure_inertie_intra_classe=inertie_intra_classe_courante;
	}
    }
  for(i=0;i<nb_individus;i++)
    {
      individus[i]->cluster=solution[i];
    }
  
  /*desallocation memoire*/
  for(i=0;i<nb_classes;i++)
    {
      free(centres[i]);
    }
  free(centres);
  free(nb_individus_par_classe);
  free(solution);
  free(points_selectionnes);
  /*fin desallocation memoire*/

  return resultat;
}

/***************************************************/
/*                                                 */
/*Procedure d'affectation des individus aux groupes*/
/*et de calcul des nouveaux centres de gravite     */
/*                                                 */
/***************************************************/

int affectation_individuskmeans(int nb_dimensions,int nb_individus,int nb_classes,
				int *nb_individus_par_classe,individu_t **individus,
				double **centres,double *inertie_intra_classe)
{
  /*declaration des variables*/
  int i,j,groupe_le_plus_proche;
  double distance_courante,distance_min;
  /*fin declaration des variables*/

  for(i=0;i<nb_classes;i++)
    {
      nb_individus_par_classe[i]=0;
    }

  for(i=0;i<nb_individus;i++)
    {
      distance_min=calcul_distance_carre(nb_dimensions,centres[0],individus[i]->valeurs_traitees);
      groupe_le_plus_proche=0;

      for(j=1;j<nb_classes;j++)
	{
	  distance_courante=calcul_distance_carre(nb_dimensions,centres[j],
						  individus[i]->valeurs_traitees);
	  if(distance_courante<distance_min)
	    {
	      distance_min=distance_courante;
	      groupe_le_plus_proche=j;
	    }
	}
      individus[i]->cluster=groupe_le_plus_proche;
      (nb_individus_par_classe[groupe_le_plus_proche])++;
    }

  for(i=0;i<nb_classes;i++)
    {
      for(j=0;j<nb_dimensions;j++)
	{
	  centres[i][j]=0;
	}
    }

  for(i=0;i<nb_individus;i++)
    {
      for(j=0;j<nb_dimensions;j++)
	{
	  centres[individus[i]->cluster][j]+=individus[i]->valeurs_traitees[j];
	}
    }
  
  for(i=0;i<nb_classes;i++)
    {
      for(j=0;j<nb_dimensions;j++)
	{
	  if(nb_individus_par_classe[i]==0)
	    {
	      return ECHEC;
	    }
	  centres[i][j]/=(double)nb_individus_par_classe[i];
	}
    }
  *inertie_intra_classe=0;
  for(i=0;i<nb_individus;i++)
    {
      for(j=0;j<nb_dimensions;j++)
	{
	  *inertie_intra_classe+=(individus[i]->valeurs_traitees[j]-centres[individus[i]->cluster][j])*
	    (individus[i]->valeurs_traitees[j]-centres[individus[i]->cluster][j]);
	}
    }  

  return SUCCES;
}

