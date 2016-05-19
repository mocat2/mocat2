#include "principal.h"
#include "tools.h"

/*************************************************************/
/*                                                           */
/*Procedure de calcul des exposants d'Holder                 */
/*Methode : on forme un hyperpave avec pour longueur         */
/*la distances entre les deux coordonnees                    */
/*les plus eloignees sur chaque axe                          */
/*                                                           */
/*************************************************************/

void calcul_exposants_Holder_points(int nb_dimensions,int nb_individus_total,
				    double **coordonnees_tous_individus,int nb_individus,
				    double **coordonnees_individus,double *exposants_Holder,
				    double *cotes_hyperpave,double *valeurs_minimum,
				    int mode_verbose,int test)
{
  /*declaration de variables*/
  int i,j,k,l,test_individu_dans_le_meme_hyperpave,nb_mesures,nb_mesures_max;
  int *individus_a_tester,nb_individus_a_tester,compteur;
  double *bornes_inf,*bornes_sup,indice_segment_axe,nb_copains;
  /*fin declaration de variables*/

  nb_mesures_max=20;

  /*allocation memoire*/
  bornes_inf=(double *)malloc(sizeof(double)*nb_dimensions);
  bornes_sup=(double *)malloc(sizeof(double)*nb_dimensions);
  individus_a_tester=(int *)malloc(sizeof(int)*nb_individus_total);
  /*fin allocation memoire*/
  
  for(i=0;i<nb_individus;i++)
    {
      nb_mesures=1;
      for(j=0;j<nb_individus_total;j++)
	{
	  individus_a_tester[j]=j;
	}
      nb_individus_a_tester=nb_individus_total;
      for(j=0;j<nb_mesures_max;j++)      
	{
	  for(k=0;k<nb_dimensions;k++)
	    {
	      indice_segment_axe=floor((double)(coordonnees_individus[i][k]-valeurs_minimum[k])/
				       ((double)cotes_hyperpave[k]/pow(2.0,(double)(j+1))));
	   if(indice_segment_axe>pow(2.0,(double)(j+1))-ZERO_LIMIT)
		{
		  indice_segment_axe-=1.0;
		}
	      bornes_inf[k]=indice_segment_axe*cotes_hyperpave[k]/pow(2.0,(double)(j+1))+
		valeurs_minimum[k];
	      bornes_sup[k]=(1.0+indice_segment_axe)*cotes_hyperpave[k]/pow(2.0,(double)(j+1))+
		valeurs_minimum[k];   
	    }
	  nb_copains=0;

	  compteur=0;

	  for(k=0;k<nb_individus_a_tester;k++)
	    {
	      test_individu_dans_le_meme_hyperpave=OUI;
	      for(l=0;l<nb_dimensions;l++)
		{
		  if((coordonnees_tous_individus[individus_a_tester[k]][l]<
		      bornes_inf[l]-ZERO_LIMIT)||
		     (coordonnees_tous_individus[individus_a_tester[k]][l]>
		      bornes_sup[l]+ZERO_LIMIT))
		    {
		      test_individu_dans_le_meme_hyperpave=NON;
		      break;
		    }
		}
	      if(test_individu_dans_le_meme_hyperpave==OUI)
		{
		  nb_copains+=1.0;
		  if(nb_copains>=2.0)
		    {
		      for(l=k;l<nb_individus_a_tester;l++)
			{
			  individus_a_tester[compteur]=individus_a_tester[l];
			  compteur++;
			}
		      break;
		    }
		  else
		    {
		      individus_a_tester[compteur]=individus_a_tester[k];
		      compteur++;
		    }
		}
	    }
	  nb_individus_a_tester=compteur;

	  if(mode_verbose==OUI)
	    {
	      printf("nb_copains : %.2f\n",nb_copains);
	    }	 

	  if(nb_copains<1.0)
	    {
	      break;
	    }
	  else if(nb_copains<2.0)
	    {
	      nb_mesures++;
	      break;
	    }
	  else
	    {
	      nb_mesures++;
	    }
	}

      /*      printf("nb_mesures : %d\n",nb_mesures);*/
      exposants_Holder[i]=(double)nb_mesures;
     
    }

  /*desallocation memoire*/
  free(bornes_inf);
  free(bornes_sup);
  free(individus_a_tester);
  /*fin desallocation memoire*/
}


