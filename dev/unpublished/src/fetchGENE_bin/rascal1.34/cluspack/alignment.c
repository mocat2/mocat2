#include "principal.h"

/*********************************************************************/
/*                                                                   */
/*Procedure de calcul des pourcentages d'identite entre les sequences*/
/*                                                                   */
/*********************************************************************/

void calcul_identites_sequences(int nb_individus,individu_t **individus,
				int longueur_alignement,char **sequences)
{
  /*declaration des variables*/
  int i,j,k,nb_identites,nb_positions_communes;
  /*fin declaration des variables*/

  for(i=0;i<nb_individus;i++)
    {
      individus[i]->valeurs_traitees[i]=1.0;
      for(j=i+1;j<nb_individus;j++)
	{
	  nb_positions_communes=0;
	  nb_identites=0;
	  for(k=0;k<longueur_alignement;k++)
	    {
	       if((sequences[i][k]!='O')&&(sequences[j][k]!='O'))
		{
		  nb_positions_communes++;
		  if(sequences[i][k]==sequences[j][k])
		    {
		      nb_identites++;
		    }
		}
	    }
	  if(nb_positions_communes==0)
	    {
	      individus[i]->valeurs_traitees[j]=0.0;
	    }
	  else
	    {
	      individus[i]->valeurs_traitees[j]=
		(double)nb_identites/(double)nb_positions_communes;
	    }
	  individus[j]->valeurs_traitees[i]=individus[i]->valeurs_traitees[j];
	}
    }

  /*printf("%d\n",nb_individus);
    for(i=0;i<nb_individus;i++)
    {   
    printf("%s ",individus[i]->nom);
    for(j=0;j<nb_individus;j++)
    {
    printf("%.2f ",individus[i]->valeurs_traitees[j]);
    }
    printf("\n");
    }  
    exit(0);*/
}

/*********************************************************************/
/*                                                                   */
/*Procedure de calcul des pourcentages d'identite entre les sequences*/
/*                                                                   */
/*********************************************************************/

void calcul_distances_sequences(int nb_individus,individu_t **individus,
				int longueur_alignement,char **sequences)
{
  /*declaration des variables*/
  int i,j,k,nb_identites,nb_positions_communes;
  /*fin declaration des variables*/

  for(i=0;i<nb_individus;i++)
    {
      individus[i]->valeurs_traitees[i]=0.0;
      for(j=i+1;j<nb_individus;j++)
	{
	  nb_positions_communes=0;
	  nb_identites=0;
	  for(k=0;k<longueur_alignement;k++)
	    {
	       if((sequences[i][k]!='O')&&(sequences[j][k]!='O'))
		{
		  nb_positions_communes++;
		  if(sequences[i][k]==sequences[j][k])
		    {
		      nb_identites++;
		    }
		}
	    }
	  if(nb_positions_communes==0)
	    {
	      individus[i]->valeurs_traitees[j]=1.0;
	    }
	  else
	    {
	      individus[i]->valeurs_traitees[j]=
		1.0-(double)nb_identites/(double)nb_positions_communes;
	    }
	  individus[j]->valeurs_traitees[i]=individus[i]->valeurs_traitees[j];
	}
    }
}
