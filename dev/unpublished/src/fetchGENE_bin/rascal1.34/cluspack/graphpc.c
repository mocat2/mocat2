#include "principal.h"
#include "divide.h"
#include "tools.h"
#include "graphpc.h"

/**************************************************************/
/*                                                            */
/*Procedure de test pour savoir si l'on doit scinder un groupe*/
/*par test d'ajustement de la distribution des similarites    */ 
/*renvoie OUI si on doit separer et NON sinon                 */
/*                                                            */
/**************************************************************/

int test_GRAPHPC(int nb_individus,int nb_dimensions,groupe_t *groupe_pere,groupe_t *groupe1,
		 groupe_t *groupe2)
{
  /*declaration de variables*/
  int i,j,compteur,on_divise;
  double *probabilites_transition_1vers2,*probabilites_transition_2vers1;
  double somme,moyenne_changement_classe1,moyenne_changement_classe2,pourcentage;
  double probabilite_attendue_changement_classe1,probabilite_attendue_changement_classe2;
  /*fin declaration de variables*/
  
  pourcentage=0.5;

  /*allocation memoire*/
  probabilites_transition_1vers2=(double *)malloc(sizeof(double)*groupe1->nb_individus);
  probabilites_transition_2vers1=(double *)malloc(sizeof(double)*groupe2->nb_individus);
  /*fin allocation memoire*/

  compteur=0;
  on_divise=NON;

  /*on commence par les transitions partant de sommets du groupe1*/
  for(i=0;i<groupe1->nb_individus;i++)
    {
      somme=0;
      for(j=0;j<groupe1->nb_individus;j++)
	{
	  if(j!=i)
	    {
	      somme+=groupe1->individus[i]->valeurs_traitees[groupe1->individus[j]->id];
	    }
	}
      probabilites_transition_1vers2[i]=0;
      for(j=0;j<groupe2->nb_individus;j++)
	{
	  probabilites_transition_1vers2[i]+=
	    groupe1->individus[i]->valeurs_traitees[groupe2->individus[j]->id];
	  somme+=groupe1->individus[i]->valeurs_traitees[groupe2->individus[j]->id];
	}

      if(somme<ZERO_LIMIT)
	{
	  on_divise=OUI;
	  break;
	}
      probabilites_transition_1vers2[i]/=somme;
    }
  if(on_divise==NON)
    {
      /*on continue par les transitions partant de sommets du groupe2*/
      for(i=0;i<groupe2->nb_individus;i++)
	{
	  somme=0;
	  for(j=0;j<groupe2->nb_individus;j++)
	    {
	      if(j!=i)
		{
		  somme+=groupe2->individus[i]->valeurs_traitees[groupe2->individus[j]->id];
		}
	    }
	  probabilites_transition_2vers1[i]=0;
	  for(j=0;j<groupe1->nb_individus;j++)
	    {
	      probabilites_transition_2vers1[i]+=
		groupe2->individus[i]->valeurs_traitees[groupe1->individus[j]->id];
	      somme+=groupe2->individus[i]->valeurs_traitees[groupe1->individus[j]->id];
	    } 

	  if(somme<ZERO_LIMIT)
	    {
	      on_divise=OUI;
	      break;
	    }
	  probabilites_transition_2vers1[i]/=somme;
	}
    }
  
  if(on_divise==NON)
    {
      /*on trie les probabilites de changement de groupe par ordre decroissant*/
      tri_rapide(probabilites_transition_1vers2,0,groupe1->nb_individus-1);
      moyenne_changement_classe1=0;
      for(i=0;i<groupe1->nb_individus;i++)
	{
	  printf("%.3f\n",probabilites_transition_1vers2[i]);
	  moyenne_changement_classe1+=probabilites_transition_1vers2[i];
	}
      moyenne_changement_classe1/=(double)groupe1->nb_individus;
      
      printf("\n**********************************\n");
      
      tri_rapide(probabilites_transition_2vers1,0,groupe2->nb_individus-1);
      moyenne_changement_classe2=0;
      for(i=0;i<groupe2->nb_individus;i++)
	{
	  printf("%.3f\n",probabilites_transition_2vers1[i]);
	  moyenne_changement_classe2+=probabilites_transition_2vers1[i];
	}
      moyenne_changement_classe2/=(double)groupe2->nb_individus;
      
      probabilite_attendue_changement_classe1=
	(double)groupe2->nb_individus/(double)(groupe_pere->nb_individus-1);
      probabilite_attendue_changement_classe2=
	(double)groupe1->nb_individus/(double)(groupe_pere->nb_individus-1);
      
      printf("moyenne1 : %f ; attendue : %f\n",moyenne_changement_classe1,
	     probabilite_attendue_changement_classe1);
      printf("moyenne2 : %f ; attendue : %f\n",moyenne_changement_classe2,
	     probabilite_attendue_changement_classe2);
      
      if((probabilite_attendue_changement_classe1-moyenne_changement_classe1>0.1)||
	 (probabilite_attendue_changement_classe2-moyenne_changement_classe2>0.1))
	{
	  on_divise=OUI;
	  printf("ON DIVISE\n");
	}
      else
	{
	  on_divise=NON;  
	  printf("ON NE DIVISE PAS\n");
	}
    }

  /*desallocation memoire*/
  free(probabilites_transition_1vers2);
  free(probabilites_transition_2vers1);
  /*fin desallocation memoire*/

  return on_divise;
}



