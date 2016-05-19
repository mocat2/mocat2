#include "principal.h"
#include "divide.h"
#include "space.h"
#include "tools.h"

/**************************************************************/
/*                                                            */
/*Procedure de test pour savoir si l'on doit scinder un groupe*/
/*par la methode DPC (Density of Points Clustering)           */ 
/*                                                            */
/**************************************************************/

int test_DPC(int nb_individus,int nb_dimensions,groupe_t *groupe_pere,
	     groupe_t *groupe1,groupe_t *groupe2,double **coordonnees_individus,
	     double *cotes_hyperpave,double *valeurs_minimum,double risque,
	     int type_donnees,int type_densite,int donnees_normalisees,int nb_simulations)
{
  /*declaration de variables*/
  FILE *fichier;
  char poubelle[200];
  int i,j,k,individu1,individu2,nombre_de_contrevenants_chez_le_pere,nb_iterations,nb_iterations2;
  int nombre_de_contrevenants_aleatoire,nombre,nb_fisher_snedecor;
  int nb_fois_hypothese_rejetee,compteur_simulations,compteur_population_reference,compteur;
  double *individu_fictif,exposant_Holder,exposant_Holder_seuil,*exposants_Holder_pere;
  double moyenne_pere,moyenne_aleatoire,moyenne;
  double *exposants_Holder_aleatoire,variance_residuelle,variance_facteur,distance_courante;
  double distance_minimum,f,P,*direction_aleatoire,signe,m,s,variance_residuelle_pere;
  double risque_encouru;
  /*fin declaration de variables*/

  nb_fois_hypothese_rejetee=0;
  risque_encouru=0;

  nb_iterations2=groupe1->nb_individus*groupe2->nb_individus/2;

  if(nb_iterations2>500)
    {
      nb_iterations2=500;
    }
  else if(nb_iterations2<2)
    {
      nb_iterations2=2;
    }

  /*allocation memoire*/
  individu_fictif=(double *)malloc(sizeof(double)*nb_dimensions);
  direction_aleatoire=(double *)malloc(sizeof(double)*nb_dimensions);
  exposants_Holder_pere=(double *)malloc(sizeof(double)*(2*groupe_pere->nb_individus+2));
  exposants_Holder_aleatoire=(double *)malloc(sizeof(double)*nb_iterations2);
  /*fin allocation memoire*/

  /*printf("taille1 : %d ; taille2 : %d\n",groupe1->nb_individus,groupe2->nb_individus);*/
  /*for(i=0;i<nb_dimensions;i++)
    {
    individu_fictif[i]=0.45;
    }
    
    calcul_exposants_Holder_points(nb_dimensions,nb_individus,coordonnees_individus,
    1,&individu_fictif,&exposant_Holder,
    cotes_hyperpave,valeurs_minimum,NON,OUI);
    for(i=0;i<nb_individus;i++)
    {
    for(j=0;j<nb_dimensions;j++)
    {
    printf("%.3f ",coordonnees_individus[i][j]);
    }
    printf("\n");
    }
    printf("exposant holder : %f\n",exposant_Holder);
    exit(0);*/

  moyenne_pere=0;
  compteur_population_reference=0;
  
  if(groupe1->nb_individus>1)
    {
      if(groupe1->nb_individus<50)
	{
	  nb_iterations=groupe1->nb_individus;
	}
      else
	{
	  nb_iterations=50;
	}
      for(i=0;i<nb_iterations;i++)
	{
	  individu1=(int)floor((double)groupe1->nb_individus*(double)rand()/(double)RAND_MAX);
	  if(individu1==groupe1->nb_individus)
	    {
	      individu1=groupe1->nb_individus-1;
	    }

	  distance_minimum=-1.0;
	  for(j=0;j<groupe1->nb_individus;j++)
	    {
	      if(j!=individu1)
		{
		  distance_courante=0;
		  for(k=0;k<nb_dimensions;k++)
		    {
		      distance_courante+=
			(groupe1->individus[individu1]->valeurs_traitees[k]-
			 groupe1->individus[j]->valeurs_traitees[k])*
			(groupe1->individus[individu1]->valeurs_traitees[k]-
			 groupe1->individus[j]->valeurs_traitees[k]);
		    }
		  if((distance_minimum<0)||(distance_courante<distance_minimum))
		    {
		      distance_minimum=distance_courante;
		      individu2=j;
		    }
		}
	    }

	  for(j=0;j<nb_dimensions;j++)
	    {
	      direction_aleatoire[j]=groupe1->individus[individu2]->valeurs_traitees[j]-
		groupe1->individus[individu1]->valeurs_traitees[j];
	    }
	 
	  compteur=0;
	  for(j=0;j<nb_dimensions;j++) 
	    { 
	      if((double)rand()/(double)RAND_MAX<0.5)
		{
		  signe=-1.0;
		  compteur++;
		}
	      else
		{
		  signe=1.0;
		}
	      if((compteur>1)&&(compteur>nb_dimensions/5)&&
		 (type_densite==DENSITE2))
		{
		  signe=1.0;
		}
	      individu_fictif[j]=groupe1->individus[individu1]->valeurs_traitees[j]+
		signe*direction_aleatoire[j];
	    }

	  if(donnees_normalisees==OUI)
	    {
	      /*on normalise l'individu fictif*/
	      m=0;
	      for(j=0;j<nb_dimensions;j++)
		{
		  m+=individu_fictif[j];
		}
	      m/=(double)nb_dimensions;

	      s=0;
	      for(j=0;j<nb_dimensions;j++)
		{
		  s+=((double)(m-individu_fictif[j]))*((double)(m-individu_fictif[j]));
		}
	      s/=(double)nb_dimensions;
	      s=sqrt((double)s);

	      for(j=0;j<nb_dimensions;j++)
		{
		  individu_fictif[j]=(individu_fictif[j]-m)/s;
		}
	    }

	  /*on verifie que l'individu fictif n'est pas hors-norme*/
	  for(j=0;j<nb_dimensions;j++)
	    {
	      if(individu_fictif[j]<valeurs_minimum[j])
		{
		  individu_fictif[j]=valeurs_minimum[j];
		}
	      else if(individu_fictif[j]>valeurs_minimum[j]+cotes_hyperpave[j])
		{
		  individu_fictif[j]=valeurs_minimum[j]+cotes_hyperpave[j];
		}
	    }
	  
	  calcul_exposants_Holder_points(nb_dimensions,nb_individus,coordonnees_individus,
					 1,&individu_fictif,&exposant_Holder,
					 cotes_hyperpave,valeurs_minimum,NON,OUI);

	  exposants_Holder_pere[compteur_population_reference]=exposant_Holder;
	  moyenne_pere+=exposant_Holder;
	  compteur_population_reference++;
	}
    }

  if(groupe2->nb_individus>1)
    {
      if(groupe2->nb_individus<50)
	{
	  nb_iterations=groupe2->nb_individus;
	}
      else
	{
	  nb_iterations=50;
	}
      for(i=0;i<nb_iterations;i++)
	{
	  individu1=(int)floor((double)groupe2->nb_individus*(double)rand()/(double)RAND_MAX);
	  if(individu1==groupe2->nb_individus)
	    {
	      individu1=groupe2->nb_individus-1;
	    }  
	  distance_minimum=-1.0;
	  for(j=0;j<groupe2->nb_individus;j++)
	    {
	      if(j!=individu1)
		{
		  distance_courante=0;
		  for(k=0;k<nb_dimensions;k++)
		    {
		      distance_courante+=
			(groupe2->individus[individu1]->valeurs_traitees[k]-
			 groupe2->individus[j]->valeurs_traitees[k])*
			(groupe2->individus[individu1]->valeurs_traitees[k]-
			 groupe2->individus[j]->valeurs_traitees[k]);
		    }
		  if((distance_minimum<0)||(distance_courante<distance_minimum))
		    {
		      distance_minimum=distance_courante;
		      individu2=j;
		    }
		}
	    }

	  for(j=0;j<nb_dimensions;j++)
	    {
	      direction_aleatoire[j]=groupe2->individus[individu2]->valeurs_traitees[j]-
		groupe2->individus[individu1]->valeurs_traitees[j];
	    }

	  compteur=0;
	  for(j=0;j<nb_dimensions;j++) 
	    { 
	      if((double)rand()/(double)RAND_MAX<0.5)
		{
		  signe=-1.0;
		  compteur++;
		}
	      else
		{
		  signe=1.0;
		}
	      if((compteur>1)&&(compteur>nb_dimensions/5)&&
		 ((type_densite==DENSITE2)))
		{
		  signe=1.0;
		}
	      individu_fictif[j]=groupe2->individus[individu1]->valeurs_traitees[j]+
		signe*direction_aleatoire[j];
	    }

	  if(donnees_normalisees==OUI)
	    {
	      /*on normalise l'individu fictif*/
	      m=0;
	      for(j=0;j<nb_dimensions;j++)
		{
		  m+=individu_fictif[j];
		}
	      m/=(double)nb_dimensions;

	      s=0;
	      for(j=0;j<nb_dimensions;j++)
		{
		  s+=((double)(m-individu_fictif[j]))*((double)(m-individu_fictif[j]));
		}
	      s/=(double)nb_dimensions;
	      s=sqrt((double)s);

	      for(j=0;j<nb_dimensions;j++)
		{
		  individu_fictif[j]=(individu_fictif[j]-m)/s;
		}
	    }

	  /*on verifie que l'individu fictif n'est pas hors-norme*/
	  for(j=0;j<nb_dimensions;j++)
	    {
	      if(individu_fictif[j]<valeurs_minimum[j])
		{
		  individu_fictif[j]=valeurs_minimum[j];
		}
	      else if(individu_fictif[j]>valeurs_minimum[j]+cotes_hyperpave[j])
		{
		  individu_fictif[j]=valeurs_minimum[j]+cotes_hyperpave[j];
		}
	    }
	  	  
	  calcul_exposants_Holder_points(nb_dimensions,nb_individus,coordonnees_individus,
					 1,&individu_fictif,&exposant_Holder,
					 cotes_hyperpave,valeurs_minimum,NON,OUI);

	  exposants_Holder_pere[compteur_population_reference]=exposant_Holder;
	  moyenne_pere+=exposant_Holder;
	  compteur_population_reference++;
	}
    }

  moyenne_pere/=(double)(compteur_population_reference);
  variance_residuelle_pere=0;
  for(i=0;i<compteur_population_reference;i++)
    { 
      variance_residuelle_pere+=((double)(exposants_Holder_pere[i]-moyenne_pere))*
	((double)(exposants_Holder_pere[i]-moyenne_pere));
    } 

  for(compteur_simulations=0;compteur_simulations<nb_simulations;compteur_simulations++)
    {
      nombre_de_contrevenants_aleatoire=0;
      moyenne_aleatoire=0;
      for(i=0;i<nb_iterations2;i++)
	{
	  /*selection aleatoire de deux individus du groupe pere*/
	  individu1=(int)floor((double)(groupe1->nb_individus-1)*(double)rand()/(double)RAND_MAX);
	  individu2=(int)floor((double)(groupe2->nb_individus-1)*(double)rand()/(double)RAND_MAX);
	
	  for(j=0;j<nb_dimensions;j++) 
	    { 
	      if((type_donnees==ALIGNEMENT)||(type_densite==DENSITE1)||(type_densite==DENSITE2))
		{
		  individu_fictif[j]=(groupe1->individus[individu1]->valeurs_traitees[j]-
				      groupe2->individus[individu2]->valeurs_traitees[j])*0.5
		    +groupe2->individus[individu2]->valeurs_traitees[j];
		}
	      else
		{
		  individu_fictif[j]=(groupe1->individus[individu1]->valeurs_traitees[j]-
				      groupe2->individus[individu2]->valeurs_traitees[j])*
		    ((double)rand()/(double)RAND_MAX)
		    +groupe2->individus[individu2]->valeurs_traitees[j];
		}
	    }
	  
	  if(donnees_normalisees==OUI)
	    {
	      /*on normalise l'individu fictif*/
	      m=0;
	      for(j=0;j<nb_dimensions;j++)
		{
		  m+=individu_fictif[j];
		}
	      m/=(double)nb_dimensions;
	      
	      s=0;
	      for(j=0;j<nb_dimensions;j++)
		{
		  s+=((double)(m-individu_fictif[j]))*((double)(m-individu_fictif[j]));
		}
	      s/=(double)nb_dimensions;
	      s=sqrt((double)s);
	      
	      for(j=0;j<nb_dimensions;j++)
		{
		  individu_fictif[j]=(individu_fictif[j]-m)/s;
		}
	    }

	  calcul_exposants_Holder_points(nb_dimensions,nb_individus,coordonnees_individus,
					 1,&individu_fictif,&exposant_Holder,
					 cotes_hyperpave,valeurs_minimum,NON,NON);
	  
	  exposants_Holder_aleatoire[i]=exposant_Holder;
	  moyenne_aleatoire+=exposant_Holder;
	}
     
      moyenne_aleatoire=moyenne_aleatoire/(double)nb_iterations2;

      /*realisation du test*/
      moyenne=(moyenne_aleatoire*(double)nb_iterations2+moyenne_pere*
	       (double)compteur_population_reference)
	/(double)(nb_iterations2+compteur_population_reference);
    
      variance_facteur=(((double)(moyenne_pere-moyenne))*((double)(moyenne_pere-moyenne))*
			(double)compteur_population_reference+((double)(moyenne_aleatoire-moyenne))*
			((double)(moyenne_aleatoire-moyenne))*(double)nb_iterations2)/
	(double)(nb_iterations2+compteur_population_reference);
     
      variance_residuelle=variance_residuelle_pere;
      for(i=0;i<nb_iterations2;i++)
	{
	  variance_residuelle+=((double)(exposants_Holder_aleatoire[i]-moyenne_aleatoire))*
	    ((double)(exposants_Holder_aleatoire[i]-moyenne_aleatoire));
	}

      variance_residuelle/=(double)(nb_iterations2+compteur_population_reference);
      

      if(variance_residuelle>ZERO_LIMIT)
	{
	  f=variance_facteur*(double)(nb_iterations2+compteur_population_reference-2)/
	    variance_residuelle;
	  
	  P=F_fisher_snedecor(1.0,(double)(nb_iterations2+compteur_population_reference-2),
			      (double)f);

	  risque_encouru+=(1.0-P); 

	  /*printf("1-P : %f ; %f ; %f\n",1.0-P,moyenne_pere,moyenne_aleatoire);*/
	    
	  if((1.0-P<risque)&&(moyenne_pere-ZERO_LIMIT>moyenne_aleatoire))
	    {  
	      nb_fois_hypothese_rejetee++;
	    } 
	}
    }

  /*desallocation memoire*/
  free(individu_fictif);
  free(direction_aleatoire);
  free(exposants_Holder_pere);
  free(exposants_Holder_aleatoire);
  /*fin desallocation memoire*/

  if(nb_fois_hypothese_rejetee>nb_simulations/2)
    {
      return OUI;
    }
  else
    {
      return NON;
    }
}

