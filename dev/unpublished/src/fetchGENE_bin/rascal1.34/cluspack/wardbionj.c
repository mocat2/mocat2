#include "principal.h"
#include "wardbionj.h"

/******************************************/
/*                                        */
/*Procedure de classification hierarchique*/
/*                                        */
/******************************************/

void ward_sur_un_arbre_bionj(int nb_sequences,noeud_t **adresses_noeuds,
			     noeud_t *noeud_fictif,int *compteur_indices_dissimilarite,
			     double *historique_indices_dissimilarite,
			     noeud_t **historique_noeuds)
{
  /*declaration des variables*/
  int i,j,k,l,groupe1,groupe2,temp_groupe,copain1,copain2,temp_copain,un_noeud_en_commun;
  int *numeros_groupes_courants,nb_groupes_courants,kcourant,lcourant;
  double **indices_dissimilarite_anciens,**indices_dissimilarite_nouveaux;
  double plus_petit_indice_dissimilarite;
  /*fin declaration des variables*/

  /*allocation memoire*/
  numeros_groupes_courants=(int *)malloc(sizeof(int)*nb_sequences);
  
  indices_dissimilarite_anciens=(double **)malloc(sizeof(double *)*nb_sequences);
  for(i=0;i<nb_sequences;i++)
    {
      indices_dissimilarite_anciens[i]=(double *)malloc(sizeof(double)*nb_sequences);
    }
  
  indices_dissimilarite_nouveaux=(double **)malloc(sizeof(double *)*nb_sequences);
  for(i=0;i<nb_sequences;i++)
    {
      indices_dissimilarite_nouveaux[i]=(double *)malloc(sizeof(double)*nb_sequences);
    }
  /*fin allocation memoire*/

  /*mise a zero des pertes d'inertie inter-classe*/
  for(i=0;i<nb_sequences;i++)
    { 
      adresses_noeuds[i]->dissimilarite=0;
    } 

  nb_groupes_courants=nb_sequences;
  
  /*initialisation des indices de dissimilarite*/
  initialisation_indices_dissimilarite(nb_sequences,indices_dissimilarite_anciens,adresses_noeuds);

  /*initialisation des numeros des groupes ;*/
  /*les numeros de groupe etant les numeros des derniers noeuds qui ont rejoint leur groupe*/
  for(i=0;i<nb_sequences;i++)
    {
      numeros_groupes_courants[i]=i;
    }

  *compteur_indices_dissimilarite=0;
  while(nb_groupes_courants>2)
    {
      /*on recherche le couple de groupes (groupe1,groupe2)*/
      /*ayant le plus petit indice de dissimilarite*/
      plus_petit_indice_dissimilarite=-1.0;
      for(i=0;i<nb_groupes_courants;i++)
	{
	  for(j=i+1;j<nb_groupes_courants;j++)
	    {
	      un_noeud_en_commun=NON;
	      for(k=0;k<3;k++)
		{
		  for(l=0;l<3;l++)
		    {
		      if((adresses_noeuds[numeros_groupes_courants[i]]->copains[k]!=NULL)&&
			 (adresses_noeuds[numeros_groupes_courants[j]]->copains[l]!=NULL)&&
			 (adresses_noeuds[numeros_groupes_courants[i]]->copains[k]==
			  adresses_noeuds[numeros_groupes_courants[j]]->copains[l]))
			{
			  un_noeud_en_commun=OUI;
			  kcourant=k;
			  lcourant=l;
			  break;
			}
		    }
		  if(un_noeud_en_commun==OUI)
		    {
		      break;
		    }
		}

	      if(((indices_dissimilarite_anciens[i][j]<=plus_petit_indice_dissimilarite)||
		  (plus_petit_indice_dissimilarite<0))&&(un_noeud_en_commun==OUI))
	    {
		  plus_petit_indice_dissimilarite=indices_dissimilarite_anciens[i][j];
		  groupe1=i;
		  groupe2=j;
		  copain1=kcourant;
		  copain2=lcourant;
		}
	    }
	}
      /*fin recherche des groupes a reunir*/

      nb_groupes_courants--;
      
      adresses_noeuds[numeros_groupes_courants[groupe1]]->copains[copain1]->dissimilarite=
	plus_petit_indice_dissimilarite;

      /*      adresses_noeuds[numeros_groupes_courants[groupe1]]->distance=
	      adresses_noeuds[numeros_groupes_courants[groupe1]]->distances[copain1];
	      adresses_noeuds[numeros_groupes_courants[groupe2]]->distance=
	      adresses_noeuds[numeros_groupes_courants[groupe2]]->distances[copain2];*/

      adresses_noeuds[numeros_groupes_courants[groupe1]]->copains[copain1]->poids=
	adresses_noeuds[numeros_groupes_courants[groupe1]]->poids+
	adresses_noeuds[numeros_groupes_courants[groupe2]]->poids;

      adresses_noeuds[numeros_groupes_courants[groupe1]]->copains[copain1]->copain1=
	adresses_noeuds[numeros_groupes_courants[groupe1]];
      adresses_noeuds[numeros_groupes_courants[groupe1]]->copains[copain1]->copain2=
	adresses_noeuds[numeros_groupes_courants[groupe2]];
       
      historique_indices_dissimilarite[*compteur_indices_dissimilarite]=
	plus_petit_indice_dissimilarite;

      historique_noeuds[*compteur_indices_dissimilarite]=
	adresses_noeuds[numeros_groupes_courants[groupe1]]->copains[copain1];

      if(groupe2<groupe1)
	{
	  temp_groupe=groupe1;
	  groupe1=groupe2;
	  groupe2=temp_groupe;

	  temp_copain=copain1;
	  copain1=copain2;
	  copain2=temp_copain;
	}

      /*mise a jour des indices de dissimilarite*/
      mise_a_jour_indices_Ward(groupe1,groupe2,adresses_noeuds,nb_groupes_courants,
			       numeros_groupes_courants,indices_dissimilarite_anciens,
			       indices_dissimilarite_nouveaux);

      /*le groupe1 happe le groupe2 qui disparait*/
      for(i=groupe2;i<nb_groupes_courants;i++)
	{
	  numeros_groupes_courants[i]=numeros_groupes_courants[i+1];
	} 
      numeros_groupes_courants[groupe1]=
	adresses_noeuds[numeros_groupes_courants[groupe1]]->copains[copain1]->numero; 

      (*compteur_indices_dissimilarite)++;
    }
  
  for(i=0;i<3;i++)
    {
      if(adresses_noeuds[numeros_groupes_courants[0]]->copains[i]==
	 adresses_noeuds[numeros_groupes_courants[1]])
	{
	  copain1=i;
	  break;
	}
    }
   for(i=0;i<3;i++)
    {
      if(adresses_noeuds[numeros_groupes_courants[1]]->copains[i]==
	 adresses_noeuds[numeros_groupes_courants[0]])
	{
	  copain2=i;
	  break;
	}
    }

   strcpy(noeud_fictif->etiquette,"noeud_fictif");
   noeud_fictif->copains[0]=NULL;
   noeud_fictif->copains[1]=adresses_noeuds[numeros_groupes_courants[0]];
   noeud_fictif->copains[2]=adresses_noeuds[numeros_groupes_courants[1]];
   noeud_fictif->copain1=adresses_noeuds[numeros_groupes_courants[0]];
   noeud_fictif->copain2=adresses_noeuds[numeros_groupes_courants[1]];
   noeud_fictif->dissimilarite=indices_dissimilarite_anciens[0][1];
   /*   noeud_fictif->distances[1]=adresses_noeuds[numeros_groupes_courants[0]]->distance;
	noeud_fictif->distances[2]=adresses_noeuds[numeros_groupes_courants[1]]->distance;*/

   adresses_noeuds[numeros_groupes_courants[0]]->copains[copain1]=noeud_fictif;
   adresses_noeuds[numeros_groupes_courants[1]]->copains[copain2]=noeud_fictif;
   
   /*   adresses_noeuds[numeros_groupes_courants[0]]->distance=
	adresses_noeuds[numeros_groupes_courants[0]]->distances[copain1]/2.0;
	adresses_noeuds[numeros_groupes_courants[1]]->distance=
	adresses_noeuds[numeros_groupes_courants[1]]->distances[copain2]/2.0;*/
   
   historique_noeuds[*compteur_indices_dissimilarite]=noeud_fictif;
   historique_indices_dissimilarite[*compteur_indices_dissimilarite]=noeud_fictif->dissimilarite;
   (*compteur_indices_dissimilarite)++;

   /*desallocation memoire*/
   free(numeros_groupes_courants);
   
   for(i=0;i<nb_sequences;i++)
     {
       free(indices_dissimilarite_anciens[i]);
     }
   free(indices_dissimilarite_anciens);
   
   for(i=0;i<nb_sequences;i++)
     {
       free(indices_dissimilarite_nouveaux[i]);
     }
   free(indices_dissimilarite_nouveaux);
   /*fin desallocation memoire*/
}

/*********************************************************/
/*                                                       */
/*Procedure d'initialisation des indices de dissimilarite*/
/*                                                       */
/*********************************************************/

void initialisation_indices_dissimilarite(int nb_individus,double **indices_dissimilarite_anciens,
					  noeud_t **adresses_noeuds)
{
  /*declaration des variables*/
  int i,j;
  /*fin declaration des variables*/
 
  /*indices de Ward*/
  for(i=0;i<nb_individus-1;i++)
    {
      for(j=i+1;j<nb_individus;j++)
	{
	  indices_dissimilarite_anciens[i][j]=
	    (adresses_noeuds[i]->poids*adresses_noeuds[j]->poids)*
	    pow(adresses_noeuds[i]->individu->valeurs_traitees[adresses_noeuds[j]->individu->id],
		2.0)/(adresses_noeuds[i]->poids+adresses_noeuds[j]->poids);
	  indices_dissimilarite_anciens[j][i]=indices_dissimilarite_anciens[i][j];
	}
    } 
}


/***************************************************/
/*                                                 */
/*Procedure de mise a jour des indices de Ward     */
/*lors de la classification hierarchique ascendante*/
/*                                                 */
/***************************************************/

void mise_a_jour_indices_Ward(int groupe1,int groupe2,noeud_t **adresses_noeuds,
			      int nb_groupes_courants,int *numeros_groupes_courants,
			      double **indices_Ward_anciens,double **indices_Ward_nouveaux)
{
  /*declaration des variables*/
  int i,j;
  /*fin declaration des variables*/
  
  /******************************************/
  /*mise a jour des indices de Ward nouveaux*/
  /******************************************/

  /*mise a jour des indices de groupe1*/
  for(i=0;i<groupe2;i++)
    {
      if(i!=groupe1)
	{
	  indices_Ward_nouveaux[groupe1][i]=
	    (indices_Ward_anciens[groupe1][i]*
	     (adresses_noeuds[numeros_groupes_courants[groupe1]]->poids+
	      adresses_noeuds[numeros_groupes_courants[i]]->poids)+
	     indices_Ward_anciens[groupe2][i]*
	     (adresses_noeuds[numeros_groupes_courants[groupe2]]->poids+
	      adresses_noeuds[numeros_groupes_courants[i]]->poids)-
	     indices_Ward_anciens[groupe1][groupe2]*
	     adresses_noeuds[numeros_groupes_courants[i]]->poids)/
	    (adresses_noeuds[numeros_groupes_courants[groupe1]]->poids+
	     adresses_noeuds[numeros_groupes_courants[groupe2]]->poids+
	     adresses_noeuds[numeros_groupes_courants[i]]->poids);
      
	  indices_Ward_nouveaux[i][groupe1]=indices_Ward_nouveaux[groupe1][i];
	}
    }
  indices_Ward_nouveaux[groupe1][groupe1]=0;

  for(i=groupe2;i<nb_groupes_courants;i++)
    {
      indices_Ward_nouveaux[groupe1][i]=
	(indices_Ward_anciens[groupe1][i+1]*
	 (adresses_noeuds[numeros_groupes_courants[groupe1]]->poids+
	  adresses_noeuds[numeros_groupes_courants[i+1]]->poids)+
	 indices_Ward_anciens[groupe2][i+1]*
	 (adresses_noeuds[numeros_groupes_courants[groupe2]]->poids+
	  adresses_noeuds[numeros_groupes_courants[i+1]]->poids)-
	 indices_Ward_anciens[groupe1][groupe2]*
	 adresses_noeuds[numeros_groupes_courants[i+1]]->poids)/
	(adresses_noeuds[numeros_groupes_courants[groupe1]]->poids+
	 adresses_noeuds[numeros_groupes_courants[groupe2]]->poids+
	 adresses_noeuds[numeros_groupes_courants[i+1]]->poids);
      
      indices_Ward_nouveaux[i][groupe1]=indices_Ward_nouveaux[groupe1][i];
    }
  /*fin mise a jour des indices de groupe1*/

  for(i=0;i<groupe2;i++)
    {
      if(i!=groupe1)
	{
	  for(j=i+1;j<groupe2;j++)
	    {
	      if(j!=groupe1)
		{
		  indices_Ward_nouveaux[i][j]=indices_Ward_anciens[i][j];
		  indices_Ward_nouveaux[j][i]=indices_Ward_nouveaux[i][j];
		}
	    }
	}
    }
  
  for(i=0;i<groupe2;i++)
    {
      if(i!=groupe1)
	{
	  for(j=groupe2;j<nb_groupes_courants;j++)
	    {
	      indices_Ward_nouveaux[i][j]=indices_Ward_anciens[i][j+1];
	      indices_Ward_nouveaux[j][i]=indices_Ward_nouveaux[i][j];
	    }
	}
    }
  
  for(i=groupe2;i<nb_groupes_courants;i++)
    {
      for(j=i+1;j<nb_groupes_courants;j++)
	{
	  indices_Ward_nouveaux[i][j]=indices_Ward_anciens[i+1][j+1];
	  indices_Ward_nouveaux[j][i]=indices_Ward_nouveaux[i][j];
	}
    }

  /**********************************************/
  /*fin mise a jour des indices de Ward nouveaux*/
  /**********************************************/

  /*****************************************/
  /*mise a jour des indices de Ward anciens*/
  /*****************************************/

  for(i=0;i<nb_groupes_courants;i++)
    {
      for(j=0;j<nb_groupes_courants;j++)
	{
	  if(i!=j)
	    {
	      indices_Ward_anciens[i][j]=indices_Ward_nouveaux[i][j];
	    }
	}
    }

  /*****************************************/
  /*mise a jour des indices de Ward anciens*/
  /*****************************************/
}

