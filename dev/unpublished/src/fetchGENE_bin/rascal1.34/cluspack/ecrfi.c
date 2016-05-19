#include "principal.h"

/**********************************************************************************/
/*                                                                                */
/*Procedure d'ecriture du fichier.dst des distances en vue du calcul de phylogenie*/
/*                                                                                */
/**********************************************************************************/

void ecriture_fichier_distances_sequences(char *nom_fichier,int nb_individus,
					  individu_t **individus)
{
  /*declaration des variables*/
  FILE *fichier;
  int i,j;
  /*fin declaration des variables*/

  fichier=fopen(nom_fichier,"w+");

  fprintf(fichier,"%d\n",nb_individus+1);
  for(i=0;i<nb_individus;i++)
    {
      fprintf(fichier,"%d ",individus[i]->id);
      for(j=0;j<nb_individus;j++)
	{
	  fprintf(fichier,"%f ",(float)individus[i]->valeurs_traitees[individus[j]->id]);
	  if((j+1)%8==0)
	    {
	      fprintf(fichier,"\n");
	    }
	}
      fprintf(fichier,"1.0 \n");
    }
  fprintf(fichier,"fictive_sequence ");
  for(j=0;j<nb_individus;j++)
    {
      fprintf(fichier,"1.0 ");
      if((j+1)%8==0)
	{
	  fprintf(fichier,"\n");
	}
    }
  fprintf(fichier,"0.0 \n");
  fclose(fichier);
}

/******************************************************************/
/*                                                                */
/*Procedure d'ecriture du fichier resultat presentant les clusters*/
/*                                                                */
/******************************************************************/

void ecriture_fichier_clusters(int nb_individus,int nb_dimensions,individu_t **individus,
			       int nb_clusters,int nb_clusters_selection,int individus_orphelins,
			       int write_coordinates,char *nom_fichier)
{
  /*declaration des variables*/
  int i,j,k,*tailles_clusters,nb_individus_total;
  individu_t ***clusters;
  FILE *file;
  /*fin declaration des variables*/

  /*allocation memoire*/
  tailles_clusters=(int *)malloc(sizeof(int)*(nb_clusters+1));
  clusters=(individu_t ***)malloc(sizeof(individu_t *)*(nb_clusters+1));
  /*fin allocation memoire*/

  nb_individus_total=nb_individus; 
  for(i=0;i<nb_individus;i++)
    {
      nb_individus_total+=individus[i]->nb_individus_similaires;
    }

  for(i=0;i<nb_clusters+1;i++)
    {
      clusters[i]=(individu_t **)malloc(sizeof(individu_t *)*nb_individus_total);
    }
  for(i=0;i<nb_clusters+1;i++)
    {
      tailles_clusters[i]=0;
    }

  for(i=0;i<nb_individus;i++)
    {
      if(individus[i]->cluster>-1)
	{
	  clusters[individus[i]->cluster][tailles_clusters[individus[i]->cluster]]=individus[i];
	  (tailles_clusters[individus[i]->cluster])++;
	  if(individus_orphelins==NON)
	    {
	      for(j=0;j<individus[i]->nb_individus_similaires;j++)
		{
		  clusters[individus[i]->cluster][tailles_clusters[individus[i]->cluster]]=
		    individus[i]->individus_similaires[j];
		  (tailles_clusters[individus[i]->cluster])++;
		}
	    }
	}
      else
	{
	  clusters[nb_clusters][tailles_clusters[nb_clusters]]=individus[i];
	  (tailles_clusters[nb_clusters])++;
	  if(individus_orphelins==NON)
	    {
	      for(j=0;j<individus[i]->nb_individus_similaires;j++)
		{
		  clusters[nb_clusters][tailles_clusters[nb_clusters]]=
		    individus[i]->individus_similaires[j];
		  (tailles_clusters[nb_clusters])++;
		}
	    }
	}
    }

  file=fopen(nom_fichier,"w");
  fprintf(file,"Number of clusters : %d\n",nb_clusters);
  
  for(i=0;i<nb_clusters;i++)
    {
      fprintf(file,"\nCluster %d ; size=%d\n",i,tailles_clusters[i]);
      for(j=0;j<tailles_clusters[i];j++)
	{
	  fprintf(file,"%s",clusters[i][j]->nom);
	  if(write_coordinates==OUI)
	    {
	      for(k=0;k<nb_dimensions;k++)
		{
		  fprintf(file,"\t%f",clusters[i][j]->valeurs_traitees[k]);
		}
	    }	  
	  if(clusters[i][j]->description!=NULL)
	    {
	      fprintf(file,"%s",clusters[i][j]->description);
	    }
	  else
	    {
	      fprintf(file,"\n");
	    }
	}
    }
  if(tailles_clusters[nb_clusters]>0)
    {
      fprintf(file,"\nunclustered points %d ; size=%d\n",i,tailles_clusters[nb_clusters]);
      for(j=0;j<tailles_clusters[nb_clusters];j++)
	{
	  fprintf(file,"%s",clusters[nb_clusters][j]->nom);
	  if(write_coordinates==OUI)
	    {
	      for(k=0;k<nb_dimensions;k++)
		{
		  fprintf(file,"\t%f",clusters[nb_clusters][j]->valeurs_traitees[k]);
		}
	    }
	  if(clusters[i][j]->description!=NULL)
	    {
	      fprintf(file,"%s",clusters[i][j]->description);
	    }
	  else
	    {
	      fprintf(file,"\n");
	    }
	}
    }
  fclose(file);
 
  /*desallocation memoire*/
  free(tailles_clusters);
  for(i=0;i<nb_clusters;i++)
    {
      free(clusters[i]);
    }
  free(clusters);
  /*fin desallocation memoire*/
}

/**********************************************************/
/*                                                        */
/*Procedure d'ecriture du fichier de l'alignement multiple*/
/*avec les sequences reordonnees suivant leurs groupes    */
/*                                                        */
/**********************************************************/

void ecriture_fichier_clusters_alignement(char *nom_fichier,int longueur_alignement,
					  int nb_individus,individu_t **individus,
					  char **sequences,int individus_orphelins,
					  int nb_clusters)
{
  /*declaration des variables*/
  int i,j,k,*tailles_clusters,nb_individus_total;
  individu_t ***clusters;
  FILE *file;
  /*fin declaration des variables*/

  /*allocation memoire*/
  tailles_clusters=(int *)malloc(sizeof(int)*(nb_clusters+1));
  clusters=(individu_t ***)malloc(sizeof(individu_t *)*(nb_clusters+1));
  /*fin allocation memoire*/

  nb_individus_total=nb_individus;
  for(i=0;i<nb_individus;i++)
    {
      nb_individus_total+=individus[i]->nb_individus_similaires;
    }
  for(i=0;i<nb_clusters+1;i++)
    {
      clusters[i]=(individu_t **)malloc(sizeof(individu_t *)*nb_individus_total);
    }
  for(i=0;i<nb_clusters+1;i++)
    {
      tailles_clusters[i]=0;
    }

  for(i=0;i<nb_individus;i++)
    {
      if(individus[i]->cluster>-1)
	{
	  clusters[individus[i]->cluster][tailles_clusters[individus[i]->cluster]]=individus[i];
	  (tailles_clusters[individus[i]->cluster])++;
	  if(individus_orphelins==NON)
	    {
	      for(j=0;j<individus[i]->nb_individus_similaires;j++)
		{
		  clusters[individus[i]->cluster][tailles_clusters[individus[i]->cluster]]=
		    individus[i]->individus_similaires[j];
		  (tailles_clusters[individus[i]->cluster])++;
		}
	    }
	}
      else
	{
	  clusters[nb_clusters][tailles_clusters[nb_clusters]]=individus[i];
	  (tailles_clusters[nb_clusters])++;
	  if(individus_orphelins==NON)
	    {
	      for(j=0;j<individus[i]->nb_individus_similaires;j++)
		{
		  clusters[nb_clusters][tailles_clusters[nb_clusters]]=
		    individus[i]->individus_similaires[j];
		  (tailles_clusters[nb_clusters])++;
		}
	    }
	}
    }

  file=fopen(nom_fichier,"w");
  for(i=0;i<nb_clusters;i++)
    {
      fprintf(file,">GROUP_%d\n",i+1);
      for(j=0;j<longueur_alignement;j++)
	{
	  fprintf(file,"-");
	}
      fprintf(file,"\n");
      for(j=0;j<tailles_clusters[i];j++)
	{
	  fprintf(file,">%s\n",clusters[i][j]->nom);
	  for(k=0;k<longueur_alignement;k++)
	    {
	      if(sequences[clusters[i][j]->id][k]=='O')
		{
		  fprintf(file,"-");
		}
	      else
		{
		  fprintf(file,"%c",sequences[clusters[i][j]->id][k]);
		}
	    }
	  fprintf(file,"\n");
	}
    }
  if(tailles_clusters[nb_clusters]>0)
    {
      fprintf(file,">UNCLUSTERED\n");
      for(j=0;j<tailles_clusters[nb_clusters];j++)
	{
	  fprintf(file,">%s\n",clusters[nb_clusters][j]->nom);
	  for(k=0;k<longueur_alignement;k++)
	    {
	      if(sequences[clusters[nb_clusters][j]->id][k]=='O')
		{
		  fprintf(file,"-");
		}
	      else
		{
		  fprintf(file,"%c",sequences[clusters[nb_clusters][j]->id][k]);
		}
	    }
	  fprintf(file,"\n");
	}
    }
  fclose(file);
 
  /*desallocation memoire*/
  free(tailles_clusters);
  for(i=0;i<nb_clusters;i++)
    {
      free(clusters[i]);
    }
  free(clusters);
  /*fin desallocation memoire*/
}
