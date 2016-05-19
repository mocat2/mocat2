#include "principal.h"
#include "kmeans.h"
#include "matrice.h"
#include "mixturemodels.h"

/*************************************************************/
/*                                                           */
/*Procedure de classification en un nombre de classes donnees*/
/*                                                           */
/*************************************************************/

void mixturemodels(int nb_individus,int nb_dimensions,individu_t **individus,
		   int *nb_clusters,int nb_iterations_max,int nb_clusters_selection)
{
  /*declaration de variables*/
  int i,j,k,nb_clusters_max,meilleur_nombre_clusters,classification_reussie;
  double densite_loi_uniforme,meilleure_vraisemblance,vraisemblance,*valeurs_minimum,*valeurs_maximum;
  groupe_t *clusters,*meilleurs_clusters;
  /*fin declaration de variables*/

  nb_clusters_max=(int)sqrt((double)nb_individus);
  if((nb_clusters_max<nb_individus/(nb_dimensions+1))&&(nb_individus/(nb_dimensions+1)<60))
    {
      nb_clusters_max=nb_individus/(nb_dimensions+1);
    }

  /*desallocation memoire*/
  valeurs_minimum=(double *)malloc(sizeof(double)*nb_dimensions);
  valeurs_maximum=(double *)malloc(sizeof(double)*nb_dimensions);
  clusters=(groupe_t *)malloc(sizeof(groupe_t)*(nb_clusters_max+1));
  for(i=0;i<nb_clusters_max+1;i++)
    {
      clusters[i].individus=(individu_t **)malloc(sizeof(individu_t *)*nb_individus);
      clusters[i].mu=(double *)malloc(sizeof(double)*nb_dimensions);
      clusters[i].sigma=(double **)malloc(sizeof(double *)*nb_dimensions);
      for(j=0;j<nb_dimensions;j++)
	{
	  clusters[i].sigma[j]=(double *)malloc(sizeof(double)*nb_dimensions);
	}
    }
  meilleurs_clusters=(groupe_t *)malloc(sizeof(groupe_t)*(nb_clusters_max+1));
  for(i=0;i<nb_clusters_max+1;i++)
    {
      meilleurs_clusters[i].individus=(individu_t **)malloc(sizeof(individu_t *)*nb_individus);
    }
  for(i=0;i<nb_individus;i++)
    {
      individus[i]->vrais=(double *)malloc(sizeof(double)*(nb_clusters_max+1));
    }
  /*fin allocation memoire*/

  for(i=0;i<nb_dimensions;i++)
    {
      valeurs_minimum[i]=individus[0]->valeurs_traitees[i];
      valeurs_maximum[i]=individus[0]->valeurs_traitees[i];
    }
  
  for(i=0;i<nb_individus;i++)
    {
      for(j=0;j<nb_dimensions;j++)
	{
	  if(valeurs_minimum[j]>individus[i]->valeurs_traitees[j])
	    {
	      valeurs_minimum[j]=individus[i]->valeurs_traitees[j];
	    }
	  if(valeurs_maximum[j]<individus[i]->valeurs_traitees[j])
	    {
	      valeurs_maximum[j]=individus[i]->valeurs_traitees[j];
	    }
	}
    }
  densite_loi_uniforme=1.0;
  for(i=0;i<nb_dimensions;i++)
    {
      densite_loi_uniforme/=3.0*(valeurs_maximum[i]-valeurs_minimum[i]);
    }
  
  if(nb_clusters_selection==FIXED)
    {
      /*classification a nombre de clusters fixe*/
      classification_reussie=kmixturemodels(nb_individus,nb_dimensions,individus,*nb_clusters,
					    clusters,nb_iterations_max,&vraisemblance,
					    densite_loi_uniforme);
    }
  else
    {
      meilleure_vraisemblance=0;
      meilleur_nombre_clusters=-1;
      
      /*classification*/
      for(i=1;i<nb_clusters_max-1;i++)
	{
	  /*classification a nombre de clusters fixe*/
	  classification_reussie=kmixturemodels(nb_individus,nb_dimensions,individus,i,
						clusters,nb_iterations_max,&vraisemblance,
						densite_loi_uniforme);
	  printf("nb_groupes : %d\n",i);
	  if(classification_reussie==OUI)
	    { 
	      if(nb_clusters_selection==AIC)
		{
		  vraisemblance-=(double)((nb_dimensions+1)*nb_dimensions/2+nb_dimensions)*(double)i+
		    (double)(i+1);
		}
	      else
		{
		  vraisemblance-=((double)((nb_dimensions+1)*nb_dimensions/2+nb_dimensions)*(double)i+
		    (double)(i+1))*log((double)nb_individus);
		}
	      printf("vraisemblance : %f\n",vraisemblance);

	      if((meilleur_nombre_clusters==-1)||(vraisemblance>meilleure_vraisemblance))
		{
		  meilleure_vraisemblance=vraisemblance;
		  
		  meilleur_nombre_clusters=i;

		  for(j=0;j<meilleur_nombre_clusters+1;j++)
		    {
		      meilleurs_clusters[j].nb_individus=clusters[j].nb_individus;
		      for(k=0;k<clusters[j].nb_individus;k++)
			{
			  meilleurs_clusters[j].individus[k]=clusters[j].individus[k];
			}
		    }
		}
	    }
	}
      for(i=0;i<meilleur_nombre_clusters+1;i++)
	{
	  for(j=0;j<meilleurs_clusters[i].nb_individus;j++)
	    {
	      meilleurs_clusters[i].individus[j]->cluster=i;
	    }
	}
      *nb_clusters=meilleur_nombre_clusters+1;
    }
  
  /*desallocation memoire*/
  free(valeurs_minimum);
  free(valeurs_maximum);
  for(i=0;i<nb_clusters_max+1;i++)
    {
      free(clusters[i].individus);
      free(clusters[i].mu);
      for(j=0;j<nb_dimensions;j++)
	{
	  free(clusters[i].sigma[j]);
	} 
      free(clusters[i].sigma);
    }
  free(clusters);
  for(i=0;i<nb_clusters_max+1;i++)
    {
      free(meilleurs_clusters[i].individus);
    }
  free(meilleurs_clusters);
  for(i=0;i<nb_individus;i++)
    {
      free(individus[i]->vrais);
    }  
  /*fin desallocation memoire*/
}


/*************************************************************/
/*                                                           */
/*Procedure de classification en un nombre de classes donnees*/
/*par modeles de melanges                                    */
/*                                                           */
/*************************************************************/

int kmixturemodels(int nb_individus,int nb_dimensions,individu_t **individus,
		   int nb_clusters,groupe_t *clusters,int nb_iterations_max,
		   double *meilleure_vraisemblance,double densite_loi_uniforme)
{
  /*declaration de variables*/
  int i,j,k,meilleur_cluster,classification_reussie=NON,compteur_iterations=0,resultat;
  double **points,vraisemblance_precedente,vraisemblance_courante,log_vrai_max,trace;
  double vraisemblance_individu,vraisemblance_bruit;
  groupe_t *meilleurs_clusters;
  /*fin declaration de variables*/

  *meilleure_vraisemblance=0;

  /*allocation memoire*/
  points=(double **)malloc(sizeof(double *)*nb_individus);
  for(i=0;i<nb_individus;i++)
    {
      points[i]=(double *)malloc(sizeof(double)*nb_dimensions);
    }
  meilleurs_clusters=(groupe_t *)malloc(sizeof(groupe_t)*(nb_clusters+1));
  for(i=0;i<nb_clusters+1;i++)
    {
      meilleurs_clusters[i].individus=(individu_t **)malloc(sizeof(individu_t *)*nb_individus);
      meilleurs_clusters[i].mu=(double *)malloc(sizeof(double)*nb_dimensions);
      meilleurs_clusters[i].sigma=(double **)malloc(sizeof(double *)*nb_dimensions);
      for(j=0;j<nb_dimensions;j++)
	{
	  meilleurs_clusters[i].sigma[j]=(double *)malloc(sizeof(double)*nb_dimensions);
	}
    }
  /*fin allocation memoire*/

  /*classification initiale par kmeans*/
  resultat=kmeans(nb_individus,nb_dimensions,individus,nb_clusters);
  
  vraisemblance_precedente=0;
  while(1)
    {
      if(resultat==ECHEC)
	{
	  break;
	}

      for(i=0;i<nb_clusters+1;i++)
	{
	  clusters[i].nb_individus=0;
	}
      for(i=0;i<nb_individus;i++)
	{
	  clusters[individus[i]->cluster].individus
	    [clusters[individus[i]->cluster].nb_individus]=individus[i];
	  (clusters[individus[i]->cluster].nb_individus)++;
	}
      for(i=0;i<nb_clusters;i++)
	{
	  if(clusters[i].nb_individus<nb_dimensions+1)
	    { 
	      compteur_iterations=nb_iterations_max+1;
	      break;
	    }
	}

      if(compteur_iterations>nb_iterations_max)
	{
	  break;
	}

      /*on maximise la vraisemblance pour chaque cluster*/
      if(clusters[nb_clusters].nb_individus==0)
	{
	  vraisemblance_bruit=0.01*densite_loi_uniforme;
	}
      else
	{
	  vraisemblance_bruit=(double)clusters[nb_clusters].nb_individus/(double)nb_individus*
	    densite_loi_uniforme;
	}      
      
      for(i=0;i<nb_clusters;i++)
	{
	  for(j=0;j<clusters[i].nb_individus;j++)
	    {
	      for(k=0;k<nb_dimensions;k++)
		{
		  points[j][k]=clusters[i].individus[j]->valeurs_traitees[k];
		}
	    }
	  
	  /*maximum de vraisemblance pour une loi*/
	  /*calcul du mu moyen*/
	  calcul_mu_moyen(clusters[i].nb_individus,points,nb_dimensions,clusters[i].mu);

	  /*estimation de la matrice de covariance*/
	  calcul_matrice_covariance(nb_dimensions,clusters[i].nb_individus,points,
				    clusters[i].sigma);

	}

      /*etape CE ; on recompose les groupes*/
      calcul_vraisemblances_points(nb_individus,nb_dimensions,individus,nb_clusters,clusters);
      
      /*calcul vraisemblance du modele*/
      vraisemblance_courante=0;
      for(i=0;i<nb_individus;i++)
	{
	  vraisemblance_individu=vraisemblance_bruit;
	  for(j=0;j<nb_clusters;j++)
	    {
	      vraisemblance_individu+=exp(individus[i]->vrais[j]);
	    }
	  vraisemblance_courante+=log(vraisemblance_individu);
	}
      
      if((vraisemblance_courante>*meilleure_vraisemblance)||(compteur_iterations==0))
	{
	  *meilleure_vraisemblance=vraisemblance_courante;
	  for(i=0;i<nb_clusters+1;i++)
	    { 
	      meilleurs_clusters[i].nb_individus=clusters[i].nb_individus;
	      for(j=0;j<clusters[i].nb_individus;j++)
		{
		  meilleurs_clusters[i].individus[j]=clusters[i].individus[j];
		}
	    }
	}

      /*printf("###########################################################\n");
	printf("vraisemblance_courante : %f\n",vraisemblance_courante);
	printf("###########################################################\n");*/
      
      
      if(fabs(vraisemblance_courante-vraisemblance_precedente)<1.0)
	{
	  break;
	}
      vraisemblance_precedente=vraisemblance_courante;

      for(i=0;i<nb_individus;i++)
	{
	  log_vrai_max=individus[i]->vrais[0];
	  meilleur_cluster=0;
	  for(j=1;j<nb_clusters;j++)
	    {
	      if(log_vrai_max<individus[i]->vrais[j])
		{
		  log_vrai_max=individus[i]->vrais[j];
		  meilleur_cluster=j;
		}
	    }
	  if(log_vrai_max<log(vraisemblance_bruit))
	    {
	      /*probablement du bruit*/
	      meilleur_cluster=nb_clusters;
	    }
	  individus[i]->cluster=meilleur_cluster;
	}
      classification_reussie=OUI;
      compteur_iterations++;
    }

  if(classification_reussie==OUI)
    {
      for(i=0;i<nb_clusters+1;i++)
	{ 
	  clusters[i].nb_individus=meilleurs_clusters[i].nb_individus;
	  for(j=0;j<meilleurs_clusters[i].nb_individus;j++)
	    {
	      clusters[i].individus[j]=meilleurs_clusters[i].individus[j];
	    }
	}
    }

  /*desallocation memoire*/
  for(i=0;i<nb_individus;i++)
    {
      free(points[i]);
    }
  free(points);
  for(i=0;i<nb_clusters+1;i++)
    {
      free(meilleurs_clusters[i].individus);
      free(meilleurs_clusters[i].mu);
      for(j=0;j<nb_dimensions;j++)
	{
	  free(meilleurs_clusters[i].sigma[j]);
	}
      free(meilleurs_clusters[i].sigma);
    }
  free(meilleurs_clusters);
  /*fin desallocation memoire*/;
  
  return classification_reussie;
}

/***************************************************************/
/*                                                             */
/*Procedure de calcul des log-vraisemblances de tous les points*/
/*par rapport a tous les clusters                              */
/*                                                             */
/***************************************************************/

void calcul_vraisemblances_points(int nb_individus,int nb_dimensions,individu_t **individus,
				  int nb_clusters,groupe_t *clusters)
{
  /*declaration de variables*/
  int i,j,k,l;
  double **sigma,norme;
  double log_determinant,*difference,commun;
  /*fin declaration de variables*/

  /*allocation memoire*/
  sigma=(double **)malloc(sizeof(double *)*nb_dimensions);
  for(i=0;i<nb_dimensions;i++)
    {
      sigma[i]=(double *)malloc(sizeof(double)*nb_dimensions);
    }
  difference=(double *)malloc(sizeof(double)*nb_dimensions);
  /*fin allocation memoire*/
  
  for(i=0;i<nb_clusters;i++)
    {
      /*on copie la matrice de covariance dans sigma*/
      for(j=0;j<nb_dimensions;j++)
	{
	  for(k=0;k<nb_dimensions;k++)
	    {
	      sigma[j][k]=clusters[i].sigma[j][k];
	    }
	}

      /*on inverse la matrice de covariance et on recupere au passage*/
      /*le log du determinant de cette meme matrice*/
      log_determinant=inversion_matrice(nb_dimensions,sigma);
      
      commun=-log_determinant*0.5-(double)nb_dimensions*0.5*log(2.0*PI)+
	log((double)clusters[i].nb_individus/(double)nb_individus);

      /*calcul des vraisemblances*/
      for(j=0;j<nb_individus;j++)
	{
	  for(k=0;k<nb_dimensions;k++)
	    {
	      difference[k]=individus[j]->valeurs_traitees[k]-clusters[i].mu[k];
	    }
	  
	  /*individus[j]->vrais[i]=-log_determinant/2.0-(double)nb_dimensions/2.0*log(2.0*PI);*/
	  norme=calcul_norme(nb_dimensions,difference,sigma);
	  individus[j]->vrais[i]=commun-norme*0.5;
	  /*individus[j]->vrais[i]+=log((double)clusters[i].nb_individus/(double)nb_individus);*/
	}
    }

  /*desallocation memoire*/
  for(i=0;i<nb_dimensions;i++)
    {
      free(sigma[i]);
    }
  free(sigma);
  free(difference);
  /*fin desallocation memoire*/
}

/**************************************************/
/*                                                */
/*Procedure de calcul du mu moyen pour initialiser*/
/*la maximisation de la vraisemblance             */
/*                                                */
/**************************************************/

void calcul_mu_moyen(int nb_points,double **coordonnees,int nb_dimensions,
		     double *mu)
{
  /*declaration de variables*/
  int i,j;
  double norme;
  /*fin declaration de variables*/

  for(j=0;j<nb_dimensions;j++)
    {
      mu[j]=0;
    }

  for(i=0;i<nb_points;i++)
    {
      for(j=0;j<nb_dimensions;j++)
	{
	  mu[j]+=coordonnees[i][j];
	}
    }
  for(j=0;j<nb_dimensions;j++)
    {
      mu[j]/=(double)nb_points;
    }
}

