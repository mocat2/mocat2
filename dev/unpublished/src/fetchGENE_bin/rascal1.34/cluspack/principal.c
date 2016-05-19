#include "principal.h"
#include "lecfi.h"
#include "filtering.h"
#include "kmeans.h"
#include "ecrfi.h"
#include "divide.h"
#include "treebionj.h"
#include "alignment.h"
#include "secatorbionj.h"
#include "treeward.h"
#include "treecut.h"
#include "treedpc.h"

/**********************************************/
/*                                            */
/*Programme principal du package de clustering*/
/*                                            */
/**********************************************/

void main(int argc,char **argv)
{
  /*declaration de variables*/
  int i,j,type_donnees,clustering_method,nb_clusters_selection,nb_clusters_selected=0;
  int type_densite,weighting,donnees_standardisees,standardisation,individus_orphelins;
  int write_coordinates,filtrage_individus,longueur_alignement,nb_individus,nb_dimensions;
  int nb_clusters,nb_iterations_max,individus_tous_identiques,nb_simulations;
  char *buf,*fichier_entree,**sequences,nom_fichier[TAILLE_NOM],*output_file;
  double distance_min,seuil_dissimilarite;
  individu_t **individus;
  noeud_t *arbre;
  /*fin declaration de variables*/

  /************************************/
  /*Traitement des parametres d'entree*/ 
  /************************************/

  /*initialisation de la generation aleatoire de nombres*/
  srand(getpid());

  /*parametres par defaut*/
  weighting=NON;
  donnees_standardisees=NON;
  standardisation=NON;
  individus_orphelins=NON;
  write_coordinates=NON;
  filtrage_individus=NON;
  output_file=NULL;
  type_densite=DENSITE1;
  nb_iterations_max=-1;
  nb_simulations=10;
  /*fin parametres par defaut*/

  fichier_entree=strdup(argv[1]);
  if(strcmp(fichier_entree,"-help")==0)
    {
      printf("Program usage : cluspack file options\n\n");
      printf("********************OPTIONS BELOW********************\n\n");
      printf("-dt=[coordinates|alignment|alignmentNonProt|distances|similarities] (dt stands for data_type)\n");
      printf("-cm=[kmeans|ward|bionj|mixturemodels] (cm stands for clustering_method)\n");
      printf("-nbc=[secator|dpc|aic|bic|number]  (nbc stands for method for computing the number \
of cluster and number is really a number like 4 etc.)\n");
      printf("[-dt1|-dt2][-standardization][-standardized_data][-wc] (dt1 stands for density1 and wc for write_coordinates)\n");
      printf("[-fd=number] (dt stands for filtering_distance)\n");
      printf("[-nbsim=nbsimulations]");
      exit(0);
    }

  for(i=2;i<argc;i++)
    {
      if(strstr(argv[i],"-dt=")!=NULL)
	{
	  buf=strstr(argv[i],"=")+1;
	  if(strcmp(buf,"coordinates")==0)
	    {
	      type_donnees=COORDONNEES;
	    }
	  else if(strcmp(buf,"alignment")==0)
	    {
	      type_donnees=ALIGNEMENT;
	    }
	   else if(strcmp(buf,"alignmentNonProt")==0)
	    {
	      type_donnees=ALIGNEMENTNONPROT;
	    }
	  else if(strcmp(buf,"distances")==0)
	    {
	      type_donnees=DISTANCES;
	    }
	  else if(strcmp(buf,"similarities")==0)
	    {
	      type_donnees=SIMILARITES;
	    }
	  else 
	    {
	      printf("bad argument for -data_types !\n");
	      exit(1);
	    }
	}
      else if(strstr(argv[i],"-cm=")!=NULL)
	{
	  buf=strstr(argv[i],"=")+1;
	  if(strcmp(buf,"mm")==0)
	    {
	      clustering_method=MIXTURE_MODEL;
	    }
	  else if(strcmp(buf,"kmeans")==0)
	    {
	      clustering_method=KMEANS;
	    }
	  else if(strcmp(buf,"ward")==0)
	    {
	      clustering_method=WARD;
	    }
	  else if(strcmp(buf,"bionj")==0)
	    {
	      clustering_method=BIONJ;
	    }
	  else if(strcmp(buf,"mixturemodels")==0)
	    {
	      clustering_method=MIXTUREMODELS;
	    }
	  else if(strcmp(buf,"normalized_cut")==0)
	    {
	      clustering_method=NORMALIZED_CUT;
	    }
	  else 
	    {
	      printf("bad argument for -clustering_method !\n");
	      exit(1);
	    }
	}
      else if(strstr(argv[i],"-nbc=")!=NULL)
	{
	  buf=strstr(argv[i],"=")+1;
	  if(strcmp(buf,"dpc")==0)
	    {
	      nb_clusters_selection=DPC;
	    }
	  else if(strcmp(buf,"secator")==0)
	    {
	      nb_clusters_selection=SECATOR;
	    }
	  else if(strcmp(buf,"graphpc")==0)
	    {
	      nb_clusters_selection=GRAPHPC;
	    }
	  else if(strcmp(buf,"aic")==0)
	    {
	      nb_clusters_selection=AIC;
	    }
	  else if(strcmp(buf,"bic")==0)
	    {
	      nb_clusters_selection=BIC;
	    }
	  else
	    {
	      nb_clusters_selected=atoi(buf);
	      nb_clusters=nb_clusters_selected;
	      nb_clusters_selection=FIXED;
	    }
	}
      else if(strstr(argv[i],"-fd=")!=NULL)
	{
	  filtrage_individus=OUI;
	  distance_min=atof(strstr(argv[i],"=")+1);
	}
      else if(strstr(argv[i],"-dt1")!=NULL)
	{
	  type_densite=DENSITE1;
	}
      else if(strstr(argv[i],"-dt2")!=NULL)
	{
	  type_densite=DENSITE2;
	}
      else if(strstr(argv[i],"-weighting")!=NULL)
	{
	  weighting=OUI;
	}
      else if(strstr(argv[i],"-standardized_data")!=NULL)
	{
	  donnees_standardisees=OUI;
	}
      else if(strstr(argv[i],"-standardization")!=NULL)
	{
	  standardisation=OUI;
	  donnees_standardisees=OUI;
	}
      else if(strstr(argv[i],"-orphan_points")!=NULL)
	{
	  individus_orphelins=OUI;
	}
      else if(strstr(argv[i],"-wc")!=NULL)
	{
	  write_coordinates=OUI;
	}
      else if(strstr(argv[i],"-output")!=NULL)
	{
	  output_file=strdup(strstr(argv[i],"=")+1);
	} 
      else if(strstr(argv[i],"-nbsim=")!=NULL)
	{
	  nb_iterations_max=atof(strstr(argv[i],"=")+1);
	  nb_simulations=atoi(strstr(argv[i],"=")+1);
	}
    }

  /****************************************/
  /*Fin traitement des parametres d'entree*/
  /****************************************/

  /*********************/
  /*                   */
  /*LECTURE DES DONNEES*/
  /*                   */
  /*********************/

  /*lecture des donnees*/
  if((type_donnees==ALIGNEMENT)||(type_donnees==ALIGNEMENTNONPROT))
    {
      /*lecture du fichier de l'alignement multiple*/
      lecture_fichier_alignement(fichier_entree,&longueur_alignement,&nb_individus,
				 &individus,&sequences,type_donnees);

      nb_dimensions=nb_individus;

      if(nb_clusters_selection==SECATOR)
	{
	  /*calcul des distances entre les sequences*/
	  calcul_distances_sequences(nb_individus,individus,longueur_alignement,sequences);
	}
      else
	{
	  /*calcul des pourcentages d'identite entre les sequences*/
	  calcul_identites_sequences(nb_individus,individus,longueur_alignement,sequences);
	}
    }
  else 
    {
      /*lecture du fichier de coordonnees, similarites ou distances*/
      lecture_fichier_coordonnees(&nb_individus,&nb_dimensions,&individus,fichier_entree);
    }
  
  /***************************/
  /*                         */
  /*PRETRAITEMENT DES DONNEES*/
  /*                         */
  /***************************/

  if(filtrage_individus==OUI)
    {
      /*filtrage des individus suivant la valeur seuil donnee et le type des donnees*/
      filtre_individus(&nb_individus,nb_dimensions,individus,type_donnees,distance_min,nb_clusters_selection);

      if((type_donnees==SIMILARITES)||(type_donnees==ALIGNEMENT))
	{
	  nb_dimensions=nb_individus;
	}
    }
  if((nb_clusters_selection==GRAPHPC)&&(type_donnees==COORDONNEES))
    {
      /*calcul des similarites a partir des coordonnees en utilisant les correlations*/
      coordonnees_to_similarites(nb_individus,nb_dimensions,individus);
    }
  if(standardisation==OUI)
    {
      /*on normalise les individus*/
      standardisation_individus(nb_individus,nb_dimensions,individus);
    }

  /************/
  /*          */
  /*CLUSTERING*/
  /*          */
  /************/

  /*clustering method : KMEANS ou NORMALIZED_CUT*/
  if((clustering_method==KMEANS)||(clustering_method==NORMALIZED_CUT))
    {
      /*the number of clusters is fixed*/
      if(nb_clusters_selection==FIXED)
	{
	  kmeans(nb_individus,nb_dimensions,individus,nb_clusters_selected);
	  nb_clusters=nb_clusters_selected;
	}
      /*the number of clusters is automatically found by DPC*/
      else if((nb_clusters_selection==DPC)||(nb_clusters_selection==GRAPHPC))
	{
	  divide(nb_individus,nb_dimensions,&nb_clusters,individus,type_donnees,
		 clustering_method,nb_clusters_selection,donnees_standardisees,type_densite);
	}
    }
  /*clustering method : WARD*/
  else if(clustering_method==WARD)
    {
      classification_ward(nb_individus,nb_dimensions,individus,type_donnees,&arbre);
      
      /*the number of clusters is fixed*/
      if(nb_clusters_selection==FIXED)
	{
	  nb_clusters_to_dissimilarity_threshold(nb_individus,nb_clusters_selected,arbre,
						 &seuil_dissimilarite);
	}
      /*the number of clusters is found by Secator*/
      else if(nb_clusters_selection==SECATOR)
	{
	  secator_tree(nb_individus,arbre,&seuil_dissimilarite,&nb_clusters);
	}

      if((nb_clusters_selection==FIXED)||(nb_clusters_selection==SECATOR))
	{
	  /*creation des groupes*/
	  decoupage_arbre_seuil_dissimilarite(nb_individus,arbre,seuil_dissimilarite,&nb_clusters);
	}
      else
	{
	  decoupage_arbre_DPC(nb_individus,nb_dimensions,individus,arbre,&nb_clusters,
			      type_donnees,nb_clusters_selection,donnees_standardisees,
			      type_densite,nb_simulations);
	}
    }
  /*clustering method : mixture models*/
  else if(clustering_method==MIXTUREMODELS)
    {
      /*nb iterations for EM*/
      if(nb_iterations_max==-1)
	{
	  nb_iterations_max=60;
	}

      if((nb_clusters_selection==FIXED)||(nb_clusters_selection==AIC)||(nb_clusters_selection==BIC))
	{
	  /*the number of clusters is fixed*/
	  mixturemodels(nb_individus,nb_dimensions,individus,
			&nb_clusters,nb_iterations_max,nb_clusters_selection);
	}
    }
  /*clustering method : bionj*/
  else if(clustering_method==BIONJ)
    {
      /*on verifie que les individus ne sont pas tous identiques*/
      individus_tous_identiques=OUI;
      for(i=0;i<nb_individus;i++)
	{
	  for(j=i+1;j<nb_individus;j++)
	    {
	      if(individus[i]->valeurs_traitees[j]>ZERO_LIMIT)
		{
		  individus_tous_identiques=NON;
		}
	    }
	}

      if(individus_tous_identiques==OUI)
	{
	  nb_clusters=1;
	  for(i=0;i<nb_individus;i++)
	    {
	      individus[i]->cluster=0;
	    }
	}
      else
	{
	  if((type_donnees==ALIGNEMENT)||(type_donnees=DISTANCES))
	    {	 
	      creation_arbre_bionj(nb_individus,individus,fichier_entree,&arbre);
	    }
	  /*the number of clusters is found by Secator*/
	  if((nb_clusters_selection==SECATOR)||(nb_clusters_selection==FIXED))
	    {
	      secator_bionj(nb_individus,arbre,&nb_clusters,weighting,
			    nb_clusters_selection,nb_clusters_selected);
	    }
	}
    }
 
  /************************/
  /*                      */
  /*ECRITURE DES RESULTATS*/
  /*                      */
  /************************/

  /*ecriture du fichier resultat presentant les clusters*/
  if(output_file==NULL)
    {
      strcpy(nom_fichier,fichier_entree);
      if(strstr(nom_fichier,".")!=NULL)
	{
	  sprintf(strstr(nom_fichier,"."),".clu");
	}
      else
	{
	  strcat(nom_fichier,".clu");
	}
    }
  else
    {
      strcpy(nom_fichier,output_file);
    }

  /*ecriture du fichier resultat*/
  ecriture_fichier_clusters(nb_individus,nb_dimensions,individus,nb_clusters,nb_clusters_selection,
			    individus_orphelins,write_coordinates,nom_fichier);

  if(type_donnees==ALIGNEMENT)
    {
      if(output_file==NULL)
	{
	  strcpy(nom_fichier,fichier_entree);
	}
      else
	{
	  strcpy(nom_fichier,output_file);
	}
      if(strstr(nom_fichier,".")!=NULL)
	{
	  sprintf(strstr(nom_fichier,"."),"2.tfa");
	}
      else
	{
	  strcat(nom_fichier,"2.tfa");
	}
    
      /*ecriture du fichier de l'alignement multiple avec les sequences*/
      /*reordonnees suivant leurs groupes*/
      ecriture_fichier_clusters_alignement(nom_fichier,longueur_alignement,nb_individus,
					   individus,sequences,individus_orphelins,nb_clusters);
    }

  /*desallocation memoire*/
  for(i=0;i<nb_individus;i++)
    {
      for(j=0;j<individus[i]->nb_individus_similaires;j++)
	{
	  free(individus[i]->individus_similaires[j]->valeurs_brutes);
	  free(individus[i]->individus_similaires[j]->valeurs_traitees);
	  free(individus[i]->individus_similaires[j]->nom);
	  free(individus[i]->individus_similaires[j]);
	}
      if(individus[i]->nb_individus_similaires>0)
	{
	  free(individus[i]->individus_similaires);
	}
      free(individus[i]->nom);
      free(individus[i]->valeurs_brutes);
      free(individus[i]->valeurs_traitees);
    }
  /*fin desallocation memoire*/

  exit(0);
}




