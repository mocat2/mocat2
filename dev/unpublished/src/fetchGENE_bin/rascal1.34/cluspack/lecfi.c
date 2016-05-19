#include "principal.h"
#include "lecfi.h"

/***********************************************************/
/*                                                         */
/*Procedure de lecture du fichier d'un arbre phylogenetique*/
/*(au format bionj)                                        */
/*                                                         */
/***********************************************************/

void lecture_arbre_phylogenetique(char *nom_fichier,char **arbre)
{
  /*declaration des variables*/
  FILE *fichier_entree;
  int nb_caracteres;
  char poubelle[TAILLE_POUBELLE];
  /*fin declaration des variables*/

  /*lecture du fichier de l'arbre*/
  fichier_entree=fopen(nom_fichier,"r");
    
  if(fichier_entree==NULL)
    {
      printf("file %s does not exist !\n",nom_fichier);
      exit(0);
    }
  
  /*allocation de memoire*/
  *arbre=(char *)malloc(sizeof(char)*TAILLE_POUBELLE);
  if(*arbre==NULL)
    {
      printf("echec de l'allocation de memoire\n");
      exit(0);
    }
  /*fin allocation de memoire*/
  
  (*arbre)[0]=0;
  nb_caracteres=0;

  while(fgets(poubelle,TAILLE_POUBELLE,fichier_entree)!=NULL)
    {
      nb_caracteres+=strlen(poubelle);

      /*reallocation de memoire*/
      (*arbre)=(char *)realloc(*arbre,sizeof(char)*(nb_caracteres+1));
      if(*arbre==NULL)
	{
	  printf("echec de l'allocation de memoire\n");
	  exit(0);
	}
      /*fin reallocation de memoire*/
      
      strcat(*arbre,poubelle);
    }
  fclose(fichier_entree);
  /*fin lecture du fichier de l'arbre*/
}

/********************************************/
/*                                          */
/*Procedure de lecture du fichier contenant */
/*les coordonnees des individus a classifier*/
/*Les individus sont stockes dans un tableau*/
/*de pointeurs sur structure                */
/*                                          */
/********************************************/

void lecture_fichier_coordonnees(int *nb_individus,int *nb_dimensions,individu_t ***individus,
				 char *nom_fichier_entree)
{
  /*declaration des variables*/
  int i,j;
  FILE *fichier;
  char nom_individu[TAILLE_NOM],ligne[TAILLE_MAX_LIGNE];
  /*fin declaration des variables*/

  fichier=fopen(nom_fichier_entree,"r");
  fscanf(fichier,"%d",nb_individus);
  fscanf(fichier,"%d",nb_dimensions);
  fgets(ligne,TAILLE_MAX_LIGNE,fichier);

  /*lecture de la ligne de description des colonnes*/
  fgets(ligne,TAILLE_MAX_LIGNE,fichier);
   
  /*allocation memoire*/
  *individus = (individu_t **)malloc(sizeof(individu_t *)*(*nb_individus));
  for(i=0;i<*nb_individus;i++)
    {
      (*individus)[i]=(individu_t *)malloc(sizeof(individu_t));
      (*individus)[i]->id=i;
      (*individus)[i]->cluster=-1;
      (*individus)[i]->nb_individus_similaires=0;
      (*individus)[i]->valeurs_brutes = (double *)malloc(sizeof(double)*(*nb_dimensions));
      (*individus)[i]->valeurs_traitees = (double *)malloc(sizeof(double)*(*nb_dimensions));
      (*individus)[i]->cluster=AUCUN_GROUPE;
    }
  /*fin allocation memoire*/
  for(i=0;i<*nb_individus;i++)
    {
      fscanf(fichier,"%s",nom_individu);
      (*individus)[i]->nom = (char *)malloc(sizeof(char)*(strlen(nom_individu)+1));
      strcpy((*individus)[i]->nom,nom_individu);
      for(j=0;j<*nb_dimensions;j++)
	{
	  fscanf(fichier,"%lf",&((*individus)[i]->valeurs_brutes[j]));
	  (*individus)[i]->valeurs_traitees[j]=(*individus)[i]->valeurs_brutes[j];
	}
      fgets(ligne,TAILLE_MAX_LIGNE,fichier);
      (*individus)[i]->description=(char *)malloc(sizeof(char)*(strlen(ligne)+10));
      strcpy((*individus)[i]->description,ligne);
    }
  fclose(fichier);
}


/**********************************************************/
/*                                                        */
/*Procedure de lecture du fichier d'un alignement multiple*/
/*                                                        */
/**********************************************************/

void lecture_fichier_alignement(char *nom_alignement,int *longueur_alignement,int *nb_individus,
				individu_t ***individus,char ***sequences,int type_donnees)
{  
  /*declaration des variables*/
  int i,type_fichier;
  FILE *fichier;
  char poubelle[TAILLE_POUBELLE];
  /*fin declaration des variables*/
  
  
  /*on regarde si on a affaire a un fichier.msf ou a un fichier.tfa*/
  fichier=fopen(nom_alignement,"r");
  type_fichier=FICHIER_TFA;
  fgets(poubelle,TAILLE_POUBELLE,fichier)!=NULL;
  if(poubelle[0]=='>')
    {
      type_fichier=FICHIER_TFA;
    }
  else
    {
      type_fichier=FICHIER_INCONNU;
      if((strstr(poubelle,"MSF:")!=NULL)&&(strstr(poubelle,"..")!=NULL))
	{
	  type_fichier=FICHIER_MSF;
	}
      while(fgets(poubelle,TAILLE_POUBELLE,fichier)!=NULL)
	{
	  if((strstr(poubelle,"MSF:")!=NULL)&&(strstr(poubelle,"..")!=NULL))
	    {
	      type_fichier=FICHIER_MSF;
	      break;
	    }
	}
    }
  fclose(fichier);
  if(type_fichier==FICHIER_INCONNU)
    {
      printf("File type must be either tfa or msf !\n");
      exit(0);
    }

  if(type_fichier==FICHIER_TFA)
    {       
      /*lecture du fichier.tfa*/
      lecture_fichier_tfa(nom_alignement,longueur_alignement,nb_individus,sequences,
			  individus,type_donnees);
    }
  else
    {
      /*lecture du fichier.msf*/
      lecture_fichier_msf(nom_alignement,longueur_alignement,nb_individus,sequences,
			  individus);      
    }
}

/***************************************/
/*                                     */
/*Procedure de lecture d'un fichier.tfa*/
/*                                     */
/***************************************/

void lecture_fichier_tfa(char *nom_alignement,int *longueur_alignement,
			 int *nb_individus,char ***sequences,individu_t ***individus,
			 int type_donnees)
{
  /*declaration des variables*/
  FILE *fichier_entree;
  int i,j,compteur,longueur_sequence_courante;
  char nom_fichier_tfa[200],poubelle[TAILLE_POUBELLE];
  /*fin declaration des variables*/

  sprintf(nom_fichier_tfa,"%s",nom_alignement);

  /*premiere lecture pour compter les sequences et calculer la longueur de l'alignement*/
  fichier_entree=fopen(nom_fichier_tfa,"r");
  
  if(fichier_entree==NULL)
    {
      printf("file %s does not exist !\n",nom_fichier_tfa);
      exit(0);
    }

  *nb_individus=0;
  longueur_sequence_courante=0;
  *longueur_alignement=0;
  while(fgets(poubelle,TAILLE_POUBELLE,fichier_entree)!=NULL)
    {
      if(poubelle[0]=='>')
	{
	  (*nb_individus)++;
	  if(longueur_sequence_courante>*longueur_alignement)
	    {
	      *longueur_alignement=longueur_sequence_courante;
	    }
	  longueur_sequence_courante=0;
	}
      else
	{
	  longueur_sequence_courante+=strlen(poubelle)-1;
	}
    }

  if(longueur_sequence_courante>*longueur_alignement)
    {
      *longueur_alignement=longueur_sequence_courante;
    }
  fclose(fichier_entree);
  /*fin premiere lecture pour compter les sequences et calculer la longueur de l'alignement*/

  /*allocation memoire*/
  *sequences=(char **)malloc(sizeof(char *)*(*nb_individus+1));
  for(i=0;i<*nb_individus+1;i++)
    {
      (*sequences)[i]=(char *)malloc(sizeof(char)*(*longueur_alignement));
    }
  *individus = (individu_t **)malloc(sizeof(individu_t *)*(*nb_individus));
  for(i=0;i<*nb_individus;i++)
    {
      (*individus)[i]=(individu_t *)malloc(sizeof(individu_t));
      (*individus)[i]->id = i;
      (*individus)[i]->cluster=-1;
      (*individus)[i]->nb_individus_similaires=0;
      (*individus)[i]->valeurs_brutes = (double *)malloc(sizeof(double)*(*nb_individus));
      (*individus)[i]->valeurs_traitees = (double *)malloc(sizeof(double)*(*nb_individus));
      (*individus)[i]->nom = (char *)malloc(sizeof(char)*TAILLE_NOM);
      (*individus)[i]->description=NULL;
    }
  /*fin allocation memoire*/

  /*seconde lecture pour lire l'alignement proprement dit*/
  fichier_entree=fopen(nom_fichier_tfa,"r");
  
  if(fichier_entree==NULL)
    {
      printf("file %s does not exist !\n",nom_fichier_tfa);
      exit(0);
    }
  
  fgets(poubelle,TAILLE_POUBELLE,fichier_entree);
  for(i=0;i<*nb_individus;i++)
    {
      sscanf(&(poubelle[1]),"%s",(*individus)[i]->nom);

      compteur=0;
      fgets(poubelle,TAILLE_POUBELLE,fichier_entree);
      while(poubelle[0]!='>')
	{
	  for(j=0;j<strlen(poubelle)-1;j++)
	    {
	      if((poubelle[j]-'A'<0)||((type_donnees==ALIGNEMENT)&&(poubelle[j]-'A'>24)))
		{
		  (*sequences)[i][compteur]='O';
		}
	      else if((type_donnees==ALIGNEMENT)&&((poubelle[j]=='B')||(poubelle[j]=='J')||(poubelle[j]=='O')||
		      (poubelle[j]=='U')||(poubelle[j]=='X')))
		{
		  (*sequences)[i][compteur]='O';
		}
	      else
		{
		  (*sequences)[i][compteur]=poubelle[j];
		}
	      compteur++;
	    }
	  if(fgets(poubelle,TAILLE_POUBELLE,fichier_entree)==NULL)
	    {
	      break;
	    }
	}
      for(j=compteur;j<*longueur_alignement;j++)
	{
	  (*sequences)[i][j]='O';
	}
    }

  fclose(fichier_entree);
  /*fin seconde lecture pour lire l'alignement proprement dit*/  

  for(i=0;i<*longueur_alignement;i++)
    {
      (*sequences)[*nb_individus][i]='P';
    }
}

/***************************************/
/*                                     */
/*Procedure de lecture d'un fichier.msf*/
/*                                     */
/***************************************/

void lecture_fichier_msf(char *nom_alignement,int *longueur_alignement,
			 int *nb_individus,char ***sequences,individu_t ***individus)
{
  /*declaration des variables*/
  FILE *fichier_entree;
  int i,j,position,compteur_sequences,position_a_gauche,nb_colonnes;
  char residu,nom_fichier_msf[200],poubelle[TAILLE_POUBELLE];
  char nom_premiere_sequence[200],mot_lu[200];
  /*fin declarations des variables*/

  /*premiere lecture du fichier pour compter le nombre de sequences*/
  sprintf(nom_fichier_msf,"%s",nom_alignement);

  fichier_entree=fopen(nom_fichier_msf,"r");

  if(fichier_entree==NULL)
    {
      printf("file %s does not exist !\n",nom_fichier_msf);
      exit(0);
    }

  fgets(poubelle,TAILLE_POUBELLE,fichier_entree);  
  while(strstr(poubelle,"Name")==NULL)
    {
      fgets(poubelle,TAILLE_POUBELLE,fichier_entree);  
    }
  sscanf(strstr(poubelle,"Name:")+5,"%s",nom_premiere_sequence);
  sscanf(strstr(poubelle,"Len:")+4,"%d",longueur_alignement);

  *nb_individus=0;
  while(strstr(poubelle,"Name")!=NULL)
    {
      (*nb_individus)++;
      fgets(poubelle,TAILLE_POUBELLE,fichier_entree);
    }
  fclose(fichier_entree);
  /*fin premiere lecture du fichier pour compter le nombre de sequences*/

  /*debut allocation de memoire*/
  *sequences=(char **)malloc(sizeof(char *)*(*nb_individus+1));
  if(*sequences==NULL)
    {
      printf("probleme d'allocation memoire\n");
      exit(0);
    }
  for(i=0;i<*nb_individus+1;i++)
    {
      (*sequences)[i]=(char *)malloc(sizeof(char)*(*longueur_alignement));
      if((*sequences)[i]==NULL)
	{
	  printf("probleme d'allocation memoire\n");
	  exit(0);
	}
    }
  *individus = (individu_t **)malloc(sizeof(individu_t *)*(*nb_individus));
  for(i=0;i<*nb_individus;i++)
    {
      (*individus)[i]=(individu_t *)malloc(sizeof(individu_t));
      (*individus)[i]->id = i;
      (*individus)[i]->cluster = -1; 
      (*individus)[i]->nb_individus_similaires=0;
      (*individus)[i]->valeurs_brutes = (double *)malloc(sizeof(double)*(*nb_individus));
      (*individus)[i]->valeurs_traitees = (double *)malloc(sizeof(double)*(*nb_individus));
      (*individus)[i]->nom = (char *)malloc(sizeof(char)*TAILLE_NOM);
      (*individus)[i]->description=NULL;
      (*individus)[i]->nb_individus_similaires=0;
     }
  /*fin allocation de memoire*/

  /*deuxieme lecture*/
  fichier_entree=fopen(nom_fichier_msf,"r");
  
  fgets(poubelle,TAILLE_POUBELLE,fichier_entree);  
  while(strstr(poubelle,"Name")==NULL)
    {
      fgets(poubelle,TAILLE_POUBELLE,fichier_entree);  
    }

  i=0;
  while(strstr(poubelle,"Name")!=NULL)
    {
      sscanf(strstr(poubelle,"Name")+5,"%s",(*individus)[i]->nom);
      i++;
      fgets(poubelle,TAILLE_POUBELLE,fichier_entree);
    }

  position_a_gauche=0;
  while(position_a_gauche<*longueur_alignement)
    {
      fscanf(fichier_entree,"%s",mot_lu);
      while(strcmp(mot_lu,nom_premiere_sequence)!=0)
	{
	  fscanf(fichier_entree,"%s",mot_lu);
	}
      position=position_a_gauche;
     
      /*  for(i=0;(i<50)&&(i+position_a_gauche<*longueur_alignement);i++)*/
      while(position<*longueur_alignement)
	{
	  fscanf(fichier_entree,"%c",&residu);
	  if(residu=='\n')
	    {
	      break;
	    }
	  else if(residu==' ')
	    {
	    }
	  else if((residu-'A'<0)||(residu-'A'>24))
	    {
	      /*O signifie gap*/
	      (*sequences)[0][position]='O';
	      position++; 
	    }
	  else if((residu=='B')||(residu=='J')||(residu=='O')||
		  (residu=='U')||(residu=='X'))
	    {
	      /*O signifie gap*/
	      (*sequences)[0][position]='O';
	      position++; 
	    }
	  else
	    {
	      (*sequences)[0][position]=residu;
	      position++;
	    }
	}
      for(compteur_sequences=1;compteur_sequences<*nb_individus;
	  compteur_sequences++)
	{
	  position=position_a_gauche;
	  fscanf(fichier_entree,"%s",poubelle);
	  /*	  for(i=0;(i<50)&&(i+position_a_gauche<*longueur_alignement);i++)*/
	  while(position<*longueur_alignement)
	    {
	      fscanf(fichier_entree,"%c",&residu);
	      if(residu=='\n')
		{
		  break;
		}
	      else if(residu==' ')
		{
		}
	      else if((residu-'A'<0)||(residu-'A'>24))
		{
		  /*O signifie gap*/
		  (*sequences)[compteur_sequences][position]='O';
		  position++; 
		}
	      else if((residu=='B')||(residu=='J')||(residu=='O')||
		      (residu=='U')||(residu=='X'))
		{
		  /*O signifie gap*/
		  (*sequences)[compteur_sequences][position]='O';
		  position++; 
		}
	      else
		{
		  (*sequences)[compteur_sequences][position]=residu;
		  position++;
		}
	    }
	}
       position_a_gauche=position;
    }

  fclose(fichier_entree);
  for(i=0;i<*longueur_alignement;i++)
    {
      (*sequences)[*nb_individus][i]='P';
    }
  /*fin deuxieme lecture*/
}
