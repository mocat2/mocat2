#include "principal.h"
#include "treebionj.h"
#include "bionj.h"

/******************************************************************************************/
/*                                                                                        */
/*Procedure construisant l'arbre bionj a partir d'un alignement de sequence dans sequences*/
/*                                                                                        */
/******************************************************************************************/

void creation_arbre_bionj(int nb_individus,individu_t **individus,
			  char *fichier_entree,noeud_t **racine)
{
  /*declaration de variables*/
  char nom_fichier_entree_bionj[TAILLE_NOM],nom_fichier_sortie_bionj[TAILLE_NOM],*arbre;
  /*fin declaration de variables*/

  strcpy(nom_fichier_entree_bionj,fichier_entree);
  if(strstr(nom_fichier_entree_bionj,".")!=NULL)
    {
      sprintf(strstr(nom_fichier_entree_bionj,"."),".dst");
    }
  else
    {
      strcat(nom_fichier_entree_bionj,".dst");
    }
  
  strcpy(nom_fichier_sortie_bionj,fichier_entree);
  if(strstr(nom_fichier_sortie_bionj,".")!=NULL)
    {
      sprintf(strstr(nom_fichier_sortie_bionj,"."),".nj");
    }
  else
    {
      strcat(nom_fichier_sortie_bionj,".nj");
    }

  /*ecriture du fichier.dst des distances en vue du calcul de phylogenie*/
  ecriture_fichier_distances_sequences(nom_fichier_entree_bionj,nb_individus,individus);

  /*calcul de l'arbre phylogenetique par bionj*/
  bionj(nom_fichier_entree_bionj,nom_fichier_sortie_bionj);

  /*lecture de l'arbre phylogenetique*/
  lecture_arbre_phylogenetique(nom_fichier_sortie_bionj,&arbre);

  /*(re)construction de l'arbre phylogenetique*/
  construction_arbre(nb_individus+1,arbre,racine);

  /*enracinement de l'arbre*/ 
  enracinement_arbre(nb_individus,racine); 

  /*association des individus aux noeuds*/
  association_individus_noeuds(nb_individus,individus,*racine); 
}

/**************************************************/
/*                                                */
/*Procedure d'association des individus aux noeuds*/
/*                                                */
/**************************************************/

void association_individus_noeuds(int nb_individus,individu_t **individus,noeud_t *noeud)
{
  /*declaration de variables*/
  int i;
  char indice[100];
  /*fin declaration de variables*/

  for(i=0;i<nb_individus;i++)
    {
      sprintf(indice,"%d",individus[i]->id);
      if(strcmp(noeud->etiquette,indice)==0)
	{
	  noeud->individu=individus[i];
	  break;
	}
    }

  if(noeud->copains[1]!=NULL)
    {
      association_individus_noeuds(nb_individus,individus,noeud->copains[1]);
      association_individus_noeuds(nb_individus,individus,noeud->copains[2]);
    }
}

/*************************************************************/
/*                                                           */
/*Procedure d'enracinement de l'arbre sur la sequence fictive*/
/*                                                           */
/*************************************************************/

void enracinement_arbre(int nb_individus,noeud_t **racine)
{
  /*declaration des variables*/
  noeud_t *noeud_sequence_fictive;
  int nb_noeuds_internes,nb_feuilles;
  /*fin declaration des variables*/

  /*procedure de recherche de la sequence fictive dans l'arbre*/
  if(cherche_noeud_sequence_fictive(nb_individus,(*racine)->copains[1],
				    &noeud_sequence_fictive)==PAS_TROUVE)
    {
      cherche_noeud_sequence_fictive(nb_individus,(*racine)->copains[2],
				     &noeud_sequence_fictive);
    }

  *racine=noeud_sequence_fictive->copains[0];

  if((*racine)->copains[0]==NULL)
    {
      if((*racine)->copains[1]==NULL)
	{
	  *racine=(*racine)->copains[2];
	}
      else
	{
	  *racine=(*racine)->copains[1];
	}
    }
  else
    {
      if((*racine)->copains[1]==noeud_sequence_fictive)
	{
	  (*racine)->copains[1]=NULL;
	}
      else
	{
	  (*racine)->copains[2]=NULL;
	}
     
      /*on restructure l'arbre en partant de la sequence fictive*/
      restructuration_arbre(*racine,NULL,0);
    }

  nb_noeuds_internes=0;
  nb_feuilles=0;

  /*renommage des noeuds internes comme le nombre de noeuds a change*/
  renommage_noeuds_internes(nb_individus,*racine,&nb_noeuds_internes,&nb_feuilles); 
}


/**********************************************************************/
/*                                                                    */
/*Procedure de renommage des noeuds internes car apres restructuration*/
/*le nombre de noeuds n'est plus le meme                              */
/*                                                                    */
/**********************************************************************/

void renommage_noeuds_internes(int nb_individus,noeud_t *noeud,int *nb_noeuds_internes,
			       int *nb_feuilles)
{
  if(noeud->copains[1]!=NULL)
    {
      noeud->numero=*nb_noeuds_internes+nb_individus;
      sprintf(noeud->etiquette,"etiquette%d",*nb_noeuds_internes+nb_individus);
      (*nb_noeuds_internes)++;
      renommage_noeuds_internes(nb_individus,noeud->copains[1],nb_noeuds_internes,nb_feuilles);
      renommage_noeuds_internes(nb_individus,noeud->copains[2],nb_noeuds_internes,nb_feuilles);
    }
  else
    {
      noeud->numero=*nb_feuilles;
      (*nb_feuilles)++;
    }
}

/****************************************************************************/
/*                                                                          */
/*Procedure de restructuration de l'arbre et de calcul de la nouvelle racine*/
/*                                                                          */
/****************************************************************************/

void restructuration_arbre(noeud_t *noeud,noeud_t *nouveau_pere,double nouvelle_distance)
{
  /*declaration des variables*/
  noeud_t *ancien_pere,*ancien_fils_gauche,*ancien_fils_droite;
  double ancienne_distance;
  /*fin declaration des variables*/

  ancien_pere=noeud->copains[0];
  ancien_fils_gauche=noeud->copains[1];
  ancien_fils_droite=noeud->copains[2];
  ancienne_distance=noeud->distances[0];
  
  if(noeud->copains[0]->copains[0]!=NULL)
    {
      if(noeud->copains[1]==nouveau_pere)
	{
	  noeud->copains[1]=ancien_pere;
	}
      else
	{
	  noeud->copains[2]=ancien_pere;
	}
      noeud->copains[0]=nouveau_pere;
      noeud->distances[0]=nouvelle_distance;
      restructuration_arbre(ancien_pere,noeud,ancienne_distance);
    }
  /*on se trouve au niveau d'un fils de l'ancienne racine*/
  /*ancien_pere pointant sur l'ancienne racine*/
  else
    { 
      if(noeud->copains[1]==nouveau_pere)
	{
	  if(ancien_pere->copains[1]==noeud)
	    {
	      noeud->copains[1]=ancien_pere->copains[2];
	      ancien_pere->copains[2]->copains[0]=noeud;
	      ancien_pere->copains[2]->distances[0]+=ancienne_distance;
	    }
	  else
	    {
	      noeud->copains[1]=ancien_pere->copains[1];
	      ancien_pere->copains[1]->copains[0]=noeud;
		  ancien_pere->copains[1]->distances[0]+=ancienne_distance;
	    }
	}
      else
	{
	  if(ancien_pere->copains[1]==noeud)
	    {
	      noeud->copains[2]=ancien_pere->copains[2];
	      ancien_pere->copains[2]->copains[0]=noeud;
	      ancien_pere->copains[2]->distances[0]+=ancienne_distance;
	      ancien_pere->copains[1]->distances[0]+=ancienne_distance;
	    }
	  else
	    {
	      noeud->copains[2]=ancien_pere->copains[1];
	      ancien_pere->copains[1]->copains[0]=noeud;
	      
	    }
	}
      
      noeud->copains[0]=nouveau_pere;
      noeud->distances[0]=nouvelle_distance;
    }
}

/*******************************************************************/
/*                                                                 */
/*Procedure d'affichage de l'arbre (essentiellement pour debuggage)*/
/*                                                                 */
/*******************************************************************/

void affichage_arbre(noeud_t *noeud)
{
  printf("%s(%d) : \n",noeud->etiquette,noeud->numero);

  if(noeud->copains[1]!=NULL)
    {
      printf("fils1 : %s(%d) \n",noeud->copains[1]->etiquette,noeud->copains[1]->numero);
      printf("fils2 : %s(%d) \n",noeud->copains[2]->etiquette,noeud->copains[2]->numero);

      affichage_arbre(noeud->copains[1]);
      affichage_arbre(noeud->copains[2]);
    }
}

/************************************************************/
/*                                                          */
/*Procedure de recherche de la sequence fictive dans l'arbre*/
/*                                                          */
/************************************************************/

int cherche_noeud_sequence_fictive(int nb_individus,noeud_t *noeud,
				   noeud_t **noeud_sequence_fictive)
{
  if(noeud->copains[1]==NULL)
    {
      if(strcmp(noeud->etiquette,"fictive_sequence")==0)
	{
	  *noeud_sequence_fictive=noeud;
	  return TROUVE;
	}
      else
	{
	  return PAS_TROUVE;
	}
    }
  else
    {
      if(cherche_noeud_sequence_fictive(nb_individus,noeud->copains[1],
					noeud_sequence_fictive)==TROUVE)
	{
	  return TROUVE;
	}
      else if(cherche_noeud_sequence_fictive(nb_individus,noeud->copains[2],
					     noeud_sequence_fictive)==TROUVE)
	{
	  return TROUVE;
	}
      else
	{
	  return PAS_TROUVE;
	}
    }
}


/*****************************************************/
/*                                                   */
/*Procedure de construction de l'arbre phylogenetique*/
/*                                                   */
/*****************************************************/

void construction_arbre(int nb_individus,char *arbre,noeud_t **racine)
{
  /*declaration des variables*/
  int i,j,nombre_noeuds,compteur_individus;
  char *position_parenthese_ouvrante,*position_parenthese_fermante;
  char *position_premiers_double_points,*position_deuxiemes_double_points;
  char *position_troisiemes_double_points,*position_premiere_virgule,*position_seconde_virgule;
  char etiquette_gauche[200],etiquette_droite[200],etiquette_milieu[200];
  char nouvelle_etiquette[200],distance_gauche[200],distance_droite[200],distance_milieu[200];
  noeud_t *noeud_gauche,*noeud_milieu,*noeud_droite,*noeud_pere,**noeuds_courants;
  /*fin declaration des variables*/

  /*allocation de memoire*/
  noeuds_courants=(noeud_t **)malloc(sizeof(noeud_t *)*nb_individus);
  /*fin allocation de memoire*/

  /*construction de l'arbre*/
  nombre_noeuds=0;
  compteur_individus=0;
  for(i=0;i<nb_individus-3;i++)
    {
      /*on cherche la premiere parenthese fermante*/
      position_parenthese_fermante=strstr(arbre,")");

      j=-1;
      while(position_parenthese_fermante[j]!='(')
	{
	  j--;
	}
      position_parenthese_ouvrante=&(position_parenthese_fermante[j]);
      position_premiers_double_points=strstr(position_parenthese_ouvrante,":");
      position_premiere_virgule=strstr(position_premiers_double_points,",");
      position_deuxiemes_double_points=strstr(position_premiere_virgule,":");

      copie_partie_string(position_parenthese_ouvrante,position_premiers_double_points,
			  etiquette_gauche);

      copie_partie_string(position_premiere_virgule,position_deuxiemes_double_points,
			  etiquette_droite);

      copie_partie_string(position_premiers_double_points,position_premiere_virgule,
			  distance_gauche);
      
      copie_partie_string(position_deuxiemes_double_points,position_parenthese_fermante,
			  distance_droite);

      if(strstr(etiquette_gauche,"etiquette")==NULL)
	{
	  /*allocation memoire*/
	  noeud_gauche=(noeud_t *)malloc(sizeof(noeud_t));
	  /*fin allocation memoire*/

	  noeud_gauche->numero=compteur_individus;compteur_individus++;
	  strcpy(noeud_gauche->etiquette,etiquette_gauche);
	  
	  noeud_gauche->copains[1]=NULL;
	  noeud_gauche->copains[2]=NULL;
	  noeud_gauche->copain1=NULL;
	  noeud_gauche->copain2=NULL;
	}
      else
	{
	  for(j=0;j<nb_individus;j++)
	    {
	      if(strcmp((noeuds_courants[j])->etiquette,etiquette_gauche)==0)
		{
		  noeud_gauche=noeuds_courants[j];
		  break;
		}
	    }
	}

      if(strstr(etiquette_droite,"etiquette")==NULL)
	{
	  /*allocation memoire*/
	  noeud_droite=(noeud_t *)malloc(sizeof(noeud_t));
	  /*fin allocation memoire*/

	  noeud_droite->numero=compteur_individus;compteur_individus++;
	  strcpy(noeud_droite->etiquette,etiquette_droite);

	  noeud_droite->copains[1]=NULL;
	  noeud_droite->copains[2]=NULL;
	  noeud_droite->copain1=NULL;
	  noeud_droite->copain2=NULL;
	}
      else
	{
	  for(j=0;j<nb_individus;j++)
	    {
	      if(strcmp((noeuds_courants[j])->etiquette,etiquette_droite)==0)
		{
		  noeud_droite=noeuds_courants[j];
		  break;
		}
	    }
	}
      
      noeud_gauche->distances[0]=atof(distance_gauche);
      noeud_droite->distances[0]=atof(distance_droite);
      
      sprintf(nouvelle_etiquette,"etiquette%d",nb_individus+nombre_noeuds);

      /*allocation memoire*/
      noeud_pere=(noeud_t *)malloc(sizeof(noeud_t));
      /*fin allocation memoire*/

      strcpy(noeud_pere->etiquette,nouvelle_etiquette);
      noeud_pere->numero=nb_individus+nombre_noeuds;
      noeud_pere->copains[1]=noeud_gauche;
      noeud_pere->copains[2]=noeud_droite;
      noeud_pere->distances[1]=atof(distance_gauche);
      noeud_pere->distances[2]=atof(distance_droite);
      noeud_gauche->copains[0]=noeud_pere;
      noeud_droite->copains[0]=noeud_pere;
      
      noeuds_courants[nombre_noeuds]=noeud_pere;

      /*on remplace les deux noeuds par un nouveau noeud dans la chaine de caracteres arbre*/
      remplace_partie_string(&arbre,position_parenthese_ouvrante,position_parenthese_fermante,
			     nouvelle_etiquette);

      nombre_noeuds++;
    }

  position_premiers_double_points=strstr(arbre,":");
  position_premiere_virgule=strstr(position_premiers_double_points,",");
  position_deuxiemes_double_points=strstr(position_premiere_virgule,":");
  position_seconde_virgule=strstr(position_deuxiemes_double_points,",");
  position_troisiemes_double_points=strstr(position_seconde_virgule,":");
  position_parenthese_fermante=strstr(arbre,")");

  copie_partie_string(arbre,position_premiers_double_points,etiquette_gauche);
  copie_partie_string(position_premiere_virgule,position_deuxiemes_double_points,etiquette_milieu);
  copie_partie_string(position_seconde_virgule,position_troisiemes_double_points,etiquette_droite);
  copie_partie_string(position_premiers_double_points,position_premiere_virgule,
		      distance_gauche);
  copie_partie_string(position_deuxiemes_double_points,position_seconde_virgule,
		      distance_milieu);
  copie_partie_string(position_troisiemes_double_points,position_parenthese_fermante,
		      distance_droite);

  if(strstr(etiquette_gauche,"etiquette")==NULL)
    {
      /*allocation memoire*/
      noeud_gauche=(noeud_t *)malloc(sizeof(noeud_t));
      /*fin allocation memoire*/
      
      noeud_gauche->numero=compteur_individus;compteur_individus++;
      strcpy(noeud_gauche->etiquette,etiquette_gauche);

      noeud_gauche->copains[1]=NULL;
      noeud_gauche->copains[2]=NULL;
      noeud_gauche->copain1=NULL;
      noeud_gauche->copain2=NULL;
    }
  else
    {
      for(j=0;j<nb_individus;j++)
	{
	  if(strcmp((noeuds_courants[j])->etiquette,etiquette_gauche)==0)
	    {
	      noeud_gauche=noeuds_courants[j];
	      break;
	    }
	}
    }

   if(strstr(etiquette_milieu,"etiquette")==NULL)
     {
       /*allocation memoire*/
       noeud_milieu=(noeud_t *)malloc(sizeof(noeud_t));
       /*fin allocation memoire*/
       
       noeud_milieu->numero=compteur_individus;compteur_individus++;
       strcpy(noeud_milieu->etiquette,etiquette_milieu);
       noeud_milieu->copains[1]=NULL;
       noeud_milieu->copains[2]=NULL;
       noeud_milieu->copain1=NULL;
       noeud_milieu->copain2=NULL;
     }
   else
     {
       for(j=0;j<nb_individus;j++)
	 {
	   if(strcmp((noeuds_courants[j])->etiquette,etiquette_milieu)==0)
	     {
	       noeud_milieu=noeuds_courants[j];
	       break;
	     }
	 }
     }
   
   if(strstr(etiquette_droite,"etiquette")==NULL)
     {
       /*allocation memoire*/
       noeud_droite=(noeud_t *)malloc(sizeof(noeud_t));
       /*fin allocation memoire*/
       
       noeud_droite->numero=compteur_individus;compteur_individus++;
       strcpy(noeud_droite->etiquette,etiquette_droite);
       noeud_droite->copains[1]=NULL;
       noeud_droite->copains[2]=NULL;
       noeud_droite->copain1=NULL;
       noeud_droite->copain2=NULL;
     }
   else
     {
       for(j=0;j<nb_individus;j++)
	 {
	   if(strcmp((noeuds_courants[j])->etiquette,etiquette_droite)==0)
	     {
	       noeud_droite=noeuds_courants[j];
	       break;
	     }
	 }
     }

   /*allocation memoire*/
   noeud_pere=(noeud_t *)malloc(sizeof(noeud_t));
   /*fin allocation memoire*/
   
   sprintf(nouvelle_etiquette,"etiquette%d",nb_individus+nombre_noeuds);

   strcpy(noeud_pere->etiquette,nouvelle_etiquette);
  
   noeud_pere->numero=nb_individus+nombre_noeuds;
   nombre_noeuds++;

   noeud_gauche->distances[0]=atof(distance_gauche);
   noeud_milieu->distances[0]=atof(distance_milieu);
   noeud_droite->distances[0]=atof(distance_droite);
   noeud_pere->distances[0]=0;
   

   noeud_pere->copains[1]=noeud_gauche;
   noeud_pere->copains[2]=noeud_milieu;
   noeud_pere->distances[1]=atof(distance_gauche);
   noeud_pere->distances[2]=atof(distance_milieu);

   noeud_gauche->copains[0]=noeud_pere;
   noeud_milieu->copains[0]=noeud_pere;

   /*allocation memoire*/
   *racine=(noeud_t *)malloc(sizeof(noeud_t));
   /*fin allocation memoire*/

   sprintf(nouvelle_etiquette,"etiquette%d",nb_individus+nombre_noeuds);

   strcpy((*racine)->etiquette,nouvelle_etiquette);
   (*racine)->numero=nb_individus+nombre_noeuds;
   (*racine)->copains[1]=noeud_pere;
   (*racine)->copains[2]=noeud_droite;
   (*racine)->copains[0]=NULL;
   /*   (*racine)->distances[0]=0;
	(*racine)->distances[1]=0;
	(*racine)->distances[2]=atof(distance_droite);*/

   noeud_droite->copains[0]=*racine;
   noeud_pere->copains[0]=*racine;

   /*desallocation memoire*/
   free(noeuds_courants);
   /*fin desallocation memoire*/
}

/*******************************************************************************/
/*                                                                             */
/*Procedure de copie d'une portion de chaine de caracteres delimitee de maniere*/
/*stricte par deux pointeurs sur char dans une autre chaine de caracteres      */
/*                                                                             */
/*******************************************************************************/

void copie_partie_string(char *debut,char *fin,char *destination)
{
  /*declaration des variables*/
  int i;
  /*fin declaration des variables*/

  i=1;
  while(&(debut[i])!=fin)
    {
      destination[i-1]=debut[i];
      i++;
    } 
  destination[i-1]=0;
} 



/******************************************************************************/
/*                                                                            */
/*Procedure de remplacement d'une partie d'une chaine de caracteres, delimitee*/
/*de maniere large par debut et fin, par une autre chaine de caracteres       */
/*                                                                            */
/******************************************************************************/

void remplace_partie_string(char **chaine,char *debut,char *fin,char *chaine_de_remplacement)
{
  /*declaration des variables*/
  char *premiere_partie,*seconde_partie;
  /*fin declaration des variables*/

  /*allocation memoire*/
  premiere_partie=(char *)malloc(sizeof(char)*(strlen(*chaine)+1));
  seconde_partie=(char *)malloc(sizeof(char)*(strlen(*chaine)+1));
  /*fin allocation memoire*/

  strcpy(premiere_partie,*chaine); 
  premiere_partie[strlen(*chaine)-strlen(debut)]=0;
  strcpy(seconde_partie,&(fin[1]));

  *chaine=realloc(*chaine,sizeof(char)*(strlen(premiere_partie)+strlen(chaine_de_remplacement)+
					strlen(seconde_partie)+1));
  (*chaine)[0]=0;
  strcat(*chaine,premiere_partie);
  strcat(*chaine,chaine_de_remplacement);
  strcat(*chaine,seconde_partie);
  
  /*desallocation memoire*/
  free(premiere_partie);
  free(seconde_partie);
  /*fin desallocation memoire*/
}

