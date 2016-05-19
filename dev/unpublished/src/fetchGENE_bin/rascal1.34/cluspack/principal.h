#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*constantes booleennes*/
#define OUI 1
#define NON 0
#define NOT_DEFINED -1.0

/*reussite d'une fonction*/
#define SUCCES 1
#define ECHEC 0
 
/*constantes des types de donnees*/
#define COORDONNEES 0
#define ALIGNEMENT 1
#define ALIGNEMENTNONPROT 2
#define DISTANCES 3
#define SIMILARITES 4

/*constantes des programmes de clustering*/
#define MIXTURE_MODEL 0
#define BIONJ 1
#define WARD 2
#define KMEANS 3
#define NORMALIZED_CUT 4
#define MIXTUREMODELS 5

/*constantes des criteres de selection du nombre de groupe*/
#define DPC 0
#define SECATOR 1
#define GRAPHPC 2
#define AIC 3
#define BIC 4
#define FIXED 5

/*constantes de DPC*/
#define DENSITE1 0
#define DENSITE2 1

/*constantes de types de fichiers*/
#define FICHIER_TFA 0
#define FICHIER_MSF 1
#define FICHIER_INCONNU 2

/*constantes liees a la gestion de la memoire*/
#define TAILLE_NOM 200
#define TAILLE_POUBELLE 40000
#define TAILLE_MAX_LIGNE 40000

/*constantes de calcul*/
#define ZERO_LIMIT 0.000001
#define VALEUR_ENORME 999999999999.0

/*constantes diverses*/
#define AUCUN_GROUPE -1

#define PI 3.1415927

typedef struct individu_t {
  int id;
  char *nom;
  double *valeurs_brutes;
  double *valeurs_traitees;
  int nb_individus_similaires;
  int cluster;
  double *vrais;
  char *description;
  struct individu_t **individus_similaires;}individu_t;

/*type d'un noeud de l'arbre phylogenetique*/
typedef struct noeud_t{
  int numero;
  int groupe;
  int qualite;
  char etiquette[TAILLE_NOM];
  double poids;
  double dissimilarite;
  double distance;
  double distances[3]; /*en particulier, distances[0] est la distance au noeud superieur*/
  /*                     quand l'arbre est roote*/
  struct noeud_t *copain1;
  struct noeud_t *copain2;
  struct noeud_t *copains[3]; /*en particulier, copains[0] est le pere du noeud*/
  /*                            quand l'arbre est roote*/
  int secable;
  individu_t *individu;
}noeud_t;

typedef struct groupe_t {
  int nb_individus;
  individu_t **individus;
  int insecable;
  struct groupe_t *precedent;
  struct groupe_t *suivant;
  struct groupe_t *fils1;
  struct groupe_t *fils2;
  struct groupe_t *pere;
  double *mu;
  double **sigma;
  int qualite;}groupe_t;

typedef struct couple {
  int i;
  int j;}couple_t;










