#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PREC 8                             /* precision of branch-lengths  */
#define PRC  12
#define LEN  50                            /* length of taxon names        */

typedef struct word
{
    char name[LEN];
    struct word *suiv;
}WORD;

typedef struct pointers
{
   WORD *head;
	WORD *tail;
}POINTERS;


void   Initialize(float **delta, FILE *input, int n, POINTERS *trees);

void   Compute_sums_Sx(float **delta, int n);

void   Best_pair(float **delta, int r, int *a, int *b, int n);

void   Finish(float **delta, int n, POINTERS *trees, FILE *output);

void   Concatenate(char chain1[LEN], int ind, POINTERS *trees, int post);

void   Print_output(int i, POINTERS *trees, FILE *output);

float Distance(int i, int j, float **delta);

float Variance(int i, int j, float **delta);

float Sum_S(int i, float **delta);

float Agglomerative_criterion(int i, int j, float **delta, int r);

float Branch_length(int a, int b, float **delta, int r);

float Reduction4(int a, float la, int b, float lb, int i, float lamda,
		 float **delta);

float Reduction10(int a, int b, int i, float lamda, float vab, float
       **delta);
float Lamda(int a, int b, float vab, float **delta, int n, int r);

float Finish_branch_length(int i, int j, int k, float **delta);

int    Emptied(int i, float **delta);

int    Symmetrize(float **delta, int n);

void bionj(char *fichier_entree, char *fichier_sortie);
