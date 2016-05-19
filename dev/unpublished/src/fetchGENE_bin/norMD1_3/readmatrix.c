#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "score.h"

static int getargs(char *line,char *args[],int max);

extern float matrix[NUMRES][NUMRES];

short gon250mt[]={
  24,
   0,   0,
   5,   0, 115,
  -3,   0, -32,  47,
   0,   0, -30,  27,  36,
 -23,   0,  -8, -45, -39,  70,
   5,   0, -20,   1,  -8, -52,  66,
  -8,   0, -13,   4,   4,  -1, -14,  60,
  -8,   0, -11, -38, -27,  10, -45, -22,  40,
  -4,   0, -28,   5,  12, -33, -11,   6, -21,  32,
 -12,   0, -15, -40, -28,  20, -44, -19,  28, -21,  40,
  -7,   0,  -9, -30, -20,  16, -35, -13,  25, -14,  28,  43,
  -3,   0, -18,  22,   9, -31,   4,  12, -28,   8, -30, -22,  38,
   3,   0, -31,  -7,  -5, -38, -16, -11, -26,  -6, -23, -24,  -9,  76,
  -2,   0, -24,   9,  17, -26, -10,  12, -19,  15, -16, -10,   7,  -2,  27,
  -6,   0, -22,  -3,   4, -32, -10,   6, -24,  27, -22, -17,   3,  -9,  15,  47,
  11,   0,   1,   5,   2, -28,   4,  -2, -18,   1, -21, -14,   9,   4,   2,  -2,  22,
   6,   0,  -5,   0,  -1, -22, -11,  -3,  -6,   1, -13,  -6,   5,   1,   0,  -2,  15,  25,
   1,   0,   0, -29, -19,   1, -33, -20,  31, -17,  18,  16, -22, -18, -15, -20, -10,   0,  34,
 -36,   0, -10, -52, -43,  36, -40,  -8, -18, -35,  -7, -10, -36, -50, -27, -16, -33, -35, -26, 142,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -22,   0,  -5, -28, -27,  51, -40,  22,  -7, -21,   0,  -2, -14, -31, -17, -18, -19, -19, -11,  41,   0,  78,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0};



void get_default_matrix(void)
{
   int gg_score = 0;
   int gr_score = 0;
   int i, j, k, ix = 0;
   int ti, tj;
   int  maxres;
   int av1,av2,av3,min, max;
   int max_aa;
   float scale=0.1;
	char *amino_acid_order = "ABCDEFGHIKLMNPQRSTVWXYZ";
   char *amino_acid_codes   =    "ABCDEFGHIJKLMNOPQRSTUVWXYZ-";
   static short  *xref, *matptr;
   short    def_aa_xref[NUMRES+1];

   char c1,c2;

   max_aa = strlen(amino_acid_codes)-2;

/*
   set up cross-reference for default matrices hard-coded in matrices.h
*/
   for (i=0;i<NUMRES;i++) def_aa_xref[i] = -1;

   maxres = 0;
   for (i=0;(c1=amino_acid_order[i]);i++)
     {
         for (j=0;(c2=amino_acid_codes[j]);j++)
          {
           if (c1 == c2)
               {
                  def_aa_xref[i] = j;
                  maxres++;
                  break;
               }
          }
     }



   matptr = gon250mt;
   xref=def_aa_xref;
/*
   default - set all scores to 0
*/
   for (i=0;i<=max_aa;i++)
      for (j=0;j<=max_aa;j++)
          matrix[i][j] = 0;

   ix = 0;
   maxres = 0;
   for (i=0;i<=max_aa;i++)
    {
      ti = xref[i];
      for (j=0;j<=i;j++)
       {
          tj = xref[j]; 
          if ((ti != -1) && (tj != -1))
            {
               k = matptr[ix];
               if (ti==tj)
                  {
                     matrix[ti][ti] = k * scale;
                     maxres++;
                  }
               else
                  {
                     matrix[ti][tj] = k * scale;
                     matrix[tj][ti] = k * scale;
                  }
               ix++;
            }
       }
    }

   --maxres;

   av1 = av2 = av3 = 0;
   for (i=0;i<=max_aa;i++)
    {
      for (j=0;j<=i;j++)
       {
           av1 += matrix[i][j];
           if (i==j)
              {
                 av2 += matrix[i][j];
              }
           else
              {
                 av3 += matrix[i][j];
              }
       }
    }

}

int readmatrix(FILE *fd)
{
   float ResComp[NUMRES][NUMRES];
   char Residues[NUMRES];
   char  c1,c2,last[NUMRES][10];
   int xref[NUMRES];
   int index[NUMRES];
   int i, j, k, ix = 0;
   int found = 0;
   int  maxres, maxcodes;
   float min;
   double f;
   char inline1[1024];
   int  numargs;
   char *args[NUMRES+4];

   while (fgets(inline1,1024,fd) != NULL)
     {
        i = strlen(inline1);
        if (inline1[i-2] == '*')
          {
            found = 1;
/*
   read residue characters.
*/
                 k = 0;
                 for (j=0;j<i-3;j++)
                    if ((inline1[j] >= 'A') && (inline1[j] <= 'Z')) 
                       Residues[k++] = inline1[j];
                 Residues[k] = '\0';
                 maxres = k;
            break;
          }
     }

   if (found == 0) 
     {
        fprintf(stderr,"Error: cannot find residues in matrix\n");
        return 0;
     }

   ix = 0;
   min = 100;
   while (fgets(inline1,1024,fd) != NULL)
     {
        if (inline1[0] == '\n') continue;
        numargs = getargs(inline1, args, maxres+1);
        if (numargs != maxres+1)
          {
             fprintf(stderr,"Error: wrong format in matrix\n");
             return 0;
          }
        for (i=0;i<maxres;i++)
          {
             f = atof(args[i]);
             ResComp[ix][i] = f;
             ResComp[i][ix] = ResComp[ix][i];
	     if (min>f) min = f;
          }
        ix++;
        if (ix==numargs) break;
     }

   if (ix != numargs)
     {
        fprintf(stderr,"Error: wrong lines in matrix\n");
        return 0;
     }

   for (i=0; i<maxres; i++)
	index[i]=Residues[i]-'A';

   for (i=0; i<26; i++)
        for (j=0; j<26; j++)
               matrix[i][j]=0;

   for (i=0; i<maxres; i++)
   	for (j=0; j<maxres; j++)
		matrix[index[i]][index[j]]=ResComp[i][j];

	return 1;

}

static int getargs(char *line,char *args[],int max)
{
 
        char    *inptr;
        int     i;
 
        inptr=line;
        for (i=0;i<=max;i++)
        {
                if ((args[i]=strtok(inptr," \t\n"))==NULL)
                        break;
                inptr=NULL;
        }
        if (i>max)
        {
                fprintf(stdout,"Too many args\n");
                return(0);
        }
 
        return(i);
}

