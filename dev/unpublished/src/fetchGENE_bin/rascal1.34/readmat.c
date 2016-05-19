#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "rascal.h"
#include "matrices.h"


/*
 *   Prototypes
 */
static Boolean commentline(char *line);
static sint read_user_matrix(char *filename,COMP_MATRIXPTR matrix,sint int_scale);
static sint read_user_matrix_series(char *filename);

/*
 *   Global variables
 */
static UserMatSeries matseries;		 /* a series of matrices from a user file */

static sint 	debug=0;

/*
   Matrices are stored as integer arrays, with the scores for residue 'A' in position 0,
   'B' in position 1, .... to 'Z' in position 25.
*/

/* reads comparison matrix scores from hard-coded arrays into COMP_MATRIX 
	dnaflag		TRUE if dna, FALSE if protein
	matptr		hard-coded, 1-dimensional array of residue-on-residue scores in matrices.h
	gapptr		hard-coded, 1-dimensional array of residue-on-gap scores in matrices.h
			(only used for bin matrices)
	neg_flag	if FALSE, add (-minimum_score) to all scores to make sure matrix
			contains only positive scores 
	scale		scale to convert float values to integer
	matrix		structure to store matrix

  returns number of residues in comparison matrix
*/
	

sint get_cl_matrix(Boolean dnaflag,sint *matptr, sint *gapptr, Boolean neg_flag, sint scale, COMP_MATRIXPTR matrix)
{
	sint i, j, k, ix = 0;
	sint ti, tj;
	sint  maxres;
	sint av3,min;
	char *order;
/*
   default - set all scores to 0
*/
	for (i=0;i<NUMRES;i++) {
          	matrix->gapscore[i] = 0;
      		for (j=0;j<NUMRES;j++)
          		matrix->score[i][j] = 0;
	}

	if(dnaflag) order=nucleic_acid_order;
	else order=amino_acid_order;

	ix = 0;
	maxres = strlen(order);
	for (i=0;i<maxres;i++) {
      		ti = tolower(order[i])-'a';
      		for (j=0;j<=i;j++) {
          		tj = tolower(order[j])-'a';
               		k = matptr[ix];
               		if (ti==tj) 
               			matrix->score[ti][ti] = k * scale;
               		else {
               			matrix->score[ti][tj] = k * scale;
               			matrix->score[tj][ti] = k * scale;
               		}
               		ix++;
       		}
    	}
/* get the residue-on-gap scores */
	if (matrix->format==1) {
		for (i=0;i<maxres;i++) {
      			ti = tolower(order[i])-'a';
			matrix->gapscore[ti] = gapptr[i] * scale;
		}
	}

	min = matrix->score[0][0];
   	for (i=0;i<NUMRES;i++) {
               	for (j=0;j<=i;j++) {
			if (matrix->score[i][j] < min) min = matrix->score[i][j];
		}
	}

   av3 = 0;
   for (i=0;i<NUMRES;i++)
    {
      for (j=0;j<=i;j++)
       {
           if (i!=j)
                 av3 += matrix->score[i][j];
       }
    }
	av3 /= ((float)(NUMRES*NUMRES-NUMRES))/2.0;
	matrix->mat_avscore = -av3;

if(debug>0) fprintf(stderr,"average mismatch %d %d\n",NUMRES,av3);

/*
   if requested, make a positive matrix - add -(lowest score) to every entry
*/
	if (neg_flag == FALSE && min < 0) {
if (debug>1) fprintf(stdout,"min %d\n",(pint)min);
   		for (i=0;i<maxres;i++) {
			if (isalpha(order[i])) {
      				ti = tolower(order[i])-'a';
                 		for (j=0;j<maxres;j++) {
      					tj = tolower(order[j])-'a';
                    			matrix->score[ti][tj] -= min;
                   		}
			}
                }
	}

	return(maxres);
}

/* reads comparison matrix scores from one of a series of matrices into COMP_MATRIX 
	userfile	file name of user input file
	pcid		percent identity of sequences
	negative	if FALSE, add (-minimum_score) to all scores to make sure matrix
			contains only positive scores 
	int_scale	scale to convert float values to integer
	matrix		structure to store matrix

  returns number of residues in comparison matrix
*/

sint get_user_matrix_series(char *userfile, double pcid, Boolean negative, sint int_scale, COMP_MATRIXPTR matrix)
{
   sint  i,j,maxres;
   static Boolean found;
   static Boolean error_given=FALSE;

   if(matseries.user_series==FALSE) {
	maxres=get_user_matrix(userfile, negative, int_scale, matrix);
	return maxres;
   }
 
   found=FALSE;
   for(i=0;i<matseries.nmat;i++)
        if(pcid>=matseries.mat[i].llimit && pcid<=matseries.mat[i].ulimit)
        {
                j=i;
                found=TRUE;
                break;
        }
   if(found==FALSE)
   {
        if(!error_given)
        warning(
"\nSeries matrix not found for sequence percent identity = %d.\n"
"(Using first matrix in series as a default.)\n"
"This alignment may not be optimal!\n"
"SUGGESTION: Check your matrix series input file and try again.",(int)pcid);
        error_given=TRUE;
        j=0;
   }

if (debug>0) fprintf(stdout,"pcid %d  matrix %d\n",(pint)pcid,(pint)j+1);

  maxres = get_user_matrix(matseries.mat[j].mat_filename, negative, int_scale, matrix);
  return(maxres);
}

/* reads comparison matrix scores from a user input file into COMP_MATRIX 
	userfile	file name of user input file
	neg_flag	if FALSE, add (-minimum_score) to all scores to make sure matrix
			contains only positive scores 
	scale		scale to convert float values to integer
	matrix		structure to store matrix

  returns number of residues in comparison matrix
*/
sint get_user_matrix(char *userfile, Boolean neg_flag, sint scale, COMP_MATRIXPTR matrix)
{
	sint gg_score = 1;
	sint gr_score = 0;
	sint i, j, k, ix = 0;
	sint ti, tj;
	sint  maxres;
	sint av1,av2,av3,min, max;

	if((maxres=read_user_matrix(userfile,matrix,scale))==0) return 0;

	min = max = matrix->score[0][0];
	for (i=0;i<NUMRES;i++)
    	for (j=0;j<=i;j++) {
        	if (matrix->score[i][j] < min) min = matrix->score[i][j];
        	if (matrix->score[i][j] > max) max = matrix->score[i][j];
      	}

   av3 = 0;
   for (i=0;i<NUMRES;i++)
    {
      for (j=0;j<=i;j++)
       {
           if (i!=j)
                 av3 += matrix->score[i][j];
       }
    }
	av3 /= ((float)(NUMRES*NUMRES-NUMRES))/2.0;
	matrix->mat_avscore = -av3;
/*
   if requested, make a positive matrix - add -(lowest score) to every entry
*/
  	if (neg_flag == FALSE) {

if (debug>1) fprintf(stdout,"min %d max %d\n",(pint)min,(pint)max);
      		if (min < 0) {
           		for (i=0;i<NUMRES;i++) {
                 		for (j=0;j<NUMRES;j++) {
                    			matrix->score[i][j] -= min;
                   		}
                	}
            	}
        }
  	return(maxres);
}

/*
   Reads a series of matrices from a user file
*/

static sint read_user_matrix_series(char *filename)
{
   FILE *fd = NULL;
   char mat_filename[FILENAMELEN];
   char inline1[1024];
   sint  maxres = 0;
   sint **matrix;
   sint nmat;
   sint n,llimit,ulimit;


   if (filename[0] == '\0')
     {
        error("comparison matrix not specified");
        return((sint)0);
     }
   if ((fd=fopen(filename,"r"))==NULL) 
     {
        error("cannot open %s", filename);
        return((sint)0);
     }

/* check the first line to see if it's a series or a single matrix */
   while (fgets(inline1,1024,fd) != NULL)
     {
        if (commentline(inline1)) continue;
	if(linetype(inline1,"CLUSTAL_SERIES"))
		matseries.user_series=TRUE;
	else
		matseries.user_series=FALSE;
        break;
     }

/* it's a single matrix */
  if(matseries.user_series == FALSE)
    {
	fclose(fd);
   	maxres=check_user_matrix(filename);
   	return(maxres);
    }

/* it's a series of matrices, find the next MATRIX line */
   nmat=0;
   matseries.nmat=0;
   while (fgets(inline1,1024,fd) != NULL)
     {
        if (commentline(inline1)) continue;
	if(linetype(inline1,"MATRIX"))
	{
		if(sscanf(inline1+6,"%d %d %s",&llimit,&ulimit,mat_filename)!=3)
		{
			error("Bad format in file %s\n",filename);
   			fclose(fd);
			return((sint)0);
		}
		if(llimit<0 || llimit > 100 || ulimit <0 || ulimit>100)
		{
			error("Bad format in file %s\n",filename);
   			fclose(fd);
			return((sint)0);
		}
		if(ulimit<=llimit)
		{
			error("in file %s: lower limit is greater than upper (%d-%d)\n",filename,llimit,ulimit);
   			fclose(fd);
			return((sint)0);
		}
   		n=check_user_matrix(mat_filename);
		if(n<=0)
		{
			error("Bad format in matrix file %s\n",mat_filename);
   			fclose(fd);
			return((sint)0);
		}
		matseries.mat[nmat].llimit=llimit;
		matseries.mat[nmat].ulimit=ulimit;
		strcpy(matseries.mat[nmat].mat_filename,mat_filename);
		nmat++;
	}
    }
   fclose(fd);
   matseries.nmat=nmat;

   maxres=n;
   return(maxres);

}

/* Reads a single matrix from a user file */

static sint read_user_matrix(char *filename,COMP_MATRIXPTR matrix,sint int_scale)
{
	double f;
	FILE *fd;
	sint  numargs,farg;
	sint i, j, k = 0;
	char codes[NUMRES];
	char inline1[1024];
	char *args[NUMRES+4];
	char c1,c2;
	sint ix1, ix = 0;
	sint  maxres = 0;
	float scale;

	if (filename[0] == '\0') {
        	error("comparison matrix not specified");
       		return((sint)0);
     	}

   	if ((fd=fopen(filename,"r"))==NULL) {
       		error("cannot open %s", filename);
       		return((sint)0);
   	}
   	maxres = 0;
	matrix->format=0;
   	while (fgets(inline1,1024,fd) != NULL) {
		if(linetype(inline1,"#bin")) {
if(debug>0) fprintf(stdout,"got a bin matrix\n");
			matrix->format=1;
		}
        	if (commentline(inline1)) continue;
		if(linetype(inline1,"CLUSTAL_SERIES")) {
       			error("in %s - single matrix expected.", filename);
			fclose(fd);
       			return((sint)0);
   		}
/*
   read residue characters.
*/
        	k = 0;
        	for (j=0;j<strlen(inline1);j++) {
             		if (isalpha((int)inline1[j])) codes[k++] = inline1[j];
             		if (k>NUMRES) {
                   		error("too many entries in matrix %s",filename);
		   		fclose(fd);
                   		return((sint)0);
                	}
          	}
        	codes[k] = '\0';
        	break;
    	}

   	if (k == 0) {
        	error("wrong format in matrix %s",filename);
  		fclose(fd);
        	return((sint)0);
     	}


/*
   get the matrix scores
*/
   	ix=0;
   	while (fgets(inline1,1024,fd) != NULL) {
        	if (inline1[0] == '\n') continue;
        	if (inline1[0] == '#' || inline1[0] == '!') break;
        	numargs = getargs(inline1, args, (int)(k+1));
        	if (numargs < k) {
             		error("wrong format in matrix %s",filename);
  	     		fclose(fd);
             		return((sint)0);
          	}
        	if (isalpha(args[0][0])) farg=1;
        	else farg=0;
   		maxres = k-farg;

/* decide whether the matrix values are float or decimal */
		scale=1.0*int_scale;
		for(i=0;i<strlen(args[farg]);i++)
			if(args[farg][i]=='.') {
/* we've found a float value */
				break;
			}

        	for (i=0;i<k;i++) {
                  	f = atof(args[i+farg]);
                  	matrix->score[tolower(codes[ix])-'a'][tolower(codes[i])-'a'] = (sint)(f*scale);
                }
                f = atof(args[k+farg]);
                matrix->gapscore[tolower(codes[ix])-'a'] = (sint)(f*scale);
       		ix++;
		if(ix==k) break;
        }
   	if (ix != k) {
        	error("wrong format in matrix %s",filename);
  		fclose(fd);
        	return((sint)0);
     	}


  fclose(fd);

  return(maxres);
}

static Boolean commentline(char *line)
{
        int i;
 
        if(line[0] == '#') return TRUE;
        for(i=0;line[i]!='\n' && line[i]!=EOS;i++) {
                if(!isspace(line[i]))
			return FALSE;
        }
        return TRUE;
}

/*
   Read a user matrix file and check for format errors.
   If no filename is passed (ie. str[0]='\0'), prompts
   user to enter a filename. 
   Returns TRUE if format is good, otherwise FALSE.
*/
Boolean check_user_matrix(char *str)
{
        sint i,maxres;
	COMP_MATRIX matrix;
	char lin2[MAXLINE+1];
	sint scale=1;

        FILE *infile;
	Boolean ret;

        if(str[0]=='\0')
                getstr("Enter name of the matrix file",lin2);
        else
                strcpy(lin2,str);

        if(*lin2 == EOS) return FALSE;

        if((infile=fopen(lin2,"r"))==NULL) {
                error("Cannot find matrix file [%s]",lin2);
                return FALSE;
        }
	fclose(infile);

	strcpy(str, lin2);

	maxres = read_user_matrix(str,&matrix,scale);
        if (maxres <= 0) ret=FALSE;
	else ret=TRUE;

	return ret;
}

/* 
   Reads a series of matrices from a user file.
   Stores information in matseries for use by
   get_user_matrix_series().
*/
Boolean check_user_matrix_series(char *str,Boolean verbose)
{
        sint maxres;
	static char lin2[MAXLINE+1];

        FILE *infile;

        if(verbose)
                getstr("Enter name of the matrix file",lin2);
        else
                strcpy(lin2,str);

        if(*lin2 == EOS) return FALSE;

        if((infile=fopen(lin2,"r"))==NULL) {
                error("Cannot find matrix file [%s]",lin2);
                return FALSE;
        }
	fclose(infile);

	strcpy(str, lin2);

	maxres = read_user_matrix_series(str);
        if (maxres <= 0) return FALSE;

	return TRUE;
}

