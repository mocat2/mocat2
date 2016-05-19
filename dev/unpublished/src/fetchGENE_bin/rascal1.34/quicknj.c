/*  Last edited: Feb  1 18:22 2002 (klh) */
/***********************************************************************
 ** FILE: quicktree.c
 ** NOTES:
 **  This program does one of three things:
 ** 1. Given an aligment, construct a distance matrix.
 ** 2. Given a distance matrix, build a tree
 ** 3. Given an aligment, build a tree
 ** Which of these modes is active is dependent upon the command options
 ** 
 ***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <stdarg.h>
#include "rascal.h"
#include "quicknj.h"

/* this routine calculates a distance amtrix and an nj tree from an alignment */
void quicktree(FILE *output,SEQ *seqs,sint nseqs)
{

  unsigned int trial;
  unsigned int thisseq,i,j;
  struct Alignment *aln, *cons_aln;
  struct DistanceMatrix *mat;
  struct ClusterGroup *group;
  struct Tree *myTree, *testTree;
  static unsigned int use_kimura = 0;

  struct Tree *(*tree_func)(struct ClusterGroup *, unsigned int);
  /* step 1: copy the alignment and the distance matrix */
  aln = (struct Alignment *) malloc_util( sizeof( struct Alignment ));
  aln->numseqs = nseqs;
  aln->seqs = (struct Sequence **) malloc_util ( nseqs * sizeof( struct Sequence *));
  for(thisseq=0;thisseq<nseqs;thisseq++) {
    aln->seqs[thisseq] = empty_Sequence();
    aln->seqs[thisseq]->name = (char *) malloc_util( (strlen( seqs[thisseq].name ) + 1) * sizeof(char));
    strcpy( aln->seqs[thisseq]->name, seqs[thisseq].name );
    aln->seqs[thisseq]->length = seqs[thisseq].len;
    aln->seqs[thisseq]->seq = (char *) malloc_util( (seqs[thisseq].len + 1) * sizeof (char ) );
    for(j=0;j<seqs[thisseq].len;j++)
      aln->seqs[thisseq]->seq[j]=seqs[thisseq].data[j];
  }
  aln->length = seqs[0].len;

  group = alignment_to_ClusterGroup( aln, TRUE );
  aln = free_Alignment( aln );
  cons_aln = consensus_aln_from_ClusterGroup( group );
  group->matrix = empty_DistanceMatrix( group->numclusters );
  calc_DistanceMatrix( group->matrix, cons_aln, FALSE, use_kimura );


  /* step 2 produce tree */

  tree_func = &neighbour_joining_buildtree;
  
  myTree = (*tree_func)( group, FALSE ); 
  
  write_newhampshire_Tree( output, myTree, FALSE );

  aln = free_Alignment( aln );
  group = free_ClusterGroup( group );
  myTree = free_Tree( myTree );

}

/* this routine calculates the nj tree from a distance matrix */
void quicknj(FILE *output,SEQ *seqs,sint nseqs,double **tmat)
{

  unsigned int trial;
  unsigned int thisseq,i,j;
  struct Alignment *aln, *cons_aln;
  struct DistanceMatrix *mat;
  struct ClusterGroup *group;
  struct Tree *myTree, *testTree;

  struct Tree *(*tree_func)(struct ClusterGroup *, unsigned int);
  /* step 1: copy the alignment and the distance matrix */
  aln = (struct Alignment *) malloc_util( sizeof( struct Alignment ));
  aln->numseqs = nseqs;
  aln->seqs = (struct Sequence **) malloc_util ( nseqs * sizeof( struct Sequence *));
  for(thisseq=0;thisseq<nseqs;thisseq++) {
    aln->seqs[thisseq] = empty_Sequence();
    aln->seqs[thisseq]->name = (char *) malloc_util( (strlen( seqs[thisseq].name ) + 1) * sizeof(char));
    strcpy( aln->seqs[thisseq]->name, seqs[thisseq].name );
    aln->seqs[thisseq]->length = seqs[thisseq].len;
    aln->seqs[thisseq]->seq = (char *) malloc_util( (seqs[thisseq].len + 1) * sizeof (char ) );
    for(j=0;j<seqs[thisseq].len;j++)
      aln->seqs[thisseq]->seq[j]=seqs[thisseq].data[j];
  }
  aln->length = seqs[0].len;

  mat = empty_DistanceMatrix( nseqs );
  for(i=0;i<nseqs;i++) {
    for(j=0;j<=i;j++) {
      mat->data[i][j]=tmat[i][j];
    }
  }

  group = alignment_to_ClusterGroup( aln, FALSE );
  group->matrix = mat;

  /* step 2 produce tree */

  tree_func = &neighbour_joining_buildtree;
  
  myTree = (*tree_func)( group, FALSE ); 
  
  write_newhampshire_Tree( output, myTree, FALSE );

  aln = free_Alignment( aln );
  group = free_ClusterGroup( group );
  myTree = free_Tree( myTree );

}


/*  Last edited: Feb  1 18:15 2002 (klh) */
/**********************************************************************
 ** FILE: distancemat.c
 ** NOTES:
 **   Functions and types for the manipulation of Distance Matrices
 **********************************************************************/


/*********************** static variables *****************************/

/* This is a table of estimated PAMs for a range of percentage-differences,
   ranging from 75% dissimilarity, going up in 0.1% steps, to 93%
   dissimilarity. For percentage dissimilarty outside this range, we either
   use Kimura's formula (for < 75) or give an arbitrarily high distance 
   (for > 93) */

static int dayhoff_pams[]={
  195,   /* 75.0% observed d; 195 PAMs estimated = 195% estimated d */
  196,   /* 75.1% observed d; 196 PAMs estimated */
                  197,    198,    199,    200,    200,    201,    202,  203,    
  204,    205,    206,    207,    208,    209,    209,    210,    211,  212,    
  213,    214,    215,    216,    217,    218,    219,    220,    221,  222,    
  223,    224,    226,    227,    228,    229,    230,    231,    232,  233,    
  234,    236,    237,    238,    239,    240,    241,    243,    244,  245,    
  246,    248,    249,    250,    /* 250 PAMs = 80.3% observed d */          
                                  252,    253,    254,    255,    257,  258,    
  260,    261,    262,    264,    265,    267,    268,    270,    271,  273,    
  274,    276,    277,    279,    281,    282,    284,    285,    287,  289,    
  291,    292,    294,    296,    298,    299,    301,    303,    305,  307,    
  309,    311,    313,    315,    317,    319,    321,    323,    325,  328,    
  330,    332,    335,    337,    339,    342,    344,    347,    349,  352,    
  354,    357,    360,    362,    365,    368,    371,    374,    377,  380,    
  383,    386,    389,    393,    396,    399,    403,    407,    410,  414,    
  418,    422,    426,    430,    434,    438,    442,    447,    451,  456,    
  461,    466,    471,    476,    482,    487,    493,    498,    504,  511,    
  517,    524,    531,    538,    545,    553,    560,    569,    577,  586,    
  595,    605,    615,    626,    637,    649,    661,    675,    688,  703,    
  719,    736,    754,    775,    796,    819,    845,    874,    907,  945,
         /* 92.9% observed; 945 PAMs */    
  988    /* 93.0% observed; 988 PAMs */
};




/*********************************************************************
 FUNCTION: calc_DistanceMatrix
 DESCRIPTION: 
   Produces a distance matrix from the given multiple alignment
 RETURNS: struct DistanceMatrix
 ARGS: 
   A DistanceMatrix to fill in
   A multiple alignment
   A boolean indicating whether or not random columns should be used
     for purposes of bootstrapping
   A boolean indicating whether the Kimura distance adjustment is to 
       be used or not.
 NOTES: 
   0. the given DistanceMatrix and Alignment should be of the same order

   1. The matrix produced is in bottom-left triangular format; don't you
   go trying to access that top-right section (I'm warning you...)

   2. At the moment, the function calculates distance based on sequence
   identity, using Kimura's function if that option is raised.

   3. If use_rand_cols is true, then the matrix is constructed using
   random sampling  of columns, for the purposes of bootstrapping. At 
   the moment, the native function 'rand' is used to do this, suitable 
   seeded by time (by the caller). This may prove unsatisfactory...

   4. Where no information is available to determine the distance 
   between two sequences, a value of twice the maximum observed 
   distance is assigned (inspiration from ISMB99 poster by Huson,
   Smith and Warnow).

 *********************************************************************/

void calc_DistanceMatrix( struct DistanceMatrix *mat,
			  struct Alignment *aln,
			  unsigned int use_rand_cols,
			  unsigned int use_kimura ) {

  /* this function will take alignment and return a distance matrix */
  /* This gives a clear separation between tree making and distances */

  unsigned int i, j, k, table_index, num_undefined_distances, mem_increment;
  unsigned int *columnlist;
  double residuecount, distance, max_observed_distance, **undefined_distances;

  columnlist = (unsigned int *) malloc_util( aln->length * sizeof(unsigned int));
  for (i=0; i < aln->length; i++) {
    if (use_rand_cols) {
      /* generate random column here */
      columnlist[i] = (unsigned int) rand() % aln->length;
    }
    else {
      columnlist[i] = i;
    }
  }

  max_observed_distance = 0.0;
  mem_increment = 10;
  undefined_distances = NULL;
  num_undefined_distances = 0;

  for( i=0; i < aln->numseqs; i++) {
    mat->data[i][i] = 0.0;
    for( j=0; j < i; j++ ) {
      residuecount = distance = 0.0;
      mat->data[i][j] = 0.0;

      for( k=0; k < aln->length; k++) {
	if ( aln->seqs[i]->seq[columnlist[k]] == '.' || 
	     aln->seqs[j]->seq[columnlist[k]] == '.' ||
	     aln->seqs[i]->seq[columnlist[k]] == '-' || 
	     aln->seqs[j]->seq[columnlist[k]] == '-' ||
	     aln->seqs[i]->seq[columnlist[k]] == ' ' ||
	     aln->seqs[j]->seq[columnlist[k]] == ' ')
	  continue;

	/* neither character is a gap, so proceed */
	residuecount += 1.0;
	if ( aln->seqs[i]->seq[columnlist[k]] != aln->seqs[j]->seq[columnlist[k]]) 
	  distance += 1.0;
      }

      /* if residue count was zero here, there must have been a gap in every position;
	 in this case, %identity is undefined so a decision must be made to consider the
	 sequences 100% identical or 100% different. I go with the standard approach of
	 assigning twoce the maximum observed distance */

      if (residuecount > 0) {
	distance = distance / residuecount;
      }
      else {
	distance = -1.0;
	if (num_undefined_distances % mem_increment == 0) {
	  if (num_undefined_distances == 0) {
	    undefined_distances = (double **) 
	      malloc_util( sizeof( double *) * mem_increment);
	  }
	  else {
	    undefined_distances = (double **) 
	      realloc_util( undefined_distances, 
			    (num_undefined_distances + mem_increment) * sizeof(double *));
	  }
	}
	undefined_distances[num_undefined_distances++] = &(mat->data[i][j]);
	/* distances should be zero, but may be slightly less due to floating point 
	   arithmetic */
      }

      /* Use Kimura's formula to convert percentage dissimilarity to a distance */
      
      if (use_kimura) {
	if ( distance < 0.75) {
	  if (distance > 0.0) 
	    distance = - log( 1.0 - distance - (distance * distance * 0.20) );
	}
	else {
	  if (distance > 0.930) {
    	    distance = 10.0;
	  }
	  else {
	    table_index = (int) ((distance*1000.0) - 750.0);
	    distance = (double) dayhoff_pams[ table_index ];
	    distance /= 100.0;
	  }
	}
      }
      if (distance > max_observed_distance) {
	max_observed_distance = distance;
      }
      mat->data[i][j] = distance;
    }
  }
  
  /* before we go, lets find those undefined distances and replace them with 
     twice the maximum observed distance */

  for (i=0; i < num_undefined_distances; i++ ) {
    *undefined_distances[i] = 2 * max_observed_distance;
  }

  columnlist = free_util( columnlist );
  if (undefined_distances != NULL) {
    undefined_distances = free_util( undefined_distances );
  }
}



/*********************************************************************
 FUNCTION: clone_DistanceMatrix
 DESCRIPTION: 
   Produces a brand new DistanceMatrix, identical to the source
 RETURNS: struct DistanceMatrix
 ARGS: 
   A source distane matrix
 NOTES: 
   1. The matrix produced is in bottom-left triangular format; don't you
   go trying to access that top-right section (I'm warning you...)
 *********************************************************************/
struct DistanceMatrix *clone_DistanceMatrix( struct DistanceMatrix *source) {
  unsigned int i,j;
  struct DistanceMatrix *dest;

  if (source != NULL) {
    dest = empty_DistanceMatrix( source->size );
    
    for( i=0; i < dest->size; i++) {
      for( j=0; j <= i; j++ ) {
	dest->data[i][j] = source->data[i][j];	
      }
    }
  }
  else {
    dest = NULL;
  }

  return dest;
}



/*********************************************************************
 FUNCTION: empty_DistanceMatrix
 DESCRIPTION: 
   Produces an empty distance matrixof the given size, uninitialised
 RETURNS: struct DistanceMatrix
 ARGS: 
   The size of the matrix to be created
 NOTES: 
   1. The matrix produced is in bottom-left triangular format; don't you
   go trying to access that top-right section (I'm warning you...)
 *********************************************************************/
struct DistanceMatrix *empty_DistanceMatrix( unsigned int size) {
  unsigned int i;
  struct DistanceMatrix *mat;

  mat = (struct DistanceMatrix *) malloc_util(sizeof(struct DistanceMatrix));
  mat->size = size;
  mat->data = (double **) malloc_util( mat->size * sizeof(double *) );

  for( i=0; i < mat->size; i++)
    mat->data[i] = (double *) malloc_util( (i+1) * sizeof(double) );
 
  return mat;
}



/*********************************************************************
 FUNCTION: free_DistanceMatrix
 DESCRIPTION: 
   Frees the memory for the given distance matrix
 RETURNS:
 ARGS: 
   struct DistanceMatrix *
 NOTES: 
 *********************************************************************/

void *free_DistanceMatrix( struct DistanceMatrix *mat ) {
  int i;

  if ( mat != NULL ) {
    if (mat->data != NULL) {
      for( i=0; i < mat->size; i++ ) {
	if (mat->data[i] != NULL)
	  mat->data[i] = free_util( mat->data[i] );
      }
      mat->data = free_util( mat->data );
    }
    mat = free_util( mat );
  }

  return mat;
}




/********************************************************************** 
 FUNCTION: index_DistanceMatrix
 DESCRIPTION: 
   indexes the given distance matrix with the given indices,
   returning the appropraite distance.
 RETURNS: distance (double)
 ARGS: 
   A distance matrix *
   row index
   column index
 NOTES: 
   This function is necessary to account for the fact that the distance 
   matrix may be implemented as a symmtrical or triangular matrix.
   It therefore abstracts the internals of the distance matrix, at the
   cost of a function call for each lookup (is this wise...?)
 **********************************************************************/

double index_DistanceMatrix( struct DistanceMatrix *mat, 
			     unsigned int i, 
			     unsigned int j) {
  if (i > j) 
    return mat->data[i][j];
  else 
    return mat->data[j][i];
}


/*********************************************************************
 FUNCTION: print_DistanceMatrix
 DESCRIPTION: 
   Prints the given distance matrix.
 RETURNS:
 ARGS: 
   struct DistanceMatrix *
 NOTES: 
   A DistanceMatrix does not exist in isolation in practice but as
   part of a Cluster (this is to maintain the tight coupling between 
   the matrix and the sequences for which it is expressing the distances). 
   Therefore, to read or write a useful distance
   matrix (for compatibility with the phylip package for example)
   use write_phylip_Cluster
 *********************************************************************/

void print_DistanceMatrix( FILE *handle, struct DistanceMatrix *mat ) {
  unsigned int row, column;

  fprintf( handle, "Size:%d\n", mat->size);
  
  for(row=0; row < mat->size; row++) {
    fprintf( handle, "%5d", row);
    for(column=0; column <= row; column++)
	fprintf( handle, "%10.5f", mat->data[row][column]);
    fprintf( handle, "\n");
  }
  fflush( handle );
}




/********************************************************************* 
 FUNCTION: read_phylip_DistanceMatrix
 DESCRIPTION: 
   This function creates a DistanceMatrix from the given input file.
   It also crates a dummy alignment (sequences with just names) and
   puts it in the given Alignment pointer
 RETURNS: struct Cluster *
 ARGS: 
   A file handle
   A pointer to an Alignment pointer
 NOTES: 
   The file is assumed to be the distance matrix file format  used
   by the phlip package:

     4
  Name_1  0.0000   0.6776   0.6786  0.2342
  Name_2  0.6776   0.0000   0.1111  0.9999
  Name_3  0.6786   0.1111   0.0000  0.4444
  Name_4  0.2342   0.9999   0.4444  0.0000
 *********************************************************************/

struct DistanceMatrix *read_phylip_DistanceMatrix( FILE *handle, struct Alignment **aln_loc) {
  struct DistanceMatrix *mat;
  unsigned int size, i, j;
  char identifier[11];
  double dist;

  /* The size of the matrix will be on the first line on its own */
  if (! fscanf( handle, "%d", &size ))
    fatal_util( "Parse error: The first line should contain the size of matrix");

  *aln_loc = (struct Alignment *) malloc_util( sizeof(struct Alignment )); 

  (*aln_loc)->numseqs = size;
  (*aln_loc)->seqs = (struct Sequence **) 
    malloc_util( size * sizeof( struct Sequence *));
  (*aln_loc)->length = 0;
  mat = empty_DistanceMatrix( size );


  for (i=0; i < size; i++) {
    /* The name should be exactly 10 chars, and the scanf should place a \0
       at the end, making 11 */
    fscanf( handle, "%s", identifier );
    /* Right; the rest of the line will consist of exactly 'size; floating
       point numbers */
    (*aln_loc)->seqs[i] = empty_Sequence();
    (*aln_loc)->seqs[i]->name = (char *) malloc_util( 11 * sizeof(char));
    strcpy( (*aln_loc)->seqs[i]->name, identifier );
    for (j=0; j  < size; j++) {
      fscanf( handle, "%lf", &dist);
      if (j <= i) 
	mat->data[i][j] = dist;
    }
  } 

  return mat;
}



/********************************************************************* 
 FUNCTION: write_phylip_DistanceMatrix
 DESCRIPTION: 
   This function takes the given DistanceMatrix and writes it to the
   given file handle in phylip format. The alignment is needed for the
   Sequence names
   format
 RETURNS: 
 ARGS: 
   A file handle
   A DistanceMatrix pointer (cluster.h)
   An Alignment pointer
 NOTES: 
   The file is written in the distance matrix file format used
   by the phlip package:

     4
  Name_1  0.0000   0.6776   0.6786  0.2342
  Name_1  0.6776   0.0000   0.1111  0.9999
  Name_1  0.6786   0.1111   0.0000  0.4444
  Name_1  0.2342   0.9999   0.4444  0.0000
*********************************************************************/

void write_phylip_DistanceMatrix( FILE *handle, 
				  struct DistanceMatrix *mat,
				  struct Alignment *align) {

  unsigned int row, column;

  fprintf( handle, "\t%d\n", align->numseqs);
  
  for(row=0; row < align->numseqs; row++) {
    fprintf( handle, "%10.10s", align->seqs[row]->name);
    for(column=0; column < align->numseqs; column++) {
      if (row > column )
	fprintf( handle, "%10.5f", mat->data[row][column]);
      else 
	fprintf( handle, "%10.5f", mat->data[column][row]);
    }
    fprintf( handle, "\n");
  }
  fflush( handle );
}

/*  Last edited: Feb  1 17:52 2002 (klh) */
/**********************************************************************
 ** FILE: util.c
 ** NOTES:
 **   This file contains general utiliy functions used throughout the 
 **   application, such as those for memory management and error
 **   messaging. I have used it as a place to put other general stuff
 **   until I have somewhere better to put it.
 **********************************************************************/


/********************************************************************* 
 FUNCTION: calloc_util
 DESCRIPTION: 
   A wrapper for the stdlib.h function calloc; exits if memory cannot
   be allocated
 RETURNS: A generic pointer to the reallocated memory
 ARGS: 
   The number of bytes to be allocated
 NOTES:
 *********************************************************************/

void *calloc_util( size_t numobjs, size_t size) {
  void *ret;
	
  if ((ret = calloc( numobjs, size )) == NULL)
    fatal_util("calloc_util: Out of memory");

  return ret;	
}



/********************************************************************* 
 FUNCTION: malloc_util
 DESCRIPTION: 
   A wrapper for the stdlib.h function malloc; exits if memory cannot
   be allocated
 RETURNS: A generic pointer to the reallocated memory
 ARGS: 
   The number of bytes to be allocated
 NOTES:
 *********************************************************************/

void *malloc_util(size_t numbytes) {
  void *ret;
	
  if ((ret = malloc( numbytes )) == NULL)
    fatal_util("malloc_util: out of memory when requesting %d bytes", numbytes);

  return ret;	
}



/********************************************************************* 
 FUNCTION: realloc_util
 DESCRIPTION: 
   A wrapper for the stdlib.h function realloc; exits if memory cannot
   be allocated
 RETURNS: A generic pointer to the reallocated memory
 ARGS:
   A pointer to the memory to be reallocated
   The number of bytes to be allocated
 NOTES:
 *********************************************************************/

void *realloc_util(void *ptr, size_t bytes) {
  void *ret = NULL;

  if (ptr == NULL)
    fatal_util("Call to realloc_util with a null pointer");
  else { 
    if ((ret = realloc(ptr, bytes)) == NULL)
      fatal_util("realloc_util: out of memory");
  }

  return ret;
}  



/********************************************************************* 
 FUNCTION: free_util
 DESCRIPTION: 
   A wrapper for the stdlib.h function free; warns if the free could 
   not be performed
 RETURNS: A (hopefully null) pointer
 ARGS:
   A pointer to the memory to be deallocated
 NOTES:
 *********************************************************************/

void *free_util( void *ptr ) {
  if (ptr == NULL)
    warning_util("Call to free_util with null pointer");
  else {
    free(ptr);
    ptr = NULL;
  }
  return ptr;
}



/********************************************************************* 
 FUNCTION: fatal_util
 DESCRIPTION: 
   Prints the given formatted error message and exits
 RETURNS:
 ARGS:
   A format string + args, c.f. printf
 NOTES:
 *********************************************************************/

void fatal_util( char *fmt, ... ) {
  va_list args;
  
  va_start( args, fmt );
  fprintf( stderr, "\nA Fatal Error occurred: ");
  vfprintf( stderr, fmt, args);
  fprintf( stderr,"\n");
  va_end( args );
  exit(1);
}



/********************************************************************* 
 FUNCTION: warning_util
 DESCRIPTION: 
   Prints the given formatted warning to stderr
 RETURNS:
 ARGS:
   A format string + args, c.f. printf
 NOTES:
 *********************************************************************/

void warning_util( char *fmt, ...) {
  va_list args;
	
  va_start( args, fmt );
  fprintf( stderr, "\nWARNING: " );
  vfprintf( stderr, fmt, args );
  fprintf( stderr, "\n" );
  va_end( args );
}



/*  Last edited: Feb  1 15:35 2002 (klh) */
/**********************************************************************
 ** FILE: align.c
 ** NOTES:
 **   Functions for the manipulation of multiple sequence
 **   alignments
**********************************************************************/

static char *comment = "#";
static char *whitespace = " \t\r\n";
static char *terminator = "//";

/**********************************************************************
 FUNCTION: free_Alignment
 DESCRIPTION: 
   Frees the memory used by the given alignment
 ARGS:
   An Alignment
 RETURNS: A null pointer
 NOTES:
**********************************************************************/
void *free_Alignment( struct Alignment *al ) {
  unsigned int i;
  
  if (al != NULL) {
    if (al->seqs != NULL) {
      for( i=0; i < al->numseqs; i++) {
	al->seqs[i] = free_Sequence( al->seqs[i] );
      }
      al->seqs = free_util( al->seqs );
    }
    al = free_util( al );
  }
  return al;
}



/**********************************************************************
 FUNCTION: read_MUL_Alignment
 DESCRIPTION: 
   Reads in a muliple alignment from the given file handle and returns
   it
 ARGS: 
   A file handle
 RETURNS: A number denoting the status:
   A pointer to the struct Alignment object created, or NULL of there
   was an error parsing the file
 NOTES: It is assumed that the aligment file is in Pfam (MUL) format.
   Garbage results should be expected if the input file is not in this 
   format

   The function allocates all the memory necessary for the 
   alignment. The caller should call free_Alignment (align.h) to 
   free this memory when the alignment is no longer needed
 **********************************************************************/

struct Alignment *read_MUL_Alignment( FILE *stream ) {
  struct Alignment *aln;
  unsigned int index = 0;
  char tempname[MAX_NAME_LENGTH];
  int c;

  aln = (struct Alignment *) malloc_util( sizeof( struct Alignment ));
  aln->numseqs = 0;
  aln->seqs = NULL;
  aln->length = 0;
    

  while (( c = fgetc(stream)) != EOF) {
    if (aln->numseqs % SEQ_BLOCKSIZE == 0) {
      /* we need to allocate some more memory. But is it a malloc or realloc ? */

      if (aln->numseqs == 0) {
      	aln->seqs = (struct Sequence **) malloc_util( SEQ_BLOCKSIZE * sizeof( struct Sequence *));
      }
      else {
	aln->seqs = (struct Sequence **) 
	  realloc_util( aln->seqs, (aln->numseqs + SEQ_BLOCKSIZE) * sizeof( struct Sequence *));
      }
    }

    index = 0;
    while (! isspace(c) ) {
      tempname[index++] = c;
      c = fgetc(stream);
    }
    tempname[index] = '\0';
    while ( isspace(c = fgetc(stream)) );


    aln->seqs[aln->numseqs] = empty_Sequence();
    aln->seqs[aln->numseqs]->name = (char *) malloc_util( MAX_NAME_LENGTH * sizeof(char));
    strcpy(aln->seqs[aln->numseqs]->name, tempname);

    index = 0;
    do {
      if (index % RES_BLOCKSIZE == 0) {
      /* we need to allocate some more memory. But is it a malloc or realloc ? */
	if (index == 0) {
	  aln->seqs[aln->numseqs]->seq = (char *) malloc_util( RES_BLOCKSIZE * sizeof(char));
	}
	else {
	  aln->seqs[aln->numseqs]->seq = (char *) 
	    realloc_util( aln->seqs[aln->numseqs]->seq, (index + RES_BLOCKSIZE) * sizeof(char));
	}	
      }
      if (! isspace(c) )
	aln->seqs[aln->numseqs]->seq[index++] = c;
      c = fgetc(stream);

    } while ( c != '\n' && c != EOF);

    /* Since we now know the length of the sequence, we can do a bit of resizing */

    if (aln->numseqs == 0 || index < aln->length) {
      /* The first sequence is the trend setter */
      aln->length = index;
    }
    /* First, resize the current sequence */
    aln->seqs[aln->numseqs]->seq = (char *)
	   realloc_util( aln->seqs[aln->numseqs]->seq, aln->length * sizeof(char));

    /* if we reduced the aligment length, The earlier sequences are too big;
       Since we will only read up to the length of the smallest sequence, this
       will only prove a problem if we run out of memory; however, if the user 
       given non-flused alignments, then they are asking for everything they get...
    */

    aln->numseqs++;

  }

  /* now we can resize the seq array, and all the actual sequences themselves; this resizing
     will only save significant memory if a non-flush alignment has been given. Is it worth
     it? Well, we have to set the length of each sequence, so may as well do it while we are
     here*/

  aln->seqs = (struct Sequence **) realloc_util(  aln->seqs, aln->numseqs * sizeof( struct Sequence *));
  for (index=0; index < aln->numseqs; index++) {
    aln->seqs[index]->length = aln->length;
    aln->seqs[index]->seq = (char *) 
      realloc_util(  aln->seqs[index]->seq, aln->length * sizeof(char));
  }

  return aln;
}



/**********************************************************************
 FUNCTION: write_MUL_Alignment
 DESCRIPTION: 
   Prints a rep. of the alignment to the given handle (in MUL format)
 ARGS:
   FILE *
   struct Alignment (align.h)
 RETURNS: struct Alignment (align.h)
 NOTES:
 **********************************************************************/

void write_MUL_Alignment( FILE *handle, struct Alignment *al ) {
  unsigned int i,j;

  for( i=0; i < al->numseqs; i++ ) {
    fprintf( handle, "%-24s ", al->seqs[i]->name);
    for( j=0; j < al->length; j++) {
      fprintf( handle, "%c",  al->seqs[i]->seq[j]);
    }
    fprintf( handle, "\n");
  }
  fflush( handle );

}



/*  Last edited: Feb  1 18:21 2002 (klh) */
/**********************************************************************
 ** FILE: sequence.c
 ** NOTES:
 **   Functions and structures for the manipulation of protein 
 **   sequences. Only minimal functionality is needed
 **********************************************************************/


/**********************************************************************
 FUNCTION: clone_Sequence
 DESCRIPTION: 
   Performs a deep copy of the given Sequence and returns it
 ARGS:
   A Sequenced pointer (sequence.h)
 RETURNS:
   A pointer to a Sequence structure
 NOTES:
 **********************************************************************/

struct Sequence *clone_Sequence( struct Sequence *source) {
  unsigned int i;
  struct Sequence *dest = NULL;

  if (source != NULL) {
    dest = empty_Sequence();

    dest->name = (char *) malloc_util( (strlen(source->name)+1) * sizeof(char));
    strcpy( dest->name, source->name );
    dest->length = source->length;
    if (source->length) 
      dest->seq = (char *) malloc_util( dest->length * sizeof( char ));
    for( i=0; i < dest->length; i++) {
      dest->seq[i] = source->seq[i];
    }
    if (source->sec_struct != NULL) {
      dest->sec_struct = (char *) malloc_util( dest->length * sizeof( char ));
      for (i=0; i < dest->length; i++) {
	dest->sec_struct[i] = source->sec_struct[i];
      }
    }
    if (source->surf_acc != NULL) {
      dest->surf_acc = (char *) malloc_util( dest->length * sizeof( char ));
      for (i=0; i < dest->length; i++) {
	dest->surf_acc[i] = source->surf_acc[i];
      }
    }
    if (source->post_prob != NULL) {
      dest->post_prob = (char *) malloc_util( dest->length * sizeof( char ));
      for (i=0; i < dest->length; i++) {
	dest->post_prob[i] = source->post_prob[i];
      }
    }
    if (source->lig_bind != NULL) {
      dest->lig_bind = (char *) malloc_util( dest->length * sizeof( char ));
      for (i=0; i < dest->length; i++) {
	dest->lig_bind[i] = source->lig_bind[i];
      }
    }
  }

  return dest;
}



/**********************************************************************
 FUNCTION: empty_Sequence
 DESCRIPTION: 
   Creates and returns an new, empty Sequence object
 ARGS: 
 RETURNS:
   A pointer to a Sequence structure
 NOTES:
 **********************************************************************/
struct Sequence *empty_Sequence( void ) {
  struct Sequence *theseq;

  theseq = (struct Sequence *) malloc_util( sizeof(struct Sequence));
  theseq->length = 0;
  theseq->name = NULL;
  theseq->seq = NULL;
  theseq->sec_struct = NULL;
  theseq->surf_acc = NULL;
  theseq->post_prob = NULL;
  theseq->lig_bind = NULL;

  return theseq;
}




/**********************************************************************
 FUNCTION: free_Sequence
 DESCRIPTION: 
   frees the memory occupied by the given Sequence reference
 ARGS: 
   A pointer to a struct Sequence
 RETURNS: 
   The NULL pointer
 NOTES:
 **********************************************************************/
void *free_Sequence( struct Sequence *theseq ) {

  if (theseq != NULL) {
    if (theseq->name != NULL) 
      theseq->name = free_util( theseq->name);
    if (theseq->seq != NULL) 
      theseq->seq = free_util( theseq->seq);
    if (theseq->sec_struct != NULL)
      theseq->sec_struct = free_util( theseq->sec_struct);
    if (theseq->surf_acc != NULL)
      theseq->surf_acc = free_util( theseq->surf_acc);
    if (theseq->post_prob != NULL)
      theseq->post_prob = free_util( theseq->post_prob);
    if (theseq->lig_bind != NULL)
      theseq->lig_bind = free_util( theseq->lig_bind);

    theseq = free_util( theseq );
  }
  return theseq;
}






/*  Last edited: Mar 18 11:55 2002 (klh) */
/**********************************************************************
 ** FILE: tree.c
 ** NOTES:
 **  Functions for the manipulation of trees
 **********************************************************************/


/********************************************************************* 
 FUNCTION: assign_nodenumbers_Tnode
 DESCRIPTION: 
   This function assigns node numbers to the internal nodes of the 
   given Tnode, starting from the given, and returns the next free 
   nodenumner 
 RETURNS: unsigned int
 ARGS: 
   struct Tnode * (tree)
   unsigned int (starting node number)
 NOTES: 
 *********************************************************************/
unsigned int assign_nodenumbers_Tnode( struct Tnode *node, unsigned int start) {

  if (node != NULL) {
    start += assign_nodenumbers_Tnode( node->left, start );
    start += assign_nodenumbers_Tnode( node->right, start );

    node->nodenumber = start++;
  }

  return start; 
}



/********************************************************************* 
 FUNCTION: clone_Tnode
 DESCRIPTION: 
   This function makes a complete copy of the tree rooted at the given
   node and returns it
 RETURNS: Tnode *
 ARGS: 
   struct Tnode *
 NOTES: 
 *********************************************************************/
struct Tnode *clone_Tnode( struct Tnode *source) {
  struct Tnode *dest = NULL;

  if (source != NULL) {
    dest = (struct Tnode *) malloc_util( sizeof( struct Tnode ) );
    
    dest->distance = source->distance;
    dest->nodenumber = source->nodenumber;
    dest->bootstrap = source->bootstrap;
    dest->clust = clone_Cluster( source->clust );
    if ( (dest->left = clone_Tnode( source->left )) != NULL)
      dest->left->parent = dest;
    if ( (dest->right = clone_Tnode( source->right )) != NULL)
      dest->right->parent = dest;
  }

  return dest; 
}



/********************************************************************* 
 FUNCTION: clone_Tree
 DESCRIPTION: 
   This function makes a complete copy of the Tree and returns it
 RETURNS: struct Tree *
 ARGS: 
   struct Tree *
 NOTES: 
 *********************************************************************/
struct Tree *clone_Tree( struct Tree *source) {
  struct Tree *dest = NULL;

  if ( source != NULL ) {
    dest = (struct Tree *) malloc_util( sizeof( struct Tree ) );
    dest->numnodes = source->numnodes;
    dest->child[0] = clone_Tnode( source->child[0] );
    dest->child[1] = clone_Tnode( source->child[1] );
    dest->child[2] = clone_Tnode( source->child[2] );
  }

  return dest;

}




/********************************************************************* 
 FUNCTION: compare_to_bootstrap_sample_Tnode
 DESCRIPTION: 
   Updates the bootstrap values of the given Tnode, according to
   the topology of the given sample Tnode
 RETURNS:
 ARGS:
   Destination Tnode
   Sample Tnode
   The number of leaf nodes in each tree
   Boolean for whether the tree is binary or not
 NOTES:
   This function assumes that the given trees have been created by 
   calling either neighbourjoin_buildtree or UPGMA_buildtree with 
   the bootstrap boolean arguement set to true. A fatal error results
   if this is not the case. This is due to the fact that the information
   needed for the tree comparisons is determined when the trees are 
   constructed
 *********************************************************************/
void compare_to_bootstrap_sample_Tnode( struct Tnode *dest, 
					struct Tnode *sample,
					unsigned int numleaves,
					unsigned int is_binary) {
  unsigned int matchcounter, i;

  if (dest != NULL) {
    if (dest->child_ids != NULL) {

      /* we have to explore every non-terminal node of sample to see
	 if the bit pattern is the same, or a mirror image. If yes, then
	 this sub-tree exists in sample. The mirror image is to take care
	 of the fact that trichotomius trees may be isomorphic but 
	 'centred' at a different node.
	 
	 For non-trichotomous (i.e. rooted binary) trees, it is erroneous
	 to allow these 'mirror image' cases. It is easy to found out what
	 sort of trees we have by examining the child fields in one of the
	 trees. If only one of them is non-null, then we have a rooted 
	 binary tree.
      */

      if (sample != NULL) {
	if (sample->child_ids != NULL) {
	  matchcounter = 0;
	  for(i=0; i < numleaves; i++) {
	    if (dest->child_ids[i] == sample->child_ids[i]) {
	      matchcounter++;
	    }
	  }
	  if ((matchcounter == numleaves) || (matchcounter == 0 && ! is_binary)) {
	    dest->bootstrap++;
	    /* printf ("Incrementing node %d...\n", dest->nodenumber ); */
	  }
	  
	  compare_to_bootstrap_sample_Tnode( dest, sample->left, numleaves, is_binary );
	  compare_to_bootstrap_sample_Tnode( dest, sample->right, numleaves, is_binary );
	}
      }
    }
  }
}




/********************************************************************* 
 FUNCTION: empty_Tree
 DESCRIPTION: 
   Creates and returns a tree with null nodes
 RETURNS: struct Tree * (trees.h)
 ARGS: 
 NOTES:
 *********************************************************************/

struct Tree *empty_Tree( void ) {
  struct Tree *ret;

  ret = (struct Tree *) malloc_util(sizeof( struct Tree ));
  ret->child[0] = NULL;
  ret->child[1] = NULL;
  ret->child[2] = NULL;
  ret->numnodes = 0;

  return ret;
}



/********************************************************************* 
 FUNCTION: free_Tnode
 DESCRIPTION: 
   This function releases the memory used by this Tnode and all of its
   children
 RETURNS: A null pointer
 ARGS: 
   struct Tnode *
 NOTES: 
 *********************************************************************/

void *free_Tnode( struct Tnode *tn ) {
  if ( tn != NULL ) {
    tn->clust = free_Cluster( tn->clust );
    tn->left = free_Tnode( tn->left );
    tn->right = free_Tnode( tn->right );
    if (tn->child_ids != NULL) {
      tn->child_ids = free_util( tn->child_ids );
    }
    tn = free_util( tn );
  }
  return tn; 

}



/********************************************************************* 
 FUNCTION: free_Tree
 DESCRIPTION: 
   This function releases the memory used br the Tnode chain in the
   given Tree
 RETURNS: A null pointer
 ARGS: 
   struct Tree
 NOTES: 
 *********************************************************************/

void *free_Tree( struct Tree *t) {
  if ( t != NULL ) {
    t->child[0] = free_Tnode( t->child[0] );
    t->child[1] = free_Tnode( t->child[1] );
    t->child[2] = free_Tnode( t->child[2] );
    t = free_util( t );
  }
  return t;
}



/********************************************************************* 
 FUNCTION: get_root_Tnode
 DESCRIPTION: This function takes a three node tree (struct Tree), and
    returns a Tnode as the 'root' of the tree. The root is inserted
    somewhat arbitraily between the three top-level nodes
 RETURNS: struct Tnode *
 ARGS: struct Tree
 NOTES: 
   The nodes of the given tree are cloned for use. This means
   that the old tree is still available and safe on return, and must
   be freed when finished with. The rooted tree returned by this 
   function must be freed by a call to free_Tnode
 *********************************************************************/

struct Tree *get_root_Tnode( struct Tree *source ) {
  struct Tnode *focal, *root, *children[3];
  struct Tree *ret;
  unsigned int rootleft, focalleft, focalright;
  double maxdist;

  /***** Method **************
     0. Clone the given tree
     1. Create the imaginary node between the three nodes in the Tree
     2. Create a root node;
  **************************/

  children[0] = clone_Tnode( source->child[0] );
  children[1] = clone_Tnode( source->child[1] );
  children[2] = clone_Tnode( source->child[2] );

  focal = new_interior_Tnode( source->numnodes );
  root = new_interior_Tnode( source->numnodes + 1 );

  /* arbitrarity choose halfway along the longest branch between
     the three nodes as the position for the root */

  maxdist = children[0]->distance;
  rootleft = 0;
  focalleft = 1;
  focalright = 2;
  if (children[1]->distance > maxdist) {
    rootleft = 1;
    focalleft = 0;
    focalright = 2;
  }   
  if (children[2]->distance > maxdist) {
    rootleft = 2;
    focalleft = 0;
    focalright = 1;
  }

  /* sort out distances; root node has zero distances */

  children[rootleft]->distance = children[rootleft]->distance * 0.5;
  focal->distance = children[rootleft]->distance;

  /* now sort out the links */

  root->right = focal;
  root->left = children[rootleft];
  focal->parent = children[rootleft]->parent = root;

  focal->left = children[focalleft];
  focal->right = children[focalright];
  children[focalleft]->parent = children[focalright] = focal;

  ret = empty_Tree();
  ret->child[0] = root;
  ret->numnodes = source->numnodes + 2;

  return ret;
}



/********************************************************************* 
 FUNCTION: new_interior_Tnode (unsigned int)
 DESCRIPTION: 
   This function handles the simple task of allocating the space
   for a new interior (non-leaf) tree node, filling it with what it 
   knows, and returning it.
 RETURNS: struct Tnode * (trees.h)
 ARGS: 
   A integer to hold the node number.
 NOTES:
   This function s for creating internal nodes, which have no string
   id but just a number identifier
 *********************************************************************/

struct Tnode *new_interior_Tnode( unsigned int label ) {
  struct Tnode *newNode;

  newNode = (struct Tnode *) malloc_util(sizeof(struct Tnode));
  newNode->left = NULL;
  newNode->right = NULL;
  newNode->parent = NULL;
  newNode->distance = 0.0;
  newNode->nodenumber = label;
  newNode->clust = NULL;
  newNode->bootstrap = 0;
  newNode->child_ids = NULL;

  return newNode;
}





/********************************************************************* 
 FUNCTION: new_leaf_Tnode(int, char *)
 DESCRIPTION: 
   This function handles the simple task of allocating the space
   for a new tree node, filling it with what it knows, and 
   returning it.
 RETURNS: struct Tnode * (trees.h)
 ARGS: 
   A integer to hold the sequence number associated with the node
   The name of the node
 NOTES: 
   This function is for creating leaf nodes, which have a name.
 *********************************************************************/

struct Tnode *new_leaf_Tnode(unsigned int label, struct Cluster *given) {
  struct Tnode *newNode;

  newNode = new_interior_Tnode( label );
  newNode->clust = given;

  return newNode;
}



/********************************************************************* 
 FUNCTION: scale_bootstraps_Tnode
 DESCRIPTION: 
   This function traverses the given Tnode, dividing the bootstrap
   values by the given number
 RETURNS: 
 ARGS: 
   Tnode
   A the number of bootstrap iterations that were performed
 NOTES:

 *********************************************************************/
void scale_bootstraps_Tnode( struct Tnode *node, unsigned int iters) {

  if (node != NULL) {
    node->bootstrap = (int) (((double) node->bootstrap / (double) iters) * 100.0);
    scale_bootstraps_Tnode( node->left, iters);
    scale_bootstraps_Tnode( node->right, iters);
  }

}


/********************************************************************* 
 FUNCTION: scale_bootstraps_Tree
 DESCRIPTION: 
   This function traverses the gieven tree, dividing the bootstrap
   values by the given number
 RETURNS: 
 ARGS: 
   Tree
   A the number of bootstrap iterations that were performed
 NOTES:

 *********************************************************************/
void scale_bootstraps_Tree( struct Tree *thetree, unsigned int iters) {

  if (thetree != NULL) {
    /* The first node will always be defined...(he says) */
    scale_bootstraps_Tnode( thetree->child[0], iters );
    if (thetree->child[1] != NULL) {
      scale_bootstraps_Tnode( thetree->child[1], iters );
      if (thetree->child[2] != NULL) {
	scale_bootstraps_Tnode( thetree->child[2], iters );
      }
    }
  }
}



/********************************************************************* 
 FUNCTION: update_bootstraps_Tree
 DESCRIPTION: 
   Updates the bootstrap values of the destination tree, according to
   the topology of the given sample tree
 RETURNS:
 ARGS:
   Destination tree
   Sample tree
   the number of leaf nodes in the tree
 NOTES:
   This function assumes that the given trees have been created by 
   calling either neighbourjoin_buildtree or UPGMA_buildtree with 
   the bootstrap boolean arguement set to true. A fatal error results
   if this is not the case

   Another thing to note is that the method uses node numbers for
   comparisons. This means that both ethe sample and reference trees
   must have been built from the same initial list of nodes, which
   is fine for bootstrapping, because the leaf nodes are stored in
   the order in which appear in the alignment
 *********************************************************************/
void update_bootstraps_Tree( struct Tree *dest, struct Tree *sample, 
				  unsigned int numleaves) {
  unsigned int is_binary,i,j;

  is_binary = ( dest->child[1] == NULL && dest->child[2] == NULL )?1:0; 

  for (i=0; i < 3; i++) {
    for (j=0; j < 3; j++) {
      update_bootstraps_Tnode( dest->child[i], 
				    sample->child[j],
				    numleaves,
				    is_binary );
    }
  }
}



/********************************************************************* 
 FUNCTION: update_bootstraps_Tnode
 DESCRIPTION: 
   Updates the bootstrap values of the given Tnode, according to
   the topology of the given sample Tnode
 RETURNS:
 ARGS:
   Destination Tnode
   Sample Tnode
   The number of leaf nodes in each tree
   Boolean for whether the tree is binary or not
 NOTES:
   This function assumes that the given trees have been created by 
   calling either neighbourjoin_buildtree or UPGMA_buildtree with 
   the bootstrap boolean arguement set to true. A fatal error results
   if this is not the case
 *********************************************************************/
void update_bootstraps_Tnode( struct Tnode *dest, 
			      struct Tnode *sample,
			      unsigned int numleaves,
			      unsigned int is_binary) {


  compare_to_bootstrap_sample_Tnode( dest, sample, numleaves, is_binary);
  if (dest != NULL) {
    update_bootstraps_Tnode( dest->left, sample, numleaves, is_binary);
    update_bootstraps_Tnode( dest->right, sample, numleaves, is_binary);
  }
}




/********************************************************************* 
 FUNCTION: write_clustering_data_Tnode
 DESCRIPTION: 
   This routine prints a text description of the clustering details
   of the given Tnode. It was written for the old implementation,
   where leaves were named "leaf_1", "leaf_2" etc, and it would
   often be the case that each leaf would contain several sequences.
   With the current implementation, this function is not used; a 
   new-hampshire output of each cluster is printed in-situ, preserving
   sequence names from the original alignment
 RETURNS:
 ARGS: 
   File handle
   TNode *
 NOTES:
 *********************************************************************/
void write_clustering_data_Tnode( FILE *handle, struct Tnode *node) {
  unsigned int i;

  if (node != NULL) {
    if (node->left == NULL && node->right == NULL && node->clust != NULL) {
      /* we have a leaf node */
      if (node->clust->clustersize == 0 || node->clust->members == NULL) 
	fatal_util("Fatal Error: encountered a leaf node with no cluster members");
      else if (node->clust->clustersize > 1) {
	fprintf( handle, "Cluster_%d:\n", node->nodenumber);
	for( i=0; i < node->clust->clustersize; i++) {
	  if ( i == 0 )
 	    fprintf( handle, "\t" );
	  else if ( i % 4 == 0 )
	    fprintf( handle, "\n\t" );
	  else 
	    fprintf( handle, ", ");
	  fprintf( handle, "%s", node->clust->members[i]->name);
	}
	fprintf( handle, "\n\n");
      }
    }
    else {
      write_clustering_data_Tnode( handle, node->left );
      write_clustering_data_Tnode( handle, node->right );
    }
  }
}



/********************************************************************* 
 FUNCTION: write_debug_Tnode
 DESCRIPTION: 
   Writes the given Tnode to the give file handle in 'debug' format
 RETURNS:
 ARGS: 
   File handle
   TNode *
   Integer offset
 NOTES:
 *********************************************************************/

void write_debug_Tnode( FILE *handle, struct Tnode *node, unsigned int offset) {
  unsigned int i,j;

  if (node != NULL) {
    /* We need to determine whether the node is a leaf or internal;
       since in this implementation internal nodes do not have names,
       it is sufficient to check the nodes name for nullness; however,
       this precludes the possibilities of internal nodes being given
       names  in the future, hence the check for internalness is made
       on the basis of the nullness of the children.
    */
    
    if ( node->left == NULL && node->right == NULL ) {
      /* this is a leaf node, so its cluster must have members; pain if not */
      if (node->clust->clustersize == 0 || node->clust->members == NULL)
	fatal_util( "Fatal Error: encountered a leaf node with no cluster info"); 
      else {
	for (i=0; i < node->clust->clustersize; i++) {
	  for (j=0; j < offset; j++) fprintf( handle, " ");
	  /* all leaves in the cluster are printed at the same offset */
	  fprintf( handle, 
		   "%d:%s:%.5f\n", 
		   node->nodenumber, 
		   node->clust->members[i]->name,
		   node->distance); 
	}
      }
    }
    else if ( node->left != NULL && node->right != NULL ) {
      for (j=0; j < offset; j++) fprintf( handle, " ");
      fprintf(handle, "Node %d:%.5f\n", node->nodenumber, node->distance);
      write_debug_Tnode( handle, node->left, offset+2);
      write_debug_Tnode( handle, node->right, offset+2);
    }
    /* else do nothing */
    
    fflush( handle );
  }
}



/********************************************************************* 
 FUNCTION: write_debug_Tree
 DESCRIPTION: 
   prints the given Tree in a format suitable for debugging
 RETURNS:
 ARGS: 
   File handle
   Tree *
 NOTES:
 *********************************************************************/

void write_debug_Tree( FILE *handle, struct Tree *thetree) {
  if (thetree != NULL) {
      write_debug_Tnode( handle, thetree->child[0], 0 );
      write_debug_Tnode( handle, thetree->child[1], 0 );
      write_debug_Tnode( handle, thetree->child[2], 0 );

      fflush( handle );
  }

}



/********************************************************************* 
 FUNCTION: write_MUL_flattened_Tnode
 DESCRIPTION: 
   Prints the given tree as a MUL format alignment, with the sequences
   in 'tree order'
 RETURNS:
 ARGS: 
   File handle
   TNode *
 NOTES:
 *********************************************************************/

void write_MUL_flattened_Tnode( FILE *handle, struct Tnode *node) {
  unsigned int i,j;

  if (node != NULL) {
    write_MUL_flattened_Tnode( handle, node->left );
    if (node->clust != NULL) {
      for( i=0; i < node->clust->clustersize; i++ ) {
	fprintf( handle, "%-24s", node->clust->members[i]->name );
	for (j=0; j < node->clust->members[i]->length; j++)
	  fprintf( handle, "%c", node->clust->members[i]->seq[j] );
	fprintf( handle, "\n");
      }
    }
    write_MUL_flattened_Tnode( handle, node->right );
  }

}



/********************************************************************* 
 FUNCTION: write_MUL_flattened_Tree
 DESCRIPTION: 
   Prints the given tree as a MUL format sequence alignment, with
   the sequences in 'tree order'
 RETURNS:
 ARGS: 
   File handle
   Tree *
 NOTES:
 *********************************************************************/

void write_MUL_flattened_Tree( FILE *handle, struct Tree *tr) {
  if (tr != NULL) {
    write_MUL_flattened_Tnode( handle, tr->child[0] );
    write_MUL_flattened_Tnode( handle, tr->child[1] );
    write_MUL_flattened_Tnode( handle, tr->child[2] );
  }
  fflush( handle );

}



/********************************************************************* 
 FUNCTION: write_newhampshire_Tnode
 DESCRIPTION: 
   prints the given Tree in 'New Hampshire' text format to the given 
   file handle
 RETURNS:
 ARGS: 
   File handle
   TNode *
   Whether or not to show bootstrap values
 NOTES:
 *********************************************************************/

void write_newhampshire_Tnode( FILE *handle, struct Tnode *node, 
			       unsigned int show_bootstraps  ) {
  if (node != NULL) {
    /* We need to determine whether the node is a leaf or internal;
       since in this implementation internal nodes do not have names,
       it is sufficient to check the nodes name for nullness; however,
       this precludes the possibilities of internal nodes being given
       names  in the future, hence the check for internalness is made
       on the basis of the nullness of the children.
    */
    
    if ( node->left == NULL && node->right == NULL) {
      /* this is a leaf node, so its cluster must have members; pain if not */
      if (node->clust->clustersize == 0 || node->clust->members == NULL)
	fatal_util( "Fatal Error: encountered a leaf node with no cluster info"); 
      else if (node->clust->clustersize == 1) 
	fprintf( handle, "%s:%.5f", node->clust->members[0]->name, node->distance );
      else {
	/* if there is more than one sequence belonging to the cluster, then this piece
	   of code will generate som internal nodes in the output tree with no bootstrap
	   values. Such is life... */
        unsigned int i;

	for (i=0; i < node->clust->clustersize - 1; i++) {
	 fprintf( handle, "(\n%s:%.5f,\n", node->clust->members[i]->name, 0.0 ); 
	}
	fprintf( handle, "%s:%.5f)\n", node->clust->members[i]->name, 0.0);
	for (i=0; i < node->clust->clustersize - 2; i++) {
	 fprintf( handle, ":%.5f)\n", 0.0 ); 
	}

	fprintf( handle, ":%.5f", node->distance);
	/* fprintf( handle, "Cluster_%d:%.5f", node->nodenumber, node->distance ); */
      }
    }
    else if ( node->left != NULL && node->right != NULL ) {
      fprintf( handle, "(\n");
      write_newhampshire_Tnode( handle, node->left, show_bootstraps );
      fprintf( handle, ",\n" );
      write_newhampshire_Tnode( handle, node->right, show_bootstraps );
      if (show_bootstraps) {
	fprintf( handle, ")\n%d:%.5f", node->bootstrap, node->distance );
      }
      else {
	fprintf( handle, ")\n:%.5f", node->distance);
      }
    }
    /* else do nothing */
  }
  fflush( handle );
}



/********************************************************************* 
 FUNCTION: write_newhampshire_Tree
 DESCRIPTION: 
   prints the given Tnode in 'New Hampshire' text format to the given 
   file handle
 RETURNS:
 ARGS: 
   File handle
   TNode *
 NOTES:
 *********************************************************************/

void write_newhampshire_Tree( FILE *handle, struct Tree *thetree,
			      unsigned int show_bootstraps) {

  /* write_newhampshire_Tnode always places parenttheses around the
     sub-tree if there is more than one sub-node, but this is not
     appropriate if there there are no other sub-trees to draw at
     this level. Hence there is a special case for a single
     sub-tree (as returned by the UPGMA method) */

  if (thetree != NULL) {
    if (thetree->child[0] != NULL) {
      if (thetree->child[1] == NULL) {
	/* draw rooted tree */

	if (thetree->child[0]->left != NULL && thetree->child[0]->right != NULL) {
	  fprintf( handle, "(\n");
	  write_newhampshire_Tnode( handle, thetree->child[0]->left, show_bootstraps );
	  fprintf( handle, ",\n");
	  write_newhampshire_Tnode( handle, thetree->child[0]->right, show_bootstraps );
	  fprintf( handle, ");\n");
	}
	else {
	  /* this is a leaf node, and leaf nodes may contain a single sequence
	     or a cluser of sequences. If this leaf contains a single sequence,
	     then we have a tree of one sequence, in which case we print an
	     error, because trees of one sequence do not make sense */
	  
	  if (thetree->child[0]->clust->clustersize == 1) 
	    fatal_util( "Cannot build a tree with a single sequence %s",
		     thetree->child[0]->clust->members[0]->name);
	  else {
	    unsigned int i;
	    
	    for (i=0; i < thetree->child[0]->clust->clustersize - 1; i++)
	      fprintf( handle, "(\n%s:%.5f,\n", thetree->child[0]->clust->members[i]->name, 0.0 ); 

	    fprintf( handle, "%s:%.5f", thetree->child[0]->clust->members[i]->name, 0.0);

	    for (i=0; i < thetree->child[0]->clust->clustersize - 2; i++) 
	      fprintf( handle, ")\n:%.5f)\n", 0.0 ); 

	    fprintf( handle, ");\n"); 
	  }
	}
      }
      else {
	fprintf( handle, "(\n");
	write_newhampshire_Tnode( handle, thetree->child[0], show_bootstraps );
	fprintf( handle, ",\n");
	write_newhampshire_Tnode( handle, thetree->child[1], show_bootstraps );
	if (thetree->child[2] != NULL) {
	  fprintf( handle, ",\n");
	  write_newhampshire_Tnode( handle, thetree->child[2], show_bootstraps );
	}
	fprintf( handle, ");\n");
      }
    }
  }
  fflush( handle );
}


/*  Last edited: Feb  1 18:18 2002 (klh) */
/**********************************************************************
 ** FILE: cluster.c
 ** NOTES:
 **  A DistanceMatrix should always be part of a Cluster 
 **  It makes no sense to have a set of pairwise distances without the 
 **  associated sequences (even if we just store their names)
 **********************************************************************/



/*********************************************************************
  FUNCTION: alignment_to_ClusterGroup
  DESCRIPTION: 
    This function returns a ClusterGroup, given an Alignment. 
    if the secind arg is true, In doing indentical sequences in the 
    alignment are merged. If bootstrapping is required, the consensus
    alignment can be extracted from the ClusterGroup using 
    get_consensus_from_ClusterGroup
  RETURNS: struct ClusterGroup
  ARGS: 
    1. A source Alignment pointer
    2. A boolean specifying whether duplicate sequences should be
       merged.
  NOTES: 
*********************************************************************/
struct ClusterGroup *alignment_to_ClusterGroup( struct Alignment *aln,
						unsigned int remove_duplicates) {
  unsigned int i, j, numclusters;
  struct ClusterGroup *group;
  struct Cluster **newclusts;

  group = empty_ClusterGroup();

  /* Need to create a cluster for every sequence in the given cluster */

  newclusts = (struct Cluster **) malloc_util( aln->numseqs * sizeof( struct Cluster *) );


  for( i=0; i < aln->numseqs; i++) {
    newclusts[i] = single_Sequence_Cluster( clone_Sequence(aln->seqs[i]));
  }
  numclusters = aln->numseqs;  

  if (remove_duplicates) {
    for( i=0; i < aln->numseqs; i++) {
      if (newclusts[i] == NULL) continue;
      for( j=i+1; j < aln->numseqs; j++) {
	if (newclusts[j] == NULL) continue;
	
	if (strncmp( aln->seqs[i]->seq, aln->seqs[j]->seq, aln->length) == 0) {
	  /* these two clusters contain the same sequence, so we can merge them */
	  newclusts[j] = merge_Cluster( newclusts[i], newclusts[j] );
	  numclusters--;
	}
      }
    }
  }


  /* newclusts will now be a sparse array containing clusters of
     identical sequences */

  group->numclusters = numclusters;
  group->clusters = (struct Cluster **) malloc_util( numclusters * sizeof( struct Cluster *) );
  for(i=0, j=0; i < aln->numseqs; i++) {
    if (newclusts[i] != NULL) {
      group->clusters[j++] = newclusts[i];
    }
  }

  newclusts = free_util( newclusts );

  return group;
}


/********************************************************************* 
 FUNCTION: clone_Cluster
 DESCRIPTION: 
   This function makes a complete copy of the given Cluster
   and returns it
 RETURNS: struct Cluster *
 ARGS: 
   struct Cluster *
 NOTES: 
*********************************************************************/
struct Cluster *clone_Cluster( struct Cluster *source) {
  unsigned int i;
  struct Cluster *dest = NULL; 
  
  if (source != NULL) {
    dest = empty_Cluster();
    dest->clustersize = source->clustersize;
    dest->members = (struct Sequence **) malloc_util( dest->clustersize
						      * sizeof( struct Sequence * ));
    for( i=0; i < source->clustersize; i++) {
      dest->members[i] = clone_Sequence( source->members[i] );
    }
    dest->consensus = clone_Sequence( source->consensus );
  }
  
  return dest;
}



/*********************************************************************
  FUNCTION: consensus_aln_from_ClusterGroup
  DESCRIPTION: 
    This function creates an alignment by taking the consensus 
    sequences from each Cluster in the given ClusterGroup
  RETURNS: struct ClusterGroup
  ARGS: 
    1. A source Alignment pointer
    2. A boolean specifying whether duplicate sequences should be
       merged.
  NOTES: 
*********************************************************************/
struct Alignment *consensus_aln_from_ClusterGroup( struct ClusterGroup *grp) {
  unsigned int i;
  struct Alignment *cons;

  cons = (struct Alignment *) malloc_util( sizeof( struct Alignment ));
  cons->numseqs = grp->numclusters;
  cons->seqs = (struct Sequence **) malloc_util( grp->numclusters * sizeof( struct Sequence *));

  for (i=0; i < grp->numclusters; i++) {
    cons->seqs[i] = clone_Sequence( grp->clusters[i]->consensus );
  }
  /* All consensus seqs will be aligned, so we can pick any one to get 
     the alignment width */
  cons->length = cons->seqs[0]->length;

  return cons;
}





/********************************************************************* 
  FUNCTION: empty_Cluster
  DESCRIPTION: 
    This function handles the simple task of allocating the space
    for a new Cluster.
  RETURNS: struct Cluster *
  ARGS: 
  NOTES: 
*********************************************************************/
struct Cluster *empty_Cluster( void ) {
  struct Cluster *newclust;
  
  newclust = (struct Cluster *) malloc_util( sizeof( struct Cluster ));
  newclust->clustersize = 0;
  newclust->members = NULL;
  newclust->consensus = NULL;
  newclust->matrix = NULL; 
  
  return newclust;
  
}



/********************************************************************* 
  FUNCTION: empty_ClusterGroup
  DESCRIPTION: 
    This function handles the simple task of allocating the space
    for a new ClusterGroup
  RETURNS: struct Cluster *
  ARGS: 
  NOTES: 
*********************************************************************/
struct ClusterGroup *empty_ClusterGroup( void ) {
  struct ClusterGroup *group;
  
  group = (struct ClusterGroup *) malloc_util( sizeof(struct ClusterGroup));
  group->numclusters = 0;
  group->clusters = NULL;
  group->matrix = NULL;
  
  return group;
}




/********************************************************************* 
  FUNCTION: free_Cluster
  DESCRIPTION: 
    This function releases the memory used by this Cluster and all of its
    members
  RETURNS: A null pointer
  ARGS: 
    struct Cluster *
  NOTES: 
    In the majority of cases, all sequences in the Cluster come from
    an alignment, and if this alignment is subsequently needed (e.g.
    for bootstrapping) then then the seqs should not be freed. To 
    prevent this, the members field should be set to null by the caller
    to prevent freeing og the alignment
*********************************************************************/
void *free_Cluster( struct Cluster *given ) {
  unsigned int i;
  
  if (given != NULL) {
    if (given->members != NULL) {
      for( i=0; i < given->clustersize; i++) {
	given->members[i] =  free_Sequence( given->members[i] );
      }
      given->members = free_util( given->members );
    }
    given->matrix = free_DistanceMatrix( given->matrix );
    given->consensus = free_Sequence( given->consensus );
    given = free_util( given );
  }
  return given;
}




/********************************************************************* 
  FUNCTION: free_ClusterGroup
  DESCRIPTION: 
    This function releases the memory used by this Cluster and all of its
    members
  RETURNS: A null pointer
  ARGS: 
    struct Cluster *
  NOTES: 
*********************************************************************/
void *free_ClusterGroup( struct ClusterGroup *given ) {
  unsigned int i;
  
  if (given != NULL) {
    if (given->clusters != NULL) {
      for( i=0; i < given->numclusters; i++ ) {
	given->clusters[i] = free_Cluster( given->clusters[i] );
      }
      given->clusters = free_util( given->clusters );
    }
    given->matrix = free_DistanceMatrix( given->matrix );
    given = free_util( given );
  }
  
  return given;
}




/********************************************************************* 
  FUNCTION: merge_Cluster
  DESCRIPTION: 
    Adds the sequences in second arg to first arg, freeing the second
    arg, returning the result of this freeing (hopefully NULL);
  RETURNS: The result of freeing the second cluster (NULL if all is well)
  ARGS: 
    Destination Cluster *, 
    Source Cluster *
  NOTES:
*********************************************************************/
void *merge_Cluster( struct Cluster *dest, struct Cluster *source) {
  unsigned int i;
  
  /* take the sequences in source and add them onto the destination list */
  
  dest->members = (struct Sequence **) 
    realloc_util( dest->members, 
		  (dest->clustersize + source->clustersize) * sizeof(struct Sequence *));
  for ( i=0; i < source->clustersize; i++) {
    dest->members[ dest->clustersize++ ] = source->members[i];
    source->members[i] = NULL;
  }
  
  /* Need to update the consensus sequence with respect to the merge.
     At this stage, I am only merging identical clusters, so the
     consensus sequence is already correct in the destination 
  */

  source = free_Cluster( source );
  
  return source;
}



/********************************************************************* 
 FUNCTION: single_Sequence_Cluster
 DESCRIPTION: 
   This function handles the simple task of allocating the space
   for a new Cluster with the single given Sequence.
 RETURNS: struct Cluster *
 ARGS: 
   A pointer to a Sequence, or NULL for an empty Cluster
 NOTES: 
*********************************************************************/

struct Cluster *single_Sequence_Cluster( struct Sequence *seq) {
  struct Cluster *newclust;

  newclust = empty_Cluster();

  if (seq != NULL) {
    newclust->clustersize = 1;
    newclust->members = (struct Sequence **) malloc_util( sizeof( struct Sequence *) );
    newclust->members[0] = seq;
    newclust->consensus = clone_Sequence( seq );
    /* Note how the consensus sequence will end up with the same name as
       the first representative of trhe cluster. This is fine, because
       when the consensus is used, it is assumed that it is not a real
       sequence so its name is meaningless 
    */
  }

  return newclust;
}




/********************************************************************* 
 FUNCTION: single_Cluster_ClusterGroup
 DESCRIPTION: 
   This function takes the given cluster and very simlpy makes a 
   single-Cluster ClusterGroup from it
 RETURNS: struct ClusterGroup *
 ARGS: 
   A pointer to a Cluster
 NOTES: 
*********************************************************************/

struct ClusterGroup *single_Cluster_ClusterGroup( struct Cluster *clust ) {
  struct ClusterGroup *group;

  group = empty_ClusterGroup();

  if (clust != NULL) {
    group->numclusters = 1;
    group->clusters = (struct Cluster **) malloc_util( sizeof(struct Cluster *));
    group->clusters[0] = clust;
  }
  
  return group;
}








/*  Last edited: Mar 18 11:52 2002 (klh) */
/**********************************************************************
 ** FILE: buildtree.c
 ** NOTES:
 **  Contains functions for building trees from distance matrices
 **  (and vice versa)
 **********************************************************************/


/**********************************************************************
 FUNCTION: export_distances_buildtree
 DESCRIPTION: 
   Returns the distance matrix induced from the given tree, i.e. by summing
   the branch paths between two nodes to obtain their distance
 ARGS: 
   A Tree
   A DistanceMatrix
 RETURNS: 
 NOTES: 
   This function does not create the memory for the distance matrix, it
   merely fills in the given matrix
 **********************************************************************/
void export_distances_buildtree( struct Tree *thetree, 
				 struct DistanceMatrix *mat) {

  struct Tnode *root;
  unsigned int *considered;

  /* first, create the imaginary node that sits between the three nodes
     of the unrooted tree */

  root = new_interior_Tnode( thetree->numnodes );

  /* now set up the connections */

  root->parent = thetree->child[0];
  root->left = thetree->child[1];
  root->right = thetree->child[2];
  root->left->parent = root;
  root->right->parent = root;
  root->parent->parent = root; /* This one seems odd but we need it */
  root->distance = root->parent->distance;

  /* A boolean array to store which nodes have been considered */

  considered = (unsigned int *) 
    malloc_util( (thetree->numnodes + 1) * sizeof( unsigned int ) );

  leaf_find_buildtree( root->parent, mat, considered, thetree->numnodes + 1);
  leaf_find_buildtree( root->left, mat, considered, thetree->numnodes + 1);
  leaf_find_buildtree( root->right, mat, considered, thetree->numnodes + 1);


  considered = free_util( considered );

  /* When freeing the temporary root we have created, we must nullify
     the parent and children. We don't want to cause a cascading delete
     of the tree */

  root->parent = NULL;
  root->left = NULL;
  root->right = NULL;
  root = free_Tnode( root );

}




/**********************************************************************
 FUNCTION: find_path_buildtree
 DESCRIPTION: 
   Given a leaf node, recursively calculates the branch-length distance
   from the node to all other leaves in the tree, and places it in the
   appropriate part of the DistanceMatrix
 ARGS: 
   unsigned int (the number of the node from which we are finding all distances)
   Tnode (the current node under consideration)
   DistanceMatrix
   Boolean array (to store which nodes have already been considered)
   unsigned int (the size of this boolean array)
 RETURNS: 
 NOTES: 
 **********************************************************************/
void find_path_buildtree( unsigned int home,
			  struct Tnode *node,
			  struct DistanceMatrix *mat,
			  double dist,
			  unsigned int *considered) {

  if (node->left == NULL && node->right == NULL) {
    if (home < node->nodenumber) {
      mat->data[node->nodenumber][home] = dist;
    }
    else {
      mat->data[home][node->nodenumber] = dist;
    }
  }
  else {
    /* we have an internal node */
    considered[node->nodenumber] = 1;

    if (node->left != NULL) {
      if (! considered[node->left->nodenumber]) {
	find_path_buildtree( home, 
			     node->left, 
			     mat, 
			     dist + node->left->distance,
			     considered );
      }
    }
    if (node->right != NULL) {
      if (! considered[node->right->nodenumber]) {
	find_path_buildtree( home, 
			     node->right, 
			     mat, 
			     dist + node->right->distance,
			     considered );
      }
    }
    if (node->parent != NULL) {
      if (! considered[node->parent->nodenumber]) {
	find_path_buildtree( home, 
			     node->parent, 
			     mat, 
			     dist + node->distance,
			     considered );
      }
    }
  }

}


/**********************************************************************
 FUNCTION: leaf_find_buildtree
 DESCRIPTION: 
   Finds the leaf nodes descended from the given interior node, and
   calculates the distance from these nodes to all other nodes, by
   way of another function call
 ARGS: 
   Tnode
   DistanceMatrix
   Boolean array (to store which nodes have already been considered)
 RETURNS: 
 NOTES: 
 **********************************************************************/
void leaf_find_buildtree( struct Tnode *node, 
			  struct DistanceMatrix *mat, 
			  unsigned int *considered,
			  unsigned int size) {

  unsigned int i;

  if (node->left != NULL || node->right != NULL) {
    if (node->left != NULL) 
      leaf_find_buildtree( node->left, mat, considered, size );
    if (node->right != NULL) 
      leaf_find_buildtree( node->right, mat, considered, size );
  }
  else {

    mat->data[node->nodenumber][node->nodenumber] = 0.0;
    for(i=0; i < size; i++) considered[i] = 0;
    considered[node->nodenumber] = 1;
    find_path_buildtree( node->nodenumber, 
			 node->parent, 
			 mat, 
			 node->distance, 
			 considered);
  }

} 




/**********************************************************************
 FUNCTION: neighbour_joining_buildtree
 DESCRIPTION: 
   Returns a phylogenetic tree of the sequences in the 
   given alignment, using Saitou and Nei's neighbour-joining 
   algorithm
 ARGS: 
    A ClusterGroup pointer (cluster.h)
    Boolean, for whether to calc information needed for later bootstrapping
 RETURNS:
    A Tree (trees.h)
 NOTES: The function allocates all the memory necessary for the tree.
   The caller should call free_tree (tree.h) to free this memory when
   the tree is no longer needed
 **********************************************************************/
struct Tree *neighbour_joining_buildtree( struct ClusterGroup *group,
					  unsigned int bootstrap) { 
  unsigned int numseqs, i, j;       /*** The current pair of nodes   ***/
  unsigned int k, m, nodecount;     /*** loop counters               ***/
  unsigned int row, column;         /*** matrix indices              ***/
  unsigned int mini = 0, minj = 0;  /*** neighbouring nodes          ***/
  unsigned int nextfreenode;        /*** incremental labels to nodes ***/ 
  unsigned int leftovers[3];        /*** three remaining nodes       ***/
  struct DistanceMatrix *mat;
  struct Tree *theTree;
  struct Tnode **nodes;             /*** starts off holding leaves    ***/
  struct Tnode *newnode;
  double fnumseqs;                  /*** divisor for sums            ***/
  double dmj, dmi, ri, minsofar, dist, dij, dist_i, dist_j, dist_k;
  double *r;                        /*** stores the r values         ***/


  /* METHOD ***********************************
     I have implemented the neightbour-joining algorithm as descibed
     by Durbin, Eddy, Krogh and Mitchison (1998 pp 171-172); see the
     text for reasoning behind variable names. I have inserted several
     changes from the basic algorithm for reasons of efficiency

     To keep track of dead nodes, we use the 'nodes' pointer array.
     Whenever we amalgamate nodes, we overwrite one of the entries
     with the new node and set the other to NULL. The NULL entries
     in the nodes array thus mark the dead nodes; this is important
     when updating the distance matrix because we must not add terms
     from nodes that are gone.

     Note the complex indexing of the distance matrix, since it is
     assumed that it is triangularised at the top-right. Two Ways
     around this: 
     1. Abstract the indexing into a distancemat method
     2. Assume a square symmetrical matrix
     Method one may have performance impact; method two is more viable
     but need to be much more cafeful with the updating of the matrix
     at the end of each iteration

     I have also used a modified version of Bill Bruno's idea to attempt
     to eliminate negative branch lengths in generated trees.

  ********************************************/



  /******* intialisation ********************/

  nodes = (struct Tnode **) malloc_util( group->numclusters * sizeof(struct Tnode *));
  for( i=0; i < group->numclusters; i++) { 
    nodes[i] = new_leaf_Tnode( i, clone_Cluster( group->clusters[i]) );
  }

  mat = group->matrix;
  numseqs = nextfreenode =  mat->size;
  
  theTree = empty_Tree();

  fnumseqs = (double) numseqs;
  

  if (numseqs > 2) {

    r = (double *) malloc_util( numseqs * sizeof( double ) );
    /* Calculate r[i] for all i */
    for( i=0; i < numseqs; i++ ) {
      ri = 0.0;
      for (k=0; k < numseqs; k++) {
	if (k > i) ri += mat->data[k][i];
	else ri += mat->data[i][k];
      }
      r[i] = ri / (fnumseqs - 2.0);
    }
    
    
    /******* main loop ************************/
    
    for (nodecount=0; nodecount < numseqs-3; nodecount++) {

      /* do the intialisation necessary for each iteration here */

      minsofar = DBL_MAX;  /* from float.h */

      /******* for each pair of matrix entries *********************/

      for( i=0; i < numseqs; i++ ) {
	if (nodes[i] == NULL) continue;
	for( j=0; j < i; j++ ) {
	  if (nodes[j] == NULL) continue;
	  
	  dist = mat->data[i][j] - (r[i] + r[j]);
	  if (dist < minsofar) {
	    minsofar = dist;
	    mini = i;
	    minj = j;
	  }
	}
      }
      
      /* printf("i = %d, j = %d\n", mini, minj); */
      
      /* we have the neighbouring i, j; lets calc distances and make the new node */

      dij = mat->data[mini][minj];
      dist_i = (dij + r[mini] - r[minj]) * 0.5;
      dist_j = dij - dist_i;

      /* Adjustment to allow for negative branch lengths */
      if (dist_i < 0.0) {
	dist_i = 0.0;
	dist_j = dij;
	if (dist_j < 0.0)
	  dist_j = 0.0;
      }
      else if (dist_j < 0.0) {
	dist_j = 0.0;
	dist_i = dij;
	if (dist_i < 0.0)
	  dist_i = 0.0;
      }

      nodes[mini]->distance = dist_i;
      /* printf("DistanceI = %f\n", nodes[mini]->distance); */
      
      nodes[minj]->distance = dist_j;
      /* printf("DistanceJ = %f\n", nodes[minj]->distance); */
      
      newnode = new_interior_Tnode( nextfreenode++ );
      newnode->left = nodes[mini];
      newnode->right = nodes[minj];
      nodes[mini]->parent = newnode;
      nodes[minj]->parent = newnode;
      nodes[mini] = newnode;
      nodes[minj] = NULL;
      if (bootstrap) { 
	/* we need to create and load the 'bit' field of child ids */
	
	newnode->child_ids = (unsigned int *) 
	  malloc_util( group->numclusters * sizeof( unsigned int ) );
	for (m=0; m < group->numclusters; m++) {
	  if ( (newnode->left->child_ids != NULL && 
		    (newnode->left->child_ids[m] || newnode->left->child_ids[m])) ||
	       (newnode->right->child_ids != NULL &&
 		    (newnode->right->child_ids[m] || newnode->right->child_ids[m])) ||
	       ((newnode->left->nodenumber == m || 
		 newnode->right->nodenumber == m ||
		 newnode->nodenumber == m))) {
	    newnode->child_ids[m] = 1;
	  }
	  else {
	    newnode->child_ids[m] = 0;
	  }
	}
      }
      
      /* now update the distance matrix; This needs hackery to make sure that the
	 indexing is correct */

      r[mini] = 0.0;  /* This is the only r[i] that requires wholesale changes */
      for( m=0; m < numseqs; m++ ) {
	if (nodes[m] == NULL) continue;
	
	if (m != mini) {
	  
	  if (m > minj) dmj = mat->data[m][minj];
	  else dmj = mat->data[minj][m];
	  
	  if (m > mini) {
	    row = m;
	    column = mini;
	  }
	  else {
	    row = mini;
	    column = m;
	  }
	  dmi = mat->data[row][column];
	  
	  /* we can actually adjust r[m] here, by using the form:
	     rm = ((rm * numseqs) - dmi - dmj + dmk) / (numseqs-1)
	  */
	  
	  /* Note: in Bill Bruno's method for negative branch elimination, then if either
	     dist_i is positive and dist_j is 0, or dist_i is zero and dist_j is positive
	     (after adjustment) then the matrix entry is formed from the distance to the
	     node in question (m) to the node with the zero branch length (whichever it was).
	     I think my code already has the same effect; this is certainly true if dij is
	     equal to dist_i + dist_j, which it should have been fixed to
	  */
	  
	  mat->data[row][column] = (dmi + dmj - dij) * 0.5;
	  r[m] = ((r[m] * (fnumseqs - 2.0)) - dmi - dmj + mat->data[row][column]) / (fnumseqs - 3.0); 
	  r[mini] += mat->data[row][column];
	}
      }
      
      fnumseqs -= 1.0;
      r[mini] /= fnumseqs - 2.0;
      
    }
    /******* end of main loop ******************/
    
    
    /* Now there are just 3 nodes left. Need to locate those three nodes */
    /* The following looks dangerous, because we only have room for three nodes in the tree;
       However, all nodes except three should be NULL, or else something is seriously wrong */
    
    for(k=0, m=0; k < numseqs; k++) {
      if (nodes[k] != NULL) { 
	theTree->child[m] = nodes[k];
	leftovers[m++] = k;
	nodes[k] = NULL; /* so that the nodes are not released when input is */
      }
    }
    

    /* Now to get rid of those negative branch lengths (see you in ~70 lines....) */
    
    dist_i = theTree->child[0]->distance =
      (mat->data[leftovers[1]][leftovers[0]] +
       mat->data[leftovers[2]][leftovers[0]] -
       mat->data[leftovers[2]][leftovers[1]]) * 0.5;
    dist_j = theTree->child[1]->distance = mat->data[leftovers[1]][leftovers[0]] - theTree->child[0]->distance;
    dist_k = theTree->child[2]->distance = mat->data[leftovers[2]][leftovers[0]] - theTree->child[0]->distance;


    if (dist_i < 0.0) {
      dist_i = 0.0;
      dist_j = mat->data[leftovers[1]][leftovers[0]];
      dist_k = mat->data[leftovers[2]][leftovers[0]];
      if (dist_j < 0.0) {
	dist_j = 0.0;
	dist_k = ( mat->data[leftovers[2]][leftovers[0]] +
		   mat->data[leftovers[2]][leftovers[1]] ) * 0.5;
	if (dist_k < 0.0) 
	  dist_k = 0.0;
      }
      else if (dist_k < 0.0) {
	dist_k = 0.0;
	dist_j = 
	  (mat->data[leftovers[1]][leftovers[0]] +
	   mat->data[leftovers[2]][leftovers[1]] ) * 0.5;
	if (dist_j < 0.0)
	  dist_j = 0.0;
      }
    }
    else if (dist_j < 0.0) {
      dist_j = 0.0;
      dist_i = mat->data[leftovers[1]][leftovers[0]];
      dist_k = mat->data[leftovers[2]][leftovers[1]];
      if (dist_i < 0.0) {
	dist_i = 0.0;
	dist_k = ( mat->data[leftovers[2]][leftovers[0]] +
		   mat->data[leftovers[2]][leftovers[1]] ) * 0.5;
	if (dist_k < 0.0) 
	  dist_k = 0.0;
      }
      else if (dist_k < 0.0) {
	dist_k = 0.0;
	dist_i = 
	  (mat->data[leftovers[1]][leftovers[0]] +
	   mat->data[leftovers[2]][leftovers[0]] ) * 0.5;
	if (dist_i < 0.0)
	  dist_i = 0.0;
      }
    }
    else if (dist_k < 0.0) {
      dist_k = 0.0;
      dist_i = mat->data[leftovers[2]][leftovers[0]];
      dist_j = mat->data[leftovers[2]][leftovers[1]];
      if (dist_i < 0.0) {
	dist_i = 0.0;
	dist_j = ( mat->data[leftovers[1]][leftovers[0]] +
		   mat->data[leftovers[2]][leftovers[1]] ) * 0.5;
	if (dist_j < 0.0) 
	  dist_j = 0.0;
      }
      else if (dist_j < 0.0) {
	dist_j = 0.0;
	dist_i = 
	  (mat->data[leftovers[1]][leftovers[0]] +
	   mat->data[leftovers[2]][leftovers[0]] ) * 0.5;
	if (dist_i < 0.0)
	  dist_i = 0.0;
      }
    }



    theTree->child[0]->distance = dist_i; 
    theTree->child[1]->distance = dist_j;
    theTree->child[2]->distance = dist_k;


    r = free_util( r );
  }
  else {
    /* deal with the trivial case of less than three leaves */

    for( i=0; i < numseqs; i++ ) {
      theTree->child[i] = nodes[i];
      nodes[i] = NULL;
    }

    if (numseqs == 2) {
      theTree->child[0]->distance = mat->data[1][0] * 0.5;
      theTree->child[1]->distance = mat->data[1][0] * 0.5;
    }
  }
 
  /* nextfreenode records how many calls to new_Tnode_tree we made; therefore, at this stage, it
     holds the total number of nodes (including leaves) in the tree */
  theTree->numnodes = nextfreenode;
  
  nodes = free_util( nodes );
  
  /* The caller of the function should free the tree when finished with it */
  
  return theTree;
}



/**********************************************************************
 FUNCTION: UPGMA_buildtree
 DESCRIPTION: 
   Returns a phylogenetic tree of the sequences in the 
   given alignment, using the Unweighted Pair-Group method based on
   Arithmentic Averages (UPGMA)
 ARGS: 
   A ClusterGroup pointer (cluster.h)
   A boolean,for whether to record information needed later for bootstrapping
 RETURNS: 
   struct Tree (trees.h)
 NOTES: The function allocates all the memory necessary for the tree.
   The caller should call free_Tnode (tree.h) to free this memory when
   the tree is no longer needed

   This algorithm produces a rooted tree; hence the returned Tree
   will have the root in child[0]; child[1] and child[2] (used to
   represent the trichotomy of unrooted trees) will be NULL 
 **********************************************************************/
struct Tree *UPGMA_buildtree(struct ClusterGroup *group,
			     unsigned int bootstrap) { 

  unsigned int numseqs, i, j;       /*** The current pair of nodes   ***/
  unsigned int m, nodecount;        /*** loop counters               ***/
  unsigned int row, column;         /*** matrix indices              ***/
  unsigned int mini = 0, minj = 0;  /*** neighbouring nodes          ***/
  unsigned int nextfreenode;        /*** incremental labels to nodes ***/ 
  unsigned int *subtreesizes;
  double *heights, minsofar, newnodeheight, dmi, dmj;
  
  struct DistanceMatrix *mat = NULL;
  struct Tree *theTree = NULL;
  struct Tnode **nodes = NULL;           /*** starts off holding leaves    ***/
  struct Tnode *newnode = NULL;

  /* METHOD ***********************************
     I have implemented the UPGMA algorithm as presented 
     by Durbin, Eddy, Krogh and Mitchison (1998 pp 171-172).

     To keep track of merged nodes, we use the 'nodes' pointer array.
     Whenever we amalgamate nodes, we overwrite one of the entries
     with the new node and set the other to NULL. The NULL entries
     in the nodes array thus mark the dead nodes; this is important
     when updating the distance matrix because we must not add terms
     from nodes that are gone.

     Note the complex indexing of the distance matrix, since it is
     assumed that it is triangularised at the top-right. Two Ways
     around this: 
     1. Abstract the indexing into a distancemat method
     2. Assume a square symmetrical matrix
     Method one may have performance impact; method two is more viable
     but need to be much more cafeful with the updating of the matrix
     at the end of each iteration

  ********************************************/



  /******* intialisation ********************/

  nodes = (struct Tnode **) malloc_util( group->numclusters * sizeof(struct Tnode *));
  heights = (double *) malloc_util( group->numclusters * sizeof(double));
  subtreesizes = (unsigned int *) malloc_util( group->numclusters * sizeof(unsigned int));
  for( i=0; i < group->numclusters; i++) { 
    nodes[i] = new_leaf_Tnode( i, group->clusters[i] );
    /* clusters[i] in group has conceptually moved into the Tnode, so... */
    group->clusters[i] = NULL;
    
    heights[i] = 0.0;
    subtreesizes[i] = 1;
  }

  mat = group->matrix;
  numseqs = nextfreenode =  mat->size;

  theTree = empty_Tree();
   
  /******* main loop ************************/

  if (numseqs == 1)
    newnode = nodes[0];
  else {
    for (nodecount=0; nodecount < numseqs-1; nodecount++) {

      /* do the intialisation necessary for each iteration here */
      
      minsofar = DBL_MAX;  /* from float.h */
      
      /******* for each pair of matrix entries *********************/
      
      for( i=0; i < numseqs; i++ ) {
	if (nodes[i] == NULL) continue;
	for( j=0; j < i; j++ ) {
	  if (nodes[j] == NULL) continue;
	  
	  if (mat->data[i][j] < minsofar) {
	    minsofar = mat->data[i][j];
	    mini = i;
	    minj = j;
	  }
	}
      }
      
      /* printf("i = %d, j = %d, ri = %f, rj = %f\n", mini, minj, ri, rj); */
      
      /* we have the neighbouring i, j; lets calc distances and make the new node */
      
      newnodeheight = mat->data[mini][minj] * 0.5;
      nodes[mini]->distance = newnodeheight - heights[mini];
      nodes[minj]->distance = newnodeheight - heights[minj];
      
      newnode = new_interior_Tnode( nextfreenode++ );
      newnode->left = nodes[mini];
      newnode->right = nodes[minj];
      nodes[mini]->parent = newnode;
      nodes[minj]->parent = newnode;
      if (bootstrap) { 
	/* we need to create and load the 'bit' field of child ids */
	
	newnode->child_ids = (unsigned int *) 
	  malloc_util( group->numclusters * sizeof( unsigned int ) );
	for (i=0; i < group->numclusters; i++) {
	  if ( (newnode->left->child_ids != NULL && 
		(newnode->left->child_ids[i] || newnode->left->child_ids[i])) ||
	       (newnode->right->child_ids != NULL && 
		(newnode->left->nodenumber == i || 
		 newnode->right->nodenumber == i ||
		 newnode->nodenumber == i))) {
	    newnode->child_ids[i] = 1; 
	  }
	else {
	  newnode->child_ids[i] = 0;
	}
	}
      }
      
      /* now update the distance matrix; This needs hackery to make sure that the
	 indexing is correct */
      
      for( m=0; m < numseqs; m++ ) {
	if (nodes[m] == NULL) continue;
	
	/*                  dmini,m*|Ci| + dminj,m*|Cj|
			    dmini,m =     ---------------------------   
			    |Ci| + |Cj|
	*/
	
	if (m > minj) dmj = mat->data[m][minj];
	else dmj = mat->data[minj][m];
	
	if (m > mini) {
	  row = m;
	  column = mini;
	}
	else {
	  row = mini;
	  column = m;
	}
	dmi = mat->data[row][column];
	
	mat->data[row][column] = 
	  (( dmi * subtreesizes[mini])+( dmj * subtreesizes[minj])) / 
	  (subtreesizes[mini] + subtreesizes[minj]);    
      }
      
      heights[mini] = newnodeheight;
      subtreesizes[mini] = subtreesizes[mini] + subtreesizes[minj] + 1; 
      
      nodes[mini] = newnode;
      nodes[minj] = NULL;
      
    }
  }
  /******* end of main loop ******************/
  
    
  theTree->child[0] = newnode;
  
  /* nextfreenode records how many calls to new_Tnode_tree we made; therefore, at this stage, it
     holds the total number of nodes (including leaves) in the tree */
  
  theTree->numnodes = nextfreenode;
  nodes = free_util( nodes );
  heights = free_util( heights );
  subtreesizes = free_util( subtreesizes );

  /* The caller of the function should free the tree when finished with it */
  
  return theTree;

}
