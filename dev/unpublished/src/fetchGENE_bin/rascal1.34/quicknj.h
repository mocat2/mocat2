/*  Last edited: Aug 27 16:31 1999 (klh) */
/**********************************************************************
 ** FILE: align.h
 ** NOTES:
 **   Functions and types for the manipulation of multiple sequence
 **   alignments
**********************************************************************/

#ifndef _ALIGN
#define _ALIGN

#define RES_BLOCKSIZE 20
#define SEQ_BLOCKSIZE 10
#define MAX_LINE_LEN  4096

/******************* structure definitions ****************************/

struct Alignment {
  unsigned int numseqs;
  unsigned int length;
  struct Sequence **seqs;
};

/***********************
 This structure represents an unrooted binary tree of three or more 
   nodes; Two node trees, although trivial, can be modelled by having
   one of the three children as null. A rooted tree can be represented
   by a Tnode only, or by having two of the children of the trichotomy
   as NULL.
*************************/

struct Tree {
  struct Tnode *child[3];
  unsigned int numnodes;
};

struct Cluster {
  unsigned int clustersize;
  struct Sequence **members;
  struct Sequence *consensus;
  struct DistanceMatrix *matrix;
};


/*
  Clusters contain groups of identical sequence. I intend to investigate
  methods where clusters contain groups of similar (not necessarily
  identical) sequences. The DistanceMatric field, although not currently
  used, will allow for the building of trees from these clusters
*/


struct ClusterGroup {
  unsigned int numclusters;
  struct Cluster **clusters;
  struct DistanceMatrix *matrix;
};

/********************** function prototypes ***************************/

/**********************************************************************
 FUNCTION: free_Alignment
 DESCRIPTION: 
   Frees the memory used by the given alignment
 ARGS:
   An Alignment
 RETURNS: A null pointer
 NOTES:
 **********************************************************************/
void *free_Alignment( struct Alignment *);


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
struct Alignment *read_MUL_Alignment( FILE * );


/**********************************************************************
 FUNCTION: read_Stockholm_Alignment
 DESCRIPTION: 
   This function fills a simple alignment structure from an
   file assumed to be in Stockholm format. At this stage, I am ignoring
   all mark-up information (all lines beginning with '#' are ignored).
   The function also allows for wrapped alignments. Note that Pfam
   alignments in MUL format will be handled correctly by this function
 ARGS:
   FILE *
 RETURNS: struct Alignment (align.h)
 NOTES:
 **********************************************************************/
struct Alignment *read_Stockholm_Alignment( FILE *);



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
void write_MUL_Alignment( FILE *, struct Alignment *);

#endif

/*  Last edited: Aug 24 14:28 1999 (klh) */
/**********************************************************************
 ** FILE: buildtree.h
 ** NOTES:
 **  Contains functions for building trees from distance matrices
 **  (and vice versa)
 **********************************************************************/

#ifndef _BUILDTREE
#define _BUILDTREE

/******************* structure definitions ****************************/


/********************** function prototypes ***************************/


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
void export_distances_buildtree( struct Tree *, struct DistanceMatrix *);


/**********************************************************************
 FUNCTION: find_path_buildtree
 DESCRIPTION: 
   Given a leaf node, recursively calculates the branch-length distance
   from the node to all other leaves in the tree, and places it in the
   appropriate part of the DistanceMatrix
 ARGS: 
   unsigned int (the node number from which we are finding all distances)
   Tnode (the current node under consideration)
   DistanceMatrix
   Boolean array (to store which nodes have already been considered)
   unsigned int (the size of this boolean array)
 RETURNS: 
 NOTES: 
 **********************************************************************/
void find_path_buildtree( unsigned int,
			  struct Tnode *,
			  struct DistanceMatrix *,
			  double,
			  unsigned int *);


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
   unsigned int (the size of this boolean array)
 RETURNS: 
 NOTES: 
 **********************************************************************/
void leaf_find_buildtree( struct Tnode *, 
			  struct DistanceMatrix *, 
			  unsigned int *,
			  unsigned int);



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
struct Tree *neighbour_joining_buildtree( struct ClusterGroup *,
					  unsigned int);


/**********************************************************************
 FUNCTION: UPGMA_buildtree
 DESCRIPTION: 
   Returns a phylogenetic tree of the sequences in the 
   given alignment, using the Unweighted Pair-Group method based on
   Arithmentic Averages (UPGMA)
 ARGS: 
   A ClusterGroup pointer (cluster.h)
   Boolean, for whether to calc information needed for later bootstrapping
 RETURNS: 
   struct Tree (trees.h)
 NOTES: The function allocates all the memory necessary for the tree.
   The caller should call free_Tnode (tree.h) to free this memory when
   the tree is no longer needed

   This algorithm produces a rooted tree; hence the returned Tree
   will have the root in child[0]; child[1] and child[2] (used to
   represent the trichotomy of unrooted trees) will be NULL 
 **********************************************************************/
struct Tree *UPGMA_buildtree( struct ClusterGroup *, unsigned int);


#endif
/*  Last edited: Aug 25 15:19 1999 (klh) */
/**********************************************************************
 ** FILE: cluster.h
 ** NOTES:
 **  A DistanceMatrix should always be part of a Cluster 
 **  It makes no sense to have a set of pairwise distances without the 
 **  associated sequences (even if we just store their names)
 **********************************************************************/

#ifndef _CLUSTER
#define _CLUSTER

/******************* structure definitions ****************************/


/********************** function prototypes ***************************/



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
struct ClusterGroup *alignment_to_ClusterGroup( struct Alignment *,unsigned int);

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
struct Cluster *clone_Cluster( struct Cluster *);

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
struct Alignment *consensus_aln_from_ClusterGroup( struct ClusterGroup *);

/********************************************************************* 
 FUNCTION: empty_Cluster
 DESCRIPTION: 
   This function handles the simple task of allocating the space
   for a new Cluster.
 RETURNS: struct Cluster *
 ARGS: 
 NOTES: 
 *********************************************************************/
struct Cluster *empty_Cluster( void );

/********************************************************************* 
 FUNCTION: empty_ClusterGroup
 DESCRIPTION: 
   This function handles the simple task of allocating the space
   for a new ClusterGroup
 RETURNS: struct Cluster *
 ARGS: 
 NOTES: 
 *********************************************************************/
struct ClusterGroup *empty_ClusterGroup( void );

/********************************************************************* 
 FUNCTION: free_Cluster
 DESCRIPTION: 
   This function releases the memory used by this Cluster and all of its
   members
 RETURNS: A null pointer
 ARGS: 
   struct Cluster *
 NOTES: 
 *********************************************************************/
void *free_Cluster( struct Cluster *);

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
void *free_ClusterGroup( struct ClusterGroup *);

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
void *merge_Cluster( struct Cluster *, struct Cluster *);

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
struct Cluster *single_Sequence_Cluster( struct Sequence *);

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
struct ClusterGroup *single_Cluster_ClusterGroup( struct Cluster *);


#endif
/*  Last edited: Aug 25 15:25 1999 (klh) */
/**********************************************************************
 ** FILE: distancemat.h
 ** NOTES:
 **   Functions and types for the manipulation of Distance Matrices
 **********************************************************************/

#ifndef _DISTANCEMAT
#define _DISTANCEMAT

/******************* structure definitions ****************************/

struct DistanceMatrix {
  double **data;
  int size;
};


/********************** function prototypes ***************************/


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
void calc_DistanceMatrix(struct DistanceMatrix *, 
			 struct Alignment *,
			 unsigned int,
			 unsigned int );

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
struct DistanceMatrix *clone_DistanceMatrix( struct DistanceMatrix *);

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
struct DistanceMatrix *empty_DistanceMatrix( unsigned int );

/*********************************************************************
 FUNCTION: free_DistanceMatrix
 DESCRIPTION: 
   Frees the memory for the given distance matrix
 RETURNS:
 ARGS: 
   struct DistanceMatrix *
 NOTES: 
 *********************************************************************/
void *free_DistanceMatrix( struct DistanceMatrix *);

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
double index_DistanceMatrix( struct DistanceMatrix *, unsigned int, unsigned int );

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
void print_DistanceMatrix( FILE *handle,  struct DistanceMatrix * );


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
struct DistanceMatrix *read_phylip_DistanceMatrix( FILE *, struct Alignment **);

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

void write_phylip_DistanceMatrix( FILE *, struct DistanceMatrix *, struct Alignment *);


#endif
/*  Last edited: Feb  1 15:25 2002 (klh) */
/**********************************************************************
 ** FILE: options.h
 ** DESCRIPTION:
 **  Rudimentary provision of command-line options
 **********************************************************************/

#ifndef _GETOPTIONS
#define _GETOPTIONS 


/************************** constants *********************************/

#define NO_ARGS    0 
#define INT_ARG    1
#define FLOAT_ARG  2
#define CHAR_ARG   3
#define STRING_ARG 4

/******************* structure definitions ****************************/

struct Option {
  char *name;         /* name of option, e.g. "-option" */
  unsigned int type;  /* for typechecking, e.g. INT_ARG     */
};


/******************* function prototypes ****************************/

/*********************************************************************
 FUNCTION: get_option
 DESCRIPTION: 
   Gets an option from the given command line
 RETURNS:
   1, if a valid option was found
   0, if no valid option was found and option parsing is therefore 
      complete
 ARGS: 
 NOTES: 
*********************************************************************/
unsigned int get_option(int, char **, struct Option *, unsigned int,
			char *, unsigned int *, char **, char **);


#endif

/*  Last edited: Nov 15 12:24 1999 (klh) */
/**********************************************************************
 ** FILE: sequence.h
 ** NOTES:
 **   Functions and structures for the manipulation of protein 
 **   sequences. Only minimal functionality is needed
 **********************************************************************/

#ifndef _SEQUENCE
#define _SEQUENCE

#define MAX_NAME_LENGTH 25

/******************* structure definitions ****************************/

struct Sequence {
  unsigned int length;
  char *name;
  char *seq;
  char *sec_struct;
  char *surf_acc;
  char *trans_mem;
  char *post_prob;
  char *lig_bind;
};

/********************** function prototypes ***************************/


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
struct Sequence *clone_Sequence( struct Sequence *);

/**********************************************************************
 FUNCTION: empty_Sequence
 DESCRIPTION: 
   Creates and returns an new, empty Sequence object
 ARGS: 
 RETURNS:
   A pointer to a Sequence structure
 NOTES:
 **********************************************************************/
struct Sequence *empty_Sequence( void );

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
void *free_Sequence( struct Sequence *);


#endif
/*  Last edited: Jun 21 12:03 2001 (klh) */
/**********************************************************************
 ** FILE: tree.h
 ** NOTES:
 **  Functions and types for the manipulation of trees
 **********************************************************************/

#ifndef _TREE
#define _TREE

#define LINE_BUFFER_SIZE 40


/******************* structure definitions ****************************/

/******************
   This structure represents a tree-node. Everything is self-explanatory,
   except that much of the tree-building code relies upon a correspondence
   between node-numbers of leaf nodes, and sequence/cluster numbers. This 
   is particulary important when re-constructing distance matrices from 
   trees, and for bootstrapping trees. This may be troublesome if a tree is
   read in from a file... 
******************/


struct Tnode {
  struct Tnode *left;
  struct Tnode *right;
  struct Tnode *parent;
  double distance;
  unsigned int nodenumber;
  struct Cluster *clust;
  unsigned int bootstrap;
  unsigned int *child_ids;
};




/********************** function prototypes ***************************/

/********************************************************************* 
 FUNCTION: assign_nodesnumbers_Tnode
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
unsigned int assign_nodenumbers_Tnode( struct Tnode *, unsigned int);


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
struct Tnode *clone_Tnode( struct Tnode *);

/********************************************************************* 
 FUNCTION: clone_Tree
 DESCRIPTION: 
   This function makes a complete copy of the Tree and returns it
 RETURNS: struct Tree *
 ARGS: 
   struct Tree *
 NOTES: 
 *********************************************************************/
struct Tree *clone_Tree( struct Tree *);


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
void compare_to_bootstrap_sample_Tnode( struct Tnode *, struct Tnode *,
					unsigned int, unsigned int);

/********************************************************************* 
 FUNCTION: empty_Tree
 DESCRIPTION: 
   Creates and returns a tree with null nodes
 RETURNS: struct Tree * (trees.h)
 ARGS: 
 NOTES:
 *********************************************************************/
struct Tree *empty_Tree( void );

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
void *free_Tnode( struct Tnode *);

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
void *free_Tree( struct Tree *);

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
   function must be freed by a call to free_Tnode_tree
 *********************************************************************/
struct Tree *get_root_Tnode( struct Tree * );

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
struct Tnode *new_interior_Tnode( unsigned int );

/********************************************************************* 
 FUNCTION: new_leaf_Tnode (int, char *)
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
struct Tnode *new_leaf_Tnode( unsigned int, struct Cluster *);

/********************************************************************* 
 FUNCTION: read_newhampshire_Tnode
 DESCRIPTION: 
   Construct a Tnode from the given file handle
 RETURNS: The total number of nodes constructed as a result of the call
 ARGS: 
   File handle (assumed in New Hampshire format)
 NOTES:
   Numbering of internal nodes is abandoned in favour of giving leaf
   nodes numbers from 0 to the number of leaves in the tree. This makes
   them handy indices into a distance matrix for example
 *********************************************************************/
unsigned int read_newhampshire_Tnode( FILE *, struct Tnode **, 
				      struct Tnode *, unsigned int);

/********************************************************************* 
 FUNCTION: read_newhampshire_Tree
 DESCRIPTION: 
   Constructs a tree from the given file handle
 RETURNS:  Tree *
 ARGS: A handle to the file (assumed in New Hamshire format)
 NOTES:
   The convention that the numbering of the interior nodes starts
   after all leaf nodes have been numbered is is violated with this
   method; the numbering is performed in a bottom-up, left-to-right 
   fashion, so that the nodes in the left subtree all have number-ids
   strictly less than those in the right subtree

 *********************************************************************/
struct Tree *read_newhampshire_Tree( FILE *);

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
void scale_bootstraps_Tnode( struct Tnode *, unsigned int);


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
void scale_bootstraps_Tree( struct Tree *, unsigned int);


/********************************************************************* 
 FUNCTION: update_bootstraps_Tree
 DESCRIPTION: 
   Updates the bootstrap values of the destination Tree, according to
   the topology of the given sample tree
 RETURNS:
 ARGS:
   Destination tree
   Sample tree
   The number of leaf nodes in the trees
 NOTES:
 *********************************************************************/
void update_bootstraps_Tree( struct Tree *, struct Tree *, unsigned int);

/********************************************************************* 
 FUNCTION: update_bootstraps_Tnode
 DESCRIPTION: 
   Updates the bootstrap values of the destination Tnode, according to
   the topology of the given sample Tnode
 RETURNS:
 ARGS:
   Destination Tnode
   Sample Tnode
   The number of leaf nodes in the tree
   Boolean for whether the trees are binary or not
 NOTES:
 *********************************************************************/
void update_bootstraps_Tnode( struct Tnode *, struct Tnode *, 
			      unsigned int, unsigned int);

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
void write_clustering_data_Tnode( FILE *, struct Tnode *);


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
void write_debug_Tnode( FILE *, struct Tnode *, unsigned int);

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
void write_debug_Tree( FILE *, struct Tree *);

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
void write_MUL_flattened_Tnode( FILE *, struct Tnode *);

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
void write_MUL_flattened_Tree( FILE *, struct Tree *);


/********************************************************************* 
 FUNCTION: write_newhampshire_Tnode
 DESCRIPTION: 
   prints the given Tree in 'New Hampshire' text format to the given 
   file handle
 RETURNS:
 ARGS: 
   File handle
   TNode *
 NOTES:
 *********************************************************************/
void write_newhampshire_Tnode( FILE *, struct Tnode *, unsigned int);

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
void write_newhampshire_Tree( FILE *, struct Tree *, unsigned int);



#endif
/*  Last edited: Feb  1 15:32 2002 (klh) */
/**********************************************************************
 ** FILE: util.h
 ** NOTES:
 **   This file contains general utiliy functions used throughout the 
 **   application, such as those for memory management and error
 **   messaging. I have used it as a place to put other general stuff
 **   until I have somewhere better to put it.
 **********************************************************************/

#ifndef _UTIL
#define _UTIL

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif


/********************** function prototypes ***************************/

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
void *calloc_util( size_t, size_t );


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
void *malloc_util( size_t );


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
void *realloc_util( void *, size_t );


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
void *free_util( void * );


/********************************************************************* 
 FUNCTION: fatal_util
 DESCRIPTION: 
   Prints the given formatted error message and exits
 RETURNS:
 ARGS:
   A format string + args, c.f. printf
 NOTES:
 *********************************************************************/
void fatal_util( char *, ... );


/********************************************************************* 
 FUNCTION: warning_util
 DESCRIPTION: 
   Prints the given formatted warning to stderr
 RETURNS:
 ARGS:
   A format string + args, c.f. printf
 NOTES:
 *********************************************************************/
void warning_util( char *, ... );


#endif
