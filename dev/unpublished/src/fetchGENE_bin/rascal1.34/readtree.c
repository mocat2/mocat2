#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>
#include "rascal.h"

/*
 *   Prototypes
 */
static void create_tree(FILE *fd,INODEPTR ptree, INODEPTR parent, INODEPTR *lptr, INODEPTR *ptrs, IN_TREEPTR itree);
static void create_node(INODEPTR pptr, INODEPTR parent);
static INODEPTR insert_node(INODEPTR pptr);
static void skip_space(FILE *fd);
static INODEPTR avail(void);
static void set_info(INODEPTR p, INODEPTR parent, sint pleaf, char *pname, float pdist);
static INODEPTR reroot(INODEPTR ptree, sint nseqs,INODEPTR *leafptr, INODEPTR *ptrs);
static INODEPTR insert_root(INODEPTR p, float diff);
static float calc_root_mean(INODEPTR root, float *maxdist,INODEPTR *leafptr);
static float calc_mean(INODEPTR nptr, float *maxdist, sint nseqs,INODEPTR *leafptr);
static void order_nodes(INODEPTR *leafptr);
static void free_tree_nodes(INODEPTR p);

static sint status;
/*
   read a phylip format tree file into a structure IN_TREE
	treefile	- the name of the phylip file
	seqs		- the sequences, used to check that the tree sequence names
			  correspond to the alignment sequence names
			  set to NULL, if you don't want to compare with the alignment
	fseq		- the first sequence in 'seqs' to be used for the name check
	nseq		- the no. of sequences in 'seqs' to be used for the name check
	itree		- the structure that will contain the tree
   returns -1 for an error, 1 for a tree with branch lengths, 0 for a tree with no lengths
*/

sint read_tree(char *treefile, SEQ *seqs, sint fseq, sint nseq,IN_TREEPTR itree)
{

  FILE *fd;
  char c;
  char name1[MAXNAMES+1], name2[MAXNAMES+1];
  sint i, j, k;
  sint ret;
  Boolean found;
  INODEPTR *lptr;		/* contains a list of pointers to the tree leaves,
				   ie. the sequences */
  INODEPTR *ptrs;		/* contains a list of pointers to all tree nodes */
  INODEPTR seq_tree;		/* a pointer to the initial root of the tree */

  if ((fd = fopen(treefile, "r")) == NULL)
    {
      error("cannot open treefile %s", treefile);
      return((sint)-1);
    }

  skip_space(fd);
  c = (char)getc(fd);
  if (c != '(')
    {
      error("Wrong format in tree file %s", treefile);
      return((sint)-1);
    }
  rewind(fd);

  itree->rooted_tree = TRUE;
  itree->distance_tree = TRUE;
  itree->nnodes = 0;
  itree->nleaves = 0;

/*
  Allocate memory for tree
*/
  ptrs = (INODEPTR *)ckalloc(3*(nseq-fseq+1) * sizeof(INODEPTR));
  lptr = (INODEPTR *)ckalloc((nseq-fseq+1) * sizeof(INODEPTR));
  itree->leafptr = (INODEPTR *)ckalloc((nseq+1) * sizeof(INODEPTR));
  
  seq_tree = avail();
  set_info(seq_tree, NULL, 0, "", (float)0.0);

  status=1;
  create_tree(fd,seq_tree,NULL,lptr,ptrs,itree);
  fclose(fd);


  if (seqs!= NULL && itree->nleaves != nseq-fseq)
     {
         error("tree not compatible with alignment\n(%d sequences in alignment and %d in tree", (pint)nseq-fseq,(pint)itree->nleaves);
         return((sint)-1);
     }

  /*if(status==0) warning("negative distances in tree - cannot calculate sequence weights");*/

  lptr[itree->nleaves]=NULL;
  ptrs[itree->nnodes]=NULL;

/*
  If the tree is unrooted, reroot the tree - ie. minimise the difference
  between the mean root->leaf distances for the left and right branches of
  the tree.
*/

  if (itree->distance_tree == FALSE)
     {
  	if (itree->rooted_tree == FALSE)
          {
       	     error("input tree is unrooted and has no distances.\nCannot align sequences");
             return((sint)-1);
          }
     }

  if (itree->rooted_tree == FALSE)
     {
        itree->root = reroot(seq_tree, itree->nleaves,lptr,ptrs);
     }
  else
     {
        itree->root = seq_tree;
     }

/*
  calculate the 'order' of each node.
*/
  order_nodes(lptr);

  if (itree->nleaves >= 2)
     {
/*
  If there are more than three sequences....
*/
/*
  assign the sequence nodes (in the same order as in the alignment file)
*/
      if(seqs!=NULL)
      for (i=fseq; i<nseq; i++)
       {
         if (strlen(seqs[i].name) > MAXNAMES)
             warning("name %s is too long for PHYLIP tree format (max %d chars)", seqs[i].name,MAXNAMES);

         for (k=0; k< strlen(seqs[i].name) && k<MAXNAMES ; k++)
           {
             c = seqs[i].name[k];
             if ((c>0x40) && (c<0x5b)) c=c | 0x20;
             if (c == ' ') c = '_';
             name2[k] = c;
           }
         name2[k]='\0';
         found = FALSE;
         for (j=0; j<itree->nleaves; j++)
           {
            for (k=0; k< strlen(lptr[j]->name) && k<MAXNAMES ; k++)
              {
                c = lptr[j]->name[k];
                if ((c>0x40) && (c<0x5b)) c=c | 0x20;
                name1[k] = c;
              }
            name1[k]='\0';
            if (strcmp(name1, name2) == 0)
              {
                itree->leafptr[i] = lptr[j];
                found = TRUE;
              }
           }
         if (found == FALSE)
           {
             error("tree not compatible with alignment:\n%s not found", name2);
             ret=0;
           }
       }
      else
	for(i=0;i<itree->nleaves;i++) 
		itree->leafptr[i] = lptr[i];

     }

   if(status==0)
	ret=0;
   else if(itree->distance_tree)
   	ret=1;
   else
   	ret=(-1);

   ptrs=ckfree((void *)ptrs);
   lptr=ckfree((void *)lptr);
   return ret;
}

/* A bit dirty, but char ch is needed here for the recursive function create_tree() */
static char ch;
static void create_tree(FILE *fd,INODEPTR ptree, INODEPTR parent, INODEPTR *lptr, INODEPTR *ptrs, IN_TREEPTR itree)
{
   INODEPTR p;
   sint i, type;
   float dist;
   char name[MAXNAMES+1];

/*
  is this a node or a leaf ?
*/
  skip_space(fd);
  ch = (char)getc(fd);
  if (ch == '(')
    {  
/*
   this must be a node....
*/
      type = NODE;
      name[0] = '\0';
      ptrs[itree->nnodes++] = ptree;

      create_node(ptree, parent);

      p = ptree->left;
      create_tree(fd, p, ptree, lptr, ptrs,itree);
           
      if ( ch == ',')
       {
          p = ptree->right;
          create_tree(fd, p, ptree, lptr, ptrs,itree);
          if ( ch == ',')
            {
               ptree = insert_node(ptree);
               ptrs[itree->nnodes++] = ptree;
               p = ptree->right;
               create_tree(fd, p, ptree, lptr, ptrs,itree);
               itree->rooted_tree = FALSE;
            }
       }

      skip_space(fd);
      ch = (char)getc(fd);
    }
/*
   ...otherwise, this is a leaf
*/
  else
    {
      type = LEAF;
      ptrs[itree->nnodes++] = lptr[itree->nleaves++] = ptree;
/*
   get the sequence name
*/
      name[0] = ch;
      ch = (char)getc(fd);
      i = 1;
      while ((ch != ':') && (ch != ',') && (ch != ')'))
        {
          if (i < MAXNAMES) name[i++] = ch;
          ch = (char)getc(fd);
        }
      name[i] = '\0';
      if (ch != ':')
         {
           itree->distance_tree = FALSE;
           dist = 0.0;
         }
    }

/*
   get the distance information
*/
  dist = 0.0;
  if (ch == ':')
     {
       skip_space(fd);
       fscanf(fd,"%f",&dist);
       skip_space(fd);
       ch = (char)getc(fd);
     }
   set_info(ptree, parent, type, name, dist);

   if(dist<0.0) status=0;
}

static void create_node(INODEPTR pptr, INODEPTR parent)
{
  INODEPTR t;

  pptr->parent = parent;
  t = avail();
  pptr->left = t;
  t = avail();
  pptr->right = t;
    
}

static INODEPTR insert_node(INODEPTR pptr)
{

   INODEPTR newnode;

   newnode = avail();
   create_node(newnode, pptr->parent);

   newnode->left = pptr;
   pptr->parent = newnode;

   set_info(newnode, pptr->parent, NODE, "", (float)0.0);

   return(newnode);
}

static void skip_space(FILE *fd)
{
  int   c;
  
  do
     c = getc(fd);
  while(isspace(c));

  ungetc(c, fd);
}

static INODEPTR avail(void)
{
   INODEPTR p;
   p = ckalloc(sizeof(INODE));
   p->left = NULL;
   p->right = NULL;
   p->parent = NULL;
   p->dist = 0.0;
   p->leaf = 0;
   p->order = 0;
   p->name[0] = '\0';
   return(p);
}

void free_tree(IN_TREEPTR itree)
{
   free_tree_nodes(itree->root);
      
   itree->leafptr=ckfree((void *)itree->leafptr);
}

static void free_tree_nodes(INODEPTR p)
{
   if (p==NULL) return;
   if (p->left != NULL)
     {
       free_tree_nodes(p->left);
     }
   if (p->right != NULL)
     {
       free_tree_nodes(p->right);
     }
   p->left = NULL;
   p->right = NULL;
   p=ckfree((void *)p);   
}

static void set_info(INODEPTR p, INODEPTR parent, sint pleaf, char *pname, float pdist)
{
   p->parent = parent;
   p->leaf = pleaf;
   p->dist = pdist;
   p->order = 0;
   strcpy(p->name, pname);
   if (p->leaf == TRUE)
     {
        p->left = NULL;
        p->right = NULL;
     }
}

/*
Let nodedist[i] be the distance from leaf i to a node.
Calculate the sum of nodedist[i] for each leaf on the left side of the node
and the sum of nodedist[i] for each leaf on the right side of the new root.
Find a new root node such that the difference between the left and right sums
is a minimum.
*/
static INODEPTR reroot(INODEPTR ptree, sint nseqs,INODEPTR *lptr, INODEPTR *ptrs)
{

   INODEPTR p, rootnode, rootptr;
   float   diff, mindiff = 0.0, mindepth = 1.0, maxdist;
   sint   i;

/*
  find the difference between the means of leaf->node
  distances on the left and on the right of each node
*/
   rootptr = ptree;
   for (i=0; (p=ptrs[i])!=NULL; i++)
     {
        if (p->parent == NULL)
           diff = calc_root_mean(p, &maxdist,lptr);
        else
           diff = calc_mean(p, &maxdist, nseqs,lptr);

        if ((diff == 0) || ((diff > 0) && (diff < 2 * p->dist)))
          {
              if ((maxdist < mindepth) || (i==0))
                 {
                    rootptr = p;
                    mindepth = maxdist;
                    mindiff = diff;
                 }
           }

     }

/*
  insert a new node as the ancestor of the node which produces the shallowest
  tree.
*/
   if (rootptr == ptree)
     {
        mindiff = rootptr->left->dist + rootptr->right->dist;
        rootptr = rootptr->right;
     }
   rootnode = insert_root(rootptr, mindiff);
  
   diff = calc_root_mean(rootnode, &maxdist,lptr);

   return(rootnode);
}

static INODEPTR insert_root(INODEPTR p, float diff)
{
   INODEPTR newp, prev, q, t;
   float dist, prevdist,td;

   newp = avail();

   t = p->parent;
   prevdist = t->dist;

   p->parent = newp;

   dist = p->dist;

   p->dist = diff / 2;
   if (p->dist < 0.0) p->dist = 0.0;
   if (p->dist > dist) p->dist = dist;

   t->dist = dist - p->dist; 

   newp->left = t;
   newp->right = p;
   newp->parent = NULL;
   newp->dist = 0.0;
   newp->leaf = NODE;

   if (t->left == p) t->left = t->parent;
   else t->right = t->parent;

   prev = t;
   q = t->parent;

   t->parent = newp;

   while (q != NULL)
     {
        if (q->left == prev)
           {
              q->left = q->parent;
              q->parent = prev;
              td = q->dist;
              q->dist = prevdist;
              prevdist = td;
              prev = q;
              q = q->left;
           }
        else
           {
              q->right = q->parent;
              q->parent = prev;
              td = q->dist;
              q->dist = prevdist;
              prevdist = td;
              prev = q;
              q = q->right;
           }
    }

/*
   remove the old root node
*/
   q = prev;
   if (q->left == NULL)
      {
         dist = q->dist;
         q = q->right;
         q->dist += dist;
         q->parent = prev->parent;
         if (prev->parent->left == prev)
            prev->parent->left = q;
         else
            prev->parent->right = q;
         prev->right = NULL;
      }
   else
      {
         dist = q->dist;
         q = q->left;
         q->dist += dist;
         q->parent = prev->parent;
         if (prev->parent->left == prev)
            prev->parent->left = q;
         else
            prev->parent->right = q;
         prev->left = NULL;
      }

   return(newp);
}

static float calc_root_mean(INODEPTR root, float *maxdist,INODEPTR *lptr)
{
   float dist , lsum = 0.0, rsum = 0.0, lmean,rmean,diff;
   INODEPTR p;
   sint i;
   sint nl, nr;
   sint direction;
/*
   for each leaf, determine whether the leaf is left or right of the root.
*/
   dist = (*maxdist) = 0;
   nl = nr = 0;
   for (i=0; (p=lptr[i])!=NULL; i++)
     {
         dist = 0.0;
         while (p->parent != root)
           {
               dist += p->dist;
               p = p->parent;
           }
         if (p == root->left) direction = LEFT;
         else direction = RIGHT;
         dist += p->dist;

         if (direction == LEFT)
           {
             lsum += dist;
             nl++;
           }
         else
           {
             rsum += dist;
             nr++;
           }
        if (dist > (*maxdist)) *maxdist = dist;
     }

   lmean = lsum / nl;
   rmean = rsum / nr;

   diff = lmean - rmean;
   return(diff);
}


static float calc_mean(INODEPTR nptr, float *maxdist, sint nseqs,INODEPTR *lptr)
{
   float dist , lsum = 0.0, rsum = 0.0, lmean,rmean,diff;
   INODEPTR p, *path2root;
   float *dist2node;
   sint depth = 0, i, n = 0;
   sint nl , nr;
   sint direction;

	path2root = (INODEPTR *)ckalloc(nseqs * sizeof(INODEPTR));
	dist2node = (float *)ckalloc(nseqs * sizeof(float));
   nl = nr = 0;
   (*maxdist) = 0;
/*
   determine all nodes between the selected node and the root;
*/
   get_path2root(nptr,path2root,dist2node,&depth);
 
/*
   for each leaf, determine whether the leaf is left or right of the node.
   (where RIGHT = descendant, LEFT = not descendant)
*/
   for (i=0; (p=lptr[i])!=NULL; i++)
     {
       if (p == nptr)
         {
            direction = RIGHT;
            dist = 0.0;
         }
       else
         {
            direction = LEFT;
/*
   find the common ancestor.
*/
	    n=common_ancestor(p,path2root,depth,&dist);
            if (path2root[n] == nptr) direction = RIGHT;
         }

         if (direction == LEFT)
           {
             lsum += dist;
             lsum += dist2node[n-1];
             nl++;
           }
         else
           {
             rsum += dist;
             nr++;
           }

        if (dist > (*maxdist)) *maxdist = dist;
     }

	dist2node=ckfree((void *)dist2node);
	path2root=ckfree((void *)path2root);
	
   lmean = lsum / nl;
   rmean = rsum / nr;
   
   diff = lmean - rmean;
   return(diff);
}

void get_path2root(INODEPTR nptr,INODEPTR *path2root,float *dist2node,sint *depth)
{
	float dist;
	INODEPTR p;

	*depth = 0;
	dist = 0.0;
	p = nptr;
	while (p != NULL) {
		path2root[*depth] = p;
		dist += p->dist;
		dist2node[*depth] = dist;
		p = p->parent;
		(*depth)++;
	}
}

/*
   find the intersection of the path from node p to the root, and the path in path2root
   depth = size of path2root
   dist = dist from the top of path2root and the intersection
*/
sint common_ancestor(INODEPTR p,INODEPTR *path2root,sint depth,float *dist)
{
	INODEPTR ancestor;
	Boolean found;
	sint n,j;

        found = FALSE;
        n = 0;
	ancestor=p;
        (*dist) = 0.0;
        while ((found == FALSE) && (ancestor->parent != NULL)) {
               for (j=0; j< depth; j++)
                 if (ancestor->parent == path2root[j]) { 
                      found = TRUE;
                      n = j;
                 }
               (*dist) += ancestor->dist;
               ancestor = ancestor->parent;
        }

	return n;
}

static void order_nodes(INODEPTR *lptr)
{
   sint i;
   INODEPTR p;

   for (i=0; (p=lptr[i])!=NULL; i++)
     {
        while (p != NULL)
          {
             p->order++;
             p = p->parent;
          }
     }
}

/*
   Calculates a matrix of pairwise sequence distances from a phylogenetic tree
	seqs	the sequences
	nseqs	number of sequences
	itree	the phylogenetic tree
   returns the distance matrix
*/
double ** pw_distances_from_tree(SEQ *seqs,sint nseqs,IN_TREEPTR itree)
{
   sint depth = 0, i,j, k, n;
   sint found;
   sint nerrs, *seq1,*seq2;
   INODEPTR p, *path2root;
   float dist;
   float *dist2node, *bad_dist;
   double **tmat;
   char err_mess[1024],err1[MAXLINE];

	tmat = (double **) ckalloc( (nseqs+1) * sizeof (double *) );
        for(i=0;i<nseqs;i++)
                tmat[i] = (double *)ckalloc( (nseqs+1) * sizeof (double) );
        for(i=0;i<nseqs;i++)
                for(j=0;j<nseqs;j++)
                        tmat[i][j]=0.0;

   path2root = (INODEPTR *)ckalloc((nseqs) * sizeof(INODEPTR));
   dist2node = (float *)ckalloc((nseqs) * sizeof(float));
   seq1 = (sint *)ckalloc((nseqs) * sizeof(sint));
   seq2 = (sint *)ckalloc((nseqs) * sizeof(sint));
   bad_dist = (float *)ckalloc((nseqs) * sizeof(float));

   if (nseqs < 2) return(NULL);
/*
   for each leaf, determine all nodes between the leaf and the root;
*/
      for (i = 0;i<nseqs; i++)
       { 
          depth = 0;
          dist = 0.0;
          p = itree->leafptr[i];
          while (p != NULL)
            {
                path2root[depth] = p;
                dist += p->dist;
                dist2node[depth] = dist;
                p = p->parent;
                depth++;
            }
 
/*
   for each pair....
*/
          for (j=0; j < i; j++)
            {
              p = itree->leafptr[j];
              dist = 0.0;
/*
   find the common ancestor.
*/
              found = FALSE;
              n = 0;
              while ((found == FALSE) && (p->parent != NULL))
                {
                    for (k=0; k< depth; k++)
                      if (p->parent == path2root[k])
                         { 
                           found = TRUE;
                           n = k;
                         }
                    dist += p->dist;
                    p = p->parent;
                }
   
              tmat[j][i]=tmat[i][j] = dist + dist2node[n-1];
            }
        }

	nerrs = 0;
        for (i=0;i<nseqs;i++)
          {
             tmat[i][i] = 0.0;
             for (j=0;j<i;j++)
               {
                  if (tmat[i][j] < 0.01) tmat[i][j] = 0.01;
                  if (tmat[i][j] > 1.0) {
                  	if (tmat[i][j] > 1.1) {
				if(nerrs<nseqs) {
                  			seq1[nerrs] = i;
                  			seq2[nerrs] = j;
                  			bad_dist[nerrs] = (float)tmat[i][j];
				}
                  		nerrs++;
                  	}
                    tmat[j][i]=tmat[i][j] = 1.0;
                  }
               }
          }
        if (nerrs>0) 
          {
             strcpy(err_mess,"The following sequences are too divergent to be aligned:\n");
             for (i=0;i<nerrs && i<5;i++)
              {
             	sprintf(err1,"           %s and %s (distance %1.3f)\n",
             	                        seqs[seq1[i]+1].name,
					seqs[seq2[i]+1].name,bad_dist[i]);
             	strcat(err_mess,err1);
              }
	     strcat(err_mess,"(All distances should be between 0.0 and 1.0)\n");
	     strcat(err_mess,"This may not be fatal but you have been warned!\n");
             strcat(err_mess,"SUGGESTION: Remove one or more problem sequences and try again");
          }

   path2root=ckfree((void *)path2root);
   dist2node=ckfree((void *)dist2node);
   seq1=ckfree((void *)seq1);
   seq2=ckfree((void *)seq2);
   bad_dist=ckfree((void *)bad_dist);

   return(tmat);
}

