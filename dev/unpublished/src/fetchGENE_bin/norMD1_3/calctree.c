#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>
#include "score.h"


/*
 *   Prototypes
 */
static void create_tree(treeptr ptree, treeptr parent);
static void create_node(treeptr pptr, treeptr parent);
static treeptr insert_node(treeptr pptr);
static void skip_space(FILE *fd);
static treeptr avail(void);
static void set_info(treeptr p, treeptr parent, int pleaf, char *pname, float pdist);
static treeptr reroot(treeptr ptree, int nseqs);
static treeptr insert_root(treeptr p, float diff);
static float calc_root_mean(treeptr root, float *maxdist);
static float calc_mean(treeptr nptr, float *maxdist, int nseqs);
static void order_nodes(void);
static int calc_weight(int leaf);
static void clear_tree_nodes(treeptr p);
static void get_path2root(treeptr  nptr,treeptr  *path2root,float *dist2node,int *depth);
static int common_ancestor(treeptr p,treeptr *path2root,int depth,float *dist);
void seq_dist(char *seqname, int nseqs);

extern void *ckfree(void *ptr);
extern void *ckalloc(size_t bytes);


/*
 *   Global variables
 */
extern char **names;
extern float *seq_score;

char ch;
FILE *fd;
treeptr *lptr;
treeptr *olptr;
treeptr *nptr;
treeptr *ptrs;
int nnodes = 0;
int ntotal = 0;
Boolean rooted_tree = TRUE;
static treeptr seq_tree,root;
static int numseq;

void calc_seq_weights(int nseqs, float *sweight)
{
  int   i;
  int   temp, sum, *weight;


/*
  If there are more than three sequences....
*/
  if (nseqs >= 2)
     {
/*
  Calculate sequence weights based on Phylip tree.
*/
      weight = (int *)calloc((nseqs+1) , sizeof(int));

      for (i=0; i<nseqs; i++)
           weight[i] = calc_weight(i);

/*
  Normalise the weights, such that the sum of the weights = nseqs
*/

         sum = 0;
         for (i=0; i<nseqs; i++)
            sum += weight[i];

         if (sum == 0)
          {
            for (i=0; i<nseqs; i++)
               weight[i] = 1;
            sum = i;
          }

         for (i=0; i<nseqs; i++)
           {
              sweight[i] = (float)(weight[i] * nseqs) / (float)sum;
           }

       free((void *)weight);

     }

   else
     {
/*
  Otherwise, use identity weights.
*/
        temp = INT_SCALE_FACTOR / nseqs;
        for (i=0; i<nseqs; i++)
           sweight[i] = temp;
     }

}

int read_tree(char *treefile, int nseqs)
{

  char c;
  char name1[MAXNAMES+1], name2[MAXNAMES+1];
  int i, j, k;
  Boolean found;

  numseq = 0;
  nnodes = 0;
  ntotal = 0;
  rooted_tree = TRUE;

#ifdef VMS
  if ((fd = fopen(treefile,"r","rat=cr","rfm=var")) == NULL)
#else
  if ((fd = fopen(treefile, "r")) == NULL)
#endif
    {
      fprintf(stderr,"cannot open %s", treefile);
      return((int)0);
    }

  skip_space(fd);
  ch = (char)getc(fd);
  if (ch != '(')
    {
      fprintf(stderr,"Wrong format in tree file %s", treefile);
      return((int)0);
    }
  rewind(fd);


/*
  Allocate memory for tree
*/
  nptr = (treeptr *)calloc(3*(nseqs+1) , sizeof(treeptr));
  ptrs = (treeptr *)calloc(3*(nseqs+1) , sizeof(treeptr));
  lptr = (treeptr *)calloc((nseqs+1) , sizeof(treeptr));
  olptr = (treeptr *)calloc((nseqs+1) , sizeof(treeptr));
  
  seq_tree = avail();
  set_info(seq_tree, NULL, 0, "", 0.0);

  create_tree(seq_tree,NULL);
  fclose(fd);


  if (numseq != nseqs)
     {
         fprintf(stderr,"tree not compatible with alignment\n(%d sequences in alignment and %d in tree", (int)nseqs,(int)numseq);
         return((int)0);
     }

/*
  If the tree is unrooted, reroot the tree - ie. minimise the difference
  between the mean root->leaf distances for the left and right branches of
  the tree.
*/

  if (rooted_tree == FALSE)
     {
        root = reroot(seq_tree, nseqs+1);
     }
  else
     {
        root = seq_tree;
     }

/*
  calculate the 'order' of each node.
*/
  order_nodes();

  if (numseq >= 2)
     {
/*
  If there are more than three sequences....
*/
/*
  assign the sequence nodes (in the same order as in the alignment file)
*/
      for (i=0; i<nseqs; i++)
       {
         if (strlen(names[i]) > MAXNAMES)
             fprintf(stderr,"name %s is too long for PHYLIP tree format (max %d chars)", names[i+1],MAXNAMES);

         for (k=0; k< strlen(names[i]) && k<MAXNAMES ; k++)
           {
             c = names[i][k];
             if ((c>0x40) && (c<0x5b)) c=c | 0x20;
             if (c == ' ') c = '_';
             name2[k] = c;
           }
         name2[k]='\0';
         found = FALSE;
         for (j=0; j<numseq; j++)
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
                olptr[i] = lptr[j];
                found = TRUE;
              }
           }
         if (found == FALSE)
           {
             fprintf(stderr,"tree not compatible with alignment:\n%s not found", name2);
             return((int)0);
           }
       }

     }
   return((int)1);
}

static void create_tree(treeptr ptree, treeptr parent)
{
   treeptr p;

   int i, type;
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
      ptrs[ntotal] = nptr[nnodes] = ptree;
      nnodes++;
      ntotal++;

      create_node(ptree, parent);

      p = ptree->left;
      create_tree(p, ptree);
           
      if ( ch == ',')
       {
          p = ptree->right;
          create_tree(p, ptree);
          if ( ch == ',')
            {
               ptree = insert_node(ptree);
               ptrs[ntotal] = nptr[nnodes] = ptree;
               nnodes++;
               ntotal++;
               p = ptree->right;
               create_tree(p, ptree);
               rooted_tree = FALSE;
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
      ptrs[ntotal++] = lptr[numseq++] = ptree;
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


}

static void create_node(treeptr pptr, treeptr parent)
{
  treeptr t;

  pptr->parent = parent;
  t = avail();
  pptr->left = t;
  t = avail();
  pptr->right = t;
    
}

static treeptr insert_node(treeptr pptr)
{

   treeptr newnode;

   newnode = avail();
   create_node(newnode, pptr->parent);

   newnode->left = pptr;
   pptr->parent = newnode;

   set_info(newnode, pptr->parent, NODE, "", 0.0);

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

static treeptr avail(void)
{
   treeptr p;
   p = calloc(1,sizeof(stree));
   p->left = NULL;
   p->right = NULL;
   p->parent = NULL;
   p->dist = 0.0;
   p->leaf = 0;
   p->order = 0;
   p->name[0] = '\0';
   return(p);
}

void clear_tree(treeptr p)
{
   clear_tree_nodes(p);
      
   free((void *)nptr);
   free((void *)ptrs);
   free((void *)lptr);
   free((void *)olptr);
}

static void clear_tree_nodes(treeptr p)
{
   if (p==NULL) p = root;
   if (p->left != NULL)
     {
       clear_tree_nodes(p->left);
     }
   if (p->right != NULL)
     {
       clear_tree_nodes(p->right);
     }
   p->left = NULL;
   p->right = NULL;
   free((void *)p);   
}

static void set_info(treeptr p, treeptr parent, int pleaf, char *pname, float pdist)
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

static treeptr reroot(treeptr ptree, int nseqs)
{

   treeptr p, rootnode, rootptr;
   float   diff, mindiff = 0.0, mindepth = 1.0, maxdist;
   int   i;
   Boolean first = TRUE;

/*
  find the difference between the means of leaf->node
  distances on the left and on the right of each node
*/
   rootptr = ptree;
   for (i=0; i<ntotal; i++)
     {
        p = ptrs[i];
        if (p->parent == NULL)
           diff = calc_root_mean(p, &maxdist);
        else
           diff = calc_mean(p, &maxdist, nseqs);

        if ((diff == 0) || ((diff > 0) && (diff < 2 * p->dist)))
          {
              if ((maxdist < mindepth) || (first == TRUE))
                 {
                    first = FALSE;
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
  
   diff = calc_root_mean(rootnode, &maxdist);

   return(rootnode);
}

static treeptr insert_root(treeptr p, float diff)
{
   treeptr newp, prev, q, t;
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

static float calc_root_mean(treeptr root, float *maxdist)
{
   float dist , lsum = 0.0, rsum = 0.0, lmean,rmean,diff;
   treeptr p;
   int i;
   int nl, nr;
   int direction;
/*
   for each leaf, determine whether the leaf is left or right of the root.
*/
   dist = (*maxdist) = 0;
   nl = nr = 0;
   for (i=0; i< numseq; i++)
     {
         p = lptr[i];
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

static int common_ancestor(treeptr p,treeptr *path2root,int depth,float *dist)
{
        treeptr ancestor;
        Boolean found;
        int n,j;

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



static void get_path2root(treeptr  nptr,treeptr  *path2root,float *dist2node,int *depth)
{
        float dist;
        treeptr  p;

        *depth = dist = 0;
        p = nptr;
        while (p != NULL) {
                path2root[*depth] = p;
                dist += p->dist;
                dist2node[*depth] = dist;
                p = p->parent;
                (*depth)++;
        }
}


static float calc_mean(treeptr nptr, float *maxdist, int nseqs)
{
   float dist , lsum = 0.0, rsum = 0.0, lmean,rmean,diff;
   treeptr p, *path2root;
   float *dist2node;
   int depth = 0, i,j , n = 0;
   int nl , nr;
   int direction, found;

	path2root = (treeptr *)calloc(nseqs , sizeof(treeptr));
	dist2node = (float *)calloc(nseqs , sizeof(float));
/*
   determine all nodes between the selected node and the root;
*/
   depth = (*maxdist) = dist = 0;
   nl = nr = 0;
   p = nptr;
   while (p != NULL)
     {
         path2root[depth] = p;
         dist += p->dist;
         dist2node[depth] = dist;
         p = p->parent;
         depth++;
     }
 
/*
   *nl = *nr = 0;
   for each leaf, determine whether the leaf is left or right of the node.
   (RIGHT = descendant, LEFT = not descendant)
*/
   for (i=0; i< numseq; i++)
     {
       p = lptr[i];
       if (p == nptr)
         {
            direction = RIGHT;
            dist = 0.0;
         }
       else
         {
         direction = LEFT;
         dist = 0.0;
/*
   find the common ancestor.
*/
         found = FALSE;
         n = 0;
         while ((found == FALSE) && (p->parent != NULL))
           {
               for (j=0; j< depth; j++)
                 if (p->parent == path2root[j])
                    { 
                      found = TRUE;
                      n = j;
                    }
               dist += p->dist;
               p = p->parent;
           }
         if (p == nptr) direction = RIGHT;
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

	free((void *)dist2node);
	free((void *)path2root);
	
   lmean = lsum / nl;
   rmean = rsum / nr;
   
   diff = lmean - rmean;
   return(diff);
}

static void order_nodes(void)
{
   int i;
   treeptr p;

   for (i=0; i<numseq; i++)
     {
        p = lptr[i];
        while (p != NULL)
          {
             p->order++;
             p = p->parent;
          }
     }
}


static int calc_weight(int leaf)
{

  treeptr p;
  float weight = 0.0;

  p = olptr[leaf];
  while (p->parent != NULL)
    {
       weight += p->dist / p->order;
       p = p->parent;
    }

  weight *= 100.0;

  return((int)weight);

}

void seq_dist(char *seqname, int nseqs)
{
	FILE *fd;
	int i,j,n;
	Boolean found;
	treeptr seqptr;
	float dist,*seqdist;
	treeptr *path2root;
	float *dist2node;
	int s,depth = 0;
 

/* find the user sequence in the tree */
	found=FALSE;
	for(i=0;i<nseqs;i++) {
		if(strcmp(seqname,lptr[i]->name)==0) {
			found=TRUE;
			s=i;
			seqptr=lptr[i];
			break;
		}
	}

	if(found==FALSE) {
		fprintf(stdout,"Error: can't find %s in tree file\n",seqname);
		exit(1);
	}

/* everything's ok - go ahead and calculate distances */
        path2root = (treeptr *)ckalloc(nseqs * sizeof(treeptr));
        dist2node = (float *)ckalloc(nseqs * sizeof(float));
        seqdist = (float *)ckalloc(nseqs * sizeof(float));
 
/* first find all the nodes between the selected sequence and the root */
	get_path2root(seqptr,path2root,dist2node,&depth);

	for(i=0;i<nseqs;i++) {
		if(i!=s) {
        		n=common_ancestor(lptr[i],path2root,depth,&dist);
			seqdist[i]=dist+dist2node[n-1];
		}
	}

        for(i=0;i<nseqs;i++) {
                if(i!=s) {
                found=FALSE;
                for(j=0;j<nseqs;j++) {
                        if(strcmp(names[j],lptr[i]->name)==0) {
                                found=TRUE;
                                seq_score[j]=(1.0-seqdist[i]);
                                break;
                        }
                }
	}
	}

}

