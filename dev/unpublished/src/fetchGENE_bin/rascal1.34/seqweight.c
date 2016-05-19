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
static sint calc_weight(IN_TREEPTR itree,sint leaf);

/*
 *   Global variables
 */
static sint debug;

sint * calc_seq_weights(sint fseq, sint nseq, IN_TREEPTR itree,Boolean no_weights)
{
  sint   i;
  sint   temp, sum, *weight;
  sint   *sweight;


  sweight = (sint *)ckalloc((nseq+1) * sizeof(sint));
/*
  If there are more than three sequences....
*/

  if ((nseq-fseq >= 2) && (itree->distance_tree == TRUE) && (no_weights==FALSE))
     {
/*
  Calculate sequence weights based on Phylip tree.
*/
      weight = (sint *)ckalloc((nseq+1) * sizeof(sint));


      for (i=fseq; i<nseq; i++)
           weight[i] = calc_weight(itree,i);

/*
  Normalise the weights, such that the sum of the weights = INT_SCALE_FACTOR
*/

         sum = 0;
         for (i=fseq; i<nseq; i++)
            sum += weight[i];

         if (sum == 0)
          {
            for (i=fseq; i<nseq; i++)
               weight[i] = 1;
            sum = i;
          }

         for (i=fseq; i<nseq; i++)
           {
              sweight[i] = (weight[i] * INT_SCALE_FACTOR) / sum;
              if (sweight[i] < 1) sweight[i] = 1;
           }

       weight=ckfree((void *)weight);

     }

   else
     {
/*
  Otherwise, use identity weights.
*/
        temp = INT_SCALE_FACTOR / nseq;
        for (i=fseq; i<nseq; i++)
           sweight[i] = temp;
     }
   return sweight;

}

static sint calc_weight(IN_TREEPTR itree,sint leaf)
{

  INODEPTR p;
  float weight = 0.0;

  p = itree->leafptr[leaf];
  while (p->parent != NULL)
    {
       weight += p->dist / p->order;
       p = p->parent;
    }

  weight *= 100.0;

  return((sint)weight);

}

