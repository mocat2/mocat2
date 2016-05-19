#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "rascal.h"

/*
 *   Prototypes
 */
/*
 *   Global variables
 */

void calc_prf2(sint **profile, SEQ *seqs, sint nseqs, sint prf_no, sint *group, sint prf_length, sint *seqweight)
{

	sint sum1, sum2;	
	sint i, d;
	sint   r;

	for (r=0; r<prf_length; r++) {
/*
   calculate sum2 = number of residues found in this column
*/
       		sum2 = 0;
       		for (i=0; i<nseqs; i++) {
			if (group[i]==prf_no) {
            			sum2 += seqweight[i];
         		}
		}
/*
   only include matrix comparison scores for those residue types found in this
   column
*/
		if (sum2 == 0) {
           		for (d=0; d<NUMRES; d++)
             			profile[r+1][d] = 0;
         	}
       		else {
           		for (d=0; d<NUMRES; d++) {
                		sum1 = 0;
       				for (i=0; i<nseqs; i++) 
					if (group[i]==prf_no) 
                  				if (d == seqs[i].data[r]-'a') sum1 += seqweight[i];
                		profile[r+1][d] = (sint)(10 * (float)sum1 / (float)sum2);
             		}
         	}
    	}
}


