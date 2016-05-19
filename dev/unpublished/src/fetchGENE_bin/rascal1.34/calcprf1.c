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
void calc_prf1(sint **profile, sint **freq, SEQ *seqs, sint nseqs, sint prf_no, sint *group, 
		sint *gaps, COMP_MATRIX matrix, sint prf_length, sint *seqweight,Boolean norm_gaps)
{

	sint d, i, res; 
	sint r, pos,ns;
	int f;
	float scale,sum2;

norm_gaps=FALSE;

	sum2 = ns = 0;
        for (i=0; i<nseqs; i++)
		if(group[i]==prf_no) {
			ns++;
       			sum2 += seqweight[i];
		}
	if(sum2<1.0) sum2=1.0;

  	for (r=0; r<prf_length; r++) {
      		for (d=0; d<NUMRES; d++)
            		freq[r][d] = 0;
      		for (i=0; i<nseqs; i++)
			if(group[i]==prf_no && r<seqs[i].len && isalpha(seqs[i].data[r]))
				freq[r][seqs[i].data[r]-'a'] += seqweight[i];
   	}

	for (pos=0; pos< prf_length; pos++) {
      		if (gaps[pos] == ns) {
			for (res=0; res<NUMRES; res++) {
                		profile[pos+1][res] = 0;
             		}
        	}
      		else {
			if(norm_gaps)
           			scale = (float)(ns-gaps[pos]) / (float)ns;
			else
				scale=1.0;
           		for (res=0; res<NUMRES; res++) {
                		f = 0;
                		for (d=0; d<NUMRES; d++)
                     			f += (freq[pos][d] * matrix.score[d][res]);
                		profile[pos+1][res] = (sint  )(((float)f / sum2)*scale);
             		}
        	}
    	}

}


