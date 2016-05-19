#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rascal.h"

#define MIN(a,b) ((a)<(b)?(a):(b))

/*
 *   Prototypes
 */
static void calc_gap_stats(SEQ *seqs, sint nseqs, sint prf_no, sint *group, sint *gaps, sint *gap_length,float *mean_gap_diff,sint prf_length,MULTGAP_OPT gap_opt,sint *seq_weight);
static void calc_p_penalties(SEQ *seqs, sint nseqs, sint prf_no, sint *group, sint n, sint *weight);
static void calc_h_penalties(SEQ *seqs, sint nseqs, sint prf_no, sint *group, sint n, sint *weight,char *hyd_residues);
static sint local_penalty(sint penalty, sint n, sint *pweight, sint *hweight);


static sint debug=0;

/*
   Calculate position-specific gap penalties for a profile 
	seqs		the sequences
	nseqs		number of sequences
	group		for each sequence, 1=in profile 1, 2=in profile 2
	gaps		the number of gaps at each position in the sequences - this is calculated here and
                        returned for use in the profile score calculation
	profile		the profile to put the gap penalties into
	dnaflag		TRUE if dna
	gap_opt		multiple alignment gap options
*/
void calc_gap_coeff(ALN mult_aln, sint prf_no, sint *group, sint *gaps, sint **profile, COMP_MATRIX matrix, sint prf_length, sint gapcoef, sint lencoef, MULTGAP_OPT gap_opt,sint *seq_weight, Boolean use_ss_motifs)
{

	sint i, j;
	sint nseqs;
	sint n;
	sint *gap_pos;
	sint *p_weight=NULL, *h_weight=NULL;
	sint *gap_length;
	float *mean_gap_diff;
	float scale;
	sint gdist;                  /* local copy of gap_dist */
	float reduced_gap = 1.0;
	Boolean nvar_pen,nhyd_pen,npref_pen; /* local copies of no_hyd_penalties, no_pref_penalties */
	SEQ *seqs;
   
	nseqs=mult_aln.nseqs;
	seqs=mult_aln.seqs;

	nvar_pen = TRUE;
       	nhyd_pen = gap_opt.no_hyd_penalties;
       	npref_pen = gap_opt.no_pref_penalties;
       	gdist = gap_opt.gap_dist;
     
	gap_length=(sint *)ckalloc((prf_length+1)*sizeof(sint));
	mean_gap_diff=(float *)ckalloc((prf_length+1)*sizeof(float));

	calc_gap_stats(seqs,nseqs,prf_no,group,gaps,gap_length,mean_gap_diff,prf_length,gap_opt,seq_weight);

	if ((!mult_aln.dnaflag) && (npref_pen == FALSE)) {
        	p_weight = (sint *) ckalloc( (prf_length+2) * sizeof (sint) );
        	calc_p_penalties(seqs, nseqs, prf_no, group, prf_length, p_weight);
     	}

	
	if ((!mult_aln.dnaflag) && (nhyd_pen == FALSE)) {
        	h_weight = (sint *) ckalloc( (prf_length+2) * sizeof (sint) );
        	calc_h_penalties(seqs, nseqs, prf_no, group, prf_length, h_weight,gap_opt.hyd_residues);
     	}

	gap_pos = (sint *) ckalloc( (prf_length+2) * sizeof (sint) );
/*
    mark the residues close to an existing gap (set gaps[i] = -ve)
*/
   if (mult_aln.dnaflag || (gdist <= 0))
     {
       for (i=0;i<prf_length;i++) gap_pos[i] = gaps[i];
     }
   else
     {
       i=0;
       while (i<prf_length)
         {
            if (gaps[i] <= 0)
              {
                 gap_pos[i] = gaps[i];
                 i++;
              }
            else 
              {
                 for (j = -gdist+1; j<0; j++)
                  {
                   if ((i+j>=0) && (i+j<prf_length) &&
                       ((gaps[i+j] == 0) || (gaps[i+j] < j))) gap_pos[i+j] = j;
                  }
                 while (gaps[i] > 0)
                    {
                       if (i>=prf_length) break;
                       gap_pos[i] = gaps[i];
                       i++;
                    }
                 for (j = 0; j<gdist; j++)
                  {
                   if (gaps[i+j] > 0) break;
                   if ((i+j>=0) && (i+j<prf_length) && 
                       ((gaps[i+j] == 0) || (gaps[i+j] < -j))) gap_pos[i+j] = -j-1;
                  }
                 i += j;
              }
         }
     }
if (debug>1)
{
fprintf(stdout,"gap open %d gap ext %d\n",(pint)gapcoef,(pint)lencoef);
fprintf(stdout,"gaps:\n");
  for(i=0;i<prf_length;i++) fprintf(stdout,"%d ", (pint)gaps[i]);
  fprintf(stdout,"\n");
fprintf(stdout,"gap_pos:\n");
  for(i=0;i<prf_length;i++) fprintf(stdout,"%d ", (pint)gap_pos[i]);
  fprintf(stdout,"\n");
}


   for (j=0;j<prf_length; j++)
     {
          
        if (gap_pos[j] <= 0)
          {
/*
    apply residue-specific and hydrophilic gap penalties.
*/
	     	if (!mult_aln.dnaflag) {
              	profile[j+1][GAPCOL] = local_penalty(gapcoef, j,
                                                   p_weight, h_weight);
              	profile[j+1][LENCOL] = lencoef;
	     	}
	     	else {
              	profile[j+1][GAPCOL] = gapcoef;
              	profile[j+1][LENCOL] = lencoef;
	     	}

/*
    increase gap penalty near to existing gaps.
*/
             if (gap_pos[j] < 0)
                {
                    profile[j+1][GAPCOL] *= 2.0+2.0*(gdist+gap_pos[j])/gdist;
                }


          }
        else
          {
/* we're in a gap .... */
/* so set the gap penalties depending on the variance of the gap length */
/*
		if(gap_length[j] == 1) {
             		profile[j+1][GAPCOL] = 0.5 * gapcoef;
             		profile[j+1][LENCOL] = lencoef;
		}
		else {
                     if (mean_gap_diff[j] < 0.001) {
                        profile[j+1][GAPCOL] = gapcoef;
                        profile[j+1][LENCOL] = lencoef;
                     }
                     else {
                        profile[j+1][GAPCOL] = gapcoef/mean_gap_diff[j];
                        profile[j+1][LENCOL] = lencoef/mean_gap_diff[j];
                     }

		}
*/
/* so set the gap penalties depending on the number of sequences that have a gap at this position */

             scale = ((float)(nseqs-(float)gaps[j])/(float)nseqs) * reduced_gap;
		if(scale<0.5) scale=0.5;
             profile[j+1][GAPCOL] = scale*gapcoef;
             profile[j+1][LENCOL] = 0.5*lencoef;
             profile[j+1][LENCOL] = lencoef;

          }

/* apply the bin matrix gap extension penalties */
	if(matrix.format==1) {
		profile[j+1][LENCOL]=0;
		n=0;
		for(i=0;i<nseqs;i++)
			if(isalpha(seqs[i].data[j])) {
				profile[j+1][LENCOL]+=matrix.gapscore[seqs[i].data[j]-'a'];
				n++;
			}
		if(n>0) profile[j+1][LENCOL]/=n;
		else profile[j+1][LENCOL]=0;

	}
/*
   make sure no penalty is zero - even for all-gap positions
*/
        if (profile[j+1][GAPCOL] <= 0) profile[j+1][GAPCOL] = 1;
        if (profile[j+1][LENCOL] <= 0) profile[j+1][LENCOL] = 0;
     }

/* set the penalties at the beginning and end of the profile */
   if(gap_opt.nendgappenalties==TRUE)
     {
        profile[0][GAPCOL] = gapcoef;
        profile[0][LENCOL] = lencoef;
     }
   else
     {
        profile[0][GAPCOL] = 0;
        profile[0][LENCOL] = 0;
     }
   if(gap_opt.cendgappenalties==TRUE)
     {
        profile[prf_length][GAPCOL] = gapcoef;
        profile[prf_length][LENCOL] = lencoef;
     }
   else
     {
        profile[prf_length][GAPCOL] = 0;
        profile[prf_length][LENCOL] = 0;
     }
if (debug>0)
{
  fprintf(stdout,"Opening penalties:\n");
  for(i=0;i<=prf_length;i++) fprintf(stdout," %d:%d ",i, (pint)profile[i][GAPCOL]);
  fprintf(stdout,"\n");
}
if (debug>0)
{
  fprintf(stdout,"Extension penalties:\n");
  for(i=0;i<=prf_length;i++) fprintf(stdout,"%d:%d ",i, (pint)profile[i][LENCOL]);
  fprintf(stdout,"\n");
}

   if ((!mult_aln.dnaflag) && (npref_pen == FALSE))
        p_weight=ckfree((void *)p_weight);

   if ((!mult_aln.dnaflag) && (nhyd_pen == FALSE))
        h_weight=ckfree((void *)h_weight);

   gap_length=ckfree((void *)gap_length);
   mean_gap_diff=ckfree((void *)mean_gap_diff);

   gap_pos=ckfree((void *)gap_pos);

}              

/*
   Calculate residue-specific penalties according to relative frequencies of
   residues found next to gaps (Pascarelli and Argos).
	aln	sequences
	n	length of sequences
	fs	first sequence in aln
	ls	last sequence in aln
	weight	array to hold penalties calculated here.
*/
static void calc_p_penalties(SEQ *seqs, sint nseqs, sint prf_no, sint *group, sint n, sint *weight)
{
	char ix;
	sint j,k;
	sint i,ns;

/* pascarella probabilities for opening a gap at specific residues */
	char   pr[] =     {'a' , 'c', 'd', 'e', 'f', 'g', 'h', 'k', 'i', 'l',
                   'm' , 'n', 'p', 'q', 'r', 's', 't', 'v', 'y', 'w'};
	sint    pas_op[] = { 87, 87,104, 69, 80,139,100,104, 68, 79,
                    71,137,126, 93,128,124,111, 75,100, 77};

	ns = 0;
	for (i=0;i<n;i++) 
      		weight[i] = 0;

	for (k=0;k<nseqs;k++) {
		if(group[k]==prf_no) {
			ns++;
			for (i=0;i<n;i++) {
           			for (j=0;j<22;j++) {
                			ix = seqs[k].data[i];
                			if (!isalpha(ix)) continue;
                			if (ix == pr[j]) {
                    				weight[i] += (180-pas_op[j]);
                    				break;
                  			}
             			}
        		}
		}
    	}
	for (i=0;i<n;i++) 
      		weight[i] /= (float)ns;

}
            
/*
   Calculate penalties based on runs of hydrophylic residues.
	aln	sequences
	n	length of sequences
	fs	first sequence in aln
	ls	last sequence in aln
	weight	array to hold penalties calculated here
	hyd_residues	residues considered to be hydrophylic.
*/
static void calc_h_penalties(SEQ *seqs, sint nseqs, sint prf_no, sint *group, sint n, sint *weight,char *hyd_residues)
{

/*
   weight[] is the length of the hydrophilic run of residues.
*/
	char ix;
	sint nh,j,k,ns;
	sint i,e,s;
	sint *hyd;

	hyd = (sint *)ckalloc((n+2) * sizeof(sint));
	nh = (sint)strlen(hyd_residues);
	for (i=0;i<n;i++)
		weight[i] = 0;
	ns=0;
	for (k=0;k<nseqs;k++) {
		if(group[k]==prf_no) {
			ns++;
       			for (i=0;i<n;i++) {
             			hyd[i] = 0;
             			for (j=0;j<nh;j++) {
                   			ix = seqs[k].data[i];
                   			if (!isalpha(ix)) continue;
                   			if (ix == tolower(hyd_residues[j])) {
                         			hyd[i] = 1;
                         			break;
                      			}
                		}
          		}
       			i = 0;
       			while (i < n) {
            			if (hyd[i] == 0) i++;
            			else {
                 			s = i;
                 			while ((hyd[i] != 0) && (i<n)) i++;
                 			e = i;
                 			if (e-s > 3)
                    				for (j=s; j<e; j++) weight[j] += 100;
              			}
         		}
    		}
	}

	for (i=0;i<n;i++)
		weight[i] /= (float)ns;

	hyd=ckfree((void *)hyd);

if (debug>1)
{
  for(i=0;i<n;i++) fprintf(stdout,"%d ", (pint)weight[i]);
  fprintf(stdout,"\n");
}

}
            
static sint local_penalty(sint penalty, sint n, sint *pweight, sint *hweight)
{

	Boolean h = FALSE;
	float gw;

	gw = 1.0;

	if (hweight != NULL) {
        	if (hweight[n] > 0) {
           		gw *= 0.5;
           		h = TRUE;
         	}
    	}
	if ((pweight != NULL) && (h==FALSE)) {
       		gw *= ((float)pweight[n]/100.0);
    	}

	gw *= penalty;
	return((sint)gw);

}

static void calc_gap_stats(SEQ *seqs, sint nseqs, sint prf_no, sint *group, sint *gaps, sint *gap_length,float *mean_gap_diff,sint prf_length,MULTGAP_OPT gap_opt,sint *seq_weight)
{
  
  int i, j, k, c;
  int is, ie, si;  
   int *insert_length;
  float mu,variance,wt_mu,wt_sum,diff;
 
  insert_length = (int *)ckalloc((nseqs+1) * sizeof(int));
   	
/* calculate sum and mean of unnormalised sequence weights */
   wt_sum = 0;
   for (i=0; i<nseqs; i++)
        wt_sum += seq_weight[i];
   if(wt_sum<1.0) wt_sum=1.0;
   wt_mu = wt_sum / nseqs;
          
/*
   calculate the number of gaps at each position
*/    
	for (j=0; j<prf_length; j++) gaps[j] = 0;
	for (i=0; i<nseqs; i++) {
		if(group[i]==prf_no) {
/*
   Include end gaps as gaps ?
*/
        		is = 0;
        		ie = prf_length;
        		if (gap_opt.use_endgaps == FALSE && gap_opt.nendgappenalties==FALSE) {
          			for (j=0; j<prf_length; j++) {
              				c = seqs[i].data[j];
              				if (!isalpha(c))
                 				is++;
              				else
                 				break;
            			}
			}
        		if (gap_opt.use_endgaps == FALSE && gap_opt.cendgappenalties==FALSE) {
          			for (j=prf_length-1; j>=0; j--) {
              				c = seqs[i].data[j];
              				if (!isalpha(c))
                 				ie--;
              				else
                 				break;
            			}
        		}

        		for (j=is; j<ie && j<seqs[i].len; j++) {
              			if (!isalpha(seqs[i].data[j])) gaps[j]++;
          		}
     		}
	}

/*
  calculate the length of the gap at each position and
  the weighted mean and variance of the insert lengths for each sequence.
*/
   j = 0;
   while (j < prf_length)
     {
        if (gaps[j] <= 0)
          {
              gap_length[j] = 0;
              mean_gap_diff[j] = 0.0;
              j++;
          }
        else
          {
             is = j;
             while ((gaps[j] != 0) && (j < prf_length)) j++;
             ie = j;
             for (i=is; i<ie; i++) gap_length[i] = ie - is;
             for (si=0;si<nseqs;si++)
               {
                  insert_length[si] = 0;
                  for (k=is;k<ie && k<seqs[si].len;k++) if (!isalpha(seqs[si].data[k]))
                  	 insert_length[si]++;
               }
             mu = 0.0;
             for (si=0;si<nseqs;si++)
               {
                  mu += seq_weight[si] * insert_length[si];
               }
             mu /= wt_sum;
             
             variance = 0.0;
             for (si=0;si<nseqs;si++)
               {
                  diff = insert_length[si] - mu;
                  variance += seq_weight[si] * diff * diff;
               }
             if (wt_sum - wt_mu <= 0.001) variance = 0;
             else
                variance = sqrt(variance/(wt_sum - wt_mu));
             for (i=is; i<ie; i++) mean_gap_diff[i] = variance;
          }
     }

	ckfree(insert_length);

}



