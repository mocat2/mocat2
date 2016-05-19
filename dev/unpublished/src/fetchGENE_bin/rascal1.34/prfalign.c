#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "rascal.h"

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

/*
 *   Prototypes
 */
static lint 	pdiff(sint A,sint B,sint i,sint j,sint go1,sint go2);
static lint	prfprfscore(sint n, sint m);
static lint 	prfseqscore(sint n, sint m);
static lint 	seqseqscore(sint n, sint m);
static sint 	gap_penalty1(sint i, sint j,sint k);
static sint 	open_penalty1(sint i, sint j);
static sint 	ext_penalty1(sint i, sint j);
static sint 	gap_penalty2(sint i, sint j,sint k);
static sint 	open_penalty2(sint i, sint j);
static sint 	ext_penalty2(sint i, sint j);
static void 	padd(sint k);
static void 	pdel(sint k);
static void 	palign(void);
static sint 	ptracepath(char *path1, char *path2);
static void add_ggaps(SEQ *seqs, sint nseqs, sint *group, int len, char *path1, char *path2);
static lint prfprfscoremotif(sint n, sint m);
static lint seqseqscoremotif(sint n, sint m);
static lint prfseqscoremotif(sint n, sint m);
static sint start_pos(char *seq,sint p);
static sint next_pos(char *seq,sint len,sint p);

/*
 *   Global variables
 */
static sint    	debug=0;

/* isprf1, isprf2 are set TRUE if we're aligning a profile. If we're
   aligning a single sequence they are set FALSE. They are used in
   the Myers and Miller iterative alignment routine to decide which
   scoring function to call. If we're aligning single sequences,
   seq1 and seq2 indicate where in SEQ *seqs the sequences are. */
static Boolean    isprf1,isprf2;
static sint       seq1,seq2;
static lint (*scorefunction)(sint n, sint m);

/* these are defined here to avoid passing them everywhere */
/* they are used by the alignment prfprf, seqseq, prfseq score functions */
static Boolean nendgappenalties;
static Boolean cendgappenalties;
static sint    	**profile1, **profile2;
static sint    	**freq1, **freq2;
static sint    	prf_length1, prf_length2;
static COMP_MATRIX matrix;
static SEQ *seqs;
static sint **motif_score;
static sint    nseqs1, nseqs2;

/* variables used by the iterative Myers and Miller alignment algorithm */
static sint 	print_ptr,last_print;
static sint 	*displ;
static sint 	*gS;
static lint 	*HH, *DD, *RR, *SS;

/* hard-coded comparison matrices, in matrices.h */
extern sint  	blosum30mt[], blosum40mt[], blosum45mt[];
extern sint  	blosum62mt2[], blosum80mt[];
extern sint  	pam20mt[], pam60mt[];
extern sint  	pam120mt[], pam160mt[], pam350mt[];
extern sint  	gon40mt[], gon80mt[];
extern sint    gon120mt[], gon160mt[], gon250mt[], gon350mt[];
extern sint  	clustalvdnamt[],swgapdnamt[];
extern sint  	idmat[];
extern sint	bin50mt[],bin40mt[],bin30mt[],bin20mt[],bin10mt[];
extern sint	bin50gap[],bin40gap[],bin30gap[],bin20gap[],bin10gap[];

/*
   Do an alignment of profile/sequence versus profile/sequence. 
	aln_type	0 = multiple alignment, 1 = profile vs profile, 2 = sequence vs profile
	mult_aln	the multiple alignment (some sequences are still aligned!)
	group		for each sequence group[i]=1 if sequence in first profile
					  group[i]=2 if sequence in second profile
					  group[i]=0 means don't align sequence this time
	tmat		the pairwise distance matrix
	use_maxid	use max identity (=100-mindist) for sequence weighting
	mult_opt	multiple alignment options

   Returns alignment score
*/

sint *sweight;

lint prfalign(sint aln_type,ALNPTR mult_aln,sint *group,double **tmat,Boolean use_maxid,MULT_OPT mult_opt,ALNCOUNT *alncount,Boolean norm_gaps)
{

	static Boolean negative;
	static sint    i, j, k, n, count = 0;
	static sint    len1, len2, is, minlen;
	static sint   se1, se2, sb1, sb2;
	static sint  maxres;
	static sint int_scale;
	static sint nlen1,nlen2;
	static sint alen1,alen2;
	static sint clen1,clen2;
	static lint    score;
	static float  scale;
	static double logmin,logdiff;
	static double pcid;
	static sint *matptr,*gapptr;
	static char   	*aln_path1, *aln_path2;
	static sint    	alignment_len;
	static sint    *gaps;
	static sint    sum1,sum2;
	static sint *ulen;
	static sint    gapcoef1,gapcoef2;
	static sint    lencoef1,lencoef2;
	static sint max_aln_length;
	static float mindist;

	nseqs1 = nseqs2 = 0;
	for (i=0;i<mult_aln->nseqs;i++) {
		if (group[i] == 1) nseqs1++;
		else if (group[i] == 2) nseqs2++;
	}

	if ((nseqs1 == 0) || (nseqs2 == 0)) 
		return(0.0);

/* if necessary, swap profile1 for profile2, so that the one with the
   most sequences is profile1 */

	/*if (aln_type==0 && nseqs2 > nseqs1) {*/
	if (nseqs2 > nseqs1) {
		for (i=0;i<mult_aln->nseqs;i++) {
			if (group[i] == 1) group[i] = 2;
			else if (group[i] == 2) group[i] = 1;
		}
		i=nseqs1;
		nseqs1=nseqs2;
		nseqs2=i;
	}

/*
  Calculate number of sequences and length for each profile 
  At the same time calculate the average length of profiles, removing gaps.
*/
	prf_length1 = prf_length2 = 0;
	len1=len2=0;
	ulen=(sint *)ckalloc(mult_aln->nseqs*sizeof(sint));
	for (i=0;i<mult_aln->nseqs;i++) {
		if (group[i] == 1) {
			if (mult_aln->seqs[i].len > prf_length1) prf_length1 = mult_aln->seqs[i].len;
			is=0;
			for (j=0; j<mult_aln->seqs[i].len; j++) 
				if(isalpha(mult_aln->seqs[i].data[j])) is++;
			ulen[i]=is;
			if(is>0) {
				if(ulen[i]>len1) len1=ulen[i];
			}
			seq1=i;
		}
		else if (group[i] == 2) {
			if (mult_aln->seqs[i].len > prf_length2) prf_length2 = mult_aln->seqs[i].len;
			is=0;
			for (j=0; j<mult_aln->seqs[i].len; j++) 
				if(isalpha(mult_aln->seqs[i].data[j])) is++;
			ulen[i]=is;
			if(is>0) {
				if(ulen[i]>len2) len2=ulen[i];
			}
			sum2+=mult_aln->seqs[i].weight;
			seq2=i;
		}
	}
	if(len1<=0 || len2<=0) return 0;
	minlen = MIN(len1,len2);

	sweight=(sint *)ckalloc(mult_aln->nseqs*sizeof(sint));
	sum1 = sum2 = 0;
use_maxid=TRUE;
	if(use_maxid) {
		for (i=0;i<mult_aln->nseqs;i++) {
			if(group[i]>0) {
                		mindist = 100.0;
                		for (j=0;j<mult_aln->nseqs;j++)
                        		if(i!=j && group[j]>0 && group[i] != group[j] && mindist>tmat[i][j]) mindist = tmat[i][j];
				sweight[i]=(100.0-100.0*mindist);
				if(group[i]==1 && ulen[i]>0) sum1+=sweight[i];
				else if(group[i]==2 && ulen[i]>0) sum2+=sweight[i];
				else sweight[i]=0;
			}
		}
	} else {
		for (i=0;i<mult_aln->nseqs;i++) {
			if(group[i]>0) {
				sweight[i]=mult_aln->seqs[i].weight;
				if(group[i]==1 && ulen[i]>0) sum1+=sweight[i];
				else if(group[i]==2 && ulen[i]>0) sum2+=sweight[i];
				else sweight[i]=0;
			}
		}
	}
	ckfree(ulen);

	if(sum1<1) sum1=1;
	if(sum2<1) sum2=1;
	for (i=0;i<mult_aln->nseqs;i++) {
		if(group[i]==1) sweight[i]=(sweight[i] * INT_SCALE_FACTOR)/sum1;
		if(group[i]==2) sweight[i]=(sweight[i] * INT_SCALE_FACTOR)/sum2;
	}
	seqs=mult_aln->seqs;

/* decide if we're aligning sequences or profiles */
	if(nseqs1==1) 
		isprf1=FALSE;
	else isprf1=TRUE;
	if(nseqs2==1) 
		isprf2=FALSE;
	else isprf2=TRUE;

/* for rascal, check for sequence fragments in the block being realigned */
	if(alncount!=NULL) {
		nlen1=nlen2=0;
		alen1=alen2=0;
		clen1=clen2=0;
                mult_opt.gap_opt->nendgappenalties=TRUE;
                mult_opt.gap_opt->cendgappenalties=TRUE;
		for (i=0;i<mult_aln->nseqs;i++) {
			if (group[i] == 1) {
				if(alncount[i].n>nlen1) nlen1=alncount[i].n;
				if(alncount[i].a>alen1) alen1=alncount[i].a;
				if(alncount[i].c>clen1) clen1=alncount[i].c;
			}
			else if (group[i] == 2) {
				if(alncount[i].n>nlen2) nlen2=alncount[i].n;
				if(alncount[i].a>alen2) alen2=alncount[i].a;
				if(alncount[i].c>clen2) clen2=alncount[i].c;
			}
		}
                if((nlen1==0 && alen1>0) || 
			(nlen2==0 && alen2>0)) {
                        mult_opt.gap_opt->nendgappenalties=FALSE;
                }
                if((clen1==0 && alen1>0) || 
			(clen2==0 && alen2>0)) {
                        mult_opt.gap_opt->cendgappenalties=FALSE;
                }
	}

if (debug>0) {
fprintf(stdout,"\nsequences profile 1:\n");
for (i=0;i<mult_aln->nseqs;i++) 
if (group[i] == 1 && sweight[i]>0) 
fprintf(stdout,"%30s %s %d %d\n",mult_aln->seqs[i].name,mult_aln->seqs[i].data,mult_aln->seqs[i].len,sweight[i]);
fprintf(stdout,"sequences profile 2:\n");
for (i=0;i<mult_aln->nseqs;i++) 
if (group[i] == 2 && sweight[i]>0) 
fprintf(stdout,"%30s %s %d %d\n",mult_aln->seqs[i].name,mult_aln->seqs[i].data,mult_aln->seqs[i].len,sweight[i]);
}

/*
   used to allocate enough memory for temporary arrays, until we know how
   long the new alignment will be.
*/
	max_aln_length = prf_length1 + prf_length2+2;
  
   
/*
   calculate the mean of the sequence pc identities between the two groups
*/
        count = 0;
        pcid = 0.0;
	negative=mult_opt.neg_matrix;
        for (i=0;i<mult_aln->nseqs;i++) {
             if (group[i] == 1)
             for (j=0;j<mult_aln->nseqs;j++)
               if (group[j] == 2) {
                       count++;
                       pcid += 1.0-tmat[i][j];
               }
        }

	pcid = (100.0*pcid)/(float)count;

if (debug > 0) fprintf(stdout,"mean tmat %3.1f\n", pcid);

/* 
   get a comparison matrix depending on the average percent identity 
*/

	matrix.format=0;

	int_scale = 100;
	if (mult_aln->dnaflag) {
		scale=1.0;
		if(mult_opt.dnamtrxname[0]=='\0')
			maxres = get_user_matrix(mult_opt.dnausermtrxname, mult_opt.neg_matrix, int_scale, &matrix);
        	else {
			if (strcmp(mult_opt.dnamtrxname, "iub") == 0) 
	    			matptr=swgapdnamt;
       			else if (strcmp(mult_opt.dnamtrxname, "clustalw") == 0) {
            			scale=0.66;
            			matptr=clustalvdnamt;
			}
            		maxres = get_cl_matrix(mult_aln->dnaflag,matptr, gapptr, mult_opt.neg_matrix, int_scale, &matrix);
        	}
            	if (maxres == 0) 
			return((sint)-1);
/* transition weights for A-G, C-T, C-U */
            	matrix.score[0][6]=mult_opt.transition_weight*matrix.score[0][0];
            	matrix.score[6][0]=mult_opt.transition_weight*matrix.score[0][0];
            	matrix.score[2][19]=mult_opt.transition_weight*matrix.score[0][0];
            	matrix.score[19][2]=mult_opt.transition_weight*matrix.score[0][0];
            	matrix.score[2][20]=mult_opt.transition_weight*matrix.score[0][0];
            	matrix.score[20][2]=mult_opt.transition_weight*matrix.score[0][0];

          	gapcoef1 = gapcoef2 = 100.0 * mult_opt.dna_gap_open *scale;
          	lencoef1 = lencoef2 = 100.0 * mult_opt.dna_gap_extend *scale;
    	}
  	else {
/* for proteins:- */

       		scale=0.75;
		if (mult_opt.mtrxname[0]=='\0')
      	   		maxres = get_user_matrix_series(mult_opt.usermtrxname, pcid, negative, int_scale, &matrix);
		else {
       			if (strcmp(mult_opt.mtrxname, "blosum") == 0) {
           			scale=0.75;
           			if (pcid > 80.0)
					matptr=blosum80mt;
           			else if (pcid > 60.0)
					matptr=blosum62mt2;
           			else if (pcid > 40.0)
					matptr=blosum45mt;
           			else if (pcid > 30.0) {
                			scale=0.5;
					matptr=blosum45mt;
             			}
           			else if (pcid > 20.0) {
                			scale=0.6;
					matptr=blosum45mt;
             			}
           			else {
                			scale=0.6;
					matptr=blosum30mt;
             			}
       			}
       			else if (strcmp(mult_opt.mtrxname, "pam") == 0) {
           			scale=0.75;
           			if (pcid > 80.0) 
					matptr=pam20mt;
           			else if (pcid > 60.0) 
					matptr=pam60mt;
           			else if (pcid > 40.0)
					matptr=pam120mt;
           			else 
					matptr=pam350mt;
       			}
       			else if (strcmp(mult_opt.mtrxname, "gonnet") == 0) {
           			int_scale /= 10;
           			if (pcid > 50.0) {
                			matptr=gon80mt;
	   				scale/=2.5;
             			}
           			else if (pcid > 30.0) {
                			matptr=gon120mt;
	   				scale/=2.5;
             			}
           			else if (pcid > 25.0) {
                			if(minlen<100) matptr=gon250mt;
      	   				else matptr=gon160mt;
	   			scale/=1.7;
             			}
           			else {
                			if(minlen<100) matptr=gon350mt;
      	   				else matptr=gon250mt;
	   			scale/=1.7;
             			}
       			}
       			else if (strcmp(mult_opt.mtrxname, "bin") == 0) {
           			int_scale /= 100;
				matrix.format=1;
       				if (pcid > 45.0) {
                			matptr=bin50mt;
                			gapptr=bin50gap;
				}
       				else if (pcid > 35.0) {
                			matptr=bin40mt;
                			gapptr=bin40gap;
				}
       				else if (pcid > 25.0) {
                			matptr=bin30mt;
                			gapptr=bin30gap;
				}
				else {
                			matptr=bin20mt;
                			gapptr=bin20gap;
       				}
			}
       			else if (strcmp(mult_opt.mtrxname, "id") == 0)
      	   			matptr=idmat;
   			maxres = get_cl_matrix(mult_aln->dnaflag,matptr, gapptr, negative, int_scale, &matrix);
		}
      		if (maxres == 0) {
           		fprintf(stdout,"Error: matrix %s not found\n", mult_opt.mtrxname);
           		return(-1);
        	}
/* 
   now we've chosen a matrix, choose gap opening and extension penalties to go with it.
*/
		if(matrix.format==1) {
       			if(pcid>0.25) scale=0.2;
			else scale=0.18;
  			if(len1==0 || len2==0) {
  				logmin=1.0;
  				logdiff=1.0;
  			}  
  			else {
				if (minlen<=10) logmin=1.0;
				else if (minlen<300) logmin = 1.0/log10((double)(minlen*0.8));
				else logmin = 1.0/log10((double)(minlen*0.25));
 				if (len2<len1)
    	 				logdiff = 1.0+2.0*log10((double)((float)len2/(float)len1));
  				else if (len1<len2)
  	   				logdiff = 1.0+2.0*log10((double)((float)len1/(float)len2));
  				else logdiff=1.0;
  			}

			if (logdiff<0.5) {
				logdiff=0.5;
			}
if(debug>0)
fprintf(stdout,"%d %d logmin %f   logdiff %f\n",
(pint)len1,(pint)len2, logmin,logdiff);
			gapcoef1 = gapcoef2 = scale * matrix.mat_avscore * mult_opt.prot_gap_open/(logdiff*logmin);
			lencoef1 = lencoef2 = 0.0;
		}
		else {
			if(len1==0 || len2==0) {
				logmin=1.0;
				logdiff=1.0;
			}  
			else {
				if (minlen<=10) logmin=1.0;
				else if (minlen<100) logmin = 1.0/log10((double)(minlen));
				else if (minlen<300) logmin = 1.0/log10((double)(minlen*0.5));
				else logmin = 1.0/log10((double)(minlen*0.25));
				if(MAX(len1,len2)<15) 
					logdiff=1.0;
				else if (len2<len1)
					logdiff = 1.0+4.0*log10((double)((float)len2/(float)len1));
				else if (len1<len2)
					logdiff = 1.0+4.0*log10((double)((float)len1/(float)len2));
				else logdiff=1.0;
			}
			if (logdiff<0.5) {
				logdiff=0.5;
			}
			if (logmin<0.5) {
				logmin=0.5;
			}
if(debug>0) fprintf(stdout,"%d %d logmin %f   logdiff %f\n",
	(pint)len1,(pint)len2, logmin,logdiff);
			if (negative) {
				gapcoef1 = gapcoef2 = 50.0 * (float)(mult_opt.prot_gap_open/logmin);
				lencoef1 = lencoef2 = 500.0 * mult_opt.prot_gap_extend;
			}
			else {
				gapcoef1 = gapcoef2 = scale * matrix.mat_avscore * mult_opt.prot_gap_open/(logdiff*logmin);
				lencoef1 = lencoef2 = 2.0 * scale * matrix.mat_avscore * mult_opt.prot_gap_extend/logdiff;
			}
		}
if(debug>1) fprintf(stdout,"\npcid %3.1f scale %3.1f\n",pcid,scale);
if (debug>0) fprintf(stdout,"matavscore %d gapopen %d gapext %d\n",matrix.mat_avscore,gapcoef1,lencoef1);
	}

/*
   calculate the local motif scores
*/
        motif_score=NULL;
        n=0;
        if(mult_aln->motifs != NULL) {
                motif_score = (sint **) ckalloc( (prf_length1+1) * sizeof (sint *) );
                for(i=0; i<prf_length1; i++)
                        motif_score[i] = (sint *) ckalloc( (prf_length2+2) * sizeof(sint) );
                for (i=0; i<prf_length1; i++)
                        for (j=0; j<prf_length2; j++)
                                motif_score[i][j]=0;

        	n=calc_motif_scores(mult_aln->seqs,mult_aln->nseqs,group,mult_aln->motifs,motif_score,prf_length1,prf_length2);
        }


if (debug>1)
{
fprintf(stdout,"Gap Open1 %d  Gap Open2 %d  Gap Extend1 %d   Gap Extend2 %d\n",
   (pint)gapcoef1,(pint)gapcoef2, (pint)lencoef1,(pint)lencoef2);
fprintf(stdout,"Matrix  %s\n", mult_opt.mtrxname);
}

	profile1 = (sint **) ckalloc( (prf_length1+2) * sizeof (sint *) );
	for(i=0; i<prf_length1+2; i++)
		profile1[i] = (sint *) ckalloc( (LENCOL+2) * sizeof(sint) );
	freq1 = (sint **) ckalloc( (prf_length1+2) * sizeof (sint *) );
	for(i=0; i<prf_length1+2; i++)
		freq1[i] = (sint *) ckalloc( (LENCOL+2) * sizeof(sint) );

	profile2 = (sint **) ckalloc( (prf_length2+2) * sizeof (sint *) );
	for(i=0; i<prf_length2+2; i++)
		profile2[i] = (sint *) ckalloc( (LENCOL+2) * sizeof(sint) );
	freq2 = (sint **) ckalloc( (prf_length2+2) * sizeof (sint *) );
	for(i=0; i<prf_length2+2; i++)
		freq2[i] = (sint *) ckalloc( (LENCOL+2) * sizeof(sint) );

/*
  calculate position-specific gap penalties for the first profile.
*/
	gaps = (sint *) ckalloc( (max_aln_length+10) * sizeof (sint) );

       	calc_gap_coeff(*mult_aln, 1, group, gaps, profile1, matrix, prf_length1, gapcoef1, lencoef1, *mult_opt.gap_opt,sweight,mult_opt.use_ss_motifs);

/*
  calculate the first profile.
*/
	if(isprf1) {
		calc_prf1(profile1, freq1, mult_aln->seqs, mult_aln->nseqs, 1, group, gaps, matrix, prf_length1, sweight,norm_gaps);
if (debug>3)
{
  for (j=0;j<NUMRES;j++)
    fprintf(stdout,"%c    ", (char)j);
 fprintf(stdout,"\n");
  for (i=0;i<prf_length1;i++)
   {
    fprintf(stdout,"%d: \n",(pint)(i+1));
    for (j=0;j<NUMRES;j++)
      fprintf(stdout,"%d ", (pint)profile1[i+1][j]);
    fprintf(stdout,"%d %d\n",(pint)profile1[i+1][GAPCOL],(pint)profile1[i+1][LENCOL]);
   }
}
	}

/*
  calculate the position-specific gap penalties for the second profile.
*/

       	calc_gap_coeff(*mult_aln, 2, group, gaps, profile2, matrix, prf_length2, gapcoef2, lencoef2, *mult_opt.gap_opt,sweight,mult_opt.use_ss_motifs);

/*
  calculate the second profile.
*/

	if(isprf2) {
		calc_prf2(profile2, mult_aln->seqs, mult_aln->nseqs, 2, group, prf_length2, sweight);
		/*calc_prf1(profile2, freq2, mult_aln->seqs, mult_aln->nseqs, 2, group, gaps, matrix, prf_length2, sweight);*/
if (debug>2)
{
  for (j=0;j<NUMRES;j++)
    fprintf(stdout,"%c    ", (char)j);
 fprintf(stdout,"\n");
  for (i=0;i<prf_length2;i++)
   {
    fprintf(stdout,"%d: \n",(pint)(i+1));
    for (j=0;j<NUMRES;j++)
      fprintf(stdout,"%d ", (pint)profile2[i+1][j]);
    fprintf(stdout,"%d %d\n",(pint)profile2[i+1][GAPCOL],(pint)profile2[i+1][LENCOL]);
   }
}
	}

	gaps=ckfree((void *)gaps);

/* use Myers and Miller to align the two profiles/sequences */

	last_print = 0;
	print_ptr = 1;

	sb1 = sb2 = 0;
	se1 = prf_length1;
	se2 = prf_length2;

	HH = (lint *) ckalloc( (max_aln_length+1) * sizeof (lint) );
	DD = (lint *) ckalloc( (max_aln_length+1) * sizeof (lint) );
	RR = (lint *) ckalloc( (max_aln_length+1) * sizeof (lint) );
	SS = (lint *) ckalloc( (max_aln_length+1) * sizeof (lint) );
	gS = (sint *) ckalloc( (max_aln_length+1) * sizeof (sint) );
	displ = (sint *) ckalloc( (max_aln_length+1) * sizeof (sint) );

	nendgappenalties=mult_opt.gap_opt->nendgappenalties;
	cendgappenalties=mult_opt.gap_opt->cendgappenalties;

        if(motif_score!=NULL) {
                if(isprf1==FALSE) scorefunction=seqseqscoremotif;
                else if(isprf2==FALSE) scorefunction=prfseqscoremotif;
                else scorefunction=prfprfscoremotif;
        }
        else {
                if(isprf1==FALSE) scorefunction=seqseqscore;
                else if(isprf2==FALSE) scorefunction=prfseqscore;
                else scorefunction=prfprfscore;
        }


	score = pdiff(sb1, sb2, se1-sb1, se2-sb2, profile1[0][GAPCOL], profile1[prf_length1][GAPCOL]);

	HH=ckfree((void *)HH);
	DD=ckfree((void *)DD);
	RR=ckfree((void *)RR);
	SS=ckfree((void *)SS);
	gS=ckfree((void *)gS);

	aln_path1 = (char *) ckalloc( (max_aln_length+1) * sizeof(char) );
	aln_path2 = (char *) ckalloc( (max_aln_length+1) * sizeof(char) );

	alignment_len = ptracepath(aln_path1,aln_path2);
  
	displ=ckfree((void *)displ);

/*
    add the new gaps into the two profiles
*/
	add_ggaps(mult_aln->seqs, mult_aln->nseqs,group,alignment_len,aln_path1,aln_path2);

	for (i=0;i<prf_length1+2;i++)
		profile1[i]=ckfree((void *)profile1[i]);
	profile1=ckfree((void *)profile1);
	for (i=0;i<prf_length1+2;i++)
		freq1[i]=ckfree((void *)freq1[i]);
	freq1=ckfree((void *)freq1);

	for (i=0;i<prf_length2+2;i++)
		profile2[i]=ckfree((void *)profile2[i]);
	profile2=ckfree((void *)profile2);
	for (i=0;i<prf_length2+2;i++)
		freq2[i]=ckfree((void *)freq2[i]);
	freq2=ckfree((void *)freq2);

        if(motif_score!=NULL) {
                for (i=0;i<prf_length1;i++)
                        motif_score[i]=ckfree((void *)motif_score[i]);
                motif_score=ckfree((void *)motif_score);
        }

	aln_path1=ckfree((void *)aln_path1);
	aln_path2=ckfree((void *)aln_path2);

	sweight=ckfree((void *)sweight);

	return(score/100);
}

/*
    add new gaps into the two profiles
	seqs	the sequences
	nseqs	number of sequences
	group		for each sequence group[i]=1 if sequence in first profile
					  group[i]=2 if sequence in second profile
					  group[i]=0 otherwise
	len	the new alignment length
	aln_path1[i]=2	residue in profile1 aligned with residue in profile2
		    =1  gap in profile1 
	aln_path2[i]=2	residue in profile2 algned with residue in profile1
		    =1  gap in profile2
*/
static void add_ggaps(SEQ *seqs, sint nseqs, sint *group, sint len, char *aln_path1, char *aln_path2)
{
	sint j;
	sint i,ix;
	char *ta;

	ta = (char *) ckalloc( (len+1) * sizeof (char) );

	for (j=0;j<nseqs;j++) {
		if(group[j]==1) {
      			ix = 0;
      			for (i=0;i<len;i++) {
           			if (aln_path1[i] == 2) {
                 			if (ix < seqs[j].len) {
                    				ta[i] = seqs[j].data[ix];
					}
                 			else {
                    				ta[i] = -3;
					}
                 			ix++;
              			}
           			else if (aln_path1[i] == 1) {
/*
   insertion in first alignment...
*/
                 			ta[i] = GAP1;
              			}
           			else {
                 			fprintf(stdout,"Error in aln_path\n");
              			}
         		}
       			ta[i] = -3;
       
       			seqs[j].data = (char *)realloc(seqs[j].data, (len+2) * sizeof (char));
       			for (i=0;i<len;i++) {
         			seqs[j].data[i] = ta[i];
			}
       			seqs[j].data[i] = -3;
       			seqs[j].len = len;
      		}
	}
	for (j=0;j<nseqs;j++) {
		if (group[j]==2) {
      			ix = 0;
      			for (i=0;i<len;i++) {
           			if (aln_path2[i] == 2) {
                 			if (ix < seqs[j].len) {
                    				ta[i] = seqs[j].data[ix];
					}
                 			else {
                    				ta[i] = -3;
					}
                 			ix++;
              			}
           			else if (aln_path2[i] == 1) {
/*
   insertion in second alignment...
*/
                 			ta[i] = GAP1;
              			}
           			else {
                 			fprintf(stdout,"Error in aln_path\n");
              			}
         		}
       			ta[i] = -3;
       
       			seqs[j].data = (char *)realloc(seqs[j].data, (len+2) * sizeof (char));
       			for (i=0;i<len;i++) {
         			seqs[j].data[i] = ta[i];
			}
       			seqs[j].data[i] = -3;
       			seqs[j].len = len;
		}
	}      
	ta=ckfree((void *)ta);


if (debug>0)
{
  char c;
  sint d;

   for (i=0;i<nseqs;i++)
     {
if(sweight[i]>0) {
fprintf(stdout,"aligned %d\n",sweight[i]);
      for (j=0;j<len;j++)
       {
        if (seqs[i].data[j] == -3) break;
        else if (!isalpha(seqs[i].data[j]))  c = '-';
        else c = seqs[i].data[j];
        fprintf(stdout,"%c", c);
       }
      fprintf(stdout,"\n\n");
     }
}
}

}                  

/* 
   Calculate a score for a match in sequence versus sequence alignment
	n	residue in sequence 1
	m	residue in sequence 2
*/
static lint seqseqscore(sint n, sint m)
{
	lint  score;

	if(islower(seqs[seq1].data[n-1]) && islower(seqs[seq2].data[m-1]))
	score=matrix.score[seqs[seq1].data[n-1]-'a'][seqs[seq2].data[m-1]-'a']*10.0;
	else score=0;
	return(score/10);
   
}

static lint seqseqscoremotif(sint n, sint m)
{
        lint  score;

        if(islower(seqs[seq1].data[n-1]) && islower(seqs[seq2].data[m-1]))
        score=matrix.score[seqs[seq1].data[n-1]-'a'][seqs[seq2].data[m-1]-'a']*10.0+motif_score[n-1][m-1];
        else score=0;
        return(score/10);

}

/* 
   Calculate a score for a match in profile versus sequence alignment
	n	position in profile 1
	m	residue in sequence 2
*/
static lint prfseqscore(sint n, sint m)
{
	lint  score;

	if(islower(seqs[seq2].data[m-1]))
	score=profile1[n][seqs[seq2].data[m-1]-'a']*10.0;
	else score=0;
   	return(score/10);
}

static lint prfseqscoremotif(sint n, sint m)
{
        lint  score;

        if(islower(seqs[seq2].data[m-1]))
        score=profile1[n][seqs[seq2].data[m-1]-'a']*10.0+motif_score[n-1][m-1];
        else score=0;
        return(score/10);
}

/* 
   Calculate a score for a match in profile versus profile alignment
	n	position in profile 1
	m	position in profile 2
*/

static lint prfprfscore(sint n, sint m)
{
	sint    ix;
	lint  score;

	score = 0.0;
   	for (ix=0; ix<NUMRES; ix++) 
       		score += (profile1[n][ix] * profile2[m][ix]);
	return(score/10);
   
}

static lint prfprfscoremotif(sint n, sint m)
{  
        sint    ix;
        lint  score;

        score = 0.0;
        for (ix=0; ix<NUMRES; ix++)  
                score += (profile1[n][ix] * profile2[m][ix]);
        score+=motif_score[n-1][m-1];
        return(score/10);

}


/*
   Use the result of the Myers and Millers alignment to calculate the
   paths through the two profiles/sequences.
   From M+M, displ[i]=0	if residue in profile1 aligned with residue in profile2
	     displ[i]=n	if n>o, insert n residues in profile 1
			if n<0, insert n residues in profile 2
	     print_ptr  points to last entry in displ[].

	aln_path1[i]=2	residue in profile1 aligned with residue in profile2
		    =1  gap in profile1 
	aln_path2[i]=2	residue in profile2 algned with residue in profile1
		    =1  gap in profile2
*/
static sint 	ptracepath(char *aln_path1, char *aln_path2)
{
    sint i,j,k,pos,to_do;

    pos = 0;

    to_do=print_ptr-1;

    for(i=1;i<=to_do;++i) {
if (debug>2) fprintf(stdout,"%d ",(pint)displ[i]);
            if(displ[i]==0) {
                    aln_path1[pos]=2;
                    aln_path2[pos]=2;
                    ++pos;
            }
            else {
                    if((k=displ[i])>0) {
                            for(j=0;j<=k-1;++j) {
                                    aln_path2[pos+j]=2;
                                    aln_path1[pos+j]=1;
                            }
                            pos += k;
                    }
                    else {
                            k = (displ[i]<0) ? displ[i] * -1 : displ[i];
                            for(j=0;j<=k-1;++j) {
                                    aln_path1[pos+j]=2;
                                    aln_path2[pos+j]=1;
                            }
                            pos += k;
                    }
            }
    }
if (debug>2) fprintf(stdout,"\n");

   return pos;

}

static void pdel(sint k)
{
        if(last_print<0)
                last_print = displ[print_ptr-1] -= k;
        else
                last_print = displ[print_ptr++] = -(k);
}

static void padd(sint k)
{

        if(last_print<0) {
                displ[print_ptr-1] = k;
                displ[print_ptr++] = last_print;
        }
        else
                last_print = displ[print_ptr++] = k;
}

static void palign(void)
{
        displ[print_ptr++] = last_print = 0;
}

/* 
   Align two profiles/sequences using Myers and Miller memory efficient
   dynamic programming.
*/

static lint pdiff(sint A,sint B,sint M,sint N,sint go1, sint go2)
{
        sint midi,midj,type;
        sint t, tl, g, h;
	sint i,j;
        lint midh;
        lint thh, f, e, s;

/* Boundary cases: M <= 1 or N == 0 */
if (debug>5) fprintf(stdout,"A %d B %d M %d N %d midi %d go1 %d go2 %d\n", 
(pint)A,(pint)B,(pint)M,(pint)N,(pint)M/2,(pint)go1,(pint)go2);

/* if sequence B is empty....                                            */

        if(N<=0)  {

/* if sequence A is not empty....                                        */

                if(M>0) {

/* delete residues A[1] to A[M]                                          */

                        pdel(M);
                }
                return(-gap_penalty1(A,B,M));
        }

/* if sequence A is empty....                                            */

        if(M<=1) {
                if(M<=0) {

/* insert residues B[1] to B[N]                                          */

                        padd(N);
                        return(-gap_penalty2(A,B,N));
                }

/* if sequence A has just one residue....                                */
                if (go1 == 0)
                        midh =  -gap_penalty1(A+1,B+1,N);
                else
                        midh =  -gap_penalty2(A+1,B,1)-gap_penalty1(A+1,B+1,N);
                midj = 0;
                for(j=1;j<=N;j++) {
                        thh = -gap_penalty1(A,B+1,j-1) + scorefunction(A+1,B+j)
                            -gap_penalty1(A+1,B+j+1,N-j);
                        if(thh>midh) {
                                midh = thh;
                                midj = j;
                        }
                }

                if(midj==0) {
                        padd(N);
                        pdel(1);
                }
                else {
                        if(midj>1) padd(midj-1);
                        palign();
                        if(midj<N) padd(N-midj);
                }
                return midh;
        }


/* Divide sequence A in half: midi */

        midi = M / 2;

/* In a forward phase, calculate all HH[j] and HH[j] */

        HH[0] = 0.0;
        t = -open_penalty1(A,B+1);
        tl = -ext_penalty1(A,B+1);
        for(j=1;j<=N;j++) {
                HH[j] = t = t+tl;
                DD[j] = t-open_penalty2(A+1,B+j);
        }

		if (go1 == 0) t = 0;
		else t = -open_penalty2(A+1,B);
        tl = -ext_penalty2(A+1,B);
        for(i=1;i<=midi;i++) {
                s = HH[0];
                HH[0] = thh = t = t+tl;
                f = t-open_penalty1(A+i,B+1);

                for(j=1;j<=N;j++) {
                	g = open_penalty1(A+i,B+j);
                	h = ext_penalty1(A+i,B+j);
                        if ((thh=thh-g-h) > (f=f-h)) f=thh;
                	g = open_penalty2(A+i,B+j);
                	h = ext_penalty2(A+i,B+j);
                        if ((thh=HH[j]-g-h) > (e=DD[j]-h)) e=thh;
                        thh = s + scorefunction(A+i, B+j);
                        if (f>thh) thh = f;
                        if (e>thh) thh = e;

                        s = HH[j];
                        HH[j] = thh;
                        DD[j] = e;

                }
        }

        DD[0]=HH[0];

/* In a reverse phase, calculate all RR[j] and SS[j] */

        RR[N]=0.0;
        tl = 0.0;
        for(j=N-1;j>=0;j--) {
                g = -open_penalty1(A+M,B+j+1);
                tl -= ext_penalty1(A+M,B+j+1);
                RR[j] = g+tl;
                SS[j] = RR[j]-open_penalty2(A+M,B+j);
                gS[j] = open_penalty2(A+M,B+j);
        }

        tl = 0.0;
        for(i=M-1;i>=midi;i--) {
                s = RR[N];
                if (go2 == 0) g = 0;
                else g = -open_penalty2(A+i+1,B+N);
                tl -= ext_penalty2(A+i+1,B+N);
                RR[N] = thh = g+tl;
                t = open_penalty1(A+i,B+N);
                f = RR[N]-t;

                for(j=N-1;j>=0;j--) {
                	g = open_penalty1(A+i,B+j+1);
                	h = ext_penalty1(A+i,B+j+1);
                        if ((thh=thh-g-h) > (f=f-h-g+t)) f=thh;
                        t = g;
                	g = open_penalty2(A+i+1,B+j);
                	h = ext_penalty2(A+i+1,B+j);
                        thh=RR[j]-g-h;
                        if (i==(M-1)) {
				 e=SS[j]-h;
			}
                        else {
				e=SS[j]-h-g+open_penalty2(A+i+2,B+j);
				gS[j] = g;
			}
                        if (thh > e) e=thh;
                        thh = s + scorefunction(A+i+1, B+j+1);
                        if (f>thh) thh = f;
                        if (e>thh) thh = e;

                        s = RR[j];
                        RR[j] = thh;
                        SS[j] = e;

                }
        }
        SS[N]=RR[N];
        gS[N] = open_penalty2(A+midi+1,B+N);

/* find midj, such that HH[j]+RR[j] or DD[j]+SS[j]+gap is the maximum */

        midh=HH[0]+RR[0];
        midj=0;
        type=1;
        for(j=0;j<=N;j++) {
                thh = HH[j] + RR[j];
                if(thh>=midh)
                        if(thh>midh || (HH[j]!=DD[j] && RR[j]==SS[j])) {
                                midh=thh;
                                midj=j;
                        }
        }

        for(j=N;j>=0;j--) {
                thh = DD[j] + SS[j] + gS[j];
                if(thh>midh) {
                        midh=thh;
                        midj=j;
                        type=2;
                }
        }

/* Conquer recursively around midpoint                                   */


        if(type==1) {             /* Type 1 gaps  */
if (debug>5) fprintf(stdout,"Type 1,1: midj %d\n",(pint)midj);
                pdiff(A,B,midi,midj,go1,1);
if (debug>5) fprintf(stdout,"Type 1,2: midj %d\n",(pint)midj);
                pdiff(A+midi,B+midj,M-midi,N-midj,1,go2);
        }
        else {
if (debug>5) fprintf(stdout,"Type 2,1: midj %d\n",(pint)midj);
                pdiff(A,B,midi-1,midj,go1, 0);
                pdel(2);
if (debug>5) fprintf(stdout,"Type 2,2: midj %d\n",(pint)midj);
                pdiff(A+midi+1,B+midj,M-midi-1,N-midj,0,go2);
        }

        return midh;       /* Return the score of the best alignment */
}

/* calculate the score for opening a gap at residues A[i] and B[j]       */

static sint open_penalty1(sint i, sint j)
{
   sint g;

   if (!nendgappenalties &&(i==0)) return(0);
   if (!cendgappenalties &&(i==prf_length1)) return(0);

   g = profile2[j][GAPCOL] + profile1[i][GAPCOL];
   return(g);
}

/* calculate the score for extending an existing gap at A[i] and B[j]    */

static sint ext_penalty1(sint i, sint j)
{
   sint h;

   if (!nendgappenalties &&(i==0)) return(0);
   if (!cendgappenalties &&(i==prf_length1)) return(0);

   h = profile2[j][LENCOL]+profile1[i][LENCOL];
   return(h);
}

/* calculate the score for a gap of length k, at residues A[i] and B[j]  */

static sint gap_penalty1(sint i, sint j, sint k)
{
   sint ix;
   sint gp;
   sint g, h = 0;

   if (k <= 0) return(0);
   if (!nendgappenalties &&(i==0)) return(0);
   if (!cendgappenalties &&(i==prf_length1)) return(0);

   g = profile2[j][GAPCOL] + profile1[i][GAPCOL];
   for (ix=0;ix<k && ix+j<prf_length2;ix++)
      h += profile2[ix+j][LENCOL]+profile1[i][LENCOL];

   gp = g + h;
   return(gp);
}
/* calculate the score for opening a gap at residues A[i] and B[j]       */

static sint open_penalty2(sint i, sint j)
{
   sint g;

   if (!nendgappenalties &&(j==0)) return(0);
   if (!cendgappenalties &&(j==prf_length2)) return(0);

   g = profile1[i][GAPCOL] + profile2[j][GAPCOL];
   return(g);
}

/* calculate the score for extending an existing gap at A[i] and B[j]    */

static sint ext_penalty2(sint i, sint j)
{
   sint h;

   if (!nendgappenalties &&(j==0)) return(0);
   if (!cendgappenalties &&(j==prf_length2)) return(0);

   h = profile1[i][LENCOL]+profile2[j][LENCOL];
   return(h);
}

/* calculate the score for a gap of length k, at residues A[i] and B[j]  */

static sint gap_penalty2(sint i, sint j, sint k)
{
   sint ix;
   sint gp;
   sint g, h = 0;

   if (k <= 0) return(0);
   if (!nendgappenalties &&(j==0)) return(0);
   if (!cendgappenalties &&(j==prf_length2)) return(0);

   g = profile1[i][GAPCOL] + profile2[j][GAPCOL];
   for (ix=0;ix<k && ix+i<prf_length1;ix++)
      h += profile1[ix+i][LENCOL]+profile2[j][LENCOL];

   gp = g + h;
   return(gp);
}

sint calc_motif_scores(SEQ *seqs, sint nseqs, sint *group, MOTIF *motifs, sint **motif_score, sint prf_length1, sint prf_length2)
{

	char c;
	sint i,j, m,n1,n2;
	sint seq1, seq2, pos1, pos2, len, weight;
	sint ix1,ix2;
   
	for (i=0; i<prf_length1; i++)
		for (j=0; j<prf_length2; j++)
			motif_score[i][j]=0;

/*
   count number of sequences in each profile 
*/
	n1=n2=0;
	for(i=0;i<nseqs;i++)
		if(group[i]==1) n1++;
		else if(group[i]==2) n2++;

/*
   Check for motifs
*/
	if(motifs==NULL) return 0;

/*     check the list of motifs */
	for(m=0;motifs[m].weight!=0;m++) {
		seq1=motifs[m].seq1-1;
		seq2=motifs[m].seq2-1;
		pos1=motifs[m].pos1;
		pos2=motifs[m].pos2;
		len=motifs[m].len;
		if(group[seq1]==1 && group[seq2]==2) {
			ix1=start_pos(seqs[seq1].data,pos1);
			ix2=start_pos(seqs[seq2].data,pos2);
			weight=motifs[m].weight;
if(debug>0) fprintf(stdout,"1. %s %s %d %d %d %d weight %d len %d\n",seqs[seq1].name,seqs[seq2].name,pos1,pos2,ix1,ix2,weight,len);
			for(i=0;i<len;i++) {
				if(ix1>=prf_length1) break;
				if(ix2>=prf_length2) break;
				if(weight>motif_score[ix1][ix2]) motif_score[ix1][ix2]=weight/2.0;
				ix1=next_pos(seqs[seq1].data,prf_length1,ix1);
				ix2=next_pos(seqs[seq2].data,prf_length2,ix2);
			}
		}
		else if(group[seq2]==1 && group[seq1]==2) {
			ix1=start_pos(seqs[seq1].data,pos1);
			ix2=start_pos(seqs[seq2].data,pos2);
			weight=motifs[m].weight;
if(debug>0) fprintf(stdout,"2. %s %s %d %d %d %d weight %d len %d\n",seqs[seq1].name,seqs[seq2].name,pos1,pos2,ix1,ix2,weight,len);
			for(i=0;i<len;i++) {
				if(ix1>=prf_length2) break;
				if(ix2>=prf_length1) break;
				if(weight>motif_score[ix2][ix1]) motif_score[ix2][ix1]=weight/2.0;
				ix1=next_pos(seqs[seq1].data,prf_length2,ix1);
				ix2=next_pos(seqs[seq2].data,prf_length1,ix2);
			}
		}
	}
if(debug>5) {
fprintf(stdout,"\n");
for(i=0;i<prf_length1;i++) {
for(j=0;j<prf_length2;j++)
fprintf(stdout,"%d ",motif_score[i][j]);
fprintf(stdout,"\n");
}
}
	return m;
}
/* index into the sequence including gaps of the pth residue */
static sint start_pos(char *seq,sint p)
{
        sint i,ix=0;

        for(i=0;seq[i]!='\0';i++) {
                if(isalpha(seq[i])) ix++;
                if(ix>=p) break;
        }

        return i;
}

static sint next_pos(char *seq,sint len,sint p)
{
        p++;
        if(p>=len) return (len-1);
        while(!isalpha(seq[p]) && p<len) p++;
        return p;
}

