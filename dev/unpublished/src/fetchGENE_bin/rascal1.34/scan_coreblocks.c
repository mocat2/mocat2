#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "rascal.h"


static void all_blocks(char *infile,sint window,sint block_cutoff);

extern sint    gon250mt[];

ALN mult_aln;
sint *seq_weight;
COMP_MATRIX matrix;
sint *is, *ie;

int main(int argc,char **argv)
{
	sint i,j,k,n,s;
	sint status;
	sint result_type;
	char c;
	char infile[FILENAMELEN+1];
	char treefile[FILENAMELEN+1];
	float meanid;
	FILE *tree;
	sint maxres,*gapptr=NULL;
	IN_TREEPTR itree;
	double dscore;
	double meanscore;
	double **tmat;
	sint window;
	sint block_cutoff;
	OPT opt;

	if(argc!=2) {
		fprintf(stdout,"Usage: %s input_aln\n",argv[0]);
		exit(1);
	}

	strcpy(infile,argv[1]);

        init_options(&opt);

/* read in the sequences */
	seq_input(infile,opt.explicit_type,FALSE,&mult_aln);
	if(mult_aln.nseqs<=0) {
		error("No sequences in %s\n",infile);
		exit(1);
	}


	window=8;

/* count pairwise residue percent identities */
        tmat = (double **) ckalloc( (mult_aln.nseqs+1) * sizeof (double *) );
        for(i=0;i<mult_aln.nseqs;i++)
                tmat[i] = (double *)ckalloc( (mult_aln.nseqs+1) * sizeof (double) );

	meanscore=0;
        for (i=0,n=0;i<mult_aln.nseqs;i++) {
                for (j=i+1;j<mult_aln.nseqs;j++) {
                        dscore = countid(mult_aln.seqs[i],mult_aln.seqs[j]);
                        tmat[j][i] = tmat[i][j] = (100.0 - dscore)/100.0;
			n++;
			meanscore+=dscore;
                }
        }
	meanscore/=(float)n;

	/*if(mult_aln.nseqs<100) block_cutoff=8;
	else if(mult_aln.nseqs<250) block_cutoff=6;
	else*/ block_cutoff=5;

	if(meanscore>50) block_cutoff=50;

/* make a tree from the percent identities (used for sequence weighting) */
	strcpy(treefile,infile);
	strcat(treefile,".ph");
        if((tree = open_explicit_file(treefile))==NULL) exit(1);

        guide_tree(tree,mult_aln.seqs,mult_aln.nseqs, tmat, QUICKNJ);
        itree=(IN_TREEPTR)ckalloc(sizeof(IN_TREE));

        status = read_tree(treefile, mult_aln.seqs, 0, mult_aln.nseqs,itree);
        for(i=0;i<mult_aln.nseqs;i++)
                ckfree(tmat[i]);
        ckfree(tmat);

        if (status < 0) exit(1);


        seq_weight = calc_seq_weights(0,mult_aln.nseqs,itree,FALSE);
        free_tree(itree);
        remove(treefile);


/* find the start and end positions of each sequence */

	is = (sint *)ckalloc((mult_aln.nseqs+1) * sizeof(sint));
	ie = (sint *)ckalloc((mult_aln.nseqs+1) * sizeof(sint));
	for(s=0;s<mult_aln.nseqs;s++) {
		is[s]=0;
                ie[s] = mult_aln.seqs[s].len;
                for (i=0; i<mult_aln.seqs[s].len; i++) {
                        c = mult_aln.seqs[s].data[i];
                        if (!isalpha(c))
                                is[s]++;
                        else
                                break;
                }
                for (i=mult_aln.seqs[s].len-1; i>=0; i--) {
                        c = mult_aln.seqs[s].data[i];
                        if (!isalpha(c))
                                ie[s]--;
                        else
                                break;
                }
	}

	matrix.format=0;
	maxres = get_cl_matrix(FALSE, gon250mt, gapptr, TRUE, 100, &matrix);

	all_blocks(infile,window,block_cutoff);
}
 
static void all_blocks(char *infile,sint window,sint block_cutoff) 
{
	sint i,j,k,k1,l,s,e,s1,e1,n,n1;
	Boolean found,in_err;
	BLOCK *blocks;
	int nseqs;
	int nblocks;
	sint first,last;
	sint minlength=2;
	float conserved;

	nseqs=mult_aln.nseqs;
	if(nseqs<=20) conserved=0.9;
	else conserved=0.5;

/* define the core blocks */
	blocks=(BLOCK *)ckalloc((mult_aln.seqs[0].len+1) * sizeof(BLOCK));
	nblocks=get_blocks(mult_aln,blocks,window,block_cutoff,conserved,minlength);
	if(nblocks<=0) {
		fprintf(stdout,"No core blocks found in alignment\n");
		return;
	}
	for(i=0;i<nblocks;i++) 
		fprintf(stdout,"GLOBALCOREBLOCK %d %d\n",blocks[i].first+1,blocks[i].last+1);


}

static void calc_median_limits(float *scores,sint n,float *ll, float *ul)
{
	sint i;
	float t,q1,q3;
	float *sortedscores;

       	if(n==0) {
               	(*ul)=(*ll)=0;
       	}
       	else {
		sortedscores=ckalloc((n+1)*sizeof(float));
		for (i=0;i<n;i++)
			sortedscores[i]=scores[i];

		sort_scores(sortedscores,0,n-1);
               	t = n/4.0 + 0.5;
               	if(t - (int)t == 0.5) {
                       	q3=(sortedscores[(int)t]+sortedscores[(int)t+1])/2.0;
                       	q1=(sortedscores[n-(int)t]+sortedscores[n-(int)t-1])/2.0;
               	}
               	else if(t - (int)t > 0.5) {
                       	q3=sortedscores[(int)t+1];
                       	q1=sortedscores[n-(int)t-1];
               	}
               	else {
                       	q3=sortedscores[(int)t];
                       	q1=sortedscores[n-(int)t];
               	}
               	if (n<4) {
			(*ul)=sortedscores[0];
			(*ll)=sortedscores[n-1];
		}
              	else {
			(*ul)=q3+(q3-q1)/2.0;
			(*ll)=q1-(q3-q1);
		}
       }
	ckfree(sortedscores);
}

static void calc_mean_limits(float *scores,sint n,float *ll, float *ul)
{
	sint i;
	float mean,sd;

	mean=sd=0.0;
       	if(n==0) {
               	(*ul)=(*ll)=0;
       	}
       	else {
		for(i=0;i<n;i++)
			mean+=scores[i];
		mean/=(float)n;
		for(i=0;i<n;i++)
			sd=(mean-scores[i])*(mean-scores[i]);
		sd/=(float)n;
		sd=(float)sqrt((double)sd);
		(*ul)=mean+sd*3.0;
		(*ul)=mean-sd;
		(*ll)=mean-sd*3.0;
        }
}

static void calc_prf(ALN mult_aln,sint *seqweight,sint firstcol,sint lastcol,GROUP group,COMP_MATRIX matrix,sint *is,sint *ie,PROF *prf)
{
	sint l,s,i;
	sint res,d;
	sint f;
	sint *weight,sum;
	Boolean *fragment;
	sint **freq;
	char c;

/* normalise the sequence weights for this group to sum to 100, so that
   the profiles for each group are directly comparable. */

	if(group.len==0) return;

        weight=(sint *)ckalloc((group.len+1)*sizeof(sint));
        fragment=(Boolean *)ckalloc((group.len+1)*sizeof(Boolean));
	sum=0.0;

	for(s=0;s<group.len;s++) {
		/*if(is[group.seqs[s]]>firstcol || ie[group.seqs[s]]<lastcol)*/
		if(is[group.seqs[s]]>lastcol || ie[group.seqs[s]]<firstcol)
                	fragment[s]=TRUE;
	}

	for(s=0;s<group.len;s++) {
		if(!fragment[s]) sum+=seqweight[group.seqs[s]];
	}

	if(sum<=0) return;
	for(s=0;s<group.len;s++) 
                weight[s]=100*((float)seqweight[group.seqs[s]]/(float)sum);

        freq = (sint **) ckalloc( (lastcol-firstcol+2) * sizeof (sint *) );
        for(i=0; i<lastcol-firstcol+2; i++)
                freq[i] = (sint *) ckalloc( (LENCOL+2) * sizeof(sint) );

	for(l=firstcol;l<=lastcol;l++) {
		for(s=0;s<group.len;s++) {
			if(!fragment[s] && isalpha(mult_aln.seqs[group.seqs[s]].data[l])) {
				i=mult_aln.seqs[group.seqs[s]].data[l]-'a';
				freq[l-firstcol][i]+=weight[s];
			}
		}
	}

	for(l=firstcol;l<=lastcol;l++) {
		for (res=0; res<NUMRES; res++) {
			f=0;
			for (d=0; d<NUMRES; d++)
				f += (freq[l-firstcol][d] * matrix.score[d][res]);
			prf->data[l-firstcol][res]=(float)f;
		}
	}

	ckfree(fragment);
	ckfree(weight);
        for(i=0; i<lastcol-firstcol+2; i++)
		ckfree(freq[i]);
	ckfree(freq);

}

