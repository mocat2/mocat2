#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "rascal.h"

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

static void all_blocks(sint window);

extern sint    gon250mt[];

GROUP *groups;
sint ngroups;
ALN mult_aln;
sint *secgroup;
sint *seq_weight;
COMP_MATRIX matrix;
sint *is, *ie;
sint maxi,maxj;
double **tmat;
BLOCK *core_blocks;
sint ncore_blocks;

int main(int argc,char **argv)
{
	sint i,j,k,n,s;
	sint status;
	char c;
	char infile[FILENAMELEN+1];
	char treefile[FILENAMELEN+1];
	FILE *tree;
	sint maxres,*gapptr=NULL;
	IN_TREEPTR itree;
	double dscore;
	sint window;
	OPT opt;

	if(argc!=2) {
		fprintf(stdout,"Usage: %s input_aln\n",argv[0]);
		exit(1);
	}

	strcpy(infile,argv[1]);
	window=8;


        init_options(&opt);

/* read in the sequences */
	seq_input(infile,opt.explicit_type,FALSE,&mult_aln);
	if(mult_aln.nseqs<=0) {
		error("No sequences in %s\n",infile);
		exit(1);
	}

/* count pairwise residue percent identities */
        tmat = (double **) ckalloc( (mult_aln.nseqs+1) * sizeof (double *) );
        for(i=0;i<mult_aln.nseqs;i++)
                tmat[i] = (double *)ckalloc( (mult_aln.nseqs+1) * sizeof (double) );

        for (i=0;i<mult_aln.nseqs;i++) {
                for (j=i+1;j<mult_aln.nseqs;j++) {
                        dscore = countid1(mult_aln.seqs[i],mult_aln.seqs[j]);
                        tmat[j][i] = tmat[i][j] = (100.0 - dscore)/100.0;
                }
        }

/* make a tree from the percent identities (used for sequence weighting) */
	/*strcpy(treefile,infile);
	strcat(treefile,".ph");*/
	strcpy(treefile,".ph");
        if((tree = open_explicit_file(treefile))==NULL) exit(1);

        guide_tree(tree,mult_aln.seqs,mult_aln.nseqs, tmat, QUICKNJ);
        itree=(IN_TREEPTR)ckalloc(sizeof(IN_TREE));

        status = read_tree(treefile, mult_aln.seqs, 0, mult_aln.nseqs,itree);

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

/* get the global core blocks */
	all_blocks(window);

        for(i=0;i<mult_aln.nseqs;i++)
                ckfree(tmat[i]);
        ckfree(tmat);

}
 
static void all_blocks(sint window)
{
	sint i;
	sint block_cutoff;
	sint minlength;

        if(mult_aln.nseqs<100) block_cutoff=10;
        else if(mult_aln.nseqs<250) block_cutoff=6;
        else block_cutoff=5;
	minlength=4;

/* define the core blocks */
	core_blocks=(BLOCK *)ckalloc((mult_aln.seqs[0].len+1) * sizeof(BLOCK));
	ncore_blocks=get_blocks(mult_aln,core_blocks,window,block_cutoff,1.8,minlength);
	if(ncore_blocks<=0) {
		fprintf(stdout,"No core blocks found in alignment\n");
		return;
	}
	fprintf(stdout,"%d %d\n",core_blocks[0].first+1,core_blocks[ncore_blocks-1].last+1);


}


