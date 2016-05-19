#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "rascal.h"

#define MIN(a,b) ((a)<(b)?(a):(b))

static void subgroup_blocks(char *infile,char *clusterfile,sint window);
static sint maxscore_block_vs_block(PROF prf1,PROF prf2);
static void all_blocks(sint window);
static void output_blerr(BLOCK **blocks,sint **blockmaxi,sint **blockmaxj,sint g1,sint b1,sint g2,sint b2,sint maxscore);

extern sint    gon250mt[];

GROUP *groups;
sint ngroups;
ALN mult_aln;
sint *secgroup;
sint *orggroup;
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
	char clusterfile[FILENAMELEN+1];
	char treefile[FILENAMELEN+1];
	FILE *tree;
	sint maxres,*gapptr=NULL;
	IN_TREEPTR itree;
	double dscore;
	sint window;
	OPT opt;

	if(argc!=3) {
		fprintf(stdout,"Usage: %s input_aln cluster_file\n",argv[0]);
		exit(1);
	}

	strcpy(infile,argv[1]);
	strcpy(clusterfile,argv[2]);
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
                        dscore = countid(mult_aln.seqs[i],mult_aln.seqs[j]);
                        tmat[j][i] = tmat[i][j] = (100.0 - dscore)/100.0;
                }
        }

/* read in the clusters */
	groups=(GROUP *)ckalloc((mult_aln.nseqs*2) * sizeof(GROUP));
	secgroup=(sint *)ckalloc((mult_aln.nseqs+1)*sizeof(sint));
	orggroup=(sint *)ckalloc((mult_aln.nseqs+1)*sizeof(sint));

	ngroups=read_secator_groups(mult_aln,tmat,clusterfile,groups,secgroup,orggroup,0.4);
	if(ngroups<=0) exit(1);

	for(i=0;i<mult_aln.nseqs;i++)
		fprintf(stdout,"SECGROUP %s %d %d\n",mult_aln.seqs[i].name,secgroup[i],orggroup[i]);

/* make a tree from the percent identities (used for sequence weighting) */
	strcpy(treefile,infile);
	strcat(treefile,".ph");
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

/* get the core blocks for each sub-group */

	subgroup_blocks(infile,clusterfile,window);
        for(i=0;i<mult_aln.nseqs;i++)
                ckfree(tmat[i]);
        ckfree(tmat);

}
 
static void subgroup_blocks(char *infile,char *clusterfile,sint window) 
{
	int i,j,k,k1,l;
	Boolean found;
	BLOCK **blocks;
	PROF **profiles1;
	PROF **profiles2;
	sint nseqs;
	sint *nblocks;
	sint first,last;
	sint **maxblockscores;	/* optimal score for aligning 2 profiles */
	sint **blockmaxi;	/* end position of optimal profile alignment */
	sint **blockmaxj;	/* end position of optimal profile alignment */
	sint *blockscoresi,*blockscoresj; /* score for aligning 2 profiles, as they are
						aligned in multiple alignment */
	sint *maxblockscoresi,*maxblockscoresj; /* maximum score for aligning 2 profiles that are overlapped
						in multiple alignment */
	sint maxscore;
	sint maxblock;
	float maxblock_factor=1.6;

	nseqs=mult_aln.nseqs;
/* calculate a profile for each cluster, in each core block */
	nblocks=(sint *)ckalloc((ngroups+2)*sizeof(sint));
	blocks=(BLOCK **)ckalloc((ngroups+2)*sizeof(BLOCK *));
	profiles1=(PROF **)ckalloc((ngroups+2)*sizeof(PROF *));
	profiles2=(PROF **)ckalloc((ngroups+2)*sizeof(PROF *));
/*
		for(j=0;j<groups[i].len;j++) {
			secgroup[groups[i].seqs[j]]=i+1;
			mult_aln.seqs[groups[i].seqs[j]].simgroup=i+1;
		}
*/

	for(i=0;i<ngroups;i++) {
		if(groups[i].len<2) continue;
fprintf(stdout,"SUBGROUP %d\n",groups[i].len);
		for(j=0;j<groups[i].len;j++) {
fprintf(stdout,"%d %s\n",groups[i].seqs[j]+1,mult_aln.seqs[groups[i].seqs[j]].name);
		}
/* define the core blocks for the sub-group i */
		blocks[i]=(BLOCK *)ckalloc((mult_aln.seqs[groups[i].seqs[0]].len+1) * sizeof(BLOCK));
		nblocks[i]=get_blocks_for_subgroup(mult_aln,tmat,groups[i],blocks[i],window,0.5);
	}

	for(i=0;i<ngroups;i++) {
		if(groups[i].len<2) continue;
for(j=0;j<nblocks[i];j++)
fprintf(stdout,"LOCALCOREBLOCK %d %d %d\n",i+1,blocks[i][j].first+1,blocks[i][j].last+1);
		profiles1[i]=(PROF *)ckalloc((nblocks[i]+1)*sizeof(PROF));
		profiles2[i]=(PROF *)ckalloc((nblocks[i]+1)*sizeof(PROF));

		for(j=0;j<nblocks[i];j++) {
			l=blocks[i][j].last-blocks[i][j].first;
        		profiles1[i][j].data = (sint **) ckalloc( (l+2) * sizeof (sint *) );
        		for(k=0; k<l+2; k++)
                		profiles1[i][j].data[k] = (sint *) ckalloc( (LENCOL+2) * sizeof(sint) );
        		profiles2[i][j].data = (sint **) ckalloc( (l+2) * sizeof (sint *) );
        		for(k=0; k<l+2; k++)
                		profiles2[i][j].data[k] = (sint *) ckalloc( (LENCOL+2) * sizeof(sint) );

			calc_blockprf1(mult_aln,seq_weight,blocks[i][j].first,blocks[i][j].last,
				groups[i],matrix,is,ie,&profiles1[i][j],1.0);
			calc_blockprf2(mult_aln,seq_weight,blocks[i][j].first,blocks[i][j].last,
				groups[i],is,ie,&profiles2[i][j]);

		}

	}

fprintf(stdout,"BLOCK VERSUS BLOCK SCORES:\n");

	for(i=0;i<ngroups;i++) {
		if(groups[i].len<2) continue;
		for(j=0;j<ngroups;j++) {
			if(i==j) continue;
			if(groups[j].len<2) continue;
/* calculate the maximum possible score for all possible blocks (without destroying existing core blocks) */
			maxblockscores=(sint **)ckalloc((nblocks[i]+1)*sizeof(sint *));
			for(k=0;k<nblocks[i];k++) 
				maxblockscores[k]=(sint *)ckalloc((nblocks[j]+1)*sizeof(sint));
			blockmaxi=(sint **)ckalloc((nblocks[i]+1)*sizeof(sint *));
			for(k=0;k<nblocks[i];k++) 
				blockmaxi[k]=(sint *)ckalloc((nblocks[j]+1)*sizeof(sint));
			blockmaxj=(sint **)ckalloc((nblocks[i]+1)*sizeof(sint *));
			for(k=0;k<nblocks[i];k++) 
				blockmaxj[k]=(sint *)ckalloc((nblocks[j]+1)*sizeof(sint));
			for(k=0;k<nblocks[i];k++) {
				for(l=0;l<nblocks[j];l++) {
					maxblockscores[k][l]=maxscore_block_vs_block(profiles1[i][k],profiles2[j][l]);
					blockmaxi[k][l]=maxi;
					blockmaxj[k][l]=maxj;
fprintf(stdout,"g %d b %d (%d), g %d b %d (%d): %d %d %d \n",i+1,k+1,blocks[i][k].first,j+1,l+1,blocks[j][l].first,maxblockscores[k][l]/100,maxi,maxj);
				}
			}

/* calculate the score for the current block alignment */
			blockscoresi=(sint *)ckalloc((nblocks[i]+1)*sizeof(sint));
			maxblockscoresi=(sint *)ckalloc((nblocks[i]+1)*sizeof(sint));
			for(k=0;k<nblocks[i];k++) {
				blockscoresi[k]=0;
				maxblockscoresi[k]=0;
				for(l=0;l<nblocks[j];l++) {
					if(overlap(blocks[i][k].first,blocks[i][k].last,blocks[j][l].first,blocks[j][l].last)>0) {
						blockscoresi[k]+=score_block_vs_block(blocks[i][k].first,blocks[j][l].first,profiles1[i][k],profiles2[j][l]);
						if(maxblockscoresi[k]<maxblockscores[k][l]) maxblockscoresi[k]=maxblockscores[k][l];
					}
				}
			}
			blockscoresj=(sint *)ckalloc((nblocks[j]+1)*sizeof(sint));
			maxblockscoresj=(sint *)ckalloc((nblocks[j]+1)*sizeof(sint));
			for(l=0;l<nblocks[j];l++) {
				blockscoresj[l]=0;
				maxblockscoresj[l]=0;
				for(k=0;k<nblocks[i];k++) {
					if(overlap(blocks[i][k].first,blocks[i][k].last,blocks[j][l].first,blocks[j][l].last)>0) {
						blockscoresj[l]+=score_block_vs_block(blocks[j][l].first,blocks[i][k].first,profiles1[j][l],profiles2[i][k]);
						if(maxblockscoresj[l]<maxblockscores[k][l]) maxblockscoresj[l]=maxblockscores[k][l];
fprintf(stdout,"maxblockscoresj %d %d %d\n",l,k,maxblockscoresj[l]);
					}
				}
			}

			for(k=0;k<nblocks[i];k++) {

/* if the block overlays a core block, don't do anything */
				found=FALSE;
				for(k1=0;k1<ncore_blocks;k1++)
					if(overlap(blocks[i][k].first,blocks[i][k].last,core_blocks[k1].first,core_blocks[k1].last)>0) {
						found=TRUE;
						break;
					}
				if(found==TRUE) continue;


/* find the best block score for each block within the limits of the core blocks */
				maxscore=0;
				maxblock=(-1);

				if(ncore_blocks==0) first=0;
				else first=core_blocks[ncore_blocks-1].last;
				for(k1=0;k1<ncore_blocks;k1++)
					if(blocks[i][k].first<core_blocks[k1].last) {
						if(k1==0) first=0;
						else first=core_blocks[k1-1].last;
						break;
					}
				if(ncore_blocks==0) last=mult_aln.seqs[0].len;
				else last=core_blocks[0].first;
				for(k1=ncore_blocks-1;k1>=0;k1--)
					if(blocks[i][k].last>core_blocks[k1].first) {
						if(k1==ncore_blocks-1) last=mult_aln.seqs[0].len;
						else last=core_blocks[k1+1].first;
						break;
					}

				for(l=0;l<nblocks[j];l++) {
					if(blocks[j][l].first>first && blocks[j][l].last<last && maxscore<maxblockscores[k][l]) {
						maxscore=maxblockscores[k][l];
						maxblock=l;
					}
				}
fprintf(stdout,"test %d %d %d %d %d %d %d\n",k+1,first,last,maxblock+1,maxscore/100,blockscoresi[k]/100,blockscoresj[maxblock]/100);
				if(maxblock==(-1)) continue;

/* if the block is aligned on the best scoring block, check the blocks are aligned properly */
				if(overlap(blocks[i][k].first,blocks[i][k].last,blocks[j][maxblock].first,blocks[j][maxblock].last)>0) {
					if(blocks[i][k].first+blockmaxi[k][maxblock]!=blocks[j][maxblock].first+blockmaxj[k][maxblock]) {
fprintf(stdout,"test1 %d %d %d %d %d %d %d\n",k+1,maxblock+1,maxscore/100,blockscoresi[k]/100,blockscoresj[maxblock]/100,maxblockscoresi[k]/100,maxblockscoresj[maxblock]/100);
						if(maxscore>blockscoresi[k]*maxblock_factor && maxscore>=maxblockscoresi[k] && maxscore>=maxblockscoresj[maxblock]) {
							output_blerr(blocks,blockmaxi,blockmaxj,i,k,j,maxblock,maxscore);
}
					}

				}
/* if the block is not aligned on the best scoring block,
 find the score for the overlapping blocks in the alignment */
				else {
fprintf(stdout,"test2 %d %d %d %d %d %d %d\n",k+1,maxblock+1,maxscore/100,blockscoresi[k]/100,blockscoresj[maxblock]/100,maxblockscoresi[k]/100,maxblockscoresj[maxblock]/100);
/* if the two sequence groups are non-overlapping, increase the score cutoff to make sure the block error is real */
					if(blockscoresi[k]==0 || blockscoresj[maxblock]==0) {
						/* if(maxscore>250000) */
						if(maxscore>150000) 
						if((maxscore>blockscoresi[k]*maxblock_factor && maxscore>blockscoresj[maxblock]*maxblock_factor) &&
					        (maxscore>maxblockscoresi[k]*maxblock_factor && maxscore>maxblockscoresj[maxblock]*maxblock_factor)) {
						output_blerr(blocks,blockmaxi,blockmaxj,i,k,j,maxblock,maxscore);
}
					}
					else if((maxscore>blockscoresi[k]*maxblock_factor && maxscore>blockscoresj[maxblock]*maxblock_factor) &&
					        (maxscore>maxblockscoresi[k]*maxblock_factor && maxscore>maxblockscoresj[maxblock]*maxblock_factor)) {
						output_blerr(blocks,blockmaxi,blockmaxj,i,k,j,maxblock,maxscore);
}
				}
			}
			for(k=0;k<nblocks[i];k++) 
				ckfree(maxblockscores[k]);
			ckfree(maxblockscores);
			ckfree(blockscoresi);
			ckfree(blockscoresj);
			ckfree(maxblockscoresi);
			ckfree(maxblockscoresj);
			for(k=0;k<nblocks[i];k++) 
				ckfree(blockmaxi[k]);
			ckfree(blockmaxi);
			for(k=0;k<nblocks[i];k++) 
				ckfree(blockmaxj[k]);
			ckfree(blockmaxj);
		}
	}

}

static void output_blerr(BLOCK **blocks,sint **blockmaxi,sint **blockmaxj,sint g1,sint b1,sint g2,sint b2,sint maxscore)
{
	sint len;

	if(blockmaxi[b1][b2]<blockmaxj[b1][b2]) len=blockmaxi[b1][b2];
	else len=blockmaxj[b1][b2];

	if(maxscore>100000) {
		fprintf(stdout,"BLOCK_ERROR ");
fprintf(stdout,"g %d %d %d g %d %d %d %d\n",g1+1,blocks[g1][b1].first+blockmaxi[b1][b2]-len,blocks[g1][b1].first+blockmaxi[b1][b2],g2+1,blocks[g2][b2].first+blockmaxj[b1][b2]-len,blocks[g2][b2].first+blockmaxj[b1][b2],maxscore/100);
	}
}

static sint maxscore_block_vs_block(PROF prf1,PROF prf2)
{
	sint i,j,s;
	sint m,n,len;
	sint maxscore;
        sint INS,match,previous,t;

	static sint    *MATCH, *DEL;
	sint gop=1000000,gep=100;

	n=prf1.len;
	m=prf2.len;
	len=(float)MIN(n,m);
	if(len==0) return 0;

        MATCH = (sint *)ckalloc((n+m+1) * sizeof(sint));
        DEL = (sint *)ckalloc((n+m+1) * sizeof(sint));

	maxscore=0;
        for (i=0;i<=m;i++) {
                MATCH[i] = 0;
                DEL[i] = -gop;
                }

        for (i=1;i<=n;i++) {
                match = previous = 0;
                INS = -gop;

                for (j=1;j<=m;j++) {
                        INS -= gep;
                        t = match - gop - gep;
                        if (INS<t) INS = t;

                        DEL[j] -= gep;
                        t = MATCH[j] - gop - gep;
                        if (DEL[j]<t) DEL[j] = t;

			s=prfscore(prf1,prf2,i-1,j-1);
                        match = previous + s;
                        if (match<INS) {
                                match = INS;
                        }
                        if (match<DEL[j]) {
                                match = DEL[j];
                        }
                        if (match<0)
                                match = 0;

                        previous = MATCH[j];
                        MATCH[j] = match;

                        if (match > maxscore) {
                                maxscore = match;
				maxi=i;
				maxj=j;
                        }
                }
        }

        MATCH=ckfree((void *)MATCH);
        DEL=ckfree((void *)DEL);

	maxscore=(float)maxscore/(float)MIN(n,m);

	return maxscore;
}

static void all_blocks(sint window) 
{
	sint i;
	sint block_cutoff;
	sint minlength;

        if(mult_aln.nseqs<100) block_cutoff=8;
        else if(mult_aln.nseqs<250) block_cutoff=6;
        else block_cutoff=5;
	minlength=2;

/* define the core blocks */
	core_blocks=(BLOCK *)ckalloc((mult_aln.seqs[0].len+1) * sizeof(BLOCK));
	ncore_blocks=get_blocks(mult_aln,core_blocks,window,block_cutoff,0.8,minlength);
	if(ncore_blocks<=0) {
		fprintf(stdout,"No core blocks found in alignment\n");
		return;
	}
	for(i=0;i<ncore_blocks;i++) 
		fprintf(stdout,"GLOBALCOREBLOCK %d %d\n",core_blocks[i].first+1,core_blocks[i].last+1);


}


