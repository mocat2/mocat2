#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "rascal.h"

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

static void subgroup_blocks(char *infile,char *clusterfile,sint window,sint query_group) ;
static void all_blocks(sint window);
static float countid2(SEQ seq1,SEQ seq2);

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
	char query_name[MAXNAMES+1];
	FILE *tree;
	sint maxres,*gapptr=NULL;
	IN_TREEPTR itree;
	double dscore;
	sint window;
	sint query_seq,query_group;
	OPT opt;
	sint *istart,*iend;

	if(argc!=4) {
		fprintf(stdout,"Usage: %s input_aln cluster_file query_seq\n",argv[0]);
		exit(1);
	}

	strcpy(infile,argv[1]);
	strcpy(clusterfile,argv[2]);
	strcpy(query_name,argv[3]);
	window=8;


        init_options(&opt);

/* read in the sequences */
	seq_input(infile,opt.explicit_type,FALSE,&mult_aln);
	if(mult_aln.nseqs<=0) {
		error("No sequences in %s\n",infile);
		exit(1);
	}

/* find the query sequence */
        query_seq=(-1);
        for(i=0;i<mult_aln.nseqs;i++) {
                if(strcasecmp(query_name,mult_aln.seqs[i].name)==0) query_seq=i;
        }
        if(query_seq==-1) {
                fprintf(stdout,"Error: the specified reference sequence (%s) was not found in the alignment\n",query_name);
                exit(1);
        }

/* find the start and end position of each sequence */
        istart = (sint *) ckalloc( (mult_aln.nseqs+1) * sizeof (sint) );
        iend = (sint *) ckalloc( (mult_aln.nseqs+1) * sizeof (sint) );
        for(i=0;i<mult_aln.nseqs;i++) {
		for(j=0;j<mult_aln.seqs[i].len;j++)
			if(isalpha(mult_aln.seqs[i].data[j])) break;
		istart[i]=j;
		for(j=mult_aln.seqs[i].len-1;j>=0;j--)
			if(isalpha(mult_aln.seqs[i].data[j])) break;
		iend[i]=j;
	}


/* count pairwise residue percent identities */
        tmat = (double **) ckalloc( (mult_aln.nseqs+1) * sizeof (double *) );
        for(i=0;i<mult_aln.nseqs;i++)
                tmat[i] = (double *)ckalloc( (mult_aln.nseqs+1) * sizeof (double) );

        for (i=0;i<mult_aln.nseqs;i++) {
                for (j=i+1;j<mult_aln.nseqs;j++) {
                        dscore = countid2(mult_aln.seqs[i],mult_aln.seqs[j]);
                        tmat[j][i] = tmat[i][j] = (100.0 - dscore)/100.0;
                }
        }

/* read in the clusters */
	groups=(GROUP *)ckalloc((mult_aln.nseqs+2) * sizeof(GROUP));
	secgroup=(sint *)ckalloc((mult_aln.nseqs+1)*sizeof(sint));
	orggroup=(sint *)ckalloc((mult_aln.nseqs+1)*sizeof(sint));

	ngroups=read_secator_groups(mult_aln,tmat,clusterfile,groups,secgroup,orggroup,0.0);
	if(ngroups<=0) exit(1);

	for(i=0;i<mult_aln.nseqs;i++)
		fprintf(stdout,"SECGROUP %s %d %d\n",mult_aln.seqs[i].name,secgroup[i],orggroup[i]);

	query_group=secgroup[query_seq];

/* make a tree from the percent identities (used for sequence weighting) */
        for (i=0;i<mult_aln.nseqs;i++) {
                for (j=i+1;j<mult_aln.nseqs;j++) {
                        dscore = countid1(mult_aln.seqs[i],mult_aln.seqs[j]);
                        tmat[j][i] = tmat[i][j] = (100.0 - dscore)/100.0;
                }
        }

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

	for(i=0;i<mult_aln.nseqs;i++)
		fprintf(stdout,"SEQWEIGHT %s %d\n",mult_aln.seqs[i].name,seq_weight[i]);

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

	subgroup_blocks(infile,clusterfile,window,query_group);
        for(i=0;i<mult_aln.nseqs;i++)
                ckfree(tmat[i]);
        ckfree(tmat);

}
 
static void subgroup_blocks(char *infile,char *clusterfile,sint window,sint query_group) 
{
	int i,j,k,k1,l,seq;
	int n,len,olen;
	Boolean found;
	BLOCK **blocks;
	PROF **profiles1;
	PROF **profiles2;
	sint nseqs;
	sint *nblocks;
	sint first,last;
	sint score;
	float meanpcid;
	sint cutoff=120000;
	sint *maxscore;

	nseqs=mult_aln.nseqs;
/* calculate a profile for each cluster, in each core block */
	nblocks=(sint *)ckalloc((ngroups+2)*sizeof(sint));
	blocks=(BLOCK **)ckalloc((ngroups+2)*sizeof(BLOCK *));
	profiles1=(PROF **)ckalloc((ngroups+2)*sizeof(PROF *));
	profiles2=(PROF **)ckalloc((ngroups+2)*sizeof(PROF *));

	for(i=0;i<ngroups;i++) {
		if(groups[i].len<2) continue;
fprintf(stdout,"SUBGROUP %d\n",groups[i].len);
		for(j=0;j<groups[i].len;j++) {
fprintf(stdout,"%d %s\n",groups[i].seqs[j]+1,mult_aln.seqs[groups[i].seqs[j]].name);
		}
/* define the core blocks for the sub-group i */
		blocks[i]=(BLOCK *)ckalloc((mult_aln.seqs[groups[i].seqs[0]].len+1) * sizeof(BLOCK));
		nblocks[i]=get_blocks_for_subgroup(mult_aln,tmat,groups[i],blocks[i],window,0.0);
	}

	for(i=0;i<ngroups;i++) {
		if(groups[i].len<2) continue;
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
			meanpcid=0;
			for(k=0,n=0;k<groups[i].len;k++)
				for(l=k;l<groups[i].len;l++) {
					meanpcid+=(1.0-tmat[groups[i].seqs[k]][groups[i].seqs[l]]);
					n++;
				}
			meanpcid/=(float)n;
fprintf(stdout,"meanid %d %d %.2f\n",i,j,meanpcid);
			calc_blockprf1(mult_aln,seq_weight,blocks[i][j].first,blocks[i][j].last,
				groups[i],matrix,is,ie,&profiles1[i][j],meanpcid);
			calc_blockprf2(mult_aln,seq_weight,blocks[i][j].first,blocks[i][j].last,
				groups[i],is,ie,&profiles2[i][j]);
		}
	}
	for(i=0;i<ngroups;i++) {
		if(groups[i].len<=1) continue;
for(j=0;j<nblocks[i];j++)
fprintf(stdout,"LOCALCOREBLOCK %d %d %d %d\n",i+1,blocks[i][j].first+1,blocks[i][j].last+1,blocks[i][j].code);
	}

fprintf(stdout,"BLOCK VERSUS BLOCK SCORES:\n");

	for(i=0;i<ngroups;i++) {
		if(groups[i].len<2) continue;
		for(j=i;j<ngroups;j++)
		if(groups[j].len>=2) {
/* calculate the score for the current block alignment */
			for(k=0;k<nblocks[i];k++) {
				for(l=0;l<nblocks[j];l++) {
					len=MAX(blocks[i][k].last-blocks[i][k].first,blocks[j][l].last-blocks[j][l].first);
					olen=overlap(blocks[i][k].first,blocks[i][k].last,blocks[j][l].first,blocks[j][l].last);
					if(olen>len/4.0 || olen>20) {
fprintf(stdout,"checking grp %d %d %d %d ",i,j,k,l);
fprintf(stdout,"len %d %d %d %d %d %d\n",blocks[i][k].first,blocks[i][k].last,blocks[j][l].first,blocks[j][l].last,len,olen);
						if(profiles1[i][k].nseqs>=profiles2[j][l].nseqs) {
							score=score_block_vs_block(blocks[i][k].first,blocks[j][l].first,profiles1[i][k],profiles2[j][l]);
fprintf(stdout,"score1 %d %d\n",score,cutoff);
						}
						else {
							score=score_block_vs_block(blocks[j][l].first,blocks[i][k].first,profiles2[j][l],profiles1[i][k]);
fprintf(stdout,"score2 %d %d\n",score,cutoff);
						}
						if(score>cutoff) {

							if(blocks[i][k].code==0 && blocks[j][l].code==0) {
								blocks[i][k].code=i+1;	
								blocks[j][l].code=i+1;	
							}
							else if(blocks[i][k].code>0)
								blocks[j][l].code=blocks[i][k].code;
							else if(blocks[j][l].code>0)
								blocks[i][k].code=blocks[j][l].code;
							first=MAX(blocks[i][k].first,blocks[j][l].first);
							last=MIN(blocks[i][k].last,blocks[j][l].last);
fprintf(stdout,"SHAREDCOREBLOCK %d %d %d %d %d\n",i+1,j+1,first+1,last+1,score);
}
					}
				}
			}
		}
	}

/* now compare the orphans to query subgroup blocks */
	for(i=0;i<ngroups;i++) {
		if(groups[i].len!=1) continue;
		seq=groups[i].seqs[0];
		blocks[i]=(BLOCK *)ckalloc((mult_aln.seqs[seq].len+1) * sizeof(BLOCK));
		maxscore=(sint *)ckalloc((mult_aln.seqs[seq].len+1) * sizeof(sint));
		nblocks[i]=0;
                for(j=0;j<ngroups;j++) {
                        if(groups[j].len<2) continue;
                                for(l=0;l<nblocks[j];l++) {
                                        score=score_sequence(mult_aln.seqs[seq].data,profiles1[j][l],blocks[j][l].first,blocks[j][l].last);
					olen=blocks[j][l].last-blocks[j][l].first;
                                        if(score>cutoff) {
/* create a new block for the orphan */
                                                blocks[i][nblocks[i]].first=blocks[j][l].first;
                                                blocks[i][nblocks[i]].last=blocks[j][l].last;
                                                blocks[i][nblocks[i]].code=j+1;
                                                maxscore[nblocks[i]]=score;
                                                nblocks[i]++;
                                        }
                                        else if (olen>40) {
/* check for a partial hit */
						score=score_partial_sequence(mult_aln.seqs[seq].data,profiles1[j][l],blocks[j][l].first,blocks[j][l].last,&first,&last);
                                        	if(score>cutoff) {
fprintf(stdout,"seq_score_partial %d %s %d %s %d %d %d %d %d\n",i,mult_aln.seqs[seq].name,j,mult_aln.seqs[groups[j].seqs[0]].name, blocks[j][l].first,blocks[j][l].last,score,blocks[j][l].first+first,blocks[j][l].first+last);
/* create a new block for the orphan */
                                                	blocks[i][nblocks[i]].first=blocks[j][l].first+first;
                                                	blocks[i][nblocks[i]].last=blocks[j][l].first+last;
                                                	blocks[i][nblocks[i]].code=j+1;
                                                	maxscore[nblocks[i]]=score;
                                                	nblocks[i]++;
                                        	}
					}
                                }
                }

		for(j=0;j<nblocks[i];j++) {
			if(blocks[i][j].first<is[groups[i].seqs[0]]) blocks[i][j].first=is[groups[i].seqs[0]];
			if(blocks[i][j].last<is[groups[i].seqs[0]]) blocks[i][j].last=is[groups[i].seqs[0]];
			if(blocks[i][j].first>ie[groups[i].seqs[0]]) blocks[i][j].first=ie[groups[i].seqs[0]];
			if(blocks[i][j].last>ie[groups[i].seqs[0]]) blocks[i][j].last=ie[groups[i].seqs[0]];
fprintf(stdout,"SHAREDCOREBLOCK %d %d %d %d %d\n",blocks[i][j].code,i+1,blocks[i][j].first+1,blocks[i][j].last+1,maxscore[j]);
		}
		ckfree(maxscore);
	}

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
	ncore_blocks=get_blocks(mult_aln,core_blocks,window,block_cutoff,0.8,minlength);
	if(ncore_blocks<=0) {
		fprintf(stdout,"No core blocks found in alignment\n");
		return;
	}
	for(i=0;i<ncore_blocks;i++) 
		fprintf(stdout,"GLOBALCOREBLOCK %d %d\n",core_blocks[i].first+1,core_blocks[i].last+1);


}

static float countid2(SEQ s1,SEQ s2)
{
   sint i,j;
   sint count,len,pcid;
   sint window=20;
   float score;
   char *seq1,*seq2;

/* remove gap positions from alignment of 2 sequences */
   seq1=(char *)ckalloc((s1.len+1)*sizeof(char));
   seq2=(char *)ckalloc((s2.len+1)*sizeof(char));
   for(i=0,len=0;i<s1.len && i<s2.len;i++) {
     if(isalpha(s1.data[i])||isalpha(s2.data[i])) {
         seq1[len]=s1.data[i];
         seq2[len]=s2.data[i];
         len++;
     }
   }

/* calculate number of identities in a fixed window length along alignment and
count number of times identity is greater than a threshold */
   count = 0;
   for (i=window;i<len-window;i++) {
       pcid=0;
       for(j=i-window;j<i+window;j++) {
         if(seq1[j]==seq2[j]) {
                   pcid++;
         }
       }
       if(pcid>0.4*(float)window*2.0) count++;
   }
   if(len-window*2.0<=0) score=0;
   else score=(float)(100.0 * (float)count/(float)(len-window*2.0));

   ckfree(seq1);
   ckfree(seq2);

   return score;
}

