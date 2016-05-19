#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "rascal.h"

#define MIN(a,b) ((a)<(b)?(a):(b))

static void subgroup_blocks(char *infile,char *clusterfile,sint window);
static void calc_median_limits(float *scores,sint n,float *ll, float *ul);

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
                fprintf(stdout,"SECGROUP %s %d %d\n",mult_aln.seqs[i].name,orggroup[i],secgroup[i]);


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

/* get the core blocks for each sub-group and locate sequence errors */

	subgroup_blocks(infile,clusterfile,window);
        for(i=0;i<mult_aln.nseqs;i++)
                ckfree(tmat[i]);
        ckfree(tmat);

/* write out the sequences */
	/*for (i=0;i<mult_aln.nseqs;i++) mult_aln.seqs[i].output_index = i;*/
	n=0;
	for (i=0;i<ngroups;i++) {
		for(j=0;j<groups[i].len;j++) {
				k=groups[i].seqs[j];
				mult_aln.seqs[k].output_index = n++;
		}
	}

}
 

static void subgroup_blocks(char *infile,char *clusterfile,sint window) 
{
	int i,j,k,l,s,n,e,s1,e1,s2;
	BLOCK **blocks;
	PROF **profiles1;
	sint nseqs;
	sint *nblocks;
	sint **seq_err;
	sint *scorethreshold;
	float ul,ll;
        float *scores;
	sint **allscores;
	/*float cutoff=0.75;*/
	float cutoff=0.75;
	Boolean in_err;

	nseqs=mult_aln.nseqs;
/* calculate a profile for each cluster, in each core block */
	nblocks=(sint *)ckalloc((ngroups+2)*sizeof(sint));
	blocks=(BLOCK **)ckalloc((ngroups+2)*sizeof(BLOCK *));
	profiles1=(PROF **)ckalloc((ngroups+2)*sizeof(PROF *));

	for(i=0;i<ngroups;i++) {
		if(groups[i].len<2) continue;
/* define the core blocks for the sub-group i */
		blocks[i]=(BLOCK *)ckalloc((mult_aln.seqs[groups[i].seqs[0]].len+1) * sizeof(BLOCK));
fprintf(stdout,"GRP %d\n",i+1);
for(j=0;j<groups[i].len;j++) fprintf(stdout,"%s\n",mult_aln.seqs[groups[i].seqs[j]].name);
		nblocks[i]=get_blocks_for_subgroup(mult_aln,tmat,groups[i],blocks[i],window,0.5);
	}

	for(i=0;i<ngroups;i++) {
		if(groups[i].len<2) continue;
		profiles1[i]=(PROF *)ckalloc((nblocks[i]+1)*sizeof(PROF));
		allscores=(sint **)ckalloc((groups[i].len+1)*sizeof(sint *));
		for(s=0;s<groups[i].len;s++)
			allscores[s]=(sint *)ckalloc((nblocks[i]+1)*sizeof(sint));
		seq_err=(sint **)ckalloc((groups[i].len+1)*sizeof(sint *));
		for(s=0;s<groups[i].len;s++)
			seq_err[s]=(sint *)ckalloc((nblocks[i]+1)*sizeof(sint));

		scorethreshold=(sint *)ckalloc((nblocks[i]+1)*sizeof(sint));
		for(j=0;j<nblocks[i];j++) {
			l=blocks[i][j].last-blocks[i][j].first;
        		profiles1[i][j].data = (sint **) ckalloc( (l+2) * sizeof (sint *) );
        		for(k=0; k<l+2; k++)
                		profiles1[i][j].data[k] = (sint *) ckalloc( (LENCOL+2) * sizeof(sint) );

			calc_blockprf1(mult_aln,seq_weight,blocks[i][j].first,blocks[i][j].last,
				groups[i],matrix,is,ie,&profiles1[i][j],1.0);

			for(s=0;s<groups[i].len;s++) {
				if(is[groups[i].seqs[s]]<blocks[i][j].first && ie[groups[i].seqs[s]]>blocks[i][j].last) {
					allscores[s][j]=score_sequence(mult_aln.seqs[groups[i].seqs[s]].data,
						profiles1[i][j],blocks[i][j].first,blocks[i][j].last);
				}
			}
		}

/* for each block, calculate a threshold score based on the median score for all sequences in the sub-group */
		for(j=0;j<nblocks[i];j++) {
			l=blocks[i][j].last-blocks[i][j].first;
			n=0;
			scores=(float *)ckalloc((groups[i].len+1)*sizeof(float));
			for(s=0;s<groups[i].len;s++) {
				if(is[groups[i].seqs[s]]<blocks[i][j].first && ie[groups[i].seqs[s]]>blocks[i][j].last) {
					scores[n]=allscores[s][j];
fprintf(stdout,"%s %d %d %d %d %f\n",mult_aln.seqs[groups[i].seqs[s]].name,i+1,j+1,blocks[i][j].first,blocks[i][j].last,scores[n]);
					n++;
				}
			}
			if(n==0) continue;
			calc_median_limits(scores,n,&ll,&ul);
			scorethreshold[j]=ll;
fprintf(stdout,"%d\n",scorethreshold[j]);
			ckfree(scores);

		}
/* for each sequence, compare scores for each block */
		for(s=0;s<groups[i].len;s++) {
			for(j=0;j<nblocks[i];j++) {

				l=blocks[i][j].last-blocks[i][j].first;
				if(is[groups[i].seqs[s]]<blocks[i][j].first && ie[groups[i].seqs[s]]>blocks[i][j].last) {
					if((scorethreshold[j]>0) && allscores[s][j]/(float)scorethreshold[j]<cutoff) {
						seq_err[s][j]=2;
					}
				}
			}

		}
/* write out regions of sequences with errors against its own subgroup */
		for(s=0;s<groups[i].len;s++) {
			in_err=FALSE;
			for(j=0;j<nblocks[i];j++) {
				if(in_err==FALSE) {
					if(seq_err[s][j]>1) {
						if(j==0) s2=0;
						else s2=blocks[i][j-1].last;
						s1=blocks[i][j].first;
						in_err=TRUE;
					}
				}
				else {
					if(seq_err[s][j]<=1) {
						e=blocks[i][j].first;
						e1=blocks[i][j-1].last;
						fprintf(stdout,"SEQ_ERROR %d %s %d %d\n",groups[i].seqs[s]+1,mult_aln.seqs[groups[i].seqs[s]].name,s2+1,e+1);
						in_err=FALSE;
					}
				}
			}
			if(in_err==TRUE) {
				e=mult_aln.seqs[groups[i].seqs[s]].len-1;
				e1=blocks[i][nblocks[i]-1].last;
				fprintf(stdout,"SEQ_ERROR %d %s %d %d\n",groups[i].seqs[s]+1,mult_aln.seqs[groups[i].seqs[s]].name,s2+1,e+1);
				in_err=FALSE;
			}
		}
		ckfree(scorethreshold);
		for(s=0;s<groups[i].len;s++)
			ckfree(allscores[s]);
		ckfree(allscores);
		for(s=0;s<groups[i].len;s++)
			ckfree(seq_err[s]);
		ckfree(seq_err);
	}
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


