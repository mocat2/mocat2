#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "rascal.h"


static void all_blocks(char *infile,char *clusterfile,sint window,sint block_cutoff,sint minlength);
static void calc_prf(ALN mult_aln,sint *seqweight,sint firstcol,sint lastcol,GROUP group,COMP_MATRIX matrix,sint *is,sint *ie,PROF *prf);
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
	double **tmat;
	sint window;
	sint block_cutoff;
	sint minlength;
	OPT opt;

	if(argc!=3) {
		fprintf(stdout,"Usage: %s input_aln cluster_file\n",argv[0]);
		exit(1);
	}

	strcpy(infile,argv[1]);
	strcpy(clusterfile,argv[2]);

        init_options(&opt);

/* read in the sequences */
	seq_input(infile,opt.explicit_type,FALSE,&mult_aln);
	if(mult_aln.nseqs<=0) {
		error("No sequences in %s\n",infile);
		exit(1);
	}

	if(mult_aln.nseqs<100) block_cutoff=8;
	else if(mult_aln.nseqs<250) block_cutoff=6;
	else block_cutoff=5;
	window=8;

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
	if(ngroups<=0) {
		ngroups=1;
		groups[0].seqs=(sint *)ckalloc((mult_aln.nseqs+1)*sizeof(sint));
		groups[0].len=mult_aln.nseqs;
		for(i=0;i<mult_aln.nseqs;i++) groups[0].seqs[i]=i;
		
	}

        for(i=0;i<mult_aln.nseqs;i++) {
		/*if(groups[secgroup[i]-1].len>1)
                fprintf(stdout,"SECGROUP %s %d\n",mult_aln.seqs[i].name,secgroup[i]);
		else
                fprintf(stdout,"SECGROUP %s 0\n",mult_aln.seqs[i].name);*/
                fprintf(stdout,"SECGROUP %s %d\n",mult_aln.seqs[i].name,secgroup[i]);
	}

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

	minlength=4;
	all_blocks(infile,clusterfile,window,block_cutoff,minlength);
}

static void all_blocks(char *infile,char *clusterfile,sint window,sint block_cutoff,sint minlength) 
{
	sint i,j,k,l,s,e,s1,e1,n;
	Boolean in_err;
	BLOCK *blocks;
	PROF **profiles;
	int nseqs;
	int nblocks;
	sint **seq_err;
	sint **scorethreshold;
	float ul,ll;
        float *scores;
	sint ***allscores;
	float cutoff=0.25;

	nseqs=mult_aln.nseqs;

/* define the core blocks */
	blocks=(BLOCK *)ckalloc((mult_aln.seqs[0].len+1) * sizeof(BLOCK));
	nblocks=get_blocks(mult_aln,blocks,window,block_cutoff,0.0,minlength);
	if(nblocks<=0) {
		fprintf(stdout,"No core blocks found in alignment\n");
		return;
	}
	/*for(i=0;i<nblocks;i++) 
		fprintf(stdout,"COREBLOCK %d %d\n",blocks[i].first+1,blocks[i].last+1);*/


/* calculate a profile for each cluster, in each core block */
	profiles=(PROF **)ckalloc((ngroups+1)*sizeof(PROF *));
	for(i=0;i<ngroups;i++) {
		profiles[i]=(PROF *)ckalloc((nblocks+1)*sizeof(PROF));
	}

	for(i=0;i<ngroups;i++) {
		if(groups[i].len<=1) continue;
		for(j=0;j<nblocks;j++) {
			l=blocks[j].last-blocks[j].first;
        		profiles[i][j].data = (sint **) ckalloc( (l+2) * sizeof (sint *) );
        		for(k=0; k<l+2; k++)
                		profiles[i][j].data[k] = (sint *) ckalloc( (LENCOL+2) * sizeof(sint) );

			calc_prf(mult_aln,seq_weight,blocks[j].first,blocks[j].last,
				groups[i],matrix,is,ie,&profiles[i][j]);
		}
	}

/* for each profile, calculate the score for each sequence in the alignment against the profile,
   and calculate a threshold score for each profile */

	allscores=(sint ***)ckalloc((ngroups+1)*sizeof(sint **));
	for(i=0;i<ngroups;i++) {
		allscores[i]=(sint **)ckalloc((nblocks+1)*sizeof(sint *));
		for(j=0;j<nblocks;j++) {
			allscores[i][j]=(sint *)ckalloc((nseqs+1)*sizeof(sint));
		}
	}

	scorethreshold=(sint **)ckalloc((ngroups+1)*sizeof(sint *));
	for(i=0;i<ngroups;i++)
		if(groups[i].len>1) scorethreshold[i]=(sint *)ckalloc((nblocks+1)*sizeof(sint));

	for(j=0;j<nblocks;j++) {
		for(k=0;k<ngroups;k++) {
			if(groups[k].len<=1) continue;
			n=0;
			for(i=0;i<nseqs;i++) {
		/*		if(is[i]<blocks[j].first && ie[i]>blocks[j].last) {*/
					allscores[k][j][i]=score_sequence(mult_aln.seqs[i].data,profiles[k][j],blocks[j].first,blocks[j].last);
				/*}*/
			}
		}
	}

/* calculate the threshold score for each profile, based on the scores for the sequences in the
corresponding group */
	for(j=0;j<nblocks;j++) {
		for(k=0;k<ngroups;k++) {
			if(groups[k].len<=1) continue;
			scores=ckalloc((mult_aln.nseqs+1)*sizeof(float));
			n=0;
			for(i=0;i<groups[k].len;i++) {
				if(is[groups[k].seqs[i]]<blocks[j].first && ie[groups[k].seqs[i]]>blocks[j].last) {
					scores[n++]=allscores[k][j][groups[k].seqs[i]];
				}
			}
			if(n==0) {
				scorethreshold[k][j]=0;
			}
			else {
				calc_median_limits(scores,n,&ll,&ul);
				scorethreshold[k][j]=ll;
			}
			ckfree(scores);
		}
	}

/* for each sequence, calculate the score for the sequence against each profile
   and report the best cluster in each block. If the sequence belonged to the
   unclustered group, simply report the best score in each block. Othewise, choose
   the original cluster, unless a significant difference in scores is observed */

	seq_err=(sint **)ckalloc((nseqs+1)*sizeof(sint *));
	for(i=0;i<nseqs;i++)
		seq_err[i]=(sint *)ckalloc((nblocks+1)*sizeof(sint));

	for(i=0;i<nseqs;i++) {
		if(secgroup[i]>0 && groups[secgroup[i]-1].len>1) continue; /* only check orphans! */
fprintf(stdout,"sequence %d %s, cluster %d\n",i+1,mult_aln.seqs[i].name,secgroup[i]);
		for(j=0;j<nblocks;j++) {
fprintf(stdout,"block %d %d %d, ",j,blocks[j].first+1,blocks[j].last+1);
		/*if(is[i]<blocks[j].first && ie[i]>blocks[j].last) {*/
			for(k=0;k<ngroups;k++) {
			if(groups[k].len>1)
fprintf(stdout,"%2d:%5d (%d)",k+1,allscores[k][j][i]/100,scorethreshold[k][j]/100);
			}
			n=0;
			for(k=0;k<ngroups;k++) {
				if(groups[k].len<=1) continue;
				if(!(k==secgroup[i]-1 && groups[k].len<=1)) {
					if((scorethreshold[k][j]>0) && (float)allscores[k][j][i]/(float)scorethreshold[k][j]>=cutoff) {
						n++;
					}
				}
			}
fprintf(stdout,"\n");
/* if n==0, then the sequence has low scores against all subgroups in this block */
			if(n==0 ) {
				seq_err[i][j]=2;
			}
		/*}*/
		}
	}
/* write out regions of sequences with errors against all subgroups */
	for(i=0;i<nseqs;i++) {
		if(secgroup[i]>0 && groups[secgroup[i]-1].len>1) continue; /* only check orphans! */
		in_err=FALSE;
		for(j=0;j<nblocks;j++) {
			if(in_err==FALSE) {
				if(seq_err[i][j]>1) {
					if(j==0) s=0;
					else s=blocks[j-1].last;
					s1=blocks[j].first;
					in_err=TRUE;
				}
			}
			else {
				if(seq_err[i][j]<=1) {
					e=blocks[j].first;
					e1=blocks[j-1].last;
					if(secgroup[i]>0 && groups[secgroup[i]-1].len>1)
					fprintf(stdout,"SEQ_ERROR %d %s %d %d\n",i+1,mult_aln.seqs[i].name,s+1,e+1);
					else
					fprintf(stdout,"ORPHAN_ERROR %d %s %d %d\n",i+1,mult_aln.seqs[i].name,s+1,e+1);
					in_err=FALSE;
				}
			}
		}
		if(in_err==TRUE) {
			e=mult_aln.seqs[i].len-1;
			e1=blocks[nblocks-1].last;
			if(secgroup[i]>0 && groups[secgroup[i]-1].len>1)
			fprintf(stdout,"SEQ_ERROR %d %s %d %d\n",i+1,mult_aln.seqs[i].name,s+1,e+1);
			else
			fprintf(stdout,"ORPHAN_ERROR %d %s %d %d\n",i+1,mult_aln.seqs[i].name,s+1,e+1);
			in_err=FALSE;
		}
	}

        for(i=0;i<mult_aln.nseqs;i++)
		ckfree(seq_err[i]);
	ckfree(seq_err);

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

static void calc_prf(ALN mult_aln,sint *seqweight,sint firstcol,sint lastcol,GROUP group,COMP_MATRIX matrix,sint *is,sint *ie,PROF *prf)
{
	sint l,s,i;
	sint res,d;
	sint f;
	sint *weight,sum;
	Boolean *fragment;
	sint **freq;

/* normalise the sequence weights for this group to sum to 100, so that
   the profiles for each group are directly comparable. */

	if(group.len==0) return;

        weight=(sint *)ckalloc((group.len+1)*sizeof(sint));
        fragment=(Boolean *)ckalloc((group.len+1)*sizeof(Boolean));
	sum=0;

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
			prf->data[l-firstcol][res]=(sint)f;
		}
	}

	ckfree(fragment);
	ckfree(weight);
        for(i=0; i<lastcol-firstcol+2; i++)
		ckfree(freq[i]);
	ckfree(freq);

}

