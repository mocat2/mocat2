#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "rascal.h"


/*
 *       Prototypes
 */
static sint malign(ALNPTR mult_aln,MULT_OPT mult_opt,sint output_order,Boolean verbose,ALNCOUNT *alncount);
static void align_sequences_to_prf(ALNPTR mult_aln,double **tmat,sint *aligned,MULT_OPT mult_opt,sint output_order,Boolean use_maxid,Boolean verbose);
static void align_seq_to_prf(ALNPTR mult_aln,double **tmat,sint *aligned,sint seq1,MULT_OPT mult_opt,sint output_order,Boolean use_maxid,Boolean verbose);


static void reset_align(sint nseqs, SEQ *seqs,Boolean reset_new, Boolean reset_all);
static void reset_prf(sint prf_no,ALNPTR mult_aln,Boolean reset_new, Boolean reset_all);
static sint renumber_groups(sint nseqs,sint *aligned_group,sint first_group);
static void calc_group_dist(sint nseqs,sint naligned_groups,sint *aligned_group,double **tmat,float **group_dist);
static void calc_group_dist1(sint nseqs,sint naligned_groups,sint *aligned_group,double **tmat,float **group_dist);



/*
 *       Global Variables
 */

static sint     debug=0;

/*
   Do a progressive multiple alignment.
	mult_aln	the multiple alignment (sequences still unaligned usually!)
	mult_opt	multiple alignment options
	output_order	=ALIGNED if sequences to be output in alignment order, 
	             	=TREEORDER if sequences to be output in input tree order, 
			=INPUT if sequences to be output in same order as input

   Returns the number of sequences aligned
*/

static sint malign(ALNPTR mult_aln,MULT_OPT mult_opt,sint output_order,Boolean verbose,ALNCOUNT *alncount)
{
	Boolean use_maxid;
	sint 	i,j,n,tx,set;
	sint	ix,t1,t2;
	sint 	status,entries;
	sint	nsets;
	sint	**sets;
	sint	naligned_groups;
	sint	*aligned_group;
	sint	*aligned;
	lint	score = 0;
	sint	*seq_weight;
	IN_TREEPTR   itree;
	float min_dist,mean_dist;
	float  **group_dist;
	double **tmat;
	Boolean *delayed,ok;
	char *dsc_ss;
	char *jnet_ss;

	if (verbose) info("Start of Multiple Alignment");

/* get the guide tree (tree should be in phylip format) */

	if (mult_aln->nseqs >= 2) {
		itree=(IN_TREEPTR)ckalloc(sizeof(IN_TREE));
       		status = read_tree(mult_aln->treename, mult_aln->seqs, (sint)0, mult_aln->nseqs,itree);
       		if (status < 0) return((sint)0);
     	}

	if(status<0) use_maxid=TRUE;
	else use_maxid=FALSE;

/* group the sequences according to their relative divergence */

        sets = (sint **) ckalloc( (mult_aln->nseqs+1) * sizeof (sint *) );
        for(i=0;i<mult_aln->nseqs;i++)
           	sets[i] = (sint *)ckalloc( (mult_aln->nseqs+1) * sizeof (sint) );

        nsets=create_sets(mult_aln->nseqs,sets,itree);
        if(verbose) info("There are %d groups",(pint)nsets);
if(debug>1)
for(set=0;set<nsets;set++) {
for (i=0;i<mult_aln->nseqs;i++) {
printf("%d ",sets[set][i]);
}
printf("\n");
}

/* calculate sequence weights according to branch lengths of the tree -
   weights are normalised to sum to 100 */

	seq_weight = calc_seq_weights((sint)0, mult_aln->nseqs, itree,mult_opt.no_weights);
	for(i=0;i<mult_aln->nseqs;i++)
		mult_aln->seqs[i].weight=seq_weight[i];

/* calculate tmat matrix containing percent distances for all pairs of sequences */

	tmat = pw_distances_from_tree(mult_aln->seqs,mult_aln->nseqs,itree);
	if (tmat == NULL) return((sint)0);

/* clear the memory used for the phylogenetic tree */

        if (mult_aln->nseqs >= 2)
             	free_tree(itree);

/* start the multiple alignments.........  */

        if(verbose) info("Aligning...");

        if (output_order==INPUT) {
        	for (i=0;i<mult_aln->nseqs;i++) 
        		mult_aln->seqs[i].output_index = i;
        }
        else {
		tx=0;
        	for (i=0;i<mult_aln->nseqs;i++) 
        		mult_aln->seqs[i].output_index = -1;
        	for(set=0;set<nsets;set++) {
			for (i=0;i<mult_aln->nseqs;i++) {
				if (sets[set][i] != 0) {
					if(mult_aln->seqs[i].output_index==(-1)) {
						mult_aln->seqs[i].output_index = tx++;
					}
				}
			}
		}
	}

/* align the closely related groups first */
        aligned_group = (sint *) ckalloc( (mult_aln->nseqs+1) * sizeof (sint) );
        delayed = (Boolean *) ckalloc( (mult_aln->nseqs+1) * sizeof (Boolean) );

        for(set=0;set<nsets;set++) {
/* decide whether to do the alignment now - if the mean distance between sequences in the
two groups to be aligned is greater than the cutoff, then don't align now */
		mean_dist=0.0;
		entries=0;
		ok=TRUE;
		for (i=0,n=0;i<mult_aln->nseqs;i++) {
			if (sets[set][i] != 0 && mult_aln->seqs[i].len>0) {
				if(delayed[i]==TRUE) {
					ok=FALSE;
					break;
				}
				entries++;
				for (j=i+1;j<mult_aln->nseqs;j++) {
					if (sets[set][j] != 0  && mult_aln->seqs[j].len>0 && sets[set][i] != sets[set][j]) {
						if(delayed[j]==TRUE) {
							ok=FALSE;
							break;
						}
						mean_dist+=tmat[i][j];
						n++;
		/*				if(mean_dist>tmat[i][j]) mean_dist=tmat[i][j];*/
					}
				}
			}
		}
		if(ok==TRUE && n>0) mean_dist=mean_dist/(float)(n);

		if ((ok==TRUE) && mean_dist < 1.0-(float)mult_opt.divergence_cutoff/100.0) {
			score = prfalign(0,mult_aln,sets[set],tmat,use_maxid,mult_opt,alncount,TRUE);
			/*if (score < 0) return(-1);*/
             		if(verbose) info("Group %d: Sequences:%4d      Score:%d",
             			(pint)set+1,(pint)entries,(pint)score);
			for (i=0;i<mult_aln->nseqs;i++) {
				if (sets[set][i] != 0) {
					aligned_group[i]=set+1;
				}
			}
		}
		else {
             		if(verbose) info("Group %d: Sequences:%4d      Delayed",(pint)set+1,(pint)entries);
			for (i=0;i<mult_aln->nseqs;i++) 
				if (sets[set][i] != 0) {
					delayed[i]=TRUE;
				}
		}
        }

/* renumber the groups and assign a group number to any orphan sequences, that haven't 
been aligned yet */
	naligned_groups=renumber_groups(mult_aln->nseqs,aligned_group,0);
if (verbose) fprintf(stdout,"done groups %d\n",naligned_groups);

        for(i=0;i<mult_aln->nseqs;i++)
		mult_aln->seqs[i].simgroup=aligned_group[i];

        for(i=0;i<mult_aln->nseqs;i++)
		if(aligned_group[i]==0) {
			aligned_group[i]=naligned_groups+1;
			naligned_groups++;
		}
if (verbose) fprintf(stdout,"done groups and added orphans %d\n",naligned_groups);
	if(naligned_groups==1) return naligned_groups;

/* now align the two closest groups of sequences together */
/* the array aligned[] contains the group number of each sequence */

	ix=0;
	if(naligned_groups>1) {
		for(i=0;i<mult_aln->nseqs;i++) mult_aln->seqs[i].output_index=(-1);
        	group_dist = (float **) ckalloc( (naligned_groups+1) * sizeof (float *) );
        	for(i=0;i<naligned_groups+1;i++)
          		group_dist[i] = (float *)ckalloc( (naligned_groups+1) * sizeof (float) );

        	aligned = (sint *) ckalloc( (mult_aln->nseqs+1) * sizeof (sint) );

/* calculate the distance between each pair of groups */
		calc_group_dist(mult_aln->nseqs,naligned_groups,aligned_group,tmat,group_dist);
	
       		for(i=0;i<mult_aln->nseqs;i++)
       			sets[0][i] = 0;
/* align the two closest groups together */
		min_dist=1000.0;
		t1=t2=0;
		for (i=0;i<naligned_groups;i++) 
			for (j=i+1;j<naligned_groups;j++) 
				if(min_dist>=group_dist[i][j]) {
					min_dist=group_dist[i][j];
					t1=i+1;
					t2=j+1;
				}
       		for(i=0;i<mult_aln->nseqs;i++)
			if(aligned_group[i]==t1) {
				sets[0][i]=1;
				if (output_order==ALIGNED) mult_aln->seqs[i].output_index = ix++;
			}
       		for(i=0;i<mult_aln->nseqs;i++)
			if(aligned_group[i]==t2) {
				sets[0][i]=2;
				aligned_group[i]=t1;
				if (output_order==ALIGNED) mult_aln->seqs[i].output_index = ix++;
			}

		for (i=0,tx=0;i<mult_aln->nseqs;i++) 
			if (sets[0][i] != 0) {
				aligned[i]=1;
				tx++;
			}
		score = prfalign(0,mult_aln,sets[0],tmat,use_maxid,mult_opt,alncount,TRUE);
		/*if (score < 0) return(-1);*/
       			if(verbose) info("Group Align %d: Sequences:%4d      Score:%d",
				(pint)tx,(pint)entries,(pint)score);
	}


/* now align the remaining groups to the first two */
	naligned_groups=renumber_groups(mult_aln->nseqs,aligned_group,t1);
	t1=1;
	while(naligned_groups>1) {
/* calculate the distance of each remaining group to the first group */
		calc_group_dist1(mult_aln->nseqs,naligned_groups,aligned_group,tmat,group_dist);
	
/* align the next closest group with the first group */
       		for(i=0;i<mult_aln->nseqs;i++)
       			sets[0][i] = 0;
		min_dist=100.0;
		t2=0;
		for (j=1;j<naligned_groups;j++) 
			if(min_dist>=group_dist[0][j]) {
				min_dist=group_dist[0][j];
				t2=j+1;
			}
if (verbose) fprintf(stdout,"Aligning Group %d %d\n",t1,t2);
       		for(i=0;i<mult_aln->nseqs;i++)
			if(aligned_group[i]==1) 
				sets[0][i]=1;
			else if(aligned_group[i]==t2) {
				sets[0][i]=2;
				aligned_group[i]=1;
				if (output_order==ALIGNED) mult_aln->seqs[i].output_index = ix++;
			}
if(debug>0)
for(i=0;i<mult_aln->nseqs;i++)
fprintf(stdout,"%s %d %d\n",mult_aln->seqs[i].name,mult_aln->seqs[i].output_index,aligned_group[i]);

		for (i=0,tx=0;i<mult_aln->nseqs;i++) 
			if (sets[0][i] != 0) {
				aligned[i]=1;
				tx++;
			}
		score = prfalign(0,mult_aln,sets[0],tmat,use_maxid,mult_opt,alncount,TRUE);
		/*if (score < 0) return(-1);*/
       			if(verbose) info("Group Align %d: Sequences:%4d      Score:%d",
				(pint)tx,(pint)entries,(pint)score);
		naligned_groups=renumber_groups(mult_aln->nseqs,aligned_group,t1);
	}


        for(i=0;i<mult_aln->nseqs;i++)
                tmat[i] = ckfree(tmat[i]);
        tmat = ckfree(tmat);

       	for (i=0;i<naligned_groups;i++)
      		group_dist[i]=ckfree((void *)group_dist[i]);
       	group_dist=ckfree(group_dist);
       	aligned=ckfree(aligned);

        for (i=0;i<naligned_groups;i++)
          	sets[i]=ckfree((void *)sets[i]);
        sets=ckfree(sets);

	return((sint)1);
}

static void calc_group_dist(sint nseqs,sint naligned_groups,sint *aligned_group,double **tmat,float **group_dist)
{
	sint i,j;
	sint **n;

        n = (sint **) ckalloc( (naligned_groups+1) * sizeof (sint *) );
        for(i=0;i<naligned_groups;i++)
           	n[i] = (sint *)ckalloc( (naligned_groups+1) * sizeof (sint) );

	for (i=0;i<naligned_groups;i++) 
		for (j=0;j<naligned_groups;j++) {
			group_dist[i][j]=1000.0;
			n[i][j]=0;
		}
	for (i=0;i<nseqs;i++) 
		for (j=0;j<nseqs;j++) {
			/*group_dist[aligned_group[i]-1][aligned_group[j]-1]+=tmat[i][j];
			n[aligned_group[i]-1][aligned_group[j]-1]++;*/
			if(group_dist[aligned_group[i]-1][aligned_group[j]-1]>tmat[i][j]) group_dist[aligned_group[i]-1][aligned_group[j]-1]=tmat[i][j];
			}
/*
	for (i=0;i<naligned_groups;i++) 
		for (j=0;j<naligned_groups;j++) {
			if(n[i][j]>0) {
				group_dist[i][j]/=(float)n[i][j];
			}
		}
*/
        for (i=0;i<naligned_groups;i++)
          	n[i]=ckfree((void *)n[i]);
        n=ckfree(n);
/*
for (i=0;i<naligned_groups;i++) {
fprintf(stdout,"\ngroup dist %d ",i+1);
for (j=0;j<naligned_groups;j++) 
fprintf(stdout,"%.2f ",group_dist[i][j]);
fprintf(stdout,"\n");
}
*/
}

static void calc_group_dist1(sint nseqs,sint naligned_groups,sint *aligned_group,double **tmat,float **group_dist)
{
	sint i,j;
	sint *n;

        n = (sint *) ckalloc( (naligned_groups+1) * sizeof (sint) );

	for (j=0;j<naligned_groups;j++) {
		group_dist[0][j]=1000.0;
		n[j]=0;
	}
	for (i=0;i<nseqs;i++) 
		if(aligned_group[i]==1) 
			for (j=0;j<nseqs;j++) {
				/*group_dist[0][aligned_group[j]-1]+=tmat[i][j];
				n[aligned_group[j]-1]++;*/
				if(group_dist[0][aligned_group[j]-1]>tmat[i][j]) group_dist[0][aligned_group[j]-1]=tmat[i][j];
			}
/*
	for (j=0;j<naligned_groups;j++) {
		if(n[j]>0) {
			group_dist[0][j]/=(float)n[j];
		}
*/
        n=ckfree(n);
/*
fprintf(stdout,"\ngroup dist ");
for (j=0;j<naligned_groups;j++) 
fprintf(stdout,"%.2f ",group_dist[0][j]);
fprintf(stdout,"\n");
*/
}

static sint renumber_groups(sint nseqs,sint *aligned_group,sint first_group)
{
	sint i,j;
	sint naligned_groups;
	sint	*tmp;

        tmp = (sint *) ckalloc( (nseqs+1) * sizeof (sint) );
        for(i=0;i<nseqs;i++) {
		tmp[i]=aligned_group[i];
		aligned_group[i]=0;
	}
	naligned_groups=0;

	if(first_group>0) {
        	for(i=0;i<nseqs;i++)
			if(tmp[i]==first_group) {
				aligned_group[i]=1;
				tmp[i]=0;
			}
		naligned_groups++;
	}

        for(i=0;i<nseqs;i++)
		if(tmp[i]!=0) {
        		for(j=i+1;j<nseqs;j++)
				if(tmp[j]==tmp[i]) {
					aligned_group[j]=naligned_groups+1;
					tmp[j]=0;
				}
			aligned_group[i]=naligned_groups+1;
			tmp[i]=0;
			naligned_groups++;
		}
        tmp=ckfree(tmp);
	return naligned_groups;
}

/*
    Align all remaining sequences to a profile. Mult_aln in fact contains all the
    sequences, aligned or not.
	mult_aln	the multiple alignment (including unaligned sequences)
	tmat		pairwise sequence distance matrix
	aligned		=1 if sequences are in profile
			=0 if sequences are unaligned
	mult_opt	multiple alignment options
	output_order	=ALIGNED if sequences to be output in alignment order, 
	             	=TREEORDER if sequences to be output in input tree order, 
			=INPUT if sequences to be output in same order as input
*/
	
static void align_sequences_to_prf(ALNPTR mult_aln,double **tmat,sint *aligned,MULT_OPT mult_opt,sint output_order,Boolean use_maxid,Boolean verbose)
{
	sint i,j,ix;
	sint iseq;
	float min,*mindist;

        ix = 0;
        for (i=0;i<mult_aln->nseqs;i++) 
		if(aligned[i]==1) 
			ix++;

/* if we haven't aligned any sequences at all yet, make sure we align the
   two most closely related sequences first */

	if(ix==0) {
        	min = 1.0;
		iseq = 0;
        	for (i=0;i<mult_aln->nseqs;i++) {
             		for (j=i+1;j<mult_aln->nseqs;j++) {
                  		if (min > tmat[i][j]) {
                     			min = tmat[i][j];
                     			iseq = i;
                  		}
              		}
          	}
        	aligned[iseq]=1;
         	if (output_order == INPUT) 
            		mult_aln->seqs[iseq].output_index = iseq;
         	else if (output_order == ALIGNED)
            		mult_aln->seqs[iseq].output_index = 0;
		ix=1;
	}

   	mindist = (float *)ckalloc( (mult_aln->nseqs+1) * sizeof (float));
/* for each unaligned sequence, find it's closest pair amongst the aligned sequences.  */
    	while (ix < mult_aln->nseqs) {
		for (i=0;i<mult_aln->nseqs;i++) {
                	if (aligned[i] == 0) {
                     		mindist[i] = 1.0;
                     		for (j=0;j<mult_aln->nseqs;j++) 
                        		if ((mindist[i] > tmat[i][j]) && (aligned[j] != 0))
                            			mindist[i] = tmat[i][j];
                  	}
              	}

/* find the most closely related sequence to those already aligned */

            	min = 1.0;
	    	iseq = 0;
            	for (i=0;i<mult_aln->nseqs;i++) {
                	if ((aligned[i] == 0) && (mindist[i] < min)) {
                     		min = mindist[i];
                     		iseq = i;
                  	}
              	}


/* align this sequence to the existing alignment */

		align_seq_to_prf(mult_aln,tmat,aligned,iseq,mult_opt,output_order,use_maxid,verbose);
            	ix++;
      	}

	mindist=ckfree((void *)mindist);
}

/* remove gaps from older alignments (code = GAP1) */
/* EXCEPT for gaps that were INPUT with the seqs which have  code = GAP2  */
static void reset_align(sint nseqs, SEQ *seqs,Boolean reset_new, Boolean reset_all)
{
        register sint sl;
        sint i,j;

	if(nseqs==0) return;

	if(!reset_new && !reset_all) return;
        for(i=0;i<nseqs;i++) {
                sl=0;
                for(j=0;j<seqs[i].len;j++) {
                        if(seqs[i].data[j] == GAP1 && (reset_new || reset_all)) continue;
                        if(seqs[i].data[j] == GAP2 && (reset_all)) continue;
                        seqs[i].data[sl++]=seqs[i].data[j];
                }
                seqs[i].len=sl;
        }


}

/* remove gaps from older alignments (code = GAP1) */
/* EXCEPT for gaps that were INPUT with the seqs which have  code = GAP2  */
static void reset_prf(sint prf_no,ALNPTR mult_aln,Boolean reset_new, Boolean reset_all)
{
        register sint sl;                    
        sint i,j,fseq,nseq,len;


	if(prf_no==1) {
		fseq=0;
		nseq=mult_aln->prf1.nseqs;
	}
	else {
		fseq=mult_aln->prf1.nseqs;
		nseq=mult_aln->prf1.nseqs+mult_aln->prf2.nseqs;
	}
	if(nseq==0) return;

        len=mult_aln->seqs[fseq].len;

        for(i=fseq;i<nseq;i++) {
                sl=0;
                for(j=0;j<len;j++) {
                        if(mult_aln->seqs[i].data[j] == GAP1 && (reset_new || reset_all)) continue;
                        if(mult_aln->seqs[i].data[j] == GAP2 && (reset_all)) continue;
                        mult_aln->seqs[i].data[sl++]=mult_aln->seqs[i].data[j];
                }
                mult_aln->seqs[i].len=sl;
        }


}

void fix_gaps(ALNPTR mult_aln)   /* fix gaps introduced in older alignments (code = GAP1) */
{
        sint i,j;

	if(mult_aln->nseqs==0) return;

        for(i=0;i<mult_aln->nseqs;i++) {
                for(j=0;j<mult_aln->seqs[i].len;j++) {
                        if(mult_aln->seqs[i].data[j] == GAP1)
                                mult_aln->seqs[i].data[j]=GAP2;
                }
        }
}





void align_from_tree(ALNPTR mult_aln,OPT opt,Boolean usemenu,Boolean from_complete_align,Boolean verbose,ALNCOUNT *alncount)
{
	char path[FILENAMELEN+1],temp[MAXLINE+1];
	sint i,count,nanchors;
	
	if(mult_aln->nseqs<=0) {
		error("No sequences in memory. Load sequences first.");
		return;
	}

        if (mult_aln->nseqs == 1) {
		error("Less than 2 sequences in memory. Load sequences first.");
		return;
	}

	count = malign(mult_aln,*opt.mult_opt,opt.alnout_opt->output_order,verbose,alncount);
	if (count <= 0) return;

	if (usemenu) fprintf(stdout,"\n\n\n");

	mult_aln->treename[0]=EOS;

}


static void align_seq_to_prf(ALNPTR mult_aln,double **tmat,sint *aligned,sint iseq,MULT_OPT mult_opt,sint output_order,Boolean use_maxid,Boolean verbose)
{
	sint j,n;
	sint *group;
	lint score = 0;

	group = (sint *)ckalloc( (mult_aln->nseqs+1) * sizeof (sint));

        for (j=0,n=0;j<mult_aln->nseqs;j++)
         	if (aligned[j] != 0) {
               		group[j] = 1;
			n++;
	}
       	group[iseq] = 2;
        aligned[iseq] = 1;

        score = prfalign(2,mult_aln,group,tmat,use_maxid,mult_opt,NULL,TRUE);
        if(verbose) info("Sequence:%d     Score:%d",(pint)iseq+1,(pint)score);
        if (output_order == INPUT) {
        	mult_aln->seqs[iseq].output_index = iseq;
        }
        else if (output_order == ALIGNED)
        	mult_aln->seqs[iseq].output_index = n;

	group=ckfree((void *)group);
}
