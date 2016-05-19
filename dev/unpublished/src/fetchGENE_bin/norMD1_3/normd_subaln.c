/* NORMD version 1.2 Sept 2005 */

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "score.h"

void calc_submd(ALNPTR mult_aln,int group,int colix);
float md_subscore(ALNPTR mult_aln,int group,int colix,float matrix[NUMRES][NUMRES],float **weight,float totweight,Boolean *fragment);
static float normalise_subscore(float score,float n,float ntot,float ntotseq);

int *use_seq;
int seqlength;
int *useqlen_array;
Boolean *fragment;
int minlen,maxlen;

Boolean verbose=FALSE;

float matrix[NUMRES][NUMRES];
float go;
float ge;
float egap;


float **seqweight;
float totweight;

float gap_openscore,gap_extscore,end_gapscore;
float **gop;
float **gep;
float **egp;

float norm_md,col,normd_rs;
float max_colscore;
int query;

int ngroups;
int *seqgroup;
int *groupseed;
int norphans;

float **pcid;


void usage(char *prog) 
{
	fprintf(stderr,"%s Version 1.2\n",prog);
	fprintf(stderr,"Usage: %s in_file outfile\n",prog);
	fprintf(stderr,"   OR: \n");
	fprintf(stderr,"Usage: %s in_file outfile matrix      | calculate norMD for alignment (-v verbose)\n",prog);
}

int main(int argc, char **argv)
{
	FILE *ofd,*ifd;
        ALN mult_aln;
        OPT opt;
	char infile[FILENAMELEN+1];
	char outfile[FILENAMELEN+1];
	int nseqs;
	int  i,j,l,n,ires;
	int err,ix,ntot;
	float min_nn,nn;
	float tmp,qpw;
	Boolean eof,found;

	if(argc!=3 && argc!=4) {
		usage(argv[0]);
		return 0;
	}

	strcpy(infile,argv[1]);
	strcpy(outfile,argv[2]);

/* open the matrix file */
	verbose=FALSE;

	if(argc==3) {
		get_default_matrix();
	}
	else {
        	if((ifd=fopen(argv[3],"r"))==NULL) {
            	fprintf(stderr,"Cannot open matrix file [%s]",argv[3]);
            	return 0;
        	}
		err=readmatrix(ifd);
		if(err<=0) {
			fprintf(stderr,"Error: bad matrix in %s\n",argv[3]);
			return 0;
		}
	}

        init_options(&opt);

        (*opt.alnout_opt).output_clustal=FALSE;
        (*opt.alnout_opt).output_relacs=TRUE;

/* read in the sequences */
        seq_input(infile,opt.explicit_type,FALSE,&mult_aln);
        if(mult_aln.nseqs<=0) {
                error("No sequences in %s\n",infile);
                exit(1);
        }
        nseqs=mult_aln.nseqs;

/* remove the gaps */
	seqlength=0;
	useqlen_array=(int *)ckalloc((nseqs+1)*sizeof(int));
	for(i=0;i<nseqs;i++) {
		if(mult_aln.seqs[i].len>seqlength) seqlength=mult_aln.seqs[i].len;
		l=0;
		for(j=0;j<mult_aln.seqs[i].len;j++)
			if(isalpha(mult_aln.seqs[i].data[j])) {
				l++;
			}
		useqlen_array[i]=l;
	}
        maxlen=0;
        for(i=0;i<nseqs;i++)
                if(useqlen_array[i]>maxlen) maxlen=useqlen_array[i];
        minlen=10000;
        for(i=0;i<nseqs;i++)
                if(useqlen_array[i]<minlen) minlen=useqlen_array[i];
	

/* calculate some simple statistics */
        pcid=(float **)ckalloc((nseqs+1)*sizeof(float *));
        for(i=0;i<nseqs;i++)
                pcid[i]=(float *)ckalloc((nseqs+1)*sizeof(float));
        for(i=0;i<nseqs;i++) {
		for(j=i+1;j<nseqs;j++) {
			pcid[j][i]=pcid[i][j]=pcidentity(mult_aln,i,j);
		}
	}

/* find the nearest neighbor for each sequence */
	min_nn=1.0;
	for(i=0;i<nseqs;i++) {
		nn=0.0;
		for(j=0;j<nseqs;j++) {
			if(i!=j && pcid[i][j]>nn) nn=pcid[i][j];
		}
		if(nn<min_nn) min_nn=nn;
	}

        seqweight=(float **)ckalloc((nseqs+1)*sizeof(float *));
        for(i=0;i<nseqs;i++)
                seqweight[i]=(float *)ckalloc((nseqs+1)*sizeof(float));
        for(i=0;i<nseqs;i++) 
		for(j=i;j<nseqs;j++) 
			seqweight[j][i]=seqweight[i][j]=1.0-pcid[i][j];


	fragment=(Boolean *)ckalloc((nseqs+1)*sizeof(Boolean));


/* calculate pairwise alignment scores using k-tuple scores */
        for(i=0,n=0;i<nseqs;i++)
                for(j=i+1;j<nseqs;j++) {
			qpw=100.0*pcid[i][j];
			if(qpw>40) {
				tmp=(float)useqlen_array[i]/(float)useqlen_array[j];
				if(tmp<0.8) fragment[i]=TRUE;
				else if(tmp>1.25) fragment[j]=TRUE;
			}
			n++;
		}

	/*if(verbose)
        for(i=0;i<nseqs;i++)
		if(fragment[i]) fprintf(stdout,"%s fragment %s\n",argv[1],names[i]);*/

	ngroups=0;
	for(i=0;i<nseqs;i++) {
		if(mult_aln.seqs[i].simgroup>ngroups) ngroups=mult_aln.seqs[i].simgroup;
	}
	seqgroup=(int *)ckalloc((nseqs+1)*sizeof(int));
	for(i=0;i<nseqs;i++) {
		seqgroup[i]=mult_aln.seqs[i].simgroup-1;
	}

/* remove any column score data that already exists in the input file */
	for(i=0;i<mult_aln.ncol_scores;i++)
		ckfree(mult_aln.col_score[i].data);

	ix=0;
	for(i=0;i<ngroups;i++) {
		ntot=0;
        	for(j=0;j<nseqs;j++) 
			if(seqgroup[j]==i) ntot++;
		if(ntot>0) {
			mult_aln.col_score[ix].data=(sint *)ckalloc((seqlength+1)*sizeof(sint));
			calc_submd(&mult_aln,i,ix);
			ix++;
		}
	}
	mult_aln.ncol_scores=ix;
		
        for(i=0;i<nseqs;i++)
                ckfree(seqweight[i]);
        ckfree(seqweight);
	ckfree(fragment);
	ckfree(useqlen_array);

/* write out the sequences */
	strcpy(opt.alnout_opt->relacs_outname,outfile);
        if(!open_alignment_output(outfile,opt.alnout_opt)) exit(1);
        create_alignment_output(mult_aln,*opt.alnout_opt);


	return 0;
}

void calc_submd(ALNPTR mult_aln,int group,int colix)
{
	int i,j,p,n,ntot;
	int nseqs;
	float **weight;

	nseqs=mult_aln->nseqs;
	ntot=0;
        for(i=0;i<nseqs;i++) 
		if(seqgroup[i]==group) ntot++;
	if(ntot==1) {
		col=1.0;
    		for(p=0; p<seqlength; p++) 
			mult_aln->col_score[colix].data[p]=100;
		return;
	}

        weight=(float **)ckalloc((nseqs+1)*sizeof(float *));
        for(i=0;i<nseqs;i++)
		if(seqgroup[i]==group)
                weight[i]=(float *)ckalloc((nseqs+1)*sizeof(float));


	totweight=0;
        for(i=0;i<nseqs;i++) 
		if(seqgroup[i]==group)
		for(j=i+1;j<nseqs;j++) 
			if(seqgroup[j]==group)
				totweight+=seqweight[i][j];

	if (totweight==0) {
        for(i=0;i<nseqs;i++) 
		if(seqgroup[i]==group)
		for(j=i+1;j<nseqs;j++) 
			if(seqgroup[j]==group)
				totweight++;
	}

        for(i=0;i<nseqs;i++) 
		if(seqgroup[i]==group)
		for(j=i+1;j<nseqs;j++) 
			if(seqgroup[j]==group)
			weight[j][i]=weight[i][j]=seqweight[i][j];

/* calculate mean column score */
	col=md_subscore(mult_aln,group,colix,matrix, weight, totweight, fragment);

        for(i=0;i<nseqs;i++) {
		if(seqgroup[i]==group)
                ckfree(weight[i]);
	}
        ckfree(weight);

}


float md_subscore(ALNPTR mult_aln,int group,int colix,float matrix[NUMRES][NUMRES],float **weight,float totweight,Boolean *fragment)
{
	char c,c1;
	int i,j,k,s,s1,p,r;
	int len,is,ie;
	int nseqs;
	float *ntot,ntotseq;
	float mean,n,score;
	float dist,diff,sum;
	float seqvector[26],seqvector1[26];
	float resdist[26][26];
	float resfreq[26][26];

	nseqs=mult_aln->nseqs;
	ntotseq=0;
	for(s=0;s<nseqs;s++)
		if(seqgroup[s]==group)
			ntotseq++;

	ntot=(float *)ckalloc((seqlength+1)*sizeof(float));
    	for(p=0; p<seqlength; p++) {
		ntot[p]=0;
	}
	
/* make ntot the number of sequences at this position, excluding fragments */
	for(s=0;s<nseqs;s++) {
		if(seqgroup[s]==group) {
		len=mult_aln->seqs[s].len;
       		is=0;
		ie=len;
		if(fragment[s]) {
       			for(k=0;k<len;k++) {
               			c = mult_aln->seqs[s].data[k];
               			if (isalpha(c)) {
                       			is=k;
                       			break;
               			}
       			}
       			for(k=len-1;k>=0;k--) {
               			c = mult_aln->seqs[s].data[k];
               			if (isalpha(c)) {
                       			ie=k;
                       			break;
               			}
       			}
		}
		for(p=is;p<=ie;p++) ntot[p]++;
		}
	}

	for(i=0;i<26;i++) {
                for (r=0;r<26; r++)
                        seqvector[r]=matrix[r][i];
                resdist[i][i]=0.0;
                for(j=i+1;j<26;j++) {
                        for (r=0;r<26; r++)
                                seqvector1[r]=matrix[r][j];
                        resdist[i][j]=0.0;
                        for(r=0;r<26;r++) {
                                diff=seqvector1[r]-seqvector[r];
                                resdist[i][j]+=diff*diff;
                        }
                        resdist[i][j]=sqrt((double)resdist[i][j]);
                        resdist[j][i]=resdist[i][j];
                }
        }
	sum=0.0;
    	for(p=0; p<seqlength; p++)
	{
		if(ntot[p]<1) continue;
    		for(i=0;i<26;i++)
    			for(j=0;j<26;j++)
				resfreq[i][j]=0.0;
/* calculate mean of seq distances */
		mean=0.0;
    		for(s=0; s<nseqs; s++) {
		if(seqgroup[s]==group) {
			if(isalpha(mult_aln->seqs[s].data[p])) {
				c=toupper(mult_aln->seqs[s].data[p]);
    				for(s1=s+1; s1<nseqs; s1++)
				{
					if(isalpha(mult_aln->seqs[s1].data[p])) {
						c1=toupper(mult_aln->seqs[s1].data[p]);
						resfreq[c-'A'][c1-'A']+=weight[s][s1];
					}
				}
			}
		}
		}
    		for(i=0;i<26;i++) {
    			for(j=0;j<26;j++) {
				mean+=resdist[i][j]*resfreq[i][j];
			}
		}
		/*mean/=(float)(ntot)*(float)(ntot-1)/2.0;*/
		mean/=totweight;

/* normalise score between 0 and 1 (where 1 is an identical column) */

		score=exp((double)(-mean)/(double)2.0);

/* normalise the score for the number of sequences with residues at this position */
		n=0;
    		for(s=0; s<nseqs; s++)
		if(seqgroup[s]==group) 
			if(isalpha(mult_aln->seqs[s].data[p]) )
				n++;
		score=normalise_subscore(score,n,ntot[p],ntotseq);
/*fprintf(stdout,"%d %f %f %f %f\n",p+1,n,ntot[p],ntotseq,score);*/
		mult_aln->col_score[colix].data[p]=100*score;
		sum+=score;
	}
	mult_aln->col_score[colix].length=p;
	mult_aln->col_score[colix].name=(char *)ckalloc(10*sizeof(char));
	mult_aln->col_score[colix].owner=(char *)ckalloc(10*sizeof(char));
	sprintf(mult_aln->col_score[colix].name,"%s%d","normd_",group+1);
	sprintf(mult_aln->col_score[colix].owner,"%s%d","group_",group+1);
	ckfree(ntot);

	return sum;
}

static float normalise_subscore(float score,float n,float ntot,float ntotseq)
{
        float ret;

        if(n==0) ret=0.0;
        else
                ret=score*(float)exp((double)(-10.0*(float)(ntot-n)/((float)(ntot))));

        return ret;

}

