/* NORMD version 1.2 sept 2005 */
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "score.h"

void calc_md(ALNPTR mult_aln);
void copy_msf_tmsf(ALN mult_aln);
void cut_block(ALN mult_aln,int is, int ie);
void score_block(ALN mult_aln,int block);

int *use_seq;
int seqlength;
int *useqlen_array;
int useqlength;
Boolean *fragment;

char **tseq_array;
int *tseqlen_array;
int tseqlength;

float matrix[NUMRES][NUMRES];
float go;
float ge;

double **tmat;
float **seqweight;
float **gop;
float **gep;

float norm_md,col;
float max_colscore;
int maxlen;
int query;

int blocklen;

int *seqgroup;
int *groupseed;
int ngroups;
int norphans;

float **pcid;

int main(int argc, char **argv)
{
	FILE *ifd;
	char infile[FILENAMELEN+1];
	ALN mult_aln;
	OPT opt;
	int nseqs;
	int  i,j,l,n,ires;
	int err,ix;
	Boolean eof,found;
	int nblocks,block,is,ie;

	if(argc!=2 && argc!=6) {
		fprintf(stderr,"%s Version 1.2\n",argv[0]);
        	fprintf(stderr,"Usage: %s aln_file \n",argv[0]);
        	fprintf(stderr,"   OR: \n");
		fprintf(stderr,"Usage: %s aln_file matrix gop gep window_length\n",argv[0]);
		return 0;
	}
	strcpy(infile,argv[1]);

        if(argc==2) {
                get_default_matrix();
                go=0.0;
                ge=0.1;
                blocklen=1;
        }
        else {

/* open the matrix file */

        	if((ifd=fopen(argv[2],"r"))==NULL) {
            		fprintf(stderr,"Cannot open matrix file [%s]",argv[2]);
            		return 0;
        	}
		err=readmatrix(ifd);
		if(err<=0) {
			fprintf(stderr,"Error: bad matrix in %s\n",argv[3]);
			exit;
		}

		go=atof(argv[3]);
		ge=atof(argv[4]);

		blocklen=atoi(argv[5]);
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

	seqlength=0;
        for(i=0;i<nseqs;i++) {
		if(mult_aln.seqs[i].len>seqlength) seqlength=mult_aln.seqs[i].len;
	}
        tseq_array=(char **)ckalloc((nseqs+1)*sizeof(char *));
        tseqlen_array=(int *)ckalloc((nseqs+1)*sizeof(int));
        for(i=0;i<nseqs;i++) {
                tseq_array[i]=(char *)ckalloc((seqlength+1)*sizeof(char));
        }

	copy_msf_tmsf(mult_aln);

	if(blocklen-blocklen/2==1) 
		nblocks=seqlength/blocklen;
	else 
		nblocks=seqlength/blocklen+1;
		
fprintf(stdout,"getting %d\n",nblocks);
	norm_md=0.0;
	for (block=0;block<blocklen/2;block++) 
        	fprintf(stdout,"%d %.3f \n",block+1,norm_md);

	for (block=blocklen/2;block<tseqlength-blocklen/2;block++) {
		is=block-blocklen/2;
		ie=is+blocklen;
		if(ie>tseqlength) is=tseqlength;

		cut_block(mult_aln,is,ie);
		score_block(mult_aln,block);

	}

	norm_md=0.0;
	for (block=tseqlength-blocklen/2;block<tseqlength;block++) 
        	fprintf(stdout,"%d %.3f \n",block+1,norm_md);
}

void score_block(ALN mult_aln,int block)
{
	int i,j,l,n;
	int nseqs;
	float tmp,qpw_id;
	Boolean found;

/* remove the gaps */
	nseqs=mult_aln.nseqs;
        useqlen_array=(int *)ckalloc((nseqs+1)*sizeof(int));
	useqlength=0;
	for(i=0;i<nseqs;i++) {
		l=0;
		for(j=0;j<mult_aln.seqs[i].len;j++)
			if(isalpha(mult_aln.seqs[i].data[j])) {
				l++;
			}
		useqlen_array[i]=l;
		if (l>useqlength) useqlength=l;
	}

	
/* calculate some simple statistics */
	maxlen=0;
        tmat=(double **)ckalloc((nseqs+2)*sizeof(double *));
        for(i=0;i<=nseqs;i++)
                tmat[i]=(double *)ckalloc((nseqs+2)*sizeof(double));
	pcid=(float **)ckalloc((nseqs+1)*sizeof(float *));
        for(i=0;i<nseqs;i++)
                pcid[i]=(float *)ckalloc((nseqs+1)*sizeof(float));
        for(i=0;i<nseqs;i++) {
                l=useqlen_array[i];
                if(l>maxlen) maxlen=l;
		tmat[i+1][i+1]=0.0;
		for(j=i+1;j<nseqs;j++) {
			pcid[j][i]=pcid[i][j]=pcidentity(mult_aln,i,j);
			tmat[i+1][j+1]=1.0-pcid[i][j];
			tmat[j+1][i+1]=tmat[i+1][j+1];
		}
	}

	use_seq=(int *)ckalloc((nseqs+1)*sizeof(int));
	for(i=0;i<nseqs;i++)
		use_seq[i]=2;

	seqgroup=(int *)ckalloc((nseqs+1)*sizeof(int));
        groupseed=(int *)ckalloc((nseqs+1)*sizeof(int));

        calc_groups(query,0.6,nseqs,pcid,seqgroup,groupseed);

	for(j=0;j<nseqs;j++) use_seq[j]=(-1);
	for(i=0;i<ngroups;i++) {
                j=groupseed[i];
                use_seq[j]=2;
        }
	ckfree(seqgroup);
	ckfree(groupseed);

        seqweight=(float **)ckalloc((nseqs+1)*sizeof(float *));
        for(i=0;i<nseqs;i++)
                seqweight[i]=(float *)ckalloc((nseqs+1)*sizeof(float));

        for(i=0;i<nseqs;i++) 
		for(j=i;j<nseqs;j++) {
			seqweight[j][i]=seqweight[i][j]=tmat[i+1][j+1];
		}

/* calculate the scores for the gaps */
        gop=(float **)ckalloc((nseqs+1)*sizeof(float *));
        for(i=0;i<nseqs;i++)
                gop[i]=(float *)ckalloc((nseqs+1)*sizeof(float));
        gep=(float **)ckalloc((nseqs+1)*sizeof(float *));
        for(i=0;i<nseqs;i++)
                gep[i]=(float *)ckalloc((nseqs+1)*sizeof(float));

       	for(i=0,n=0;i<nseqs;i++) 
              	for(j=i+1;j<nseqs;j++)
			score_gaps(mult_aln,i,j);


	fragment=(Boolean *)ckalloc((nseqs+1)*sizeof(Boolean));
        for(i=0,n=0;i<nseqs;i++)
                for(j=i+1;j<nseqs;j++) {
                        qpw_id=100.0*pcid[i][j];
                        if(qpw_id>40) {
                                tmp=(float)useqlen_array[i]/(float)useqlen_array[j];
                                if(tmp<0.8) fragment[i]=TRUE;
                                else if(tmp>1.25) fragment[j]=TRUE;
                        }
                        n++;
                }

	calc_md(&mult_aln);

        fprintf(stdout,"%d %.3f \n",block+1,norm_md);

        for(i=0;i<nseqs;i++)
                ckfree(gop[i]);
        ckfree(gop);
        for(i=0;i<nseqs;i++)
                ckfree(gep[i]);
        ckfree(gep);

        for(i=0;i<=nseqs;i++)
                ckfree(tmat[i]);
        ckfree(tmat);
        for(i=0;i<nseqs;i++)
                ckfree(pcid[i]);
        ckfree(pcid);

	ckfree(use_seq);
	ckfree(fragment);
        for(i=0;i<nseqs;i++)
                ckfree(seqweight[i]);
        ckfree(seqweight);

	ckfree(useqlen_array);

}

void copy_msf_tmsf(ALN mult_aln) 
{
	int i,j;

	for(i=0;i<mult_aln.nseqs;i++) {
		for(j=0;j<seqlength;j++) 
			tseq_array[i][j]=mult_aln.seqs[i].data[j];
		tseqlen_array[i]=mult_aln.seqs[i].len;
	}
	tseqlength=seqlength;
}

void cut_block(ALN mult_aln,int is, int ie)
{
	int i,j,len;

	len=ie-is;
	for(i=0;i<mult_aln.nseqs;i++) {
		for(j=is;j<ie;j++) 
			mult_aln.seqs[i].data[j-is]=tseq_array[i][j];
		mult_aln.seqs[i].len=len;
	}
	seqlength=len;


}

void calc_md(ALNPTR mult_aln)
{
	int i,j,n,n1,ntot;
	int nseqs;
	float id,meanid;
	float totweight;
	float gap_openscore,gap_extscore;
	float ll;
	float medianqpw;
	float tl,q1,q3;
	float tmp;
	float **weight;
	float *qpw;

	nseqs=mult_aln->nseqs;
	ntot=0;
        for(i=0;i<nseqs;i++) 
		if(use_seq[i]>0) ntot++;
	if(ntot==1) {
		norm_md=1.0;
		return;
	}

        weight=(float **)ckalloc((nseqs+1)*sizeof(float *));
        for(i=0;i<nseqs;i++)
                weight[i]=(float *)ckalloc((nseqs+1)*sizeof(float));

	totweight=0;
        for(i=0;i<nseqs;i++) 
		if(use_seq[i]>0) 
		for(j=i+1;j<nseqs;j++) 
			if(use_seq[j]>0) 
				totweight+=seqweight[i][j];

	if (totweight==0) {
        for(i=0;i<nseqs;i++) 
		if(use_seq[i]>0)
		for(j=i+1;j<nseqs;j++) 
			if(use_seq[j]>0)
				totweight++;
	}

        for(i=0;i<nseqs;i++) 
		if(use_seq[i]>0)
		for(j=i+1;j<nseqs;j++) 
			if(use_seq[j]>0) {
			weight[j][i]=weight[i][j]=seqweight[i][j];
			}

/* calculate pairwise alignment scores using k-tuple scores */
	qpw=(float *)ckalloc((nseqs*nseqs+1)*sizeof(float));
        for(i=0,n=0;i<nseqs;i++)
		if(use_seq[i]>0)
                for(j=i+1;j<nseqs;j++) 
			if(use_seq[j]>0) {
                        qpw[n]=show_pair(*mult_aln,i,j);
                        n++;
                }


	gap_openscore=0.0;
	gap_extscore=0.0;
       	for(i=0,n=0;i<nseqs;i++) 
		if(use_seq[i]>0) {
              		for(j=i+1;j<nseqs;j++) {
				if(use_seq[j]>0) {
				gap_openscore+=gop[i][j]*weight[i][j];
				gap_extscore+=gep[i][j]*weight[i][j];
				n++;
				}
			}
		}
        gap_openscore/=(float)totweight;
        gap_extscore/=(float)totweight;


/* sort the pairwise k-tuple scores into ascending order */
        sort_scores(qpw,0,n-1);

/* find the lower quartile range of the pairwise k-tuple scores */
        if(n == 0)
                medianqpw = 0;
        else if(n % 2 == 0)
                medianqpw=(qpw[n/2-1]+qpw[n/2])/2.0;
        else
                medianqpw=qpw[n/2];

        if(n==0)
                ll=0;
        else {
                tl = (float)n/4.0 + 0.5;
                if(tl - (int)tl == 0.5) {
                        q3=(qpw[(int)tl]+qpw[(int)tl+1])/2.0;
                        q1=(qpw[n-(int)tl]+qpw[n-(int)tl-1])/2.0;
                }
                else if(tl - (int)tl > 0.5) {
                        q3=qpw[(int)tl+1];
                        q1=qpw[n-(int)tl-1];
                }
                else {
                        q3=qpw[(int)tl];
                        q1=qpw[n-(int)tl];
                }
/* use lower quartile range if more than 10 sequences, otherwise use mean */
        }
	if(q1<10) q1=10;
	q1=1.14*q1-6.4;

/* calculate mean column score */
	max_colscore=max_score(nseqs,maxlen,fragment);
	col=md_score(mult_aln,0,matrix,weight,totweight,fragment);
	if(blocklen==1) norm_md=col;
	else {
        	if(max_colscore>0.0) {
			norm_md=100.0*(col-gap_openscore*go-gap_extscore*ge)/(max_colscore*q1);
		}
        	else norm_md=0.0;
		if (norm_md<0) norm_md=0;
	}

        for(i=0;i<nseqs;i++)
                ckfree(weight[i]);
        ckfree(weight);
        ckfree(qpw);


}
