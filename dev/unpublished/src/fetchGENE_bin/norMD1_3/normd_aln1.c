/* NORMD version 1.2 Sept 2005 */

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "score.h"

void calc_md(ALNPTR mult_aln,int colix);
float set_mdcutoff(int ntot,float q1);

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

float tl,q1,q3;
float medianqpw;
float *qpw;
float **qpw_id;

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
	fprintf(stderr,"Usage: %s in_file outfile matrix gop gep endgap [-v]     | calculate norMD for alignment (-v verbose)\n",prog);
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
	float tmp;
	Boolean eof,found;

	if(argc!=3 && argc!=7 && argc!=8) {
		usage(argv[0]);
		return 0;
	}

	strcpy(infile,argv[1]);
	strcpy(outfile,argv[2]);

/* open the matrix file */
	verbose=FALSE;

	if(argc==3) {
		get_default_matrix();
		go=0.0;
		ge=0.1;
		egap=0.0;
	}
	else {
		if(argc==8) verbose=TRUE;
	

        	if((ifd=fopen(argv[3],"r"))==NULL) {
            	fprintf(stderr,"Cannot open matrix file [%s]",argv[3]);
            	return 0;
        	}
		err=readmatrix(ifd);
		if(err<=0) {
			fprintf(stderr,"Error: bad matrix in %s\n",argv[3]);
			return 0;
		}

		go=atof(argv[4]);
		ge=atof(argv[5]);
		egap=atof(argv[6]);
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
	
/* remove any column score data that already exists in the input file */
	/*if (mult_aln.ncol_scores==1)
		ckfree(mult_aln.col_score[0].data);*/
	ix=mult_aln.ncol_scores;
	mult_aln.col_score[ix].data=(sint *)ckalloc((seqlength+1)*sizeof(sint));
	mult_aln.ncol_scores=ix+1;

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
		for(j=i;j<nseqs;j++) {
			seqweight[j][i]=seqweight[i][j]=1.0-pcid[i][j];
		}


	fragment=(Boolean *)ckalloc((nseqs+1)*sizeof(Boolean));


/* calculate pairwise alignment scores using k-tuple scores */
        qpw_id=(float **)ckalloc((nseqs+1)*sizeof(float *));
        for(i=0;i<nseqs;i++)
                qpw_id[i]=(float *)ckalloc((nseqs+1)*sizeof(float));
        for(i=0,n=0;i<nseqs;i++)
                for(j=i+1;j<nseqs;j++) {
			qpw_id[i][j]=100.0*pcid[i][j];
                        if(qpw_id[i][j]<60) {
				qpw_id[i][j]=show_pair(mult_aln,i,j);
			}
			if(qpw_id[i][j]>40) {
				tmp=(float)useqlen_array[i]/(float)useqlen_array[j];
				if(tmp<0.8) fragment[i]=TRUE;
				else if(tmp>1.25) fragment[j]=TRUE;
			}
			n++;
		}

	/*if(verbose)
        for(i=0;i<nseqs;i++)
		if(fragment[i]) fprintf(stdout,"%s fragment %s\n",argv[1],names[i]);*/

/* calculate sequence groups and keep first sequence in each group for processing */
        use_seq=(int *)ckalloc((nseqs+1)*sizeof(int));
	for(i=0;i<nseqs;i++)
		use_seq[i]=2;

	query=0;
	seqgroup=(int *)ckalloc((nseqs+1)*sizeof(int));
	groupseed=(int *)ckalloc((nseqs+1)*sizeof(int));
        calc_groups(query,0.7,nseqs,pcid,seqgroup,groupseed);

	if(ngroups<=0) {
		fprintf(stderr,"Error: problem with sequence grouping\n");
		exit(1);
	}
        for(j=0;j<nseqs;j++) use_seq[j]=(-1);
        for(i=0;i<ngroups;i++) {
		j=groupseed[i];
		use_seq[j]=2;
        }


        for(i=0;i<nseqs;i++)
                ckfree(pcid[i]);
        ckfree(pcid);
	ckfree(groupseed);
	ckfree(seqgroup);


	qpw=(float *)ckalloc((ngroups*ngroups+1)*sizeof(float));

        for(i=0,n=0;i<nseqs;i++)
		if(use_seq[i]>1) 
                for(j=i+1;j<nseqs;j++)
			if(use_seq[j]>1)
                        qpw[n++]=qpw_id[i][j];
        for(i=0;i<nseqs;i++)
                ckfree(qpw_id[i]);
        ckfree(qpw_id);

/* sort the pairwise k-tuple scores into ascending order */
        sort_scores(qpw,0,n-1);

/* calculate the scores for the gaps */
        gop=(float **)ckalloc((nseqs+1)*sizeof(float *));
        for(i=0;i<nseqs;i++)
                gop[i]=(float *)ckalloc((nseqs+1)*sizeof(float));
        gep=(float **)ckalloc((nseqs+1)*sizeof(float *));
        for(i=0;i<nseqs;i++)
                gep[i]=(float *)ckalloc((nseqs+1)*sizeof(float));
        egp=(float **)ckalloc((nseqs+1)*sizeof(float *));
        for(i=0;i<nseqs;i++)
                egp[i]=(float *)ckalloc((nseqs+1)*sizeof(float));

       	for(i=0;i<nseqs;i++) 
		if(use_seq[i]>1)
              	for(j=i+1;j<nseqs;j++)
			if(use_seq[j]>1)
			score_gaps(mult_aln,i,j);

	calc_md(&mult_aln,ix);
	for(i=0,ntot=0;i<nseqs;i++)
		if(use_seq[i]>1) {
			ntot++;
		}
	tmp=set_mdcutoff(ntot,q1);
	normd_rs/=tmp;
	tmp=1.0;
	norm_md/=tmp;

	mult_aln.alnscore=norm_md;
	mult_aln.validalnscore=TRUE;

	if(!verbose) {
        	fprintf(stdout,"%s\t%.3f\n",argv[1],norm_md);
        	/*fprintf(stdout,"%.3f\n",norm_md);*/
	} else {
		/*fprintf(stdout,"%s %.3f %.3f %.3f %.3f %.3f %.3f %d %d %d\n",
		argv[1],norm_md,normd_rs,col,max_colscore,gap_extscore*0.1,q1,nseqs,minlen,maxlen);*/
		fprintf(stdout,"FILE  %s\n",argv[1]);
        	fprintf(stdout,"norMD    %.3f\n",norm_md);
        	/*fprintf(stdout,"norMD_of %.3f\n",norm_md);
        	fprintf(stdout,"norMD_rs %.3f\n",normd_rs);*/
        	fprintf(stdout,"NSEQS    %d\n",nseqs);
        	fprintf(stdout,"MD       %.3f\n",col);
        	fprintf(stdout,"maxMD    %.3f\n",max_colscore);
        	fprintf(stdout,"GOP      %.3f\n",gap_openscore);
        	fprintf(stdout,"GEP      %.3f\n",gap_extscore);
        	fprintf(stdout,"LQR      %.3f\n",q1);
	}

        for(i=0;i<nseqs;i++)
                ckfree(gop[i]);
        ckfree(gop);
        for(i=0;i<nseqs;i++)
                ckfree(gep[i]);
        ckfree(gep);
        for(i=0;i<nseqs;i++)
                ckfree(egp[i]);
        ckfree(egp);

        for(i=0;i<nseqs;i++)
                ckfree(seqweight[i]);
        ckfree(seqweight);
	ckfree(qpw);
	ckfree(use_seq);
	ckfree(fragment);
	ckfree(useqlen_array);

/* write out the sequences */
	strcpy(opt.alnout_opt->relacs_outname,outfile);
        if(!open_alignment_output(outfile,opt.alnout_opt)) exit(1);
        create_alignment_output(mult_aln,*opt.alnout_opt);


	return 0;
}

float set_mdcutoff(int ntot,float q1)
{
        float n,tmp;

        if(ntot>15) n=15;
        else n=ntot;
        tmp=(300.0)/((n*n)*4.0+100.0);
        return tmp;
}


void calc_md(ALNPTR mult_aln,int colix)
{
	int i,j,n,ntot;
	int nseqs;
	float **weight;

	nseqs=mult_aln->nseqs;
	ntot=0;
        for(i=0;i<nseqs;i++) 
		if(use_seq[i]>1) ntot++;
	if(ntot==1) {
		col=1.0;
		max_colscore=1.0;
		norm_md=1.0;
		normd_rs=100.0;
		return;
	}
/*fprintf(stdout,"nseqs %d ntot %d\n",nseqs,ntot);*/

        weight=(float **)ckalloc((nseqs+1)*sizeof(float *));
        for(i=0;i<nseqs;i++)
                weight[i]=(float *)ckalloc((nseqs+1)*sizeof(float));


	totweight=0;
        for(i=0;i<nseqs;i++) 
		if(use_seq[i]>1) 
		for(j=i+1;j<nseqs;j++) 
			if(use_seq[j]>1) 
				totweight+=seqweight[i][j];

	if (totweight==0) {
        for(i=0;i<nseqs;i++) 
		if(use_seq[i]>1)
		for(j=i+1;j<nseqs;j++) 
			if(use_seq[j]>1)
				totweight++;
	}

        for(i=0;i<nseqs;i++) 
		if(use_seq[i]>1)
		for(j=i+1;j<nseqs;j++) 
			if(use_seq[j]>1) {
			weight[j][i]=weight[i][j]=seqweight[i][j];
			}

	gap_openscore=0.0;
	gap_extscore=0.0;
	end_gapscore=0.0;
       	for(i=0,n=0;i<nseqs;i++) 
		if(use_seq[i]>1) {
              		for(j=i+1;j<nseqs;j++) {
				if(use_seq[j]>1) {
				gap_openscore+=gop[i][j]*weight[i][j];
				gap_extscore+=gep[i][j]*weight[i][j];
				end_gapscore+=egp[i][j]*weight[i][j];
				n++;
				}
			}
		}
        gap_openscore/=(float)totweight;
        gap_extscore/=(float)totweight;
        end_gapscore/=(float)totweight;

/* calculate pairwise alignment scores using k-tuple scores */
        for(i=0,n=0;i<nseqs;i++)
		if(use_seq[i]>1)
                	for(j=i+1;j<nseqs;j++) 
                        	if(use_seq[j]>1) 
                        		n++;


/* find the lower quartile range of the pairwise k-tuple scores */
        if(n == 0)
                medianqpw = 0;
        else if(n % 2 == 0)
                medianqpw=(qpw[n/2-1]+qpw[n/2])/2.0;
        else
                medianqpw=qpw[n/2];

        if(n==0)
                q1=0;
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
        }
	q1=(1.14*q1-6.4)/100.0;

/* calculate mean column score */
	max_colscore=max_score(nseqs,maxlen, fragment);
	col=md_score(mult_aln,colix,matrix, weight, totweight, fragment);
	strcpy(mult_aln->col_score[mult_aln->ncol_scores-1].name,"group_all");

       	if(max_colscore>0) {
		norm_md=(col-gap_openscore*go-gap_extscore*ge-end_gapscore*egap)/(max_colscore*q1);
		normd_rs=100.0*norm_md*q1;
	}
       	else {
		norm_md=0.0;
		normd_rs=0.0;
	}

        for(i=0;i<nseqs;i++) {
                ckfree(weight[i]);
	}
        ckfree(weight);

}
