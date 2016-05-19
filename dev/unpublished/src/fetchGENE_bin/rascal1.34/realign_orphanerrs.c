#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "rascal.h"

#define MIN(a,b) ((a)<(b)?(a):(b))

typedef struct {
        sint seq;
        char name[MAXNAMES+1];
	sint first_col;
	sint last_col;
	sint first_res;
	sint last_res;
} SEQERR,*SEQERR_PTR;

static sint count_seqerrs(char *filename);
static void read_seqerrs(char *filename,ALNPTR mult_aln,SEQERR *seqerrs);
static void realign_seq(ALNPTR mult_aln,double **tmat,sint iseq,sint iseqlen,sint *secgroup,sint first_res,sint last_res,OPT opt);
static void align_seq_to_groups(ALNPTR mult_aln,double **tmat,sint *secgroup,sint iseq,MULT_OPT mult_opt,ALNCOUNT *alncount);
static void remove_gaps(ALNPTR mult_aln,sint *group,sint iseq);
static float countpcid(SEQ seq1,SEQ seq2);

int main(int argc,char **argv)
{
	int i,j,s;
	char infile[FILENAMELEN+1];
	char errfile[FILENAMELEN+1];
	char outfile[FILENAMELEN+1];
	ALN mult_aln;
	OPT opt;
	SEQERR *seqerrs;
	sint nseqerrs;
	sint nseqs,len,status;
	sint *secgroup;
	sint *orggroup;
	double dscore,**tmat;
	sint *is, *ie;
	char c;



	if(argc!=4) {
		fprintf(stdout,"Usage: %s input_aln err_file output_aln\n",argv[0]);
		exit(1);
	}
	strcpy(infile,argv[1]);
	strcpy(errfile,argv[2]);
	strcpy(outfile,argv[3]);

        init_options(&opt);

	(*opt.alnout_opt).output_clustal=FALSE;
	(*opt.alnout_opt).output_gcg=TRUE;

/* read in the sequences */
	seq_input(infile,opt.explicit_type,FALSE,&mult_aln);
	if(mult_aln.nseqs<=0) {
		error("No sequences in %s\n",infile);
		exit(1);
	}
	nseqs=mult_aln.nseqs;

/* read in the groups */
	secgroup=(sint *)ckalloc((nseqs+1)*sizeof(sint));
	orggroup=(sint *)ckalloc((nseqs+1)*sizeof(sint));
	status=get_groups(errfile,&mult_aln,secgroup,orggroup);
	ckfree(orggroup);
	if(status<=0) {
	strcpy((*opt.alnout_opt).gcg_outname, outfile);
	for (i=0;i<mult_aln.nseqs;i++) mult_aln.seqs[i].output_index = i;

	if(!open_alignment_output(infile,opt.alnout_opt)) exit(1);
        create_alignment_output(mult_aln,*opt.alnout_opt);
		exit(1);
	}

        for (i=0;i<mult_aln.nseqs;i++) 
		mult_aln.seqs[i].simgroup=secgroup[i];

	nseqerrs=count_seqerrs(errfile);
	if(nseqs <=6 || nseqerrs<=0) {
		fprintf(stdout,"No orphan errors in file\n");
	strcpy((*opt.alnout_opt).gcg_outname, outfile);
	for (i=0;i<mult_aln.nseqs;i++) mult_aln.seqs[i].output_index = i;

	if(!open_alignment_output(infile,opt.alnout_opt)) exit(1);
        create_alignment_output(mult_aln,*opt.alnout_opt);
		exit(1);
	}
	seqerrs=(SEQERR *)ckalloc((nseqerrs+1) * sizeof(SEQERR));
	read_seqerrs(errfile,&mult_aln,seqerrs);


/* count pairwise residue percent identities */
        tmat = (double **) ckalloc( (mult_aln.nseqs+1) * sizeof (double *) );
        for(i=0;i<mult_aln.nseqs;i++)
                tmat[i] = (double *)ckalloc( (mult_aln.nseqs+1) * sizeof (double) );

        for (i=0;i<mult_aln.nseqs;i++) {
                for (j=i+1;j<mult_aln.nseqs;j++) {
			dscore = countpcid(mult_aln.seqs[i],mult_aln.seqs[j]);
                        tmat[j][i] = tmat[i][j] = (100.0 - dscore)/100.0;
                }
        }

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

	/*opt.mult_opt->gap_opt->no_hyd_penalties = TRUE;
	opt.mult_opt->gap_opt->no_pref_penalties = TRUE;*/
	opt.mult_opt->gap_opt->use_endgaps = TRUE;
	opt.mult_opt->gap_opt->nendgappenalties=FALSE;
	opt.mult_opt->gap_opt->cendgappenalties=FALSE;
	opt.mult_opt->gap_opt->gap_dist=0;
	opt.mult_opt->prot_gap_open=7.0;
	opt.mult_opt->prot_gap_extend=0.2;

	for(i=0;i<nseqerrs;i++) {
		len=seqerrs[i].last_res-seqerrs[i].first_res+1;
fprintf(stdout,"Realigning orphan %s : %d %d (%d %d) len %d\n",seqerrs[i].name,seqerrs[i].first_col,seqerrs[i].last_col,seqerrs[i].first_res,seqerrs[i].last_res,len);
		realign_seq(&mult_aln,tmat,seqerrs[i].seq,ie[seqerrs[i].seq],secgroup,seqerrs[i].first_res,seqerrs[i].last_res,opt);
	}

/* write out the sequences */
	strcpy((*opt.alnout_opt).gcg_outname, outfile);
	for (i=0;i<mult_aln.nseqs;i++) mult_aln.seqs[i].output_index = i;

fprintf(stdout,"writing file %s\n",outfile);
	if(!open_alignment_output(infile,opt.alnout_opt)) exit(1);
        create_alignment_output(mult_aln,*opt.alnout_opt);
}


static sint count_seqerrs(char *filename)
{
	FILE *fin;
	char line[MAXLINE+1];
	int nseqerrs;

        if((fin=fopen(filename,"r"))==NULL) {
                error("Could not open file %s",filename);
                return 0;    
        }

	nseqerrs=0;
	while(fgets(line,MAXLINE+1,fin)) {
		if(keyword(line,"ORPHAN_ERROR")) {
			nseqerrs++;
		}
	}
	
	fclose(fin);
	return nseqerrs;
}

static void read_seqerrs(char *filename,ALNPTR mult_aln,SEQERR *seqerrs)
{
	FILE *fin;
	char line[MAXLINE+1];
	int nseqerrs;
	int iseq,f,l,fr,lr;
	char sname[MAXLINE+1];
	char tmp[MAXLINE+1];

        if((fin=fopen(filename,"r"))==NULL) {
                error("Could not open file %s",filename);
                return;    
        }

	nseqerrs=0;
	while(fgets(line,MAXLINE+1,fin)) {
		if(keyword(line,"ORPHAN_ERROR")) {
			sscanf(line,"%s %d %s %d %d \n",tmp,&iseq,sname,&f,&l);
			if(iseq>=1 && iseq<=mult_aln->nseqs && f<l && l<=mult_aln->seqs[iseq-1].len) {
				seqerrs[nseqerrs].seq=iseq-1;
				strcpy(seqerrs[nseqerrs].name,sname);
				seqerrs[nseqerrs].first_col=f-1;
				seqerrs[nseqerrs].last_col=l-1;
				col2pos(mult_aln->seqs[iseq-1].data,f-1,l-1,&fr,&lr);
				seqerrs[nseqerrs].first_res=fr;
				seqerrs[nseqerrs].last_res=lr;
				nseqerrs++;
			}
		}
	}
	
	fclose(fin);
}

static void realign_seq(ALNPTR mult_aln,double **tmat,sint iseq,sint iseqlen,sint *secgroup,sint first_res,sint last_res,OPT opt)
{
	sint i,j;
	sint save_order,length,length1,length2;
	sint maxlen;
	sint first,last;
	ALNCOUNT *alncount;
	sint *tmplen_array;
	char **tmp_array;
	sint *newlen_array;
	char **new_array;

	pos2col(mult_aln->seqs[iseq].data,first_res,last_res,&first,&last);
	if(first>=iseqlen) return;

/* save the alignment into a temporary area */

	maxlen=0;
	for (i=0;i<mult_aln->nseqs;i++)
		if(mult_aln->seqs[i].len>maxlen) maxlen=mult_aln->seqs[i].len;
	if(last>=maxlen) last=maxlen-1;

	tmplen_array=(sint *)ckalloc((mult_aln->nseqs+2) * sizeof(sint));
	tmp_array=(char **)ckalloc((mult_aln->nseqs+2) * sizeof(char *));
	for (i=0;i<mult_aln->nseqs;i++)
	{
		tmplen_array[i]=maxlen;
		tmp_array[i]=(char *)ckalloc((maxlen+2) * sizeof(char));
		for(j=0;j<mult_aln->seqs[i].len;j++)
			tmp_array[i][j]=mult_aln->seqs[i].data[j];
		for(j=mult_aln->seqs[i].len;j<maxlen;j++)
			tmp_array[i][j]=GAP2;
	}

/* copy the selected residue range to the clustal alignment arrays */

	length=last-first+1;
	for (i=0;i<mult_aln->nseqs;i++)
	{
		mult_aln->seqs[i].len=length;
		/*realloc_seq(&mult_aln->seqs[i],length);*/
		ckfree(mult_aln->seqs[i].data);
		mult_aln->seqs[i].data = (char *)ckalloc((length+2) * sizeof (char));
		for(j=first;j<=last;j++)
			mult_aln->seqs[i].data[j-first]=tmp_array[i][j];
		mult_aln->seqs[i].data[j-first]=EOS;
	}

/* for each sequence, calculate the number of residues to be aligned this time, the number of
residues left on n-ter and on c-ter. these will be used in prfalign to set appropriate end gap
penalties */
        alncount=(ALNCOUNT *)ckalloc((mult_aln->nseqs+1) * sizeof(ALNCOUNT));
        for (i=0;i<mult_aln->nseqs;i++) {
                alncount[i].n=alncount[i].a=alncount[i].c=0;
                for(j=0;j<first;j++)
                        if(isalpha(tmp_array[i][j])) alncount[i].n++;
                for(j=first;j<=last;j++)
                        if(isalpha(tmp_array[i][j])) alncount[i].a++;
                for(j=last+1;j<tmplen_array[i];j++)
                        if(isalpha(tmp_array[i][j])) alncount[i].c++;
        }

/* temporarily set the output order to be the same as the input */
	save_order=(*opt.alnout_opt).output_order;
	(*opt.alnout_opt).output_order=INPUT;

/* align the residue range */
        align_seq_to_groups(mult_aln,tmat,secgroup,iseq,*opt.mult_opt,alncount);
	ckfree(alncount);

	(*opt.alnout_opt).output_order=save_order;

/* save the new alignment into another temporary area */
	newlen_array=(sint *)ckalloc((mult_aln->nseqs+2) * sizeof(sint));
	new_array=(char **)ckalloc((mult_aln->nseqs+2) * sizeof(char *));
	for (i=0;i<mult_aln->nseqs;i++)
	{
		newlen_array[i]=mult_aln->seqs[i].len;
		new_array[i]=(char *)ckalloc((mult_aln->seqs[i].len+2) * sizeof(char));
		for(j=0;j<mult_aln->seqs[i].len;j++)
			new_array[i][j]=mult_aln->seqs[i].data[j];
	}

/* paste the realigned range back into the alignment */
	length1=length2=0;
	for (i=0;i<mult_aln->nseqs;i++)
	{
		length1=tmplen_array[i]-length+newlen_array[i];
		length2=newlen_array[i];
		mult_aln->seqs[i].len=length1;
		/*realloc_seq(&mult_aln->seqs[i],length1);*/
		ckfree(mult_aln->seqs[i].data);
		mult_aln->seqs[i].data = (char *)ckalloc((length1+2) * sizeof (char));
		for(j=0;j<first;j++)
			mult_aln->seqs[i].data[j]=tmp_array[i][j];
		for(j=first;j<first+length2;j++)
			mult_aln->seqs[i].data[j]=new_array[i][j-first];
		for(j=first+length2;j<length1;j++)
			mult_aln->seqs[i].data[j]=tmp_array[i][last+j-first-length2+1];
	}
	ckfree(tmplen_array);
	for(i=0;i<mult_aln->nseqs;i++)
		ckfree(tmp_array[i]);
	ckfree(tmp_array);
	ckfree(newlen_array);
	for(i=0;i<mult_aln->nseqs;i++)
		ckfree(new_array[i]);
	ckfree(new_array);

}

static void align_seq_to_groups(ALNPTR mult_aln,double **tmat,sint *secgroup,sint iseq,MULT_OPT mult_opt,ALNCOUNT *alncount)
{
        sint i;
        sint *group;
        lint score = 0;
	Boolean use_maxid=FALSE;
	Boolean norm_gaps=TRUE;

        group = (sint *)ckalloc( (mult_aln->nseqs+1) * sizeof (sint));

        for (i=0;i<mult_aln->nseqs;i++) {
                group[i] = 1;
                if (i==iseq || secgroup[i]>0) 
			mult_aln->seqs[i].weight=(100.0-tmat[i][iseq]);
		else
			mult_aln->seqs[i].weight=0;
	}
        group[iseq] = 2;

/* remove positions that contain just gaps */
	remove_gaps(mult_aln,group,iseq);

        score = prfalign(2,mult_aln,group,tmat,use_maxid,mult_opt,alncount,norm_gaps);

        group=ckfree((void *)group);
}

static void remove_gaps(ALNPTR mult_aln,sint *group,sint iseq)
{
        int i,j,k,l,n,ngaps,maxlen;

	maxlen=n=0;
        for (j=0;j<mult_aln->nseqs;j++) {
                if(group[j]!=2) n++;
		if(mult_aln->seqs[j].len>maxlen) maxlen=mult_aln->seqs[j].len;
	}

/* first remove gap columns from the alignment */
        for (i=0,l=0;i<maxlen;i++) {
                ngaps=0;
                for (j=0;j<mult_aln->nseqs;j++) {
                        if(group[j]!=2) {
				if(l>=mult_aln->seqs[j].len)
					 ngaps++;
				else
					if (!isalpha(mult_aln->seqs[j].data[l])) ngaps++;
			}
		}
                if (ngaps==n) {
                        for (j=0;j<mult_aln->nseqs;j++) {
                        	if(group[j]!=2) {
                                	for(k=l+1;k<=mult_aln->seqs[j].len;k++)
                                        	mult_aln->seqs[j].data[k-1]=mult_aln->seqs[j].data[k];
                                	mult_aln->seqs[j].len--;
				}
                        }
                }
                else l++;
        }
/* then remove all gaps from the sequence to be realigned */
        for (i=0,l=0;i<maxlen;i++) {
		if (!isalpha(mult_aln->seqs[iseq].data[l])) {
                       	for(k=l+1;k<=mult_aln->seqs[iseq].len;k++)
                               	mult_aln->seqs[iseq].data[k-1]=mult_aln->seqs[iseq].data[k];
                       	mult_aln->seqs[iseq].len--;
                }
                else l++;
        }
}

static float countpcid(SEQ seq1,SEQ seq2)
{
   char c1,c2;
   sint i;
   sint count,total,len;
   float score;

   len=MIN(seq1.reslen,seq2.reslen);

   count = total = 0;
   for (i=0;i<seq1.len && i<seq2.len;i++) {
     c1 = seq1.data[i];
     c2 = seq2.data[i];
     if (isalpha(c1) && isalpha(c2)) {
       total++;
       if (c1 == c2) count++;
     }

   }

   if(total==0) score=0;
   else
   score = 100.0 * (float)count / (float)len;
   return(score);

}

