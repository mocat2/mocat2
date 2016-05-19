#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "rascal.h"

#define MIN(a,b) ((a)<(b)?(a):(b))

typedef struct {
        sint first_col; 
        sint last_col; 
	sint nseqs;
} CBLOCK,*CBLOCKPTR;

static sint count_blocks(char *filename);
static void read_blocks(char *filename,ALN mult_aln,CBLOCK *blocks);
static void realign_block(ALNPTR mult_aln,sint first,sint last,OPT opt);
static void align_block(ALNPTR mult_aln,OPT opt,ALNCOUNT *alncount);
static float countpcid(SEQ seq1,SEQ seq2);
char infile[FILENAMELEN+1];

int main(int argc,char **argv)
{
	FILE *tree;
	int i,j,k,n;
	char errfile[FILENAMELEN+1];
	char outfile[FILENAMELEN+1];
	ALN mult_aln;
	OPT opt;
	CBLOCK *blocks;
	int nblocks;
	int nblocks1;
	int nseqs;
	sint len;
	double dscore,**tmat;
	Boolean realign_ends=TRUE;



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

/* read in the global core blocks */
	nblocks=count_blocks(errfile);
	if(nblocks<0) {
		exit(1);
	}
	if(nblocks<=0) {
		fprintf(stdout,"No errors in block file\n");
		strcpy((*opt.alnout_opt).gcg_outname, outfile);
		for (i=0;i<mult_aln.nseqs;i++) mult_aln.seqs[i].output_index = i;

		for (i=0;i<mult_aln.nseqs;i++) mult_aln.seqs[i].simgroup = (-1);

		if(!open_alignment_output(infile,opt.alnout_opt)) exit(1);
        	create_alignment_output(mult_aln,*opt.alnout_opt);
		exit(1);
	}

	blocks=(CBLOCK *)ckalloc((nblocks+1) * sizeof(CBLOCK));
	if(nblocks<=0) {
		blocks[0].first_col=0;
		blocks[0].last_col=mult_aln.seqs[0].len-1;
		opt.mult_opt->neg_matrix=FALSE;
	}
	else {
	read_blocks(errfile,mult_aln,blocks);
	}

	nblocks1=0;
	for(i=0;i<nblocks;i++) {
		n=0;
		for(j=0;j<mult_aln.nseqs;j++) {
			len=0;
			for(k=blocks[i].first_col;k<blocks[i].last_col;k++)
				if(isalpha(mult_aln.seqs[j].data[k])) len++;
			if(len>0) n++;
		}
		blocks[i].nseqs=n;
		if(n>mult_aln.nseqs*0.5) nblocks1++;
	}
	if(nblocks1<=0) {
		fprintf(stdout,"No core blocks in file\n");
		strcpy((*opt.alnout_opt).gcg_outname, outfile);
		for (i=0;i<mult_aln.nseqs;i++) mult_aln.seqs[i].output_index = i;

		for (i=0;i<mult_aln.nseqs;i++) mult_aln.seqs[i].simgroup = (-1);

		if(!open_alignment_output(infile,opt.alnout_opt)) exit(1);
        	create_alignment_output(mult_aln,*opt.alnout_opt);
		exit(1);
	}

/* make a tree from the percent identities in the full alignment */

        tmat=(double **)ckalloc((nseqs+2)*sizeof(double *));
        for(i=0;i<=nseqs;i++)
                tmat[i]=(double *)ckalloc((nseqs+2)*sizeof(double));
        for(i=0;i<nseqs;i++) {
                tmat[i+1][i+1]=0.0;
                for(j=i+1;j<nseqs;j++) {
			dscore = countpcid(mult_aln.seqs[i],mult_aln.seqs[j]);
			tmat[j][i] = tmat[i][j] = (100.0 - dscore)/100.0;
                }
        }

	strcpy(mult_aln.treename,infile);
	strcat(mult_aln.treename,".ph");
        if((tree=fopen(mult_aln.treename,"w"))==NULL) {
            fprintf(stdout,"\nError: Cannot write tree file [%s]",mult_aln.treename);
            exit(1);
        }
	guide_tree(tree,mult_aln.seqs,nseqs,tmat,QUICKNJ);

/* realign errors */

	opt.mult_opt->gap_opt->use_endgaps = FALSE;
	/*opt.mult_opt->divergence_cutoff=100.0;*/

	if(realign_ends) {
fprintf(stdout,"Realigning column %d-%d \n",blocks[nblocks-1].last_col+1,mult_aln.seqs[0].len);
	opt.mult_opt->gap_opt->nendgappenalties=TRUE;
	opt.mult_opt->gap_opt->cendgappenalties=FALSE;
	opt.mult_opt->prot_gap_open=10.0;
	opt.mult_opt->divergence_cutoff=25;
	realign_block(&mult_aln,blocks[nblocks-1].last_col+1,mult_aln.seqs[0].len-1,opt);
	}

	opt.mult_opt->gap_opt->nendgappenalties=TRUE;
	opt.mult_opt->gap_opt->cendgappenalties=TRUE;
	for(i=nblocks-2;i>=0;i--) {
fprintf(stdout,"Realigning column %d-%d \n",blocks[i].last_col+1,blocks[i+1].first_col-1);
		realign_block(&mult_aln,blocks[i].last_col+1,blocks[i+1].first_col-1,opt);
	}
	if(realign_ends && blocks[0].first_col>0) {
fprintf(stdout,"Realigning column 0-%d \n",blocks[0].first_col-1);
	opt.mult_opt->prot_gap_open=10.0;
	opt.mult_opt->divergence_cutoff=25;
	opt.mult_opt->gap_opt->nendgappenalties=FALSE;
	opt.mult_opt->gap_opt->cendgappenalties=TRUE;
	realign_block(&mult_aln,0,blocks[0].first_col-1,opt);
	}


/* write out the sequences */
	strcpy((*opt.alnout_opt).gcg_outname, outfile);
	for (i=0;i<mult_aln.nseqs;i++) mult_aln.seqs[i].output_index = i;

	for (i=0;i<mult_aln.nseqs;i++) mult_aln.seqs[i].simgroup = (-1);

	if(!open_alignment_output(infile,opt.alnout_opt)) exit(1);
        create_alignment_output(mult_aln,*opt.alnout_opt);
	fprintf(stdout,"\n");

	remove(mult_aln.treename);
}


static sint count_blocks(char *filename)
{
	FILE *fin;
	char line[MAXLINE+1];
	int nblocks;

        if((fin=fopen(filename,"r"))==NULL) {
                error("Could not open block file %s",filename);
                return -1;    
        }

	nblocks=0;
	while(fgets(line,MAXLINE+1,fin)) {
		if(keyword(line,"GLOBALCOREBLOCK")) {
			nblocks++;
		}
	}
	
	fclose(fin);
	return nblocks;
}

static void read_blocks(char *filename,ALN mult_aln,CBLOCK *blocks)
{
	FILE *fin;
	char line[MAXLINE+1];
	int nblocks;
	int f,l;
	char tmp[MAXLINE+1];

        if((fin=fopen(filename,"r"))==NULL) {
                error("Could not open block file %s",filename);
                return;    
        }
	nblocks=0;
	while(fgets(line,MAXLINE+1,fin)) {
		if(keyword(line,"GLOBALCOREBLOCK")) {
			sscanf(line,"%s %d %d \n",tmp,&f,&l);
			if(f<l ) {
				blocks[nblocks].first_col=f-1;
				blocks[nblocks].last_col=l-1;
				nblocks++;
			}
		}
	}
	
	fclose(fin);
}

static void realign_block(ALNPTR mult_aln,sint first,sint last,OPT opt)
{
	sint i,j,k;
	sint save_order,length,length1,length2;
	sint maxlen;
	ALNCOUNT *alncount;
	sint *tmplen_array;
	sint tmp_alnlen;
	char **tmp_array;
	sint *newlen_array;
        char **new_array;

        length=last-first+2;
        if(length<=1) return;

/* save the alignment into a temporary area */

        maxlen=0;
        for (i=0;i<mult_aln->nseqs;i++)
                if(mult_aln->seqs[i].len>maxlen) maxlen=mult_aln->seqs[i].len;

        tmplen_array=(sint *)ckalloc((mult_aln->nseqs+2) * sizeof(sint));
        tmp_array=(char **)ckalloc((mult_aln->nseqs+2) * sizeof(char *));
        for (i=0;i<mult_aln->nseqs;i++)
        {
                tmplen_array[i]=mult_aln->seqs[i].len;
                tmp_array[i]=(char *)ckalloc((maxlen+2) * sizeof(char));
                for(j=0;j<mult_aln->seqs[i].len;j++)
                        tmp_array[i][j]=mult_aln->seqs[i].data[j];
                for(j=mult_aln->seqs[i].len;j<maxlen;j++)
                        tmp_array[i][j]=GAP2;
        }
	tmp_alnlen=maxlen;

/* copy the selected residue range to the clustal alignment arrays, removing any gaps */

        for (i=0;i<mult_aln->nseqs;i++) {
                mult_aln->seqs[i].len=length;
                /*realloc_seq(&mult_aln->seqs[i],length);*/
                ckfree(mult_aln->seqs[i].data);
                mult_aln->seqs[i].data = (char *)ckalloc((length+2) * sizeof (char));
                for(j=first,k=0;j<=last;j++)
                        if(isalpha(tmp_array[i][j])) mult_aln->seqs[i].data[k++]=tmp_array[i][j];
                mult_aln->seqs[i].data[k]=EOS;
                mult_aln->seqs[i].len=k;
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
        align_block(mult_aln,opt,alncount);
	ckfree(alncount);

        (*opt.alnout_opt).output_order=save_order;

/* save the new alignment into another temporary area */
        newlen_array=(sint *)ckalloc((mult_aln->nseqs+2) * sizeof(sint));
        new_array=(char **)ckalloc((mult_aln->nseqs+2) * sizeof(char *));
        maxlen=0;
        for (i=0;i<mult_aln->nseqs;i++)
                if(mult_aln->seqs[i].len>maxlen) maxlen=mult_aln->seqs[i].len;
        for (i=0;i<mult_aln->nseqs;i++)
        {
                new_array[i]=(char *)ckalloc((maxlen+2) * sizeof(char));
                for(j=0;j<mult_aln->seqs[i].len;j++)
                        new_array[i][j]=mult_aln->seqs[i].data[j];
                for(j=mult_aln->seqs[i].len;j<maxlen;j++)
                        new_array[i][j]='.';
                newlen_array[i]=maxlen;
        }

/* paste the realigned range back into the alignment */
/* tmp_array contains the complete sequences of the old alignment, new_array contains the realigned regions only,
   length is the old length of the region to be realigned */
        length1=length2=0;
        for (i=0;i<mult_aln->nseqs;i++)
        {
                length1=tmp_alnlen-length+newlen_array[i]+1;
                length2=newlen_array[i];
                mult_aln->seqs[i].len=length1;
                /*realloc_seq(&mult_aln->seqs[i],length1);*/
                ckfree(mult_aln->seqs[i].data);
                mult_aln->seqs[i].data = (char *)ckalloc((length1+2) * sizeof (char));
                for(j=0;j<first;j++)
                        mult_aln->seqs[i].data[j]=tmp_array[i][j];
/* copy the new alignment */
			for(j=first;j<first+length2;j++) 
				mult_aln->seqs[i].data[j]=new_array[i][j-first];
/* copy the c-ter sequences */
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

static void align_block(ALNPTR mult_aln,OPT opt,ALNCOUNT *alncount)
{
        sint i,n;
	sint meanlen,minlen,maxlen,len;

	meanlen=0;
	maxlen=0;
	minlen=alncount[0].a;
        for(i=0,n=0;i<mult_aln->nseqs;i++) {
		len=alncount[i].a;
		if(len>0) n++;
		meanlen+=len;
		if(len>0 && minlen>len) minlen=len;
		if(maxlen<len) maxlen=len;
	}
	if(n>0)
	meanlen/=(float)n;
fprintf(stdout,"meanlen %d\n",meanlen);

	if(meanlen<=7) {
		opt.mult_opt->gap_opt->no_hyd_penalties = TRUE;
		opt.mult_opt->gap_opt->no_pref_penalties = TRUE;
		opt.mult_opt->prot_gap_open=10.0;
		opt.mult_opt->divergence_cutoff=100;
	}
	else if(meanlen<=15) {
		opt.mult_opt->gap_opt->no_hyd_penalties = FALSE;
		opt.mult_opt->gap_opt->no_pref_penalties = FALSE;
		opt.mult_opt->prot_gap_open=15.0;
		opt.mult_opt->divergence_cutoff=40;
	}
	else if(meanlen<=20) {
		opt.mult_opt->gap_opt->no_hyd_penalties = FALSE;
		opt.mult_opt->gap_opt->no_pref_penalties = FALSE;
		opt.mult_opt->prot_gap_open=10.0;
		opt.mult_opt->divergence_cutoff=35;
	}
	else if(meanlen<=40) {
		opt.mult_opt->gap_opt->no_hyd_penalties = FALSE;
		opt.mult_opt->gap_opt->no_pref_penalties = FALSE;
		opt.mult_opt->prot_gap_open=10.0;
		opt.mult_opt->divergence_cutoff=30;
	}
	else if(meanlen<=60) {
		opt.mult_opt->gap_opt->no_hyd_penalties = FALSE;
		opt.mult_opt->gap_opt->no_pref_penalties = FALSE;
		opt.mult_opt->prot_gap_open=8.0;
		opt.mult_opt->divergence_cutoff=30;
	}
	else {
		opt.mult_opt->gap_opt->no_hyd_penalties = FALSE;
		opt.mult_opt->gap_opt->no_pref_penalties = FALSE;
		opt.mult_opt->prot_gap_open=7.0;
		opt.mult_opt->divergence_cutoff=30;
	}
	strcpy(mult_aln->treename,infile);
	strcat(mult_aln->treename,".ph");
 	align_from_tree(mult_aln,opt,FALSE,TRUE,FALSE,alncount);

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

