#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "rascal.h"


typedef struct {
	sint g1;
	sint f1;
	sint l1;
	sint g2;
	sint f2;
	sint l2;
	sint score;
	sint newblock;
} BLOCKERR,*BLOCKERR_PTR;

typedef struct {
	sint n;
	sint len;
	sint *g;
	sint *fres;
	sint *fcol;
}NEWBLOCK,*NEWBLOCK_PTR;

static sint count_blockerrs(char *filename,sint low_cutoff,sint high_cutoff);
static void read_blockerrs(char *filename,ALNPTR mult_aln,BLOCKERR *blockerrs,sint low_cutoff,sint high_cutoff);
static sint check_block_conflict(BLOCKERR *blockerrs, sint n,BLOCKERR block);
static Boolean block_overlap(BLOCKERR block1,BLOCKERR block2);
static void swap_errs(BLOCKERR *blockerrs,int s1, int s2);
static void sort_errs(BLOCKERR *blockerrs,int f,int l);
static void realign_blocks(ALNPTR mult_aln,sint *secgroup,NEWBLOCK newblock);
static void align_blocks(ALNPTR mult_aln,sint *secgroup,NEWBLOCK newblock,sint first,sint last);
static void c2p(char *seq,sint cstart,sint *pstart);
static void p2c(char *seq,sint seqlen,sint pstart,sint *cstart);
static void remove_gaps(ALNPTR mult_aln);

sint ngroups;
sint *groupseq;

int main(int argc,char **argv)
{
	int i,j,k,n;
	char infile[FILENAMELEN+1];
	char errfile[FILENAMELEN+1];
	char outfile[FILENAMELEN+1];
	ALN mult_aln;
	OPT opt;
	Boolean found;
	sint found1,found2;
	BLOCKERR *blockerrs;
	BLOCKERR *blockwarns;
	NEWBLOCK *newblocks;
	sint nblockerrs,nblockwarns;
	sint iseq,nseqs,status;
	sint *secgroup;
	sint *orggroup;
	sint nnewblocks;
	sint block_errcutoff=1600;
	sint block_warncutoff=1000;


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
	ngroups=get_groups(errfile,&mult_aln,secgroup,orggroup);
	ckfree(orggroup);
	if(ngroups<=0) {
/* write out the sequences */
		strcpy((*opt.alnout_opt).gcg_outname, outfile);
		if(!open_alignment_output(infile,opt.alnout_opt)) exit(1);
       		for (i=0;i<mult_aln.nseqs;i++) 
			mult_aln.seqs[i].simgroup=secgroup[i];
        	create_alignment_output(mult_aln,*opt.alnout_opt);
			exit(1);
	}

	groupseq=(sint *)ckalloc((ngroups+1)*sizeof(sint));
	for(j=0;j<ngroups;j++) {
		for(i=0;i<mult_aln.nseqs;i++) {
			if(secgroup[i]==j+1) {
				groupseq[j]=i;
				break;
			}
		}
	}

	strcpy((*opt.alnout_opt).gcg_outname, outfile);
	for (i=0;i<mult_aln.nseqs;i++) mult_aln.seqs[i].output_index = i;

/* read in the block errors and warnings */
	nblockerrs=count_blockerrs(errfile,block_errcutoff,100000);
	if(nblockerrs<=0) {
		fprintf(stdout,"No block errors in file\n");
        	for (i=0;i<mult_aln.nseqs;i++) 
			mult_aln.seqs[i].simgroup=secgroup[i];
		if(!open_alignment_output(infile,opt.alnout_opt)) exit(1);
        	create_alignment_output(mult_aln,*opt.alnout_opt);
		exit(1);
	}

	nblockwarns=count_blockerrs(errfile,block_warncutoff,block_errcutoff);

	blockerrs=(BLOCKERR *)ckalloc((nblockerrs+nblockwarns+1) * sizeof(BLOCKERR));
	read_blockerrs(errfile,&mult_aln,blockerrs,block_errcutoff,100000);

/* sort the errors into descending score order */
	sort_errs(blockerrs,0,nblockerrs-1);

/* assign the block errors to overlapping groups */
	blockerrs[0].newblock=nnewblocks=1;
	i=0;
	for(i=1;i<nblockerrs;i++) {
		status=check_block_conflict(blockerrs,i,blockerrs[i]);
		if(status<=0) {
			continue;
		}
		found=FALSE;
		for(j=0;j<nblockerrs;j++) {
			if(blockerrs[j].newblock>0 && block_overlap(blockerrs[i],blockerrs[j])) {
				blockerrs[i].newblock=blockerrs[j].newblock;
				found=TRUE;
				break;
			}
		}
		if(found==FALSE) {
			blockerrs[i].newblock=++nnewblocks;
		}
	}

	if(nblockwarns>0) {
		blockwarns=(BLOCKERR *)ckalloc((nblockwarns+1) * sizeof(BLOCKERR));
		read_blockerrs(errfile,&mult_aln,blockwarns,block_warncutoff,block_errcutoff);
/* sort the warnings into descending score order */
		sort_errs(blockwarns,0,nblockwarns-1);
	}

/* check whether the block warning overlaps an existing block error */
	n=0;
	for(i=0;i<nblockwarns;i++) {
		status=check_block_conflict(blockerrs,nblockerrs+n,blockwarns[i]);
		if(status<=0) continue;
		for(j=0;j<nblockerrs;j++) {
			if(block_overlap(blockwarns[i],blockerrs[j])) {
				blockerrs[nblockerrs+n].g1=blockwarns[i].g1;
				blockerrs[nblockerrs+n].f1=blockwarns[i].f1;
				blockerrs[nblockerrs+n].l1=blockwarns[i].l1;
				blockerrs[nblockerrs+n].g2=blockwarns[i].g2;
				blockerrs[nblockerrs+n].f2=blockwarns[i].f2;
				blockerrs[nblockerrs+n].l2=blockwarns[i].l2;
				blockerrs[nblockerrs+n].score=blockwarns[i].score;
				blockerrs[nblockerrs+n].newblock=blockerrs[j].newblock;
				n++;
				break;
			}
		}
	}

	nblockerrs+=n;

/* sort the pairwise block errors into new block regions */
	newblocks=(NEWBLOCK *)ckalloc((nnewblocks+1) * sizeof(NEWBLOCK));
	for(i=0;i<nnewblocks;i++) {
		newblocks[i].g=(sint *)ckalloc((ngroups+1)*sizeof(sint));
		newblocks[i].fres=(sint *)ckalloc((ngroups+1)*sizeof(sint));
		newblocks[i].fcol=(sint *)ckalloc((ngroups+1)*sizeof(sint));
		n=0;
		for(j=0;j<nblockerrs;j++) {
			if(blockerrs[j].newblock==i+1) {
				if(n==0) {
					newblocks[i].g[n]=blockerrs[j].g1;
					newblocks[i].fcol[n]=blockerrs[j].f1;
					iseq=groupseq[newblocks[i].g[n]-1];
					c2p(mult_aln.seqs[iseq].data,newblocks[i].fcol[n],&newblocks[i].fres[n]);
					newblocks[i].len=blockerrs[j].l1-blockerrs[j].f1+1;
					n++;
					newblocks[i].g[n]=blockerrs[j].g2;
					newblocks[i].fcol[n]=blockerrs[j].f2;
					iseq=groupseq[newblocks[i].g[n]-1];
					c2p(mult_aln.seqs[iseq].data,newblocks[i].fcol[n],&newblocks[i].fres[n]);
					n++;
				}
				else {
					found1=found2=(-1);
					for(k=0;k<n;k++) {
						if(blockerrs[j].g1==newblocks[i].g[k]) found1=k;
						if(blockerrs[j].g2==newblocks[i].g[k]) found2=k;
					}
					if(found1>=0 && found2==(-1)) {
						newblocks[i].g[n]=blockerrs[j].g2;
						newblocks[i].fcol[n]=blockerrs[j].f2+newblocks[i].fcol[found1]-blockerrs[j].f1;
						iseq=groupseq[newblocks[i].g[n]-1];
						c2p(mult_aln.seqs[iseq].data,newblocks[i].fcol[n],&newblocks[i].fres[n]);
						n++;
					}
					else if (found2>=0 && found1==(-1)) {
						newblocks[i].g[n]=blockerrs[j].g1;
						newblocks[i].fcol[n]=blockerrs[j].f1+newblocks[i].fcol[found2]-blockerrs[j].f2;
						iseq=groupseq[newblocks[i].g[n]-1];
						c2p(mult_aln.seqs[iseq].data,newblocks[i].fcol[n],&newblocks[i].fres[n]);
						n++;
					}
				}
			}
		}
		newblocks[i].n=n;
	}

	for(i=0;i<nnewblocks;i++) {
		realign_blocks(&mult_aln,secgroup,newblocks[i]);
	}

/* write out the sequences */
	strcpy((*opt.alnout_opt).gcg_outname, outfile);
	for (i=0;i<mult_aln.nseqs;i++) {
		mult_aln.seqs[i].output_index = i;
		for(j=0;j<mult_aln.seqs[i].len;j++)
			if(mult_aln.seqs[i].data[j]=='_') mult_aln.seqs[i].data[j]='-';
	}

	remove_gaps(&mult_aln);

	if(!open_alignment_output(infile,opt.alnout_opt)) exit(1);
       	for (i=0;i<mult_aln.nseqs;i++) 
		mult_aln.seqs[i].simgroup=secgroup[i];
        create_alignment_output(mult_aln,*opt.alnout_opt);
}

/* check that this block doesn't conflict with any higher scoring blocks */
static sint check_block_conflict(BLOCKERR *blockerrs, sint n,BLOCKERR block)
{
	sint i;
	Boolean found1,found2;

	found1=found2=FALSE;
	for (i=0;i<n;i++) {
		if(blockerrs[i].newblock==0) continue;
		if((blockerrs[i].g1==block.g1 && overlap(blockerrs[i].f1,blockerrs[i].l1,block.f1,block.l1)>0) ||
	   		(blockerrs[i].g2==block.g1 && overlap(blockerrs[i].f2,blockerrs[i].l2,block.f1,block.l1)>0)) {
			found1=TRUE;
		}
	   	if((blockerrs[i].g1==block.g2 && overlap(blockerrs[i].f1,blockerrs[i].l1,block.f2,block.l2)>0) ||
	   		(blockerrs[i].g2==block.g2 && overlap(blockerrs[i].f2,blockerrs[i].l2,block.f2,block.l2)>0)) {
			found2=TRUE;
		}
		if(found1 && found2) return 0;
	}
	return 1;
}

static sint count_blockerrs(char *filename,sint low_cutoff,sint high_cutoff)
{
	FILE *fin;
	char line[MAXLINE+1];
	int nblockerrs;
	char tmp[MAXLINE+1];
	sint g1,f1,l1,g2,f2,l2,score;

        if((fin=fopen(filename,"r"))==NULL) {
                error("Could not open file %s",filename);
                return 0;    
        }

	nblockerrs=0;
	while(fgets(line,MAXLINE+1,fin)) {
		if(keyword(line,"BLOCK_ERROR")) {
			sscanf(line,"%s %s %d %d %d %s %d %d %d %d \n",tmp,tmp,&g1,&f1,&l1,tmp,&g2,&f2,&l2,&score);
			if(score>low_cutoff && score<high_cutoff) nblockerrs++;
		}
	}
	
	fclose(fin);
	return nblockerrs;
}

static void read_blockerrs(char *filename,ALNPTR mult_aln,BLOCKERR *blockerrs,sint low_cutoff,sint high_cutoff)
{
	FILE *fin;
	char line[MAXLINE+1];
	int nblockerrs;
	char tmp[MAXLINE+1];
	sint g1,f1,l1,g2,f2,l2,score;

        if((fin=fopen(filename,"r"))==NULL) {
                error("Could not open file %s",filename);
                return;    
        }

	nblockerrs=0;
	while(fgets(line,MAXLINE+1,fin)) {
		if(keyword(line,"BLOCK_ERROR")) {
			sscanf(line,"%s %s %d %d %d %s %d %d %d %d \n",tmp,tmp,&g1,&f1,&l1,tmp,&g2,&f2,&l2,&score);
			if(score>low_cutoff && score<high_cutoff && f1<l1 && l1<=mult_aln->seqs[0].len && f2<l2 && l2<=mult_aln->seqs[0].len) {
				blockerrs[nblockerrs].g1=g1;
				blockerrs[nblockerrs].f1=f1;
				blockerrs[nblockerrs].l1=l1;
				blockerrs[nblockerrs].g2=g2;
				blockerrs[nblockerrs].f2=f2;
				blockerrs[nblockerrs].l2=l2;
				blockerrs[nblockerrs].score=score;
				nblockerrs++;
			}
		}
	}
	
	fclose(fin);
}

static void align_blocks(ALNPTR mult_aln,sint *secgroup,NEWBLOCK newblock,sint first,sint last)
{
	sint i,j,k,g;
	sint maxn,maxc;
	sint *block_fcol,*block_lcol;
	sint new_length,maxlen;
	sint *tmplen_array;
	char **tmp_array;

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

/* find the start and end positions of the block in each group of sequences */
	block_fcol=(sint *)ckalloc((ngroups+1)*sizeof(sint));
	block_lcol=(sint *)ckalloc((ngroups+1)*sizeof(sint));
	for(i=0;i<ngroups;i++) block_fcol[i]=(-1);
	for(i=0;i<newblock.n;i++) {
		block_fcol[newblock.g[i]-1]=newblock.fcol[i]-first;
		block_lcol[newblock.g[i]-1]=newblock.fcol[i]+newblock.len-first-1;
	}

	maxn=maxc=0;
	for(i=0;i<ngroups;i++) {
		if(block_fcol[i]>maxn) maxn=block_fcol[i];
		if(last-first-block_lcol[i]>maxc) maxc=last-first-block_lcol[i];
	}
	new_length=maxn+newblock.len+maxc;

	for(i=0;i<mult_aln->nseqs;i++) {
                mult_aln->seqs[i].len=new_length;
                ckfree(mult_aln->seqs[i].data);
                mult_aln->seqs[i].data = (char *)ckalloc((new_length+2) * sizeof (char));
		g=secgroup[i]-1;
		for(j=0;j<block_fcol[g];j++)
                        mult_aln->seqs[i].data[j]=tmp_array[i][j];
		for(;j<maxn;j++)
			mult_aln->seqs[i].data[j]='_';
		for(j=maxn;j<maxn+newblock.len;j++)
			mult_aln->seqs[i].data[j]=tmp_array[i][block_fcol[g]+j-maxn];
		for(k=block_lcol[g]+1,j=maxn+newblock.len;k<tmplen_array[i];k++) {
			mult_aln->seqs[i].data[j]=tmp_array[i][k];
			j++;
		}
		for(;j<new_length;j++)
			mult_aln->seqs[i].data[j]='_';
                mult_aln->seqs[i].data[j]=EOS;
			
	}

	ckfree(block_fcol);
	ckfree(block_lcol);
}

static void realign_blocks(ALNPTR mult_aln,sint *secgroup,NEWBLOCK newblock)
{
	sint i,j,k,n,n1;
	sint first,last;
	sint maxlen,iseq;
	sint length,length1,length2;
	sint *tmplen_array;
	char **tmp_array;
	sint *newlen_array;
        char **new_array;
	sint *fcol;
	sint *nfcol;
	sint consensus_fcol;
	Boolean found;

	for(i=0;i<newblock.n;i++) {
		iseq=groupseq[newblock.g[i]-1];
                p2c(mult_aln->seqs[iseq].data,mult_aln->seqs[iseq].len,newblock.fres[i],&newblock.fcol[i]);               
	}

	fcol=(sint *)ckalloc((ngroups+1)*sizeof(sint));
	nfcol=(sint *)ckalloc((ngroups+1)*sizeof(sint));
	for(j=0;j<ngroups;j++) {
		fcol[j]=0;
		nfcol[j]=0;
	}
	n1=0;
	for(i=0;i<newblock.n;i++) {
		fcol[i]=newblock.fcol[i];
		nfcol[i]++;
		found=FALSE;
		for(k=0;k<n1;k++)
			if(newblock.fcol[i]==fcol[k]) {
				found=TRUE;
				nfcol[k]++;
				break;
			}
		if(found==FALSE) {
			fcol[n1]=newblock.fcol[i];
			nfcol[n1]=1;
			n1++;
		}
	}
/* decide to correct the block, if (i) all blocks have signficant error scores or (ii) there is a consensus start column */
	if(newblock.n!=ngroups) {
		n=0;
		consensus_fcol=(-1);
		for(j=0;j<n1;j++) {
			if(nfcol[j]>n) {
				n=nfcol[j];
				consensus_fcol=fcol[j];
			}
		}
		if(n<=1) {
			return;
		}

/* use the consensus start column for blocks that do not have significant errors */
		n=newblock.n;
		for(j=0;j<ngroups;j++) {
			found=FALSE;
			for(k=0;k<n;k++)
				if(newblock.g[k]==j+1) {
					found=TRUE;
					break;
				}
			if (found==FALSE) {
				newblock.g[n]=j+1;
				newblock.fcol[n]=consensus_fcol;
				iseq=groupseq[j];
				c2p(mult_aln->seqs[iseq].data,newblock.fcol[n],&newblock.fres[n]);
				n++;
			}
		}
		newblock.n=n;
	}
	ckfree(fcol);
	ckfree(nfcol);

/* find the limits of the region to be realigned */
	first=mult_aln->seqs[0].len;
	last=0;
	for(i=0;i<newblock.n;i++) {
		if(newblock.fcol[i]<first) first=newblock.fcol[i];
		if(newblock.fcol[i]>last) last=newblock.fcol[i];
fprintf(stdout,"Realigning group %d: columns %d-%d\n",i+1,newblock.fcol[i],newblock.fcol[i]+newblock.len);
	}
	last+=newblock.len-1;
        maxlen=0;
        for (i=0;i<mult_aln->nseqs;i++)
                if(mult_aln->seqs[i].len>maxlen) maxlen=mult_aln->seqs[i].len;
	if(last>maxlen-1) last=maxlen-1;
	length=last-first+1;
fprintf(stdout,"Realigning %d-%d (%d)\n",first,last,length);
	if(length<=1 || length>1000) return;

/* save the alignment into a temporary area */


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
/* copy the selected residue range to the clustal alignment arrays */

        for (i=0;i<mult_aln->nseqs;i++) {
                mult_aln->seqs[i].len=length;
                ckfree(mult_aln->seqs[i].data);
                mult_aln->seqs[i].data = (char *)ckalloc((length+2) * sizeof (char));
                for(j=first;j<=last;j++)
                        mult_aln->seqs[i].data[j-first]=tmp_array[i][j];
                mult_aln->seqs[i].data[j-first]=EOS;
        }

/* align the residue range */
        align_blocks(mult_aln,secgroup,newblock,first,last);

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

static Boolean block_overlap(BLOCKERR block1,BLOCKERR block2)
{
        Boolean ret=FALSE;

	if((block1.g1==block2.g1 && overlap(block1.f1,block1.l1,block2.f1,block2.l1)>0) ||
	   (block1.g1==block2.g2 && overlap(block1.f1,block1.l1,block2.f2,block2.l2)>0) ||
	   (block1.g2==block2.g1 && overlap(block1.f2,block1.l2,block2.f1,block2.l1)>0) ||
	   (block1.g2==block2.g2 && overlap(block1.f2,block1.l2,block2.f2,block2.l2)>0)) {
		ret=TRUE;
	}
	return ret;
}


static void sort_errs(BLOCKERR *blockerrs,int f,int l)
{
        int i,last;

        if(f>=l) return;

        swap_errs(blockerrs,f,(f+l)/2);
        last=f;
        for(i=f+1;i<=l;i++)
        {
                if(blockerrs[i].score>blockerrs[f].score)
                        swap_errs(blockerrs,++last,i);
        }
        swap_errs(blockerrs,f,last);
        sort_errs(blockerrs,f,last-1);
        sort_errs(blockerrs,last+1,l);

}

static void swap_errs(BLOCKERR *blockerrs,int s1, int s2)
{
	sint g1,f1,l1,g2,f2,l2,score,newblock;

	g1=blockerrs[s1].g1;
	f1=blockerrs[s1].f1;
	l1=blockerrs[s1].l1;
	g2=blockerrs[s1].g2;
	f2=blockerrs[s1].f2;
	l2=blockerrs[s1].l2;
	score=blockerrs[s1].score;
	newblock=blockerrs[s1].newblock;
        blockerrs[s1].g1=blockerrs[s2].g1;
        blockerrs[s1].f1=blockerrs[s2].f1;
        blockerrs[s1].l1=blockerrs[s2].l1;
        blockerrs[s1].g2=blockerrs[s2].g2;
        blockerrs[s1].f2=blockerrs[s2].f2;
        blockerrs[s1].l2=blockerrs[s2].l2;
        blockerrs[s1].score=blockerrs[s2].score;
        blockerrs[s1].newblock=blockerrs[s2].newblock;
        blockerrs[s2].g1=g1;
        blockerrs[s2].f1=f1;
        blockerrs[s2].l1=l1;
        blockerrs[s2].g2=g2;
        blockerrs[s2].f2=f2;
        blockerrs[s2].l2=l2;
        blockerrs[s2].score=score;
        blockerrs[s2].newblock=newblock;
}

static void p2c(char *seq,sint seqlen,sint pstart,sint *cstart)
{
        int i,ix;

        ix=(-1);
        if(pstart<0)
        {
                (*cstart)=0;
                return;
        }
        for(i=0;i<seqlen;i++)
        {
                if(seq[i]!='_') ix++;
                if(ix==pstart) break;
        }
        (*cstart)=i;

}

static void c2p(char *seq,sint cstart,sint *pstart)
{
        int i,ix;

        ix=0;
        if(cstart<0)
        {
                (*pstart)=-1;
                return;
        }
        for(i=0;i<strlen(seq);i++)
        {
                if(i==cstart) break;
                if(seq[i]!='_') ix++;
        }
        (*pstart)=ix;

        if(*pstart<=0) (*pstart)=0;
}

static void remove_gaps(ALNPTR mult_aln)
{
        int i,j,k,ngaps;

        for (i=0;i<mult_aln->seqs[0].len;)
        {
                ngaps=0;
                for (j=0;j<mult_aln->nseqs;j++)
                        if(!isalpha(mult_aln->seqs[j].data[i])) ngaps++;
                if (ngaps==mult_aln->nseqs)
                {
                        for (j=0;j<mult_aln->nseqs;j++)
                        {
                                for(k=i+1;k<=mult_aln->seqs[j].len;k++)
                                        mult_aln->seqs[j].data[k-1]=mult_aln->seqs[j].data[k];
                                mult_aln->seqs[j].len--;
                        }
                        if(mult_aln->seqs[0].len<=0) break;
                }
                else i++;
        }
}

