/********* Sequence input routines for CLUSTAL W *******************/
/* DES was here.  FEB. 1994 */
/* Now reads PILEUP/MSF and CLUSTAL alignment files */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "rascal.h"	

#define MIN(a,b) ((a)<(b)?(a):(b))



/*
*	Prototypes
*/
static sint get_seq_with_index(SEQ *seqs,sint nseqs,sint ii);
static void xmlprinttext(FILE *fd,char *el,char *content);

/*
 *	Global variables
 */
static sint debug;


Boolean open_alignment_output(char *filename,ALNOUT_OPTPTR alnout_opt)
{

	if(!alnout_opt->output_clustal && !alnout_opt->output_nbrf && !alnout_opt->output_gcg &&
		 !alnout_opt->output_phylip && !alnout_opt->output_gde && !alnout_opt->output_nexus && !alnout_opt->output_gscope && !alnout_opt->output_relacs && !alnout_opt->output_rsf && !alnout_opt->output_tfa) {
                error("You must select an alignment output format");
                return FALSE;
        }

	if(alnout_opt->output_clustal) 
		if((alnout_opt->clustal_outfile=open_output_file(filename, "CLUSTAL", "aln", alnout_opt->clustal_outname))==NULL) return FALSE;

	if(alnout_opt->output_nbrf) 
		if((alnout_opt->nbrf_outfile=open_output_file(filename, "NBRF", "pir", alnout_opt->nbrf_outname))==NULL) return FALSE;

	if(alnout_opt->output_gcg) 
		if((alnout_opt->gcg_outfile=open_output_file(filename, "GCG", "msf", alnout_opt->gcg_outname))==NULL) return FALSE;

	if(alnout_opt->output_phylip) 
		if((alnout_opt->phylip_outfile=open_output_file(filename, "PHYLIP", "phy", alnout_opt->phylip_outname))==NULL) return FALSE;

	if(alnout_opt->output_gde) 
		if((alnout_opt->gde_outfile=open_output_file(filename, "GDE", "gde", alnout_opt->gde_outname))==NULL) return FALSE;

	if(alnout_opt->output_tfa) 
		if((alnout_opt->tfa_outfile=open_output_file(filename, "FASTA", "tfa", alnout_opt->tfa_outname))==NULL) return FALSE;

	if(alnout_opt->output_nexus) 
	if(alnout_opt->output_nexus) 
		if((alnout_opt->nexus_outfile=open_output_file(filename, "NEXUS", "nxs", alnout_opt->nexus_outname))==NULL) return FALSE;

	if(alnout_opt->output_gscope) 
		if((alnout_opt->gscope_outfile=open_output_file(filename, "GSCOPE", "msf", alnout_opt->gscope_outname))==NULL) return FALSE;

	if(alnout_opt->output_relacs) 
		if((alnout_opt->relacs_outfile=open_output_file(filename, "RELACS", "xml", alnout_opt->relacs_outname))==NULL) return FALSE;

	if(alnout_opt->output_rsf) {
		if((alnout_opt->rsf_outfile=open_output_file(filename, "RSF", "rsf", alnout_opt->rsf_outname))==NULL) return FALSE;
	}

	return TRUE;
}

/* gets a filename from the user and opens the file 
returns the file handle (the filename is written into char *file_name) */

FILE *  open_output_file(char *in_name, char *prompt, char *file_extension, char *out_name)
{
	char temp[FILENAMELEN+1];
	char path[FILENAMELEN+1];
	char local_prompt[MAXLINE];
	FILE * file_handle;

/* if the output filename is already specified, just open the file and return */
	if (out_name[0]!=EOS) {
		file_handle = open_explicit_file(out_name);
		return file_handle;
	}
}

/*
	writes the selected sequences to an output file.
	Only those sequences with output_index>=0 will be included in the file.
	Sequences will be written in order of output_index.
	eg. the sequence with output_index=3 will be the third sequence written to the file.

	mult_aln	the sequences
	alnout_opt	output options, including file format etc.
*/

void create_alignment_output(ALN mult_aln,ALNOUT_OPT alnout_opt)
{
	sint i,length,nseqs;

	length=nseqs=0;
        for (i=0;i<mult_aln.nseqs;i++)
                if (mult_aln.seqs[i].output_index>=0) {
			nseqs++;
			if(length < mult_aln.seqs[i].len) length = mult_aln.seqs[i].len;
		}

	if(alnout_opt.output_clustal) {
		clustal_out(alnout_opt.clustal_outfile, 0, length, alnout_opt.seq_numbers,mult_aln);
		fclose(alnout_opt.clustal_outfile);
		info("CLUSTAL-Alignment file created  [%s]",alnout_opt.clustal_outname);
	}
	if(alnout_opt.output_nbrf)  {
		nbrf_out(alnout_opt.nbrf_outfile, 0, length, mult_aln);
		fclose(alnout_opt.nbrf_outfile);
		info("NBRF/PIR-Alignment file created [%s]",alnout_opt.nbrf_outname);
	}
	if(alnout_opt.output_gcg)  {
		gcg_out(alnout_opt.gcg_outfile, 0, length, mult_aln);
		fclose(alnout_opt.gcg_outfile);
		info("GCG-Alignment file created      [%s]",alnout_opt.gcg_outname);
	}
	if(alnout_opt.output_phylip)  {
		phylip_out(alnout_opt.phylip_outfile, 0, length, mult_aln);
		fclose(alnout_opt.phylip_outfile);
		info("PHYLIP-Alignment file created   [%s]",alnout_opt.phylip_outname);
	}
	if(alnout_opt.output_gde)  {
		gde_out(alnout_opt.gde_outfile, 0, length, alnout_opt.lowercase,mult_aln);
		fclose(alnout_opt.gde_outfile);
		info("GDE-Alignment file created      [%s]",alnout_opt.gde_outname);
	}
	if(alnout_opt.output_tfa)  {
		tfa_out(alnout_opt.tfa_outfile, 0, length, mult_aln);
		fclose(alnout_opt.tfa_outfile);
		info("FASTA-Alignment file created      [%s]",alnout_opt.tfa_outname);
	}
	if(alnout_opt.output_nexus)  {
		nexus_out(alnout_opt.nexus_outfile, 0, length, mult_aln);
		fclose(alnout_opt.nexus_outfile);
		info("NEXUS-Alignment file created    [%s]",alnout_opt.nexus_outname);
	}
	if(alnout_opt.output_gscope)  {
		gscope_out(alnout_opt.gscope_outfile, 0, length, mult_aln);
		fclose(alnout_opt.gscope_outfile);
		info("Gscope-Alignment file created      [%s]",alnout_opt.gscope_outname);
	}
	if(alnout_opt.output_relacs)  {
		relacs_out(alnout_opt.relacs_outfile, 0, length, mult_aln, alnout_opt.relacs_outname);
		fclose(alnout_opt.relacs_outfile);
		info("RELACS-Alignment file created      [%s]",alnout_opt.relacs_outname);
	}
	if(alnout_opt.output_rsf)  {
		rsf_out(alnout_opt.rsf_outfile, 0, length, mult_aln, alnout_opt);
		fclose(alnout_opt.rsf_outfile);
		info("RSF-Alignment file created      [%s]",alnout_opt.rsf_outname);
	}
}

void clustal_out(FILE *clusout, sint fres, sint len,Boolean seq_numbers,ALN mult_aln)
{
	static char *seq1;
	static sint *seq_no;
	static sint *print_seq_no;
	char 	    temp[MAXLINE];
	char c;
	char val;
	sint ii,lv1,catident1[NUMRES],catident2[NUMRES],ident,chunks;
	sint i,j,k,l;
	sint pos,ptr,nseqs;
	sint line_length;
	sint max_aln_length,max_names;

/*These are all the positively scoring groups that occur in the Gonnet Pam250
matrix. There are strong and weak groups, defined as strong score >0.5 and
weak score =<0.5. Strong matching columns to be assigned ':' and weak matches
assigned '.' in the clustal output format.
*/

char *res_cat1[] = {
                "sta",
                "neqk",
                "nhqk",
                "ndeq",
                "qhrk",
                "milv",
                "milf",
                "hy",
                "fyw",
                NULL };

char *res_cat2[] = {
                "csa",
                "atv",
                "sag",
                "stnk",
                "stpa",
                "sgnd",
                "sndeqk",
                "ndeqhk",
                "neqhrk",
                "fvlim",
                "hfy",
                NULL };


	nseqs=0;
        for (i=0;i<mult_aln.nseqs;i++)
                if (mult_aln.seqs[i].output_index>=0) 
			nseqs++;

	seq_no = (sint *)ckalloc((nseqs+1) * sizeof(sint));
	print_seq_no = (sint *)ckalloc((nseqs+1) * sizeof(sint));
	max_aln_length=max_names=0;
	for (i=0;i<mult_aln.nseqs;i++)
        {
		if(mult_aln.seqs[i].output_index>=0) {
			 if(mult_aln.seqs[i].len>max_aln_length) max_aln_length=mult_aln.seqs[i].len;
			if(strlen(mult_aln.seqs[i].name)>max_names) max_names=strlen(mult_aln.seqs[i].name);
			print_seq_no[i] = seq_no[i] = 0;
                	for(j=0;j<fres;j++) {
				val = mult_aln.seqs[i].data[j];
				if(isalpha(val)) seq_no[i]++;
			}
		}
        }


	seq1 = (char *)ckalloc((max_aln_length+1) * sizeof(char));

	fprintf(clusout,"CLUSTAL multiple sequence alignment\n\n");

/* decide the line length for this alignment - maximum is LINELENGTH */
	line_length=PAGEWIDTH-max_names;
	line_length=line_length-line_length % 10; /* round to a multiple of 10*/
	if (line_length > LINELENGTH) line_length=LINELENGTH;

	chunks = len/line_length;
	if(len % line_length != 0)
		++chunks;
		
	for(lv1=0;lv1<chunks;lv1++) {
		pos = lv1*line_length;
		ptr = (len<pos+line_length) ? len : pos+line_length;

		fprintf(clusout,"\n");
		
		for(ii=0;ii<nseqs;ii++) {
			i=get_seq_with_index(mult_aln.seqs,mult_aln.nseqs,ii);
			print_seq_no[i] = 0;
			for(j=pos;j<ptr;j++) {
				if (j+fres<mult_aln.seqs[i].len)
					val = mult_aln.seqs[i].data[j+fres];
				else val = -3;
				if(val == EOS) break;
				else if(!isalpha(val))
                                        seq1[j]='-';
				else {
					seq1[j]=toupper(val);
					seq_no[i]++;
					print_seq_no[i]=1;
				}
			}
			for(;j<ptr;j++) seq1[j]='-';
			strncpy(temp,&seq1[pos],ptr-pos);
			temp[ptr-pos]=EOS;
			fprintf(clusout,"%-*s %s",max_names+5,mult_aln.seqs[i].name,temp);
			if (seq_numbers && print_seq_no[i])
					fprintf(clusout," %d",seq_no[i]);
			fprintf(clusout,"\n");
		}
	
		for(i=pos;i<ptr;i++) {
			seq1[i]=' ';
			ident=0;
			for(j=0;res_cat1[j]!=NULL;j++) catident1[j] = 0;
			for(j=0;res_cat2[j]!=NULL;j++) catident2[j] = 0;
			for(j=0;j<mult_aln.nseqs;j++) {
				if(isalpha(mult_aln.seqs[0].data[i])) {
					if(mult_aln.seqs[0].data[i] == mult_aln.seqs[j].data[i])
					++ident;
					for(k=0;res_cat1[k]!=NULL;k++) {
					        for(l=0;(c=res_cat1[k][l]);l++) {
							if (mult_aln.seqs[j].data[i]==c)
							{
								catident1[k]++;
								break;
							}
						}
					}
					for(k=0;res_cat2[k]!=NULL;k++) {
					        for(l=0;(c=res_cat2[k][l]);l++) {
							if (mult_aln.seqs[j].data[i]==c)
							{
								catident2[k]++;
								break;
							}
						}
					}
				}
			}
			if(ident==mult_aln.nseqs)
				seq1[i]='*';
			else if (!mult_aln.dnaflag) {
				for(k=0;res_cat1[k]!=NULL;k++) {
					if (catident1[k]==mult_aln.nseqs) {
						seq1[i]=':';
						break;
					}
				}
				if(seq1[i]==' ')
					for(k=0;res_cat2[k]!=NULL;k++) {
						if (catident2[k]==mult_aln.nseqs) {
							seq1[i]='.';
							break;
						}
					}
			}
		}
		strncpy(temp,&seq1[pos],ptr-pos);
		temp[ptr-pos]=EOS;
		for(k=0;k<max_names+6;k++) fprintf(clusout," ");
		fprintf(clusout,"%s\n",temp);
	}
		

	seq1=ckfree((void *)seq1);
	
}


static sint get_seq_with_index(SEQ *seqs,sint nseqs,sint ii)
{
	sint i;

	for(i=0;i<nseqs;i++) {
		if(seqs[i].output_index==ii) break;
	}
	if(i==nseqs) i=(-1);

	return i;
}



void gcg_out(FILE *gcgout, sint fres, sint len, ALN mult_aln)
{
	char *seq, residue;
	char val;
	sint *all_checks;
	sint i,ii,chunks,block,nseqs;
	sint j,pos1,pos2;	
	sint max_aln_length,max_names;
	long grand_checksum;

	nseqs=0;
        for (i=0;i<mult_aln.nseqs;i++)
                if (mult_aln.seqs[i].output_index>=0) 
			nseqs++;

	max_aln_length=max_names=0;
	for(i=0; i<mult_aln.nseqs; i++) {
		if(mult_aln.seqs[i].output_index>=0 && mult_aln.seqs[i].len>max_aln_length) max_aln_length=mult_aln.seqs[i].len;
		if(mult_aln.seqs[i].output_index>=0 && strlen(mult_aln.seqs[i].name)>max_names) max_names=strlen(mult_aln.seqs[i].name);
	}

	seq = (char *)ckalloc((max_aln_length+1) * sizeof(char));
	all_checks = (sint *)ckalloc((mult_aln.nseqs+1) * sizeof(sint));
	
	for(i=0; i<nseqs; i++) {
		if(mult_aln.seqs[i].output_index>=0) {
			for(j=fres; j<fres+len; j++) {
				val = mult_aln.seqs[i].data[j];
				if(val == EOS) break;
				else if(!isalpha(val))
					residue = '.';
                        	else {
					residue = val;
				}
				seq[j-fres] = toupper(residue);
			}
/* pad any short sequences with gaps, to make all sequences the same length */
			for(; j<fres+len; j++) 
				seq[j-fres] = '.';
			all_checks[i] = SeqGCGCheckSum(seq, (int)len);
		}
	}	

	grand_checksum = 0;
	for(i=0; i<mult_aln.nseqs; i++) if(mult_aln.seqs[i].output_index>=0) grand_checksum += all_checks[i];
	grand_checksum = grand_checksum % 10000;
        fprintf(gcgout,"PileUp\n\n");
	fprintf(gcgout,"\n\n   MSF:%5d  Type: ",(pint)len);
	if(mult_aln.dnaflag)
		fprintf(gcgout,"N");
	else
		fprintf(gcgout,"P");
	fprintf(gcgout,"    Check:%6ld   .. \n\n", (long)grand_checksum);

	for(ii=0; ii<nseqs; ii++)  {
		i = get_seq_with_index(mult_aln.seqs,mult_aln.nseqs,ii);
		if(i<0) continue;
		fprintf(gcgout,
			" Name: %s oo  Len:%5d  Check:%6ld  Weight:  %.1f\n",
			mult_aln.seqs[i].name,(pint)len,(long)all_checks[i],(float)mult_aln.seqs[i].weight*100.0/(float)INT_SCALE_FACTOR);
        }
	fprintf(gcgout,"\n//\n");  

	chunks = len/GCG_LINELENGTH;
	if(len % GCG_LINELENGTH != 0) ++chunks;

	for(block=0; block<chunks; block++) {
		fprintf(gcgout,"\n\n");
		pos1 = (block * GCG_LINELENGTH);
		pos2 = (len<pos1+GCG_LINELENGTH)? len : pos1+GCG_LINELENGTH;
		for(ii=0; ii<nseqs; ii++) {
			i = get_seq_with_index(mult_aln.seqs,mult_aln.nseqs,ii);
			if(i<0) continue;
			fprintf(gcgout,"\n%-*s ",max_names+5,mult_aln.seqs[i].name);
			for(j=pos1; j<pos2; j++) {
/*
    JULIE -
    check for sint sequences - pad out with '.' characters to end of alignment
*/
				if (j+fres<mult_aln.seqs[i].len)
					val = mult_aln.seqs[i].data[j+fres];
				else val = -3;
				if(val == EOS)
					residue = '.';
				else if(!isalpha(val))
					residue = '.';
				else {
					residue = toupper(val);
				}
				fprintf(gcgout,"%c",residue);
				if((j+1) % 10 == 0) fprintf(gcgout," ");
			}
		}
	}

	seq=ckfree((void *)seq);
	all_checks=ckfree((void *)all_checks);
	
	fprintf(gcgout,"\n\n");
}




void gscope_out(FILE *gcgout, sint fres, sint len, ALN mult_aln)
{
	char *seq, residue;
	char val;
	sint *all_checks;
	sint i,ii,chunks,block;
	sint j,pos1,pos2,nseqs;	
	sint max_aln_length,max_names;
	long grand_checksum;

	nseqs=0;
        for (i=0;i<mult_aln.nseqs;i++)
                if (mult_aln.seqs[i].output_index>=0) 
			nseqs++;

	max_aln_length=max_names=0;
	for(i=0; i<mult_aln.nseqs; i++) {
		if(mult_aln.seqs[i].output_index>=0 && mult_aln.seqs[i].len>max_aln_length) max_aln_length=mult_aln.seqs[i].len;
		if(mult_aln.seqs[i].output_index>=0 && strlen(mult_aln.seqs[i].name)>max_names) max_names=strlen(mult_aln.seqs[i].name);
	}

	seq = (char *)ckalloc((max_aln_length+1) * sizeof(char));
	all_checks = (sint *)ckalloc((mult_aln.nseqs+1) * sizeof(sint));
	
	for(i=0; i<mult_aln.nseqs; i++) {
		if(mult_aln.seqs[i].output_index>=0) {
			for(j=fres; j<fres+len; j++) {
				val = mult_aln.seqs[i].data[j];
				if(val == EOS) break;
				else if(!isalpha(val))
					residue = '.';
                        	else {
					residue = val;
				}
				seq[j-fres] = toupper(residue);
			}
/* pad any short sequences with gaps, to make all sequences the same length */
			for(; j<fres+len; j++) 
				seq[j-fres] = '.';
			all_checks[i] = SeqGCGCheckSum(seq, (int)len);
		}
	}	

	grand_checksum = 0;
	for(i=0; i<mult_aln.nseqs; i++) if(mult_aln.seqs[i].output_index>=0) grand_checksum += all_checks[i];
	grand_checksum = grand_checksum % 10000;
        fprintf(gcgout,"PileUp\n");
        fprintf(gcgout,"ANCHORS %d\n\n",(pint)mult_aln.nanchors);
	fprintf(gcgout,"\n\n   MSF:%5d  Type: ",(pint)len);
	if(mult_aln.dnaflag)
		fprintf(gcgout,"N");
	else
		fprintf(gcgout,"P");
	fprintf(gcgout,"    Check:%6ld   .. \n\n", (long)grand_checksum);
	for(ii=0; ii<nseqs; ii++)  {
		i = get_seq_with_index(mult_aln.seqs,mult_aln.nseqs,ii);
		if(i<0) continue;
		fprintf(gcgout,
			" Name: %s oo  Len:%5d  Check:%6ld  Weight:  %.1f\n",
			mult_aln.seqs[i].name,(pint)len,(long)all_checks[i],(float)mult_aln.seqs[i].weight*100.0/(float)INT_SCALE_FACTOR);
        }
	fprintf(gcgout,"\n//\n");  

	chunks = len/GSCOPE_LINELENGTH;
	if(len % GSCOPE_LINELENGTH != 0) ++chunks;

	for(block=0; block<chunks; block++) {
		fprintf(gcgout,"\n\n");
		pos1 = (block * GSCOPE_LINELENGTH);
		pos2 = (len<pos1+GSCOPE_LINELENGTH)? len : pos1+GSCOPE_LINELENGTH;
		for(ii=0; ii<nseqs; ii++) {
			i = get_seq_with_index(mult_aln.seqs,mult_aln.nseqs,ii);
			if(i<0) continue;
			fprintf(gcgout,"\n%-*s ",max_names+5,mult_aln.seqs[i].name);
			for(j=pos1; j<pos2; j++) {
/*
    JULIE -
    check for sint sequences - pad out with '.' characters to end of alignment
*/
				if (j+fres<mult_aln.seqs[i].len)
					val = mult_aln.seqs[i].data[j+fres];
				else val = -3;
				if(val == EOS)
					residue = '.';
				else if(!isalpha(val))
					residue = '.';
				else {
					residue = toupper(val);
				}
				fprintf(gcgout,"%c",residue);
				if((j+1) % 10 == 0) fprintf(gcgout," ");
			}
		}
	}

	seq=ckfree((void *)seq);
	all_checks=ckfree((void *)all_checks);
	
	fprintf(gcgout,"\n\n");
}



void nexus_out(FILE *nxsout, sint fres, sint len, ALN mult_aln)
{
	char residue;
	char val;
	sint i,ii,chunks,block;	
	sint j,pos1,pos2,nseqs;	
	
	nseqs=0;
        for (i=0;i<mult_aln.nseqs;i++)
                if (mult_aln.seqs[i].output_index>=0) 
			nseqs++;

	chunks = len/GCG_LINELENGTH;
	if(len % GCG_LINELENGTH != 0) ++chunks;

        fprintf(nxsout,"#NEXUS\n");

        fprintf(nxsout,"BEGIN DATA;\n");
        fprintf(nxsout,"dimensions ntax=%d nchar=%d;\n",(pint)nseqs,(pint)len);
        fprintf(nxsout,"format missing=?\n");
        fprintf(nxsout,"symbols=\"");
        for(i=0;i<strlen(mult_aln.alphabet)-1;i++)
                fprintf(nxsout,"%c",mult_aln.alphabet[i]);
        fprintf(nxsout,"\"\n");
        fprintf(nxsout,"interleave datatype=");
        fprintf(nxsout, mult_aln.dnaflag ? "DNA " : "PROTEIN ");
        fprintf(nxsout,"gap= -;\n");
        fprintf(nxsout,"\nmatrix");



	for(block=0; block<chunks; block++) {
		pos1 = (block * GCG_LINELENGTH);
		pos2 = (len<pos1+GCG_LINELENGTH)? len : pos1+GCG_LINELENGTH;
		for(ii=0; ii<nseqs; ii++)  {
			i = get_seq_with_index(mult_aln.seqs,mult_aln.nseqs,ii);
			if(i<0) continue;
			fprintf(nxsout,"\n%-10s ",mult_aln.seqs[i].name);
			for(j=pos1; j<pos2; j++) {
				if (j+fres<mult_aln.seqs[i].len)
					val = mult_aln.seqs[i].data[j+fres];
				else val = -3;
				if(val == EOS)
					break;
				else if(!isalpha(val))
					residue = '-';
				else {
					residue = toupper(val);
				}
				fprintf(nxsout,"%c",residue);
			}
		}
		fprintf(nxsout,"\n");
	}

        fprintf(nxsout,";\nend;\n");

}


void phylip_out(FILE *phyout, sint fres, sint len, ALN mult_aln)
{
	char residue;
	char val;
	sint i,ii,chunks,block,nseqs;	
	sint j,pos1,pos2;	
	sint name_len;
	Boolean warn;
	char **snames;
	
	nseqs=0;
        for (i=0;i<mult_aln.nseqs;i++)
                if (mult_aln.seqs[i].output_index>=0) 
			nseqs++;

	snames=(char **)ckalloc((mult_aln.nseqs+1)*sizeof(char *));
	name_len=0;
	for(i=0; i<mult_aln.nseqs; i++)  {
		if(mult_aln.seqs[i].output_index>=0) {
			snames[i]=(char *)ckalloc((11)*sizeof(char));
			ii=strlen(mult_aln.seqs[i].name);
			strncpy(snames[i],mult_aln.seqs[i].name,10);
			if(name_len<ii) name_len=ii;
		}
	}
	if(name_len>10) {
		warn=FALSE;
		for(i=0; i<mult_aln.nseqs; i++)  {
			if(mult_aln.seqs[i].output_index>=0) {
                		for(j=i+1;j<mult_aln.nseqs;j++) {
					if(mult_aln.seqs[j].output_index>=0) {
                        			if (strcmp(snames[i],snames[j]) == 0) 
							warn=TRUE;
                			}
				}
			}
        	}
		if(warn)
			warning("Truncating sequence names to 10 characters for PHYLIP output.\n"
			"Names in the PHYLIP format file are NOT unambiguous.");
		else
			warning("Truncating sequence names to 10 characters for PHYLIP output.");
	}


	chunks = len/GCG_LINELENGTH;
	if(len % GCG_LINELENGTH != 0) ++chunks;

	fprintf(phyout,"%6d %6d",(pint)nseqs,(pint)len);

	for(block=0; block<chunks; block++) {
		pos1 = (block * GCG_LINELENGTH);
		pos2 = (len<pos1+GCG_LINELENGTH)? len : pos1+GCG_LINELENGTH;
		for(ii=0; ii<nseqs; ii++)  {
			i = get_seq_with_index(mult_aln.seqs,mult_aln.nseqs,ii);
			if(i<0) continue;
			if(block == 0)  {
				fprintf(phyout,"\n%-10s ",snames[i]);
			}
			else
				fprintf(phyout,"\n           ");
			for(j=pos1; j<pos2; j++) {
				if (j+fres<mult_aln.seqs[i].len)
					val = mult_aln.seqs[i].data[j+fres];
				else val = -3;
				if(val == EOS)
					break;
				else if(!isalpha(val))
					residue = '-';
				else {
					residue = toupper(val);
				}
				fprintf(phyout,"%c",residue);
				if((j+1) % 10 == 0) fprintf(phyout," ");
			}
		}
		fprintf(phyout,"\n");
	}

	for(i=0;i<mult_aln.nseqs;i++)
		if(mult_aln.seqs[i].output_index>=0) 
		ckfree(snames[i]);
	ckfree(snames);
}



void fasta_out(FILE *out, sint fres, sint len, ALN mult_aln)
{
	char *seq, residue;
	char val;
	sint i,ii;
	sint j,slen,nseqs;	
	sint line_length;
	sint max_aln_length,max_names;

	nseqs=0;
	max_aln_length=max_names=0;
	for(i=0; i<mult_aln.nseqs; i++)  {
                if (mult_aln.seqs[i].output_index>=0) {
			if(mult_aln.seqs[i].len>max_aln_length) max_aln_length=mult_aln.seqs[i].len;
			if(strlen(mult_aln.seqs[i].name)>max_names) max_names=strlen(mult_aln.seqs[i].name);
			nseqs++;
		}
	}

	seq = (char *)ckalloc((max_aln_length+1) * sizeof(char));
	
/* decide the line length for this alignment - maximum is LINELENGTH */
	line_length=PAGEWIDTH-max_names;
	line_length=line_length-line_length % 10; /* round to a multiple of 10*/
	if (line_length > LINELENGTH) line_length=LINELENGTH;

	for(ii=0; ii<nseqs; ii++) {
		i = get_seq_with_index(mult_aln.seqs,mult_aln.nseqs,ii);
		if(i<0) continue;
		fprintf(out, ">%s\n%s\n", mult_aln.seqs[i].name, mult_aln.seqs[i].title);
		slen = 0;
		for(j=fres; j<fres+len; j++) {
			val = mult_aln.seqs[i].data[j];
			if(j>=mult_aln.seqs[i].len)
				residue = '-';
			else if(!isalpha(val))
				residue = '-';
			else {
				residue = val;
			}
			seq[j-fres] = toupper(residue);
			slen++;
		}
		for(j=0; j<slen; j++) {
			fprintf(out,"%c",seq[j]);
			if(((j+1) % line_length == 0) || (j == slen-1)) 
				fprintf(out,"\n");
		}
		fprintf(out,"\n");
	}	

	seq=ckfree((void *)seq);
}




void nbrf_out(FILE *nbout, sint fres, sint len, ALN mult_aln)
{
	char *seq, residue;
	char val;
	sint i,ii,nseqs;
	sint j,slen;	
	sint line_length;
	sint max_aln_length,max_names;

	nseqs=0;
	max_aln_length=max_names=0;
	for(i=0; i<mult_aln.nseqs; i++)  {
                if (mult_aln.seqs[i].output_index>=0) {
			if(mult_aln.seqs[i].len>max_aln_length) max_aln_length=mult_aln.seqs[i].len;
			if(strlen(mult_aln.seqs[i].name)>max_names) max_names=strlen(mult_aln.seqs[i].name);
			nseqs++;
		}
	}

	seq = (char *)ckalloc((max_aln_length+1) * sizeof(char));
	
/* decide the line length for this alignment - maximum is LINELENGTH */
	line_length=PAGEWIDTH-max_names;
	line_length=line_length-line_length % 10; /* round to a multiple of 10*/
	if (line_length > LINELENGTH) line_length=LINELENGTH;

	for(ii=0; ii<nseqs; ii++) {
		i = get_seq_with_index(mult_aln.seqs,mult_aln.nseqs,ii);
		if(i<0) continue;
		fprintf(nbout, mult_aln.dnaflag ? ">DL;" : ">P1;");
		fprintf(nbout, "%s\n%s\n", mult_aln.seqs[i].name, mult_aln.seqs[i].title);
		slen = 0;
		for(j=fres; j<fres+len; j++) {
			val = mult_aln.seqs[i].data[j];
			if(val == EOS)
				break;
			else if(!isalpha(val))
				residue = '-';
			else {
				residue = val;
			}
			seq[j-fres] = toupper(residue);
			slen++;
		}
		for(j=0; j<slen; j++) {
			fprintf(nbout,"%c",seq[j]);
			if(((j+1) % line_length == 0) || (j == slen-1)) 
				fprintf(nbout,"\n");
		}
		fprintf(nbout,"*\n");
	}	

	seq=ckfree((void *)seq);
}


void rsf_out(FILE *rsfout, sint fres, sint len, ALN mult_aln, ALNOUT_OPT alnout_opt)
{
	char *seq, residue;
	char val;
	sint n;
	sint i,ii,nseqs;
	sint j,k,slen,cstart,cend;	
	sint start,query_start,query_end;
	sint line_length=60;
	sint max_aln_length,max_names;
	sint color;
	char shape[20];
	char type[20];
	sint nrep,*reptype;

	nseqs=0;
	max_aln_length=max_names=0;
	for(i=0; i<mult_aln.nseqs; i++)  {
                if (mult_aln.seqs[i].output_index>=0) {
			if(mult_aln.seqs[i].len>max_aln_length) max_aln_length=mult_aln.seqs[i].len;
			if(strlen(mult_aln.seqs[i].name)>max_names) max_names=strlen(mult_aln.seqs[i].name);
                        nseqs++;
		}
	}

	seq = (char *)ckalloc((max_aln_length+1) * sizeof(char));

	fprintf(rsfout,"!!RICH_SEQUENCE 1.0\n..\n");

/* set n to the number of sequence groups */
	n=0;
	for(i=0; i<nseqs; i++) 
		if(mult_aln.seqs[i].simgroup>n) n=mult_aln.seqs[i].simgroup;
	n++;

	for(ii=0; ii<nseqs; ii++) {
		i = get_seq_with_index(mult_aln.seqs,mult_aln.nseqs,ii);
		if(i<0) continue;
		fprintf(rsfout, "{\nname ");
		if(alnout_opt.output_names==1)
			fprintf(rsfout, "%s\n", mult_aln.seqs[i].nid);
		else if(alnout_opt.output_names==2)
			fprintf(rsfout, "%s\n", mult_aln.seqs[i].access);
		else
			fprintf(rsfout, "%s\n", mult_aln.seqs[i].name);
		fprintf(rsfout, "descrip %s\n", mult_aln.seqs[i].title);
		fprintf(rsfout, "creator %s %s\n", mult_aln.seqs[i].org,mult_aln.seqs[i].lifedomain);
		fprintf(rsfout, "type    ");
		fprintf(rsfout, mult_aln.dnaflag ? "DNA" : "PROTEIN");
		fprintf(rsfout, "\n");
		if(mult_aln.seqs[i].simgroup>0) 
			fprintf(rsfout, "group %d\n",mult_aln.seqs[i].simgroup);
		else {
			fprintf(rsfout, "group %d\n",n);
		}
		if(mult_aln.seqs[i].sense!=0) 
			fprintf(rsfout, "strand %d\n",mult_aln.seqs[i].sense);

		for(j=0;j<MAXFTTYPE;j++) {
			for(k=0;k<mult_aln.ft[i].nentries[j];k++) {
				if(mult_aln.ft[i].data[j][k].start<0) continue;
				pos2col(mult_aln.seqs[i].data,mult_aln.ft[i].data[j][k].start,mult_aln.ft[i].data[j][k].end,&cstart,&cend);
				color=mult_aln.ft[i].data[j][k].color;
				if(color>=0) {
				if(color>12) color=12;
				fprintf(rsfout, "feature %d %d %d %s %s %s\n",cstart+1,cend+1,
										color,
										"square",
										"solid",
										mult_aln.ft[i].data[j][k].type);
				fprintf(rsfout,"  FT   %-10s%5d  %5d  %s\n",
								mult_aln.ft[i].data[j][k].type,
								mult_aln.ft[i].data[j][k].start+1,
								mult_aln.ft[i].data[j][k].end+1,
								mult_aln.ft[i].data[j][k].name);
				}
			}
		}

		slen = 0;
		for(j=fres; j<fres+len; j++) {
			val = mult_aln.seqs[i].data[j];
			if(val==EOS)
				break;
			else if(!isalpha(val))
				residue = '.';
			else {
				residue = val;
			}
			seq[j-fres] = (char)toupper((int)residue);
			slen++;
		}
		fprintf(rsfout,"sequence\n  ");
		for(j=0; j<slen; j++) {
			fprintf(rsfout,"%c",seq[j]);
			if( (j == slen-1)) 
				fprintf(rsfout,"\n");
			else if(( j!=0 && (j+1) % line_length == 0)) 
				fprintf(rsfout,"\n  ");
		}
		fprintf(rsfout,"}\n");
	}	

	seq=ckfree((void *)seq);
}

void tfa_out(FILE *tfaout, sint fres, sint len,ALN mult_aln)
{
	char *seq, residue;
	char val;
	sint i,ii;
	sint j,slen;	
	sint line_length;
	sint max_aln_length,max_names,nseqs;

	nseqs=0;
	max_aln_length=max_names=0;
	for(i=0; i<mult_aln.nseqs; i++)  {
		if(mult_aln.seqs[i].output_index>=0) {
			if(mult_aln.seqs[i].len>max_aln_length) max_aln_length=mult_aln.seqs[i].len;
			if(strlen(mult_aln.seqs[i].name)>max_names) max_names=strlen(mult_aln.seqs[i].name);
			nseqs++;
		}
	}

	seq = (char *)ckalloc((max_aln_length+1) * sizeof(char));

/* decide the line length for this alignment - maximum is LINELENGTH */
	line_length=PAGEWIDTH-max_names;
	line_length=line_length-line_length % 10; /* round to a multiple of 10*/
	if (line_length > LINELENGTH) line_length=LINELENGTH;

	for(ii=0; ii<nseqs; ii++) {
		i = get_seq_with_index(mult_aln.seqs,mult_aln.nseqs,ii);
		if(i<0) continue;
		fprintf(tfaout, ">");
		fprintf(tfaout, "%s\n", mult_aln.seqs[i].name);
		slen = 0;
		for(j=fres; j<fres+len; j++) {
			val = mult_aln.seqs[i].data[j];
			if(val==EOS)
				break;
			else if(!isalpha(val))
				residue = '-';
			else {
				residue = val;
			}
			seq[j-fres] = (char)toupper((int)residue);
			slen++;
		}
		for(j=0; j<slen; j++) {
			fprintf(tfaout,"%c",seq[j]);
			if(( j!=0 && (j+1) % line_length == 0) || (j == slen-1)) 
				fprintf(tfaout,"\n");
		}
	}	

	seq=ckfree((void *)seq);
}


void gde_out(FILE *gdeout, sint fres, sint len, Boolean lowercase,ALN mult_aln)
{
	char *seq, residue;
	char val;
	sint i,ii;
	sint j,slen;	
	sint line_length;
	sint max_aln_length,max_names,nseqs;

	nseqs=0;
	max_aln_length=max_names=0;
	for(i=0; i<mult_aln.nseqs; i++)  {
		if(mult_aln.seqs[i].output_index>=0) {
			if(mult_aln.seqs[i].len>max_aln_length) max_aln_length=mult_aln.seqs[i].len;
			if(strlen(mult_aln.seqs[i].name)>max_names) max_names=strlen(mult_aln.seqs[i].name);
			nseqs++;
		}
	}

	seq = (char *)ckalloc((max_aln_length+1) * sizeof(char));

/* decide the line length for this alignment - maximum is LINELENGTH */
	line_length=PAGEWIDTH-max_names;
	line_length=line_length-line_length % 10; /* round to a multiple of 10*/
	if (line_length > LINELENGTH) line_length=LINELENGTH;

	for(ii=0; ii<nseqs; ii++) {
		i = get_seq_with_index(mult_aln.seqs,mult_aln.nseqs,ii);
		if(i<0) continue;
		fprintf(gdeout, mult_aln.dnaflag ? "#" : "%%");
		fprintf(gdeout, "%s\n", mult_aln.seqs[i].name);
		slen = 0;
		for(j=fres; j<fres+len; j++) {
			val = mult_aln.seqs[i].data[j];
			if(val==EOS)
				break;
			else if(!isalpha(val))
				residue = '-';
			else {
				residue = val;
			}
			if (lowercase)
				seq[j-fres] = (char)tolower((int)residue);
			else
				seq[j-fres] = (char)toupper((int)residue);
			slen++;
		}
		for(j=0; j<slen; j++) {
			fprintf(gdeout,"%c",seq[j]);
			if(( j!=0 && (j+1) % line_length == 0) || (j == slen-1)) 
				fprintf(gdeout,"\n");
		}
	}	

	seq=ckfree((void *)seq);
}



void relacs_out(FILE *rout, sint fres, sint len, ALN mult_aln, char *aln_name)
{
	char *seq, residue;
	char *aln_root;
	char val;
	sint i,ii,k;
	sint cpos1;	
	sint j,slen;	
	sint line_length;
	sint max_aln_length,max_names;
	sint nseqs,nentries;

	fprintf(rout,"<?xml version=\"1.0\"?>\n");
 	fprintf(rout,"<!DOCTYPE macsim SYSTEM \"http://www-bio3d-igbmc.u-strasbg.fr/macsim.dtd\">\n");
 	fprintf(rout,"<macsim>\n");
 	fprintf(rout,"<alignment>\n");
	aln_root=(char *)ckalloc((strlen(aln_name)+1)*sizeof(char));
	strcpy(aln_root,aln_name);
	for(i=strlen(aln_root)-1;i>=0;i--)
		if(aln_root[i]=='.') {
			aln_root[i]='\0';
			break;
		}
 	fprintf(rout,"<aln-name>%s</aln-name>\n",aln_root);
	if(mult_aln.validalnscore)
                fprintf(rout,"<aln-score>%.3f</aln-score>\n",mult_aln.alnscore);


	nseqs=0;
	max_aln_length=max_names=0;
	for(i=0; i<mult_aln.nseqs; i++)  {
		if(mult_aln.seqs[i].output_index>=0) {
			if(mult_aln.seqs[i].len>max_aln_length) max_aln_length=mult_aln.seqs[i].len;
			if(strlen(mult_aln.seqs[i].name)>max_names) max_names=strlen(mult_aln.seqs[i].name);
			nseqs++;
		}
	}

	seq = (char *)ckalloc((max_aln_length+1) * sizeof(char));

/* decide the line length for this alignment - maximum is LINELENGTH */
	line_length=PAGEWIDTH-max_names;
	line_length=line_length-line_length % 10; /* round to a multiple of 10*/
	if (line_length > LINELENGTH) line_length=LINELENGTH;

	for(ii=0; ii<nseqs; ii++) {
		i = get_seq_with_index(mult_aln.seqs,mult_aln.nseqs,ii);
		if(i<0) continue;
		fprintf(rout, "<sequence seq-type=");
		fprintf(rout, mult_aln.dnaflag ? "\"DNA\"" : "\"Protein\"");
		fprintf(rout, ">\n");
		fprintf(rout, "<seq-name>%s</seq-name>\n", mult_aln.seqs[i].name);
		fprintf(rout, "<seq-info>\n");
		xmlprinttext(rout,"accession",mult_aln.seqs[i].access);
		xmlprinttext(rout,"nid",mult_aln.seqs[i].nid);
		xmlprinttext(rout,"definition",mult_aln.seqs[i].title);
		xmlprinttext(rout,"organism",mult_aln.seqs[i].org);
		xmlprinttext(rout,"lifedomain",mult_aln.seqs[i].lifedomain);
                if(mult_aln.seqs[i].ntaxons>0) fprintf(rout,"<lineage>\n");
                for(j=0;j<mult_aln.seqs[i].ntaxons;j++) {
                        fprintf(rout, "<taxon>%s</taxon>\n",
                                                        mult_aln.seqs[i].taxon[j]);
                }
                if(mult_aln.seqs[i].ntaxons>0) fprintf(rout,"</lineage>\n");
		xmlprinttext(rout,"taxid",mult_aln.seqs[i].taxid);
		if(mult_aln.seqs[i].ec[0]>0)
		fprintf(rout, "<ec>%d.%d.%d.%d</ec>\n", mult_aln.seqs[i].ec[0],mult_aln.seqs[i].ec[1],mult_aln.seqs[i].ec[2],mult_aln.seqs[i].ec[3]);
		if(mult_aln.seqs[i].hydrophobicity!=0.0)
		fprintf(rout, "<hydrophobicity>%.2f</hydrophobicity>\n",mult_aln.seqs[i].hydrophobicity);
		fprintf(rout, "<group>%d</group>\n",mult_aln.seqs[i].simgroup);
		if(mult_aln.seqs[i].fragment==TRUE) fprintf(rout, "<fragment status=\"Yes\"/>\n");
		if(mult_aln.seqs[i].sense!=0) fprintf(rout, "<sense>%d</sense>\n",mult_aln.seqs[i].sense);

		nentries=0;
		for(j=0;j<MAXFTTYPE;j++)
			nentries+=mult_aln.ft[i].nentries[j];
		if(nentries>0) fprintf(rout,"<ftable>\n");
		for(j=0;j<MAXFTTYPE;j++) {
			for(k=0;k<mult_aln.ft[i].nentries[j];k++) {
				if(mult_aln.ft[i].data[j][k].start<0) continue;
				if((strcmp(mult_aln.ft[i].data[j][k].name,"STRAND")==0) && mult_aln.ft[i].data[j][k].color!=8) continue;
				if((strcmp(mult_aln.ft[i].data[j][k].name,"HELIX")==0) && mult_aln.ft[i].data[j][k].color!=0) continue;
				fprintf(rout, "<fitem><ftype>%s</ftype><fstart>%d</fstart><fstop>%d</fstop><fcolor>%d</fcolor><fscore>%.2f</fscore><fnote>%s</fnote></fitem>\n",
							mult_aln.ft[i].data[j][k].type,
							mult_aln.ft[i].data[j][k].start+1,
							mult_aln.ft[i].data[j][k].end+1,
							mult_aln.ft[i].data[j][k].color,
							mult_aln.ft[i].data[j][k].score,
							mult_aln.ft[i].data[j][k].name);
			}
		}
		if(nentries>0) fprintf(rout,"</ftable>\n");

                if(mult_aln.seqs[i].accessibility!=NULL) {
                        fprintf(rout, "<surface-accessibility>\n");
                        fprintf(rout, "<suracc-name>suracc</suracc-name>\n");
                        fprintf(rout, "<suracc-type>float</suracc-type>\n");
                        fprintf(rout, "<suracc-data>");
                        for(j=0;j<mult_aln.seqs[i].naccessibility;j++)
                                fprintf(rout, " %.2f",mult_aln.seqs[i].accessibility[j]);
                        fprintf(rout, " </suracc-data>\n");
                        fprintf(rout, "</surface-accessibility>\n");
                }

                if(mult_aln.seqs[i].rcontacts!=NULL) {
                        for(j=0;j<mult_aln.seqs[i].len;j++) {
                                if(mult_aln.seqs[i].rcontacts[j].n>0) {
                                        fprintf(rout, "<residue-contact-list>\n");
                                        col2pos1(mult_aln.seqs[i].data,j,&cpos1);
                                        fprintf(rout, "<contact-residue1>%d</contact-residue1>\n",cpos1+1);
                                        for(k=0;k<mult_aln.seqs[i].rcontacts[j].n;k++) {
                                                col2pos1(mult_aln.seqs[i].data,mult_aln.seqs[i].rcontacts[j].res[k],&cpos1);
                                                fprintf(rout, "<contact-residue2>%d</contact-residue2>\n",cpos1+1);
                                        }
                                        fprintf(rout, "</residue-contact-list>\n");
                                }
                        }
                }

		if(mult_aln.go[i].ngorefs>0) fprintf(rout,"<goxreflist>\n");
		for(j=0;j<mult_aln.go[i].ngorefs;j++) {
			fprintf(rout, "<goxref><goid>%s</goid><goclass>%c</goclass>",
							mult_aln.go[i].goref[j].id,
							mult_aln.go[i].goref[j].class);
			if(mult_aln.go[i].goref[j].desc!=NULL) fprintf(rout, "<godesc>%s</godesc>", mult_aln.go[i].goref[j].desc);
			if(mult_aln.go[i].goref[j].evidence!=NULL) fprintf(rout, "<goevidence>%s</goevidence>", mult_aln.go[i].goref[j].evidence);
			fprintf(rout, "</goxref>\n");
		}
		if(mult_aln.go[i].ngorefs>0) fprintf(rout,"</goxreflist>\n");

		if(mult_aln.seqs[i].nkeywords>0) fprintf(rout,"<keywordlist>\n");
		for(j=0;j<mult_aln.seqs[i].nkeywords;j++) {
			fprintf(rout, "<keyword>%s</keyword>\n",
							mult_aln.seqs[i].keyword[j]);
		}
		if(mult_aln.seqs[i].nkeywords>0) fprintf(rout,"</keywordlist>\n");

		fprintf(rout, "<length>%d</length>\n",mult_aln.seqs[i].len);
		fprintf(rout, "<weight>%d</weight>\n",mult_aln.seqs[i].weight);
		if(mult_aln.seqs[i].simgroup>0) 
			fprintf(rout, "<group>%d</group>\n",mult_aln.seqs[i].simgroup);
		fprintf(rout, "</seq-info>\n");
		slen = 0;
		for(j=fres; j<fres+len; j++) {
			val = mult_aln.seqs[i].data[j];
			if(val==EOS)
				break;
			else if(!isalpha(val))
				residue = '-';
			else {
				residue = val;
			}
			seq[j-fres] = residue;
			slen++;
		}
		fprintf(rout, "<seq-data>");
		for(j=0; j<slen; j++) {
			fprintf(rout,"%c",seq[j]);
			/*if(( j!=0 && (j+1) % line_length == 0) || (j == slen-1)) 
				fprintf(rout,"\n");*/
		}
				fprintf(rout,"\n");
		fprintf(rout, "</seq-data></sequence>\n");
	}	

	for(i=0;i<mult_aln.ncol_scores;i++) {
		fprintf(rout, "<column-score>\n");
		fprintf(rout, "<colsco-name>%s</colsco-name>\n",mult_aln.col_score[i].name);
		fprintf(rout, "<colsco-owner>%s</colsco-owner>\n",mult_aln.col_score[i].owner);
		fprintf(rout, "<colsco-type>int</colsco-type>\n");
		fprintf(rout, "<colsco-data>\n");
		for(j=fres; j<mult_aln.col_score[i].length && j<fres+len;j++)
			fprintf(rout,"%d ",mult_aln.col_score[i].data[j]);
		fprintf(rout, "</colsco-data>\n");
		fprintf(rout, "</column-score>\n");
	}


	seq=ckfree((void *)seq);
 	fprintf(rout,"</alignment>\n");
 	fprintf(rout,"</macsim>\n");
}

static void xmlprinttext(FILE *fd,char *el,char *content)
{
	if(content!=NULL) 
		if(content[0]!='\0') fprintf(fd, "<%s>%s</%s>\n",el,content,el);
}
