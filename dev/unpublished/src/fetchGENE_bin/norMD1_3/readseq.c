/********* Sequence input routines for CLUSTAL W *******************/
/* DES was here.  FEB. 1994 */
/* Now reads PILEUP/MSF and CLUSTAL alignment files */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "score.h"	

#define MIN(a,b) ((a)<(b)?(a):(b))



/*
*	Prototypes
*/

static char * get_seq(char *sname,sint *len,char *nid,char *access,char *tit,sint *seq_type,char *org,char *strand,Boolean dnaflag,FTPTR ft);
static char * get_clustal_seq(char *sname,sint *len,sint seqno);
static char * get_msf_seq(char *sname,sint *len,sint seqno);
static void check_infile(sint *noseqs,Boolean *dnaflag);
static void p_encode(char *seq, char *naseq, sint l);
static void n_encode(char *seq, char *naseq, sint l);
static sint res_index(char *codes,char c);
static Boolean check_dnaflag(char *seq, sint len);
static sint count_clustal_seqs(void);
static sint count_pir_seqs(void);
static sint count_msf_seqs(void);
static sint count_rsf_seqs(void);
static Boolean cl_blankline(char *line);


/*
 *	Global variables
 */
FILE *fin;
static sint debug;

static sint seqFormat;
static char *formatNames[] = {"unknown","EMBL/Swiss-Prot","PIR",
			      "Pearson","GDE","Clustal","Pileup/MSF","RSF","USER"};

static char * get_msf_seq(char *sname,sint *len,sint seqno)
/* read the seqno_th. sequence from a PILEUP multiple alignment file */
{
	static char line[MAXLINE+1];
	char *seq = NULL;
	sint i,j,k;
	unsigned char c;

	fseek(fin,0,0); 		/* start at the beginning */

	*len=0;				/* initialise length to zero */
        for(i=0;;i++) {
		if(fgets(line,MAXLINE+1,fin)==NULL) return NULL; /* read the title*/
		if(linetype(line,"//") ) break;		    /* lines...ignore*/
	}

	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(!blankline(line)) {

			for(i=0;i<seqno;i++) fgets(line,MAXLINE+1,fin);
                        for(j=0;j<=strlen(line);j++) if(line[j] != ' ') break;
			for(k=j;k<=strlen(line);k++) if(line[k] == ' ') break;
			strncpy(sname,line+j,MIN(MAXNAMES,k-j)); 
			sname[MIN(MAXNAMES,k-j)]=EOS;
			rtrim(sname);
                       	blank_to_(sname);

			if(seq==NULL)
				seq=(char *)ckalloc((MAXLINE+2)*sizeof(char));
			else
				seq=(char *)ckrealloc(seq,((*len)+MAXLINE+2)*sizeof(char));
			for(i=k;i<MAXLINE;i++) {
				c=line[i];
				if(c == '.' || c == '~' ) c = GAP2;
				if(c == '*') c = 'X';
				if(c == '\n' || c == EOS) break; /* EOL */
				if(isalpha(c)||c==GAP2) seq[(*len)++]=c;
			}

			for(i=0;;i++) {
				if(fgets(line,MAXLINE+1,fin)==NULL) return seq;
				if(blankline(line)) break;
			}
		}
	}
	return seq;
}

static Boolean cl_blankline(char *line)
{
	int i;

	if (line[0] == '!') return TRUE;
	
	for(i=0;line[i]!='\n' && line[i]!=EOS;i++) {
		if( isdigit(line[i]) ||
		    isspace(line[i]) ||
		    (line[i] == '*') ||
		    (line[i] == ':') ||
                    (line[i] == '.')) 
			;
		else
			return FALSE;
	}
	return TRUE;
}

static char * get_clustal_seq(char *sname,sint *len,sint seqno)
/* read the seqno_th. sequence from a clustal multiple alignment file */
{
	static char line[MAXLINE+1];
	static char tseq[MAXLINE+1];
	char *seq = NULL;
	sint i,j;
	unsigned char c;

	fseek(fin,0,0); 		/* start at the beginning */

	*len=0;				/* initialise length to zero */
	fgets(line,MAXLINE+1,fin);	/* read the title line...ignore it */

	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(!cl_blankline(line)) {

			for(i=0;i<seqno;i++) fgets(line,MAXLINE+1,fin);
			for(j=0;j<=strlen(line);j++) if(line[j] != ' ') break;

			sscanf(line,"%s%s",sname,tseq);
			for(j=0;j<MAXNAMES;j++) if(sname[j] == ' ') break;
			sname[j]=EOS;
			rtrim(sname);
                       	blank_to_(sname);

			if(seq==NULL)
				seq=(char *)ckalloc((MAXLINE+2)*sizeof(char));
			else
				seq=(char *)ckrealloc(seq,((*len)+MAXLINE+2)*sizeof(char));
			for(i=0;i<MAXLINE;i++) {
				c=tseq[i];
				/*if(c == '\n' || c == EOS) break;*/ /* EOL */
				if(isspace(c) || c == EOS) break; /* EOL */
				if(isalpha(c)||c==GAP2) seq[(*len)++]=c;
			}

			for(i=0;;i++) {
				if(fgets(line,MAXLINE+1,fin)==NULL) return seq;
				if(cl_blankline(line)) break;
			}
		}
	}

	return seq;
}

static char * get_seq(char *sname,sint *len,char *nid,char *access,char *tit,sint *seq_type,char *org,char *strand,Boolean dnaflag,FTPTR ft)
{
	static char line[MAXLINE+1];
	static char ft_type[12];
	static char ft_name[102];
	sint ft_start,ft_end,ft_pstart,ft_pend;
	char *seq = NULL;
	sint i, j, l, nft, ti, oi, offset = 0;
	sint i1,i2;
        unsigned char c=EOS;
	Boolean got_seq=FALSE;
	Boolean found;
	Boolean old_balibase=FALSE;
	char str1[MAXLINE+1],str2[MAXLINE+1],feature[MAXLINE+1];
	sint  color,start_pos, end_pos;
	sint type=0;

	ti=nft=0;
	switch(seqFormat) {

/************************************/
		case EMBLSWISS:
			while( !linetype(line,"ID") )
				if (fgets(line,MAXLINE+1,fin) == NULL) return NULL;
			
                        for(i=5;i<=strlen(line);i++)  /* DES */
				if(line[i] != ' ') break;
			strncpy(sname,line+i,MAXNAMES); /* remember entryname */
                	for(i=0;i<=strlen(sname);i++)
                        	if(sname[i] == ' ') {
                                	sname[i]=EOS;
                                	break;
                        	}

			sname[MAXNAMES]=EOS;
			rtrim(sname);
                        blank_to_(sname);
			strcpy(nid,sname);

			while( !linetype(line,"AC") )
				if (fgets(line,MAXLINE+1,fin) == NULL) return NULL;
			
                        for(i=5;i<=strlen(line);i++)  /* DES */
				if(line[i] != ' ') break;
			strncpy(access,line+i,MAXNAMES); /* remember access number */
                	for(i=0;i<=strlen(access);i++)
                        	if(access[i] == ' ') {
                                	access[i]=EOS;
                                	break;
                        	}

			access[MAXNAMES]=EOS;
			rtrim(access);
                        blank_to_(access);

						
			while( !linetype(line,"SQ") ) {
				if (linetype(line,"DE") && ti<MAXTITLES) {
                			for(i=5;i<strlen(line)-1;i++) {
						tit[ti++]=line[i];	
						if(ti>=MAXTITLES) break;
					}
				}
				if (linetype(line,"OS") && oi<MAXORGANISMS) {
                			for(i=5;i<strlen(line)-1;i++) {
						org[oi++]=line[i];	
						if(oi>=MAXORGANISMS) break;
					}
				}
				else if (linetype(line,"FT") ) {
                			for(i=5,j=0;i<13;i++) 
						if(!isspace(line[i])) ft_type[j++]=line[i];
					ft_type[j]='\0';
					ft_name[0]='\0';
					strcpy(ft_name,&line[34]);
					if(check_ft_type(ft_type,ft_name,&type)==1) {
						nft=ft->nentries[type];
						ft_start=atoi(&line[13]);
						ft_end=atoi(&line[20]);
						if(ft_start>0 && ft_end>0) {
							if(nft>=MAXFT) {
								warning("too many features for %s\n",sname);
							}
							else {
								alloc_ft_entry(&ft->data[type][nft]);
								strcpy(ft->data[type][nft].type,ft_type);
								strcpy(ft->data[type][nft].name,ft_name);
								l=strlen(ft->data[type][nft].name);
								for(i=33;i<strlen(line)-1;i++)
									ft->data[type][nft].name[l++]=line[i];
								ft->data[type][nft].start=ft_start-1;
								ft->data[type][nft].end=ft_end-1;
								nft++;
								ft->nentries[type]=nft;
							}
						}
					}
				}
				if (fgets(line,MAXLINE+1,fin) == NULL) return NULL;
			}
			tit[ti]='\0';
				
			*len=0;
			while(fgets(line,MAXLINE+1,fin)) {
				/* if(got_seq && blankline(line)) break; */
				if( linetype(line,"ID") ) break;
 				if( strlen(line) > 2 && line[strlen(line)-2]=='.' && line[strlen(line)-3]=='.' ) 
					continue;
				if(seq==NULL)
					seq=(char *)ckalloc((MAXLINE+2)*sizeof(char));
				else
					seq=(char *)ckrealloc(seq,((*len)+MAXLINE+2)*sizeof(char));
				for(i=0;i<MAXLINE;i++) {
					c=line[i];
				if(c == '\n' || c == EOS || c == '/')
					break;			/* EOL */
				if(isalpha(c)||c==GAP2) {
					got_seq=TRUE;
					seq[(*len)++]=c;
				}
				}
			if(c == '/') break;
			}
		break;
		
/************************************/
		case PIR:
			while(*line != '>')
				fgets(line,MAXLINE+1,fin);			
                        for(i=4;i<=strlen(line);i++)  /* DES */
				if(line[i] != ' ') break;
			strncpy(sname,line+i,MAXNAMES); /* remember entryname */
			sname[MAXNAMES]=EOS;
			rtrim(sname);
                        blank_to_(sname);
			strcpy(nid,sname);
			strcpy(access,sname);

			fgets(line,MAXLINE+1,fin);
			strncpy(tit,line,MAXTITLES);
			tit[MAXTITLES]=EOS;
			i=strlen(tit);
			if(tit[i-1]=='\n') tit[i-1]=EOS;
			
			*len=0;
			while(fgets(line,MAXLINE+1,fin)) {
				if(seq==NULL)
					seq=(char *)ckalloc((MAXLINE+2)*sizeof(char));
				else
					seq=(char *)ckrealloc(seq,((*len)+MAXLINE+2)*sizeof(char));
				for(i=0;i<MAXLINE;i++) {
					c=line[i];
				if(c == '\n' || c == EOS || c == '*')
					break;			/* EOL */
			
				if(isalpha(c)||c==GAP2) seq[(*len)++]=c;
				}
			if(c == '*') break;
			}
		break;
/***********************************************/
		case PEARSON:
			while(*line != '>')
				fgets(line,MAXLINE+1,fin);
			
                        for(i=1;i<=strlen(line);i++)  /* DES */
				if(line[i] != ' ') break;
			strncpy(sname,line+i,MAXNAMES); /* remember entryname */
                        for(i=1;i<=strlen(sname);i++)  /* DES */
				if(sname[i] == ' ') break;
			sname[i]=EOS;
			rtrim(sname);
                        blank_to_(sname);
			strcpy(nid,sname);
			strcpy(access,sname);

			*tit=EOS;
			
			*len=0;
			while(fgets(line,MAXLINE+1,fin)) {
				if(seq==NULL)
					seq=(char *)ckalloc((MAXLINE+2)*sizeof(char));
				else
					seq=(char *)ckrealloc(seq,((*len)+MAXLINE+2)*sizeof(char));
				for(i=0;i<MAXLINE;i++) {
					c=line[i];
				if(c == '\n' || c == EOS || c == '>')
					break;			/* EOL */
			
				if(isalpha(c)||c==GAP2) seq[(*len)++]=c;
			}
			if(c == '>') break;
			}
		break;
/**********************************************/
		case GDE:
			if (dnaflag) {
				while(*line != '#')
					fgets(line,MAXLINE+1,fin);
			}
			else {
				while(*line != '%')
					fgets(line,MAXLINE+1,fin);
			}
			
			for (i=1;i<=MAXNAMES;i++) {
				if (line[i] == '(' || line[i] == '\n')
                                    break;
				sname[i-1] = line[i];
			}
			i--;
			sname[i]=EOS;
			if (sname[i-1] == '(') sscanf(&line[i],"%d",&offset);
			else offset = 0;
			for(i--;i > 0;i--) 
				if(isspace(sname[i])) {
					sname[i]=EOS;	
				}
				else break;		
                        blank_to_(sname);
			strcpy(nid,sname);
			strcpy(access,sname);

			*tit=EOS;
			*org=EOS;
			
			*len=0;
			for (i=0;i<offset;i++) seq[(*len)++] = GAP2;
			while(fgets(line,MAXLINE+1,fin)) {
			if(*line == '%' || *line == '#' || *line == '"') break;
				if(seq==NULL)
					seq=(char *)ckalloc((MAXLINE+2)*sizeof(char));
				else
					seq=(char *)ckrealloc(seq,((*len)+MAXLINE+2)*sizeof(char));
				for(i=0;i<MAXLINE;i++) {
					c=line[i];
				if(c == '\n' || c == EOS) 
					break;			/* EOL */
			
				if(isalpha(c)||c==GAP2) seq[(*len)++]=c;
				}
			}
		break;
/***********************************************/
		case RSF:
			while(*line != '{')
				if (fgets(line,MAXLINE+1,fin) == NULL) return NULL;
			
			while( !keyword(line,"name") )
				if (fgets(line,MAXLINE+1,fin) == NULL) return NULL;
			
                        for(i=5;i<=strlen(line);i++)  /* DES */
				if(line[i] != ' ') break;
			strncpy(sname,line+i,MAXNAMES); /* remember entryname */
                	for(i=0;i<=strlen(sname);i++)
                        	if(sname[i] == ' ') {
                                	sname[i]=EOS;
                                	break;
                        	}

			sname[MAXNAMES]=EOS;
			rtrim(sname);
                        blank_to_(sname);
			strcpy(nid,sname);
			strcpy(access,sname);

						
			while( !keyword(line,"sequence") ) {
				if (keyword(line,"descrip") ) {
					for(i=8;i<strlen(line)-1;i++)
						tit[i-8]=line[i];
					tit[i-8]='\0';
				}
				else if (keyword(line,"strand") ) {
					for(i=8;i<strlen(line)-1;i++)
						strand[i-8]=line[i];
					strand[i-8]='\0';
				}
				else if (keyword(line,"creator") ) {
					for(i=8;i<strlen(line)-1;i++)
						org[i-8]=line[i];
					org[i-8]='\0';
				}
				else if (keyword(line,"feature") ) {
					if (sscanf(&line[8],"%d%d%d%s%s%s",&ft_start,&ft_end,&color,str1,str2,feature) == 6) {
/* for old balibase files */
						if(keyword(str1,"circle")) {
							old_balibase=TRUE;
                                        		strcpy(ft_type,"STRUCT");
                                        		strcpy(ft_name,"HELIX");
							type=STRUCT;
                                        		nft=ft->nentries[type];
                                        		alloc_ft_entry(&ft->data[type][nft]);
                                        		strcpy(ft->data[type][nft].type,ft_type);
                                        		strcpy(ft->data[type][nft].name,ft_name);
                                        		ft->data[type][nft].color=3;
                                        		ft->data[type][nft].start=ft_start-1;
                                        		ft->data[type][nft].end=ft_end-1;
                                        		nft++;
                                        		ft->nentries[type]=nft;
                                		}
                                		else if(keyword(str1,"right_arrow")) {
							old_balibase=TRUE;
                                        		strcpy(ft_type,"STRUCT");
                                        		strcpy(ft_name,"STRAND");
							type=STRUCT;
                                        		nft=ft->nentries[type];
                                        		alloc_ft_entry(&ft->data[type][nft]);
                                        		strcpy(ft->data[type][nft].type,ft_type);
                                        		strcpy(ft->data[type][nft].name,ft_name);
                                        		ft->data[type][nft].start=ft_start-1;
                                        		ft->data[type][nft].end=ft_end-1;
                                        		ft->data[type][nft].color=0;
                                        		nft++;
                                        		ft->nentries[type]=nft;
						}
						ft_name[0]='\0';
						if(check_ft_type(feature,ft_name,&type)==1) {
							nft=ft->nentries[type];
							if(ft_start>0 && ft_end>0) {
								if(nft>=MAXFT) {
									warning("too many features for %s\n",sname);
								}
								else {
/* get the line starting "FT" */

									if (fgets(line,MAXLINE+1,fin) == NULL) return NULL;
									for(i=31;line[i]!='\n' && line[i]!='\0';i++) 
										ft_name[i-31]=line[i];
									ft_name[i-31]='\0';
									alloc_ft_entry(&ft->data[type][nft]);
									strcpy(ft->data[type][nft].type,feature);
									strcpy(ft->data[type][nft].name,ft_name);
/* remember that the start and end entries here refer to alignment columns not residue positions, so need to
convert later when we've got the sequence */
									ft->data[type][nft].start=ft_start-1;
									ft->data[type][nft].end=ft_end-1;
									ft->data[type][nft].color=color;
									nft++;
									ft->nentries[type]=nft;
								}
							}
						}
					}
				}
				if (fgets(line,MAXLINE+1,fin) == NULL) return NULL;
			}
	
			*len=0;
			while(fgets(line,MAXLINE+1,fin)) {
				if(seq==NULL)
					seq=(char *)ckalloc((MAXLINE+2)*sizeof(char));
				else
					seq=(char *)ckrealloc(seq,((*len)+MAXLINE+2)*sizeof(char));
				for(i=0;i<MAXLINE;i++) {
					c=line[i];
				if(c == EOS || c == '}')
					break;			/* EOL */
				if( c=='.') c=GAP2;
				if(isalpha(c)||c==GAP2) seq[(*len)++]=c;
				}
				if(c == '}') break;
			}
		break;
/***********************************************/
	}
	
	seq[*len+1]=EOS;

/* convert column positions to residues */
	if(seqFormat==RSF) {
		for(i=0;i<MAXFTTYPE;i++) {
			for(j=0;j<ft->nentries[i];j++) {
				col2pos(seq,ft->data[i][j].start,ft->data[i][j].end,&ft_pstart,&ft_pend);
                		ft->data[i][j].start=ft_pstart;
                		ft->data[i][j].end=ft_pend;
			}
		}
	}

		
	(*seq_type)=type;
	return seq;
}

/*first_seq is the #no. of the first seq. to read */
sint readseqs(char *filename,sint prf_no,sint first_seq,sint explicit_type,ALNPTR mult_aln)
{
	static char *seq1,sname1[MAXNAMES+1],title[MAXTITLES+1],organism[MAXORGANISMS+1];
	char strand[5];
	char nid[MAXNAMES+1],access[MAXNAMES+1];
	char tstr[MAXTITLES+1];
	char tmp[MAXTITLES+1];
	sint i,j,k;
	sint no_seqs;
	sint seq_type;
	static sint l1;
	static sint max_aln_length;
	static Boolean dnaflag1;
	Boolean found;
	
	if(*filename == EOS) return -1;
	if((fin=fopen(filename,"r"))==NULL) {
		error("Could not open sequence file %s",filename);
		return -1;      /* DES -1 => file not found */
	}
	no_seqs=0;

	check_infile(&no_seqs,&mult_aln->dnaflag);
	/*info("Sequence format is %s",formatNames[seqFormat]);*/

/* DES DEBUG 
	fprintf(stdout,"\n\n File name = %s\n\n",filename);
*/
	if(no_seqs == 0)
		return 0;       /* return the number of seqs. (zero here)*/

/* if this is a multiple alignment, or profile 1 - free any memory used
by previous alignments, then allocate memory for the new alignment */
	if(first_seq == 0) {
		free_aln(mult_aln);
		alloc_aln(no_seqs,mult_aln);
	}
/* otherwise, this is a profile 2, and we need to reallocate the arrays,
leaving the data for profile 1 intact */
	else {
		realloc_aln(first_seq,no_seqs,mult_aln);
	}
	if (seqFormat == RELACS) {
	}
	else {
		for(i=first_seq;i<first_seq+no_seqs;i++) {    /* get the seqs now*/
                	mult_aln->seqs[i].output_index=i;
			seq_type=0;
			if(seqFormat == CLUSTAL) {
				seq1=get_clustal_seq(sname1,&l1,i-first_seq);
				strcpy(nid,sname1);
				strcpy(access,sname1);
			}
			else if(seqFormat == MSF){
			    	seq1=get_msf_seq(sname1,&l1,i-first_seq);
				strcpy(nid,sname1);
				strcpy(access,sname1);
			}
			else
				seq1=get_seq(sname1,&l1,nid,access,title,&seq_type,organism,strand,mult_aln->dnaflag,&mult_aln->ft[i]);
	
			if(seq1==NULL) break;
/* JULIE */
/*  Set max length of dynamically allocated arrays in prfalign.c */
			if (l1 > max_aln_length) max_aln_length = l1;
			mult_aln->seqs[i].len=l1;                   /* store the length */
			strcpy(mult_aln->seqs[i].name,sname1);              /*    "   "  name   */
			strcpy(mult_aln->seqs[i].nid,nid);              /*    "   "  identity  */
			strcpy(mult_aln->seqs[i].access,access);              /*    "   "  access  */
			strcpy(mult_aln->seqs[i].title,title);              /*    "   "  title  */
			strcpy(mult_aln->seqs[i].org,organism);              /*    "   "  organism  */
			mult_aln->seqs[i].sense=atoi(strand);
			mult_aln->seqs[i].type=seq_type;

			if(explicit_type==0) {
				dnaflag1 = check_dnaflag(seq1,l1); /* check DNA/Prot */
		        	if(i == 1) mult_aln->dnaflag = dnaflag1;
			}			/* type decided by first seq*/
			else if(explicit_type==1)
				dnaflag1 = TRUE;
			else
				dnaflag1 = FALSE;

			alloc_seq(&mult_aln->seqs[i],l1);

			for(j=0;j<l1;j++) {
				if(isalpha(seq1[j])) mult_aln->seqs[i].data[j]=tolower(seq1[j]);
				else mult_aln->seqs[i].data[j]=GAP2;
			}
			if(seq1!=NULL) seq1=ckfree(seq1);

			mult_aln->seqs[i].reslen=0;
			for(j=0;j<mult_aln->seqs[i].len;j++)
				if(isalpha(mult_aln->seqs[i].data[j])) mult_aln->seqs[i].reslen++;
		}
	}

	max_aln_length *= 2;
/*
   JULIE
   check sequence names are all different - otherwise phylip tree is 
   confused.
	for(i=0;i<first_seq+no_seqs;i++) {
		for(j=i+1;j<first_seq+no_seqs;j++) {
			if (strncmp(mult_aln->seqs[i].name,mult_aln->seqs[j].name,MAXNAMES) == 0) {
				error("Multiple sequences found with same name, %s (first %d chars are significant) seq1:%d %s seq2: %d %s", mult_aln->seqs[i].name,MAXNAMES,i,mult_aln->seqs[i].name,j,mult_aln->seqs[j].name);
				return 0;
			}
		}
	}
*/
	for(i=first_seq;i<first_seq+no_seqs;i++)
	{
		if(mult_aln->seqs[i].len>max_aln_length)
			max_aln_length=mult_aln->seqs[i].len;
	}
	
	fclose(fin);
			
	return no_seqs;    /* return the number of seqs. read in this call */
}

static Boolean check_dnaflag(char *seq, sint slen)
/* check if DNA or Protein
   The decision is based on counting all A,C,G,T,U or N. 
   If >= 85% of all characters (except -) are as above => DNA  */
{
	sint i, c, nresidues, nbases;
	float ratio;
	char *dna_codes="ACGTUN";
	
	nresidues = nbases = 0;	
	for(i=0; i < slen; i++) {
		if(seq[i] != GAP2) {
			nresidues++;
			if(seq[i] == 'N')
				nbases++;
			else {
				c = res_index(dna_codes, seq[i]);
				if(c >= 0)
					nbases++;
			}
		}
	}
	if( (nbases == 0) || (nresidues == 0) ) return FALSE;
	ratio = (float)nbases/(float)nresidues;
/* DES 	fprintf(stdout,"\n nbases = %d, nresidues = %d, ratio = %f\n",
		(pint)nbases,(pint)nresidues,(pint)ratio); */
	if(ratio >= 0.85) 
		return TRUE;
	else
		return FALSE;
}



static void check_infile(sint *nseqs,Boolean *dnaflag)
{
	char line[MAXLINE+1];
	sint i;	

	*nseqs=0;
	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(!blankline(line)) 
			break;
	}

	for(i=strlen(line)-1;i>=0;i--)
		if(isgraph(line[i])) break;
	line[i+1]=EOS;
        
	for(i=0;i<=6;i++) line[i] = toupper(line[i]);

	if( linetype(line,"ID") ) {					/* EMBL/Swiss-Prot format ? */
		seqFormat=EMBLSWISS;
		(*nseqs)++;
	}
        else if( linetype(line,"CLUSTAL") ) {
		seqFormat=CLUSTAL;
	}
 	else if( linetype(line,"PILEUP") ) {
		seqFormat = MSF;
	}
 	else if( linetype(line,"!!AA_MULTIPLE_ALIGNMENT") ) {
		seqFormat = MSF;
		*dnaflag = FALSE;
	}
 	else if( linetype(line,"!!NA_MULTIPLE_ALIGNMENT") ) {
		seqFormat = MSF;
		*dnaflag = TRUE;
	}
 	else if( strstr(line,"MSF") && line[strlen(line)-1]=='.' &&
                                 line[strlen(line)-2]=='.' ) {
		seqFormat = MSF;
	}
 	else if( linetype(line,"!!RICH_SEQUENCE") ) {
		seqFormat = RSF;
	}
 	else if( linetype(line,"<?XML") ) {
		seqFormat = RELACS;
	}
	else if(*line == '>') {						/* no */
		seqFormat=(line[3] == ';')?PIR:PEARSON; /* distinguish PIR and Pearson */
		(*nseqs)++;
	}
	else if((*line == '"') || (*line == '%') || (*line == '#')) {
		seqFormat=GDE; /* GDE format */
		if (*line == '%') {
                        (*nseqs)++;
			*dnaflag = FALSE;
		}
		else if (*line == '#') {
			(*nseqs)++;
			*dnaflag = TRUE;
		}
	}
	else {
		seqFormat=UNKNOWN;
		return;
	}

	while(fgets(line,MAXLINE+1,fin) != NULL) {
		switch(seqFormat) {
			case EMBLSWISS:
				if( linetype(line,"ID") )
					(*nseqs)++;
				break;
			case PIR:
				*nseqs = count_pir_seqs();
				fseek(fin,0,0);
				return;
			case PEARSON:
                                if( *line == '>' )
                                        (*nseqs)++;
                                break;
			case GDE:
				if(( *line == '%' ) && ( *dnaflag == FALSE))
					(*nseqs)++;
				else if (( *line == '#') && ( *dnaflag == TRUE))
					(*nseqs)++;
				break;
			case CLUSTAL:
				*nseqs = count_clustal_seqs();
/* DES */ 			/* fprintf(stdout,"\nnseqs = %d\n",(pint)*nseqs); */
				fseek(fin,0,0);
				return;
			case MSF:
				*nseqs = count_msf_seqs();
				fseek(fin,0,0);
				return;
			case RSF:
				fseek(fin,0,0);
				*nseqs = count_rsf_seqs();
				fseek(fin,0,0);
				return;
			case RELACS:
				fseek(fin,0,0);
				fseek(fin,0,0);
				return;
			case USER:
			default:
				break;
		}
	}
	fseek(fin,0,0);
}


static sint count_pir_seqs(void)
/* count the number of sequences in a pir alignment file */
{
	char line[MAXLINE+1],c;
	sint  nseqs, i;
	Boolean seq_ok;

	seq_ok = FALSE;
	while (fgets(line,MAXLINE+1,fin) != NULL) { /* Look for end of first seq */
		if(*line == '>') break;
		for(i=0;seq_ok == FALSE;i++) {
			c=line[i];
			if(c == '*') {
				seq_ok = TRUE;	/* ok - end of sequence found */
				break;
			}			/* EOL */
			if(c == '\n' || c == EOS)
				break;			/* EOL */
		}
		if (seq_ok == TRUE)
			break;
	}
	if (seq_ok == FALSE) {
		error("PIR format sequence end marker '*'\nmissing for one or more sequences.");
		return (sint)0;	/* funny format*/
	}
	
	
	nseqs = 1;
	
	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(*line == '>') {		/* Look for start of next seq */
			seq_ok = FALSE;
			while (fgets(line,MAXLINE+1,fin) != NULL) { /* Look for end of seq */
				if(*line == '>') {
					error("PIR format sequence end marker '*' missing for one or more sequences.");
					return (sint)0;	/* funny format*/
				}
				for(i=0;seq_ok == FALSE;i++) {
					c=line[i];
					if(c == '*') {
						seq_ok = TRUE;	/* ok - sequence found */
						break;
					}			/* EOL */
					if(c == '\n' || c == EOS)
						break;			/* EOL */
				}
				if (seq_ok == TRUE) {
					nseqs++;
					break;
				}
			}
		}
	}
	return (sint)nseqs;
}


static sint count_clustal_seqs(void)
/* count the number of sequences in a clustal alignment file */
{
	char line[MAXLINE+1];
	sint  nseqs;

	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(!cl_blankline(line)) break;		/* Look for next non- */
	}						/* blank line */
	nseqs = 1;

	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(cl_blankline(line)) return nseqs;
		nseqs++;
	}

	return (sint)0;	/* if you got to here-funny format/no seqs.*/
}

static sint count_msf_seqs(void)
{
/* count the number of sequences in a PILEUP alignment file */

	char line[MAXLINE+1];
	sint  nseqs;

	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(linetype(line,"//")) break;
	}

	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(!blankline(line)) break;		/* Look for next non- */
	}						/* blank line */
	nseqs = 1;

	while (fgets(line,MAXLINE+1,fin) != NULL) {
		if(blankline(line)) return nseqs;
		nseqs++;
	}

	return (sint)0;	/* if you got to here-funny format/no seqs.*/
}

static sint count_rsf_seqs(void)
{
/* count the number of sequences in a GCG RSF alignment file */

	char line[MAXLINE+1];
	sint  nseqs;

	nseqs = 0;
/* skip the comments */
	while (fgets(line,MAXLINE+1,fin) != NULL) {
 		if(line[strlen(line)-2]=='.' &&
                                 line[strlen(line)-3]=='.')
			break;
	}

	while (fgets(line,MAXLINE+1,fin) != NULL) {
                if( *line == '{' )
                      nseqs++;
	}
	return (sint)nseqs;
}

static sint res_index(char *t,char c)
{
	register sint i;
	
	for(i=0;t[i] && t[i] != c;i++)
		;
	if(t[i]) return(i);
	else return -1;
}

sint seq_input(char *filename,sint explicit_type,Boolean append,ALNPTR mult_aln)
{
        sint i;
	sint local_nseqs;
	sint max_names;


        mult_aln->dnaflag=FALSE;
        strcpy(mult_aln->filename,"");
        mult_aln->prf1.nseqs=0;
        strcpy(mult_aln->prf1.filename,"");
        mult_aln->prf2.nseqs=0;
        strcpy(mult_aln->prf2.filename,"");
       if (append)
          	local_nseqs = readseqs(filename,(sint)0,mult_aln->nseqs,explicit_type,mult_aln);
       else {
        	mult_aln->nseqs=0;
          	local_nseqs = readseqs(filename,(sint)0,(sint)0,explicit_type,mult_aln);  
	}
       if(local_nseqs < 0)               /* file could not be opened */
           { 
		return local_nseqs;
           }
       else if(local_nseqs == 0)         /* no sequences */
           {
	       error("No sequences in file!  Bad format?");
               return local_nseqs;
           }
       else 
           {
		strcpy(mult_aln->filename,filename);
		if(append) mult_aln->nseqs+=local_nseqs;
		else mult_aln->nseqs=local_nseqs;
		/*info("Sequences assumed to be %s",
			mult_aln->dnaflag?"DNA":"PROTEIN");
		info("\n\n");*/
		max_names=0;
               	for(i=0; i<mult_aln->nseqs; i++) {
			if(strlen(mult_aln->seqs[i].name)>max_names)
				max_names=strlen(mult_aln->seqs[i].name);
		}
		if(max_names<10) max_names=10;
               	/*for(i=0; i<mult_aln->nseqs; i++) {
                       	info("Sequence %d: %-*s   %6.d %s",
                       	(pint)i+1,max_names,mult_aln->seqs[i].name,(pint)mult_aln->seqs[i].len,mult_aln->dnaflag?"bp":"aa");
               	}*/
	   }
	return local_nseqs;	
}


sint profile_input(char *filename,sint prf_no,sint explicit_type,ALNPTR mult_aln)   /* read a profile   */
{                                           /* prf_no is 1 or 2  */
        sint local_nseqs, i;
	sint max_names;
	
        if(prf_no == 2 && mult_aln->prf1.nseqs<=0) 
           {
             error("You must read in profile number 1 first");
             return 0;
           }

        mult_aln->dnaflag=FALSE;
        mult_aln->nseqs=0;
        strcpy(mult_aln->filename,"");
    if(prf_no == 1)     /* for the 1st profile */
      {
        mult_aln->prf1.nseqs=0;
        strcpy(mult_aln->prf1.filename,"");
        mult_aln->prf2.nseqs=0;
        strcpy(mult_aln->prf2.filename,"");
       local_nseqs = readseqs(filename,(sint)1,(sint)0,explicit_type,mult_aln); /* (1) means 1st seq to be read = no. 1 */
       if(local_nseqs == 0)         /* no sequences  */
           {
	       error("No sequences in file!  Bad format?");
           }
       else if (local_nseqs > 0)
           { 				/* success; found some seqs. */
		strcpy(mult_aln->filename,filename);
		strcpy(mult_aln->prf1.filename,filename);
                mult_aln->nseqs = mult_aln->prf1.nseqs = local_nseqs;
		info("No. of seqs=%d",(pint)mult_aln->nseqs);
	   }
      }
    else
      {			        /* first seq to be read = profile1_nseqs + 1 */
        mult_aln->prf2.nseqs=0;
        strcpy(mult_aln->prf2.filename,"");
       local_nseqs = readseqs(filename,(sint)2,mult_aln->prf1.nseqs,explicit_type,mult_aln); 
       if(local_nseqs == 0)         /* no sequences */
           {
	       error("No sequences in file!  Bad format?");
           }
       else if(local_nseqs > 0)
           {
		info("No. of seqs in profile=%d",(pint)local_nseqs);
		strcpy(mult_aln->filename,filename);
		strcpy(mult_aln->prf2.filename,filename);
                mult_aln->nseqs = mult_aln->prf1.nseqs + local_nseqs;
		mult_aln->prf2.nseqs = local_nseqs;
                info("Total no. of seqs     =%d",(pint)mult_aln->nseqs);
	   }

      }

	if(local_nseqs<=0) return local_nseqs;
	
	/*info("Sequences assumed to be %s",
		mult_aln->dnaflag?"DNA":"PROTEIN");
	info("\n\n");*/
	max_names=0;
        for(i=(mult_aln->prf2.nseqs==0)?0:mult_aln->prf1.nseqs; i<mult_aln->nseqs; i++) {
		if(strlen(mult_aln->seqs[i].name)>max_names)
			max_names=strlen(mult_aln->seqs[i].name);
	}
	if(max_names<10) max_names=10;
        /*for(i=(mult_aln->prf2.nseqs==0)?0:mult_aln->prf1.nseqs; i<mult_aln->nseqs; i++) {
               	info("Sequence %d: %-*s   %6.d %s",
                  	(pint)i+1,max_names,mult_aln->seqs[i].name,(pint)mult_aln->seqs[i].len,mult_aln->dnaflag?"bp":"aa");
	}*/

	return local_nseqs;
}


