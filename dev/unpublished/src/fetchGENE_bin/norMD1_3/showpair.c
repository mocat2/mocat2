#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "score.h"	

static void make_p_ptrs(int *tptr, int *pl, char *seq1, int l);
static void put_frag(int fs, int v1, int v2, int flen);
static int frag_rel_pos(int a1, int b1, int a2, int b2);
static void des_quick_sort(int *array1, int *array2, int array_size);
static void pair_align(char *seq, int l1, int l2);

extern void *ckfree(void *ptr);
extern void *ckalloc(size_t bytes);



/*
*	 Global variables
*/

static int 	next;
static int 	curr_frag,maxsf,vatend;
static int 	**accum;
static int 	*diag_index;
static char 	*slopes;

int ktup=1,window=5,wind_gap=3,signif=5;    		      /* Pairwise aln. params */
int *displ;
int *zza, *zzb, *zzc, *zzd;
int max_aln_length;

Boolean percent=TRUE;


static void make_p_ptrs(int *tptr,int *pl,char *seq,int l)
{
	static int a[10];
	int i,j,limit,code,flag;
	char residue;
	
	for (i=1;i<=ktup;i++)
           a[i] = (int) pow((double)(27),(double)(i-1));

	limit = (int) pow((double)(27),(double)ktup);
	for(i=1;i<=limit;++i)
		pl[i]=0;
	for(i=1;i<=l;++i)
		tptr[i]=0;
	
	for(i=1;i<=(l-ktup+1);++i) {
		code=0;
		flag=FALSE;
		for(j=1;j<=ktup;++j) {
			residue = seq[i+j-1];
			if((residue<0) || (residue > 26)){
				flag=TRUE;
				break;
			}
			code += ((residue) * a[j]);
		}
		if(flag)
			continue;
		++code;
		if(pl[code]!=0)
			tptr[i]=pl[code];
		pl[code]=i;
	}
}


static void put_frag(int fs,int v1,int v2,int flen)
{
	int end;
	accum[0][curr_frag]=fs;
	accum[1][curr_frag]=v1;
	accum[2][curr_frag]=v2;
	accum[3][curr_frag]=flen;
	
	if(!maxsf) {
		maxsf=1;
		accum[4][curr_frag]=0;
		return;
	}
	
        if(fs >= accum[0][maxsf]) {
		accum[4][curr_frag]=maxsf;
		maxsf=curr_frag;
		return;
	}
	else {
		next=maxsf;
		while(TRUE) {
			end=next;
			next=accum[4][next];
			if(fs>=accum[0][next])
				break;
		}
		accum[4][curr_frag]=next;
		accum[4][end]=curr_frag;
	}
}


static int frag_rel_pos(int a1,int b1,int a2,int b2)
{
	int ret;
	
	ret=FALSE;
	if(a1-b1==a2-b2) {
		if(a2<a1)
			ret=TRUE;
	}
	else {
		if(a2+ktup-1<a1 && b2+ktup-1<b1)
			ret=TRUE;
	}
	return ret;
}


static void des_quick_sort(int *array1, int *array2, int array_size)
/*  */
/* Quicksort routine, adapted from chapter 4, page 115 of software tools */
/* by Kernighan and Plauger, (1986) */
/* Sort the elements of array1 and sort the */
/* elements of array2 accordingly */
/*  */
{
	int temp1, temp2;
	int p, pivlin;
	int i, j;
	int lst[50], ust[50];       /* the maximum no. of elements must be*/
								/* < log(base2) of 50 */

	lst[1] = 1;
	ust[1] = array_size-1;
	p = 1;

	while(p > 0) {
		if(lst[p] >= ust[p])
			p--;
		else {
			i = lst[p] - 1;
			j = ust[p];
			pivlin = array1[j];
			while(i < j) {
				for(i=i+1; array1[i] < pivlin; i++)
					;
				for(j=j-1; j > i; j--)
					if(array1[j] <= pivlin) break;
				if(i < j) {
					temp1     = array1[i];
					array1[i] = array1[j];
					array1[j] = temp1;
					
					temp2     = array2[i];
					array2[i] = array2[j];
					array2[j] = temp2;
				}
			}
			
			j = ust[p];

			temp1     = array1[i];
			array1[i] = array1[j];
			array1[j] = temp1;

			temp2     = array2[i];
			array2[i] = array2[j];
			array2[j] = temp2;

			if(i-lst[p] < ust[p] - i) {
				lst[p+1] = lst[p];
				ust[p+1] = i - 1;
				lst[p]   = i + 1;
			}
			else {
				lst[p+1] = i + 1;
				ust[p+1] = ust[p];
				ust[p]   = i - 1;
			}
			p = p + 1;
		}
	}
	return;

}





static void pair_align(char *seq,int l1,int l2)
{
	int pot[8],i,j,l,m,flag,limit,pos,tl1,vn1,vn2,flen,osptr,fs;
	int tv1,tv2,encrypt,subt1,subt2,rmndr;
	char residue;
	
	for (i=1;i<=ktup;i++)
       		pot[i] = (int) pow((double)(27),(double)(i-1));
	limit = (int) pow((double)(27),(double)ktup);
	
	tl1 = (l1+l2)-1;
	
	for(i=1;i<=tl1;++i) {
		slopes[i]=displ[i]=0;
		diag_index[i] = i;
	}
	

/* increment diagonal score for each k_tuple match */

	for(i=1;i<=limit;++i) {
		vn1=zzc[i];
		while(TRUE) {
			if(!vn1) break;
			vn2=zzd[i];
			while(vn2 != 0) {
				osptr=vn1-vn2+l2;
				++displ[osptr];
				vn2=zzb[vn2];
			}
			vn1=zza[vn1];
		}
	}

/* choose the top SIGNIF diagonals */

	des_quick_sort(displ, diag_index, tl1);

	j = tl1 - signif + 1;
	if(j < 1) j = 1;
 
/* flag all diagonals within WINDOW of a top diagonal */

	for(i=tl1; i>=j; i--) 
		if(displ[i] > 0) {
			pos = diag_index[i];
			l = (1  >pos-window) ? 1   : pos-window;
			m = (tl1<pos+window) ? tl1 : pos+window;
			for(; l <= m; l++) 
				slopes[l] = 1;
		}

	for(i=1; i<=tl1; i++)  displ[i] = 0;

	
	curr_frag=maxsf=0;
	
	for(i=1;i<=(l1-ktup+1);++i) {
		encrypt=flag=0;
		for(j=1;j<=ktup;++j) {
			residue = seq[i+j-1];
			if((residue<0) || (residue>26)) {
				flag=TRUE;
				break;
			}
			encrypt += ((residue)*pot[j]);
		}
		if(flag) continue;
		++encrypt;
	
		vn2=zzd[encrypt];
	
		flag=FALSE;
		while(TRUE) {
			if(!vn2) {
				flag=TRUE;
				break;
			}
			osptr=i-vn2+l2;
			if(slopes[osptr]!=1) {
				vn2=zzb[vn2];
				continue;
			}
			flen=0;
			fs=ktup;
			next=maxsf;		
		
		/*
		* A-loop
		*/
		
			while(TRUE) {
				if(!next) {
					++curr_frag;
					if(curr_frag>=2*max_aln_length) {
						vatend=1;
						return;
					}
					displ[osptr]=curr_frag;
					put_frag(fs,i,vn2,flen);
				}
				else {
					tv1=accum[1][next];
					tv2=accum[2][next];
					if(frag_rel_pos(i,vn2,tv1,tv2)) {
						if(i-vn2==accum[1][next]-accum[2][next]) {
							if(i>accum[1][next]+(ktup-1))
								fs=accum[0][next]+ktup;
							else {
								rmndr=i-accum[1][next];
								fs=accum[0][next]+rmndr;
							}
							flen=next;
							next=0;
							continue;
						}
						else {
							if(displ[osptr]==0)
								subt1=ktup;
							else {
								if(i>accum[1][displ[osptr]]+(ktup-1))
									subt1=accum[0][displ[osptr]]+ktup;
								else {
									rmndr=i-accum[1][displ[osptr]];
									subt1=accum[0][displ[osptr]]+rmndr;
								}
							}
							subt2=accum[0][next]-wind_gap+ktup;
							if(subt2>subt1) {
								flen=next;
								fs=subt2;
							}
							else {
								flen=displ[osptr];
								fs=subt1;
							}
							next=0;
							continue;
						}
					}
					else {
						next=accum[4][next];
						continue;
					}
				}
				break;
			}
		/*
		* End of Aloop
		*/
		
			vn2=zzb[vn2];
		}
	}
	vatend=0;
}		 

float show_pair(ALN mult_aln,sint i,sint j)
{
	int k,dsr;
	char *seq1,*seq2;
	int len1,len2;
	float calc_score;
	
	max_aln_length=0;
	for(k=0;k<mult_aln.nseqs;++k) 
		if(max_aln_length<mult_aln.seqs[k].len) max_aln_length=mult_aln.seqs[k].len;
	max_aln_length*=2;
	if(max_aln_length<30) max_aln_length=30;

	seq1 = (char *)ckalloc( (2*max_aln_length+1) *sizeof (char) );
	seq2 = (char *)ckalloc( (2*max_aln_length+1) *sizeof (char) );

	accum = (int **)ckalloc( 5*sizeof (int *) );
	for (k=0;k<5;k++)
		accum[k] = (int *) ckalloc((2*max_aln_length+1) * sizeof (int) );

	displ      = (int *) ckalloc( (2*max_aln_length +1) * sizeof (int) );
	slopes     = (char *)ckalloc( (2*max_aln_length +1) * sizeof (char));
	diag_index = (int *) ckalloc( (2*max_aln_length +1) * sizeof (int) );

	zza = (int *)ckalloc( (max_aln_length+1) * sizeof (int) );
	zzb = (int *)ckalloc( (max_aln_length+1) * sizeof (int) );

	zzc = (int *)ckalloc( (max_aln_length+1) * sizeof (int) );
	zzd = (int *)ckalloc( (max_aln_length+1) * sizeof (int) );

		len1=1;
		for(k=0;k<mult_aln.seqs[i].len;k++)
			if(isalpha(mult_aln.seqs[i].data[k])) seq1[len1++]=toupper(mult_aln.seqs[i].data[k])-'A';
		make_p_ptrs(zza,zzc,seq1,len1);
			len2=1;
			for(k=0;k<mult_aln.seqs[j].len;k++) 
				if(isalpha(mult_aln.seqs[j].data[k])) seq2[len2++]=toupper(mult_aln.seqs[j].data[k])-'A';
			make_p_ptrs(zzb,zzd,seq2,len2);
			pair_align(seq1,len1,len2);
			if(!maxsf)
				calc_score=0.0;
			else {
				calc_score=accum[0][maxsf];
				dsr=(len1<len2) ? len1 : len2;
				calc_score = (calc_score/(float)dsr) * 100.0;

			}
/*
			if(calc_score>0.1) 
				printf("Sequences (%d:%d) Aligned. Score: %f\n",
               			i,j,calc_score);
			else
				printf("Sequences (%d:%d) Not Aligned\n",
						i,j);
*/
	for (k=0;k<5;k++)
	   accum[k]=ckfree((void *)accum[k]);
	accum=ckfree((void *)accum);

	displ=ckfree((void *)displ);
	slopes=ckfree((void *)slopes);
	diag_index=ckfree((void *)diag_index);

	zza=ckfree((void *)zza);
	zzb=ckfree((void *)zzb);
	zzc=ckfree((void *)zzc);
	zzd=ckfree((void *)zzd);

	seq1=ckfree((void *)seq1);
	seq2=ckfree((void *)seq2);

	return calc_score;
}

