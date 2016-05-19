#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "rascal.h"

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

static void calc_normd(ALN mult_aln,sint *winscore,sint score_win,float conserved);
static void calc_normd_for_subgroup(ALN mult_aln,double **tmat,GROUP group,sint *winscore,sint score_win,float conserved);

int read_secator_groups(ALN mult_aln,double **tmat,char *filename,GROUP *groups,sint *secgroup,sint *orggroup,float icutoff)
{
	FILE *fin;
	char line[MAXLINE+1];
	sint i,j,k,g,l,s,ngroups;
	sint nseqs,n;
	sint len;
	sint *orphans,norphans;
	sint *tmpgrp;
	sint maxgroup;
	float meanid,maxid;
	float cutoff;
	char sname[MAXLINE+1];
	Boolean found,got_orphans;

	nseqs=mult_aln.nseqs;
        if((fin=fopen(filename,"r"))==NULL) {
                error("Could not open cluster file %s",filename);
                return 0;    
        }


	if(fgets(line,MAXLINE+1,fin)==NULL) return NULL;
	for(i=0;i<strlen(line);i++) 
		if(line[i]==':') break;
	sscanf(&line[i+1],"%d",&g);
	ngroups=g+1;
	if(ngroups<=0 || ngroups>nseqs+1) {
		error("Problem in cluster file %s", filename);
		return 0;
	}

	g=(-1);
	got_orphans=FALSE;
	while(fgets(line,MAXLINE+1,fin)) {
		if(keyword(line,"Cluster")) {
			g++;
			for(i=0;i<strlen(line);i++) 
                		if(line[i]=='=') break;
			sscanf(&line[i+1],"%d",&l);
			if(l<=0) {
				error("No sequences in cluster %d",g);
				return 0;
			}
			groups[g].len=l;
/* in case we want to add some orphans to the group, allocate more memory */
        		groups[g].seqs=(sint *)ckalloc((nseqs+1) * sizeof(sint));
			s=0;
		}
		else if(keyword(line,"unclustered")) {
			g++;
			for(i=0;i<strlen(line);i++) 
                		if(line[i]=='=') break;
			sscanf(&line[i+1],"%d",&l);
			if(l<0) {
				error("No sequences in cluster %d",g);
				return 0;
			}
			groups[g].len=l;
        		groups[g].seqs=(sint *)ckalloc((nseqs+1) * sizeof(sint));
			s=0;
			got_orphans=TRUE;
		}
		else {
			if(!isspace(line[0])) {
				sscanf(line,"%s",sname);
				found=FALSE;
				for(i=0;i<nseqs;i++) {
					if(strcmp(sname,mult_aln.seqs[i].name)==0) {
						groups[g].seqs[s++]=i;
						found=TRUE;
						break;
					}
				}
				/*if(found==FALSE) {
					error("sequence %s in cluster file is not in alignment",sname);
					return 0;
				}*/
			}
		}
	}

	if(got_orphans==FALSE) {
		groups[ngroups-1].len=0;
		groups[ngroups-1].seqs=(sint *)ckalloc((nseqs+1) * sizeof(sint));
	}

	for(i=0;i<ngroups;i++)
		for(j=0;j<groups[i].len;j++)
			orggroup[groups[i].seqs[j]]=i+1;

	if(ngroups<=0) return 0;

/* now we've got the groups defined by secator - check that all the groups are consistent */

/* for the clustered sequences.... */
/* check for orphan sequences in the group - remove any sequences with less than 30% mean identity
   with the other sequences in the group */
        for (i=0;i<ngroups-1;i++) {
		if(groups[i].len==1) {
			groups[ngroups-1].seqs[groups[ngroups-1].len]=groups[i].seqs[0];
			groups[ngroups-1].len++;
			groups[i].seqs[0]=0;
			groups[i].len=0;
		}
		orphans=(sint *)ckalloc((groups[i].len+1)*sizeof(sint));
		for(j=0;j<groups[i].len;j++) {
			meanid=0;
			n=0;
			for(k=0;k<groups[i].len;k++) {
				if(k!=j) {
					meanid+=tmat[groups[i].seqs[j]][groups[i].seqs[k]];
					n++;
				}
			}
			if(n>0) meanid=1.0-meanid/(float)n;
/*			if(icutoff==0.0) {
				len=MAX(mult_aln.seqs[groups[i].seqs[j]].reslen,mult_aln.seqs[groups[i].seqs[k]].reslen);
				cutoff=((float)len-60.0)/(float)len;
        			if(cutoff<0.6) cutoff=0.6;
			}
			else cutoff=icutoff;
			if(meanid<cutoff) {*/
			if(meanid<0.3) {
				orphans[j]=TRUE;
			}
			else orphans[j]=FALSE;
		}

		tmpgrp=(sint *)ckalloc((groups[i].len+1)*sizeof(sint));
		for(j=0;j<groups[i].len;j++) {
			tmpgrp[j]=groups[i].seqs[j];
			groups[i].seqs[j]=0;
		}
		len=groups[i].len;
		groups[i].len=0;
		for(j=0;j<len;j++) {
			if(orphans[j]==TRUE) {
				groups[ngroups-1].seqs[groups[ngroups-1].len]=tmpgrp[j];
				groups[ngroups-1].len++;
			}
			else {
				groups[i].seqs[groups[i].len]=tmpgrp[j];
				groups[i].len++;
				secgroup[tmpgrp[j]]=i+1;
                		mult_aln.seqs[tmpgrp[j]].simgroup=i+1;
			}
		}
		ckfree(orphans);
		ckfree(tmpgrp);
/* check the group again, once the orphans are removed */
		if(groups[i].len==1) {
			groups[ngroups-1].seqs[groups[ngroups-1].len]=groups[i].seqs[0];
			groups[ngroups-1].len++;
			groups[i].seqs[0]=0;
			groups[i].len=0;
		}
		else {
			meanid=0;
			n=0;
			for(j=0;j<groups[i].len;j++) {
				for(k=j+1;k<groups[i].len;k++) {
					meanid+=tmat[groups[i].seqs[j]][groups[i].seqs[k]];
					n++;
				}
			}
			if(n>0) meanid=1.0-meanid/(float)n;
		}
	}

/* for the unclustered sequences.... */
	if(groups[ngroups-1].len>0) {
		ngroups--;
		orphans=(sint *)ckalloc((groups[ngroups].len+1)*sizeof(sint));
		for(j=0;j<groups[ngroups].len;j++) {
			secgroup[groups[ngroups].seqs[j]]=0;
			orphans[j]=groups[ngroups].seqs[j];
			norphans=groups[ngroups].len;
		}
		ckfree(groups[ngroups].seqs);
		groups[ngroups].len=0;

/* check whether the unclustered sequences can now fit into one of the groups */
        	for (i=0;i<norphans;i++) {
			if(secgroup[orphans[i]]==0) {
				maxid=0;
				maxgroup=(-1);
        			for (j=0;j<ngroups;j++) {
					meanid=0.0;
					n=0;
					if(groups[j].len>1) {
						for(k=0;k<groups[j].len;k++) {
							meanid+=tmat[orphans[i]][groups[j].seqs[k]];
							n++;
						}
					}
					if(n>0) meanid=1.0-meanid/(float)n;
					if(meanid>maxid) {
						maxid=meanid;
						maxgroup=j;
					}
				}
				if(icutoff==0.0 && maxgroup>=0) {
					len=MAX(mult_aln.seqs[orphans[i]].reslen,mult_aln.seqs[groups[maxgroup].seqs[0]].reslen);
					cutoff=((float)len-60.0)/(float)len;
        				if(cutoff<0.6) cutoff=0.6;
				}
				else cutoff=icutoff;
				if(maxid>cutoff) {
					n=groups[maxgroup].len;
					groups[maxgroup].seqs[n]=orphans[i];
					groups[maxgroup].len=n+1;
					secgroup[orphans[i]]=maxgroup+1;
fprintf(stdout,"adding %s %d\n",mult_aln.seqs[orphans[i]].name,secgroup[orphans[i]]);
                			mult_aln.seqs[orphans[i]].simgroup=maxgroup+1;
				}
				else {
/* check whether 2 orphans can make a new group */
        				for (j=0;j<norphans;j++) {
						if(icutoff==0.0) {
							len=MAX(mult_aln.seqs[orphans[i]].reslen,mult_aln.seqs[orphans[j]].reslen);
							cutoff=((float)len-60.0)/(float)len;
        						if(cutoff<0.6) cutoff=0.6;
						}
						else cutoff=icutoff;
						if(j!=i && secgroup[orphans[i]]==0 && tmat[orphans[i]][orphans[j]]<(1.0-cutoff)) {
        						groups[ngroups].seqs=(sint *)ckalloc((mult_aln.nseqs+1) * sizeof(sint));
							groups[ngroups].seqs[0]=orphans[i];
							groups[ngroups].seqs[1]=orphans[j];
							groups[ngroups].len=2;
							secgroup[orphans[i]]=ngroups+1;
                					mult_aln.seqs[orphans[i]].simgroup=ngroups+1;
							secgroup[orphans[j]]=ngroups+1;
                					mult_aln.seqs[orphans[j]].simgroup=ngroups+1;
fprintf(stdout,"grouping %s %s %d %d\n",mult_aln.seqs[orphans[i]].name,mult_aln.seqs[orphans[j]].name,secgroup[orphans[i]],secgroup[orphans[j]]);
							ngroups++;
							break;
						}
					}
				}
			}
		}
		for (i=0;i<norphans;i++) {
			if(secgroup[orphans[i]]==0) {
				groups[ngroups].seqs=(sint *)ckalloc((mult_aln.nseqs+1) * sizeof(sint));
				groups[ngroups].seqs[0]=orphans[i];
				groups[ngroups].len=1;
				secgroup[orphans[i]]=ngroups+1;
				mult_aln.seqs[orphans[i]].simgroup=ngroups+1;
				ngroups++;
			}
		}

	}

	return ngroups;
}

int get_blocks(ALN mult_aln,BLOCK *blocks,sint window,sint block_cutoff,float conserved,sint minlength)
{
        int i,nblocks;
        sint *colscore;
        Boolean inregion;
        int first,last;

        colscore = (sint *) ckalloc( (mult_aln.seqs[0].len+1) * sizeof (sint) );
        calc_normd(mult_aln,colscore,window,conserved);

        inregion=FALSE;
        nblocks=0;
        for(i=0;i<mult_aln.seqs[0].len;i++) {
                if(colscore[i]>=block_cutoff) {
                        if(inregion==FALSE) {
                                first=i;
                                inregion=TRUE;
                        }
                }
                else {
                        if(inregion==TRUE) {
                                last=i-1;
                                if(last>=first+minlength) {
                                        blocks[nblocks].first=first;
                                        blocks[nblocks].last=last;
                                        nblocks++;
                                }
                                inregion=FALSE;
                        }
                }
        }
        if(inregion==TRUE) {
                last=mult_aln.seqs[0].len;
                if(last>=first+minlength) {
                        blocks[nblocks].first=first;
                        blocks[nblocks].last=last;
                        nblocks++;
                }
        }

        ckfree(colscore);
        return nblocks;

}


static void calc_normd(ALN mult_aln,sint *winscore,sint score_win,float conserved)
{
	char c;
	sint n,i,j,k,s,s1,p,r,r1;
	sint is,ie,len,nseqs;
	sint *colscore,last_score,next_score;
	sint *ns;
	sint half_win;
	float mean,meanid;
	float t,q1,q3,ul;
	float seqdist,diff;
	float **resdist;
	float *ntot,tmp;
	double dscore;
	sint maxres;
	sint *seqvector;
	sint *seqvector1;
        COMP_MATRIX matrix;
	sint *matptr,*gapptr=NULL;
	extern sint   swgapdnamt[],gon250mt[];

	matrix.format=0;

	if(mult_aln.dnaflag)
		matptr=swgapdnamt;
	else 
		matptr=gon250mt;
	maxres = get_cl_matrix(mult_aln.dnaflag,matptr,gapptr,FALSE,100,&matrix);
	if (maxres == 0)
	{
		error("matrix not found for column score");
		return;
	}

        colscore = (int *) ckalloc( (mult_aln.seqs[0].len+1) * sizeof (int) );
	meanid=n=0;
	nseqs=0;
        for (i=0;i<mult_aln.nseqs;i++) {
			nseqs++;
               		for (j=i+1;j<mult_aln.nseqs;j++) {
                       		dscore = countid(mult_aln.seqs[i],mult_aln.seqs[j]);
				meanid+=dscore;
				n++;
               		}
        }
	meanid/=100.0*(float)n;

	ntot = (float *) ckalloc( (mult_aln.seqs[0].len+2) * sizeof(float) );
	ns = (sint *) ckalloc( (mult_aln.seqs[0].len+2) * sizeof(sint) );
        for(p=0; p<mult_aln.seqs[0].len; p++) {
                ntot[p]=0;
        }
/* make ntot the number of sequences at this position, excluding fragments */
        for(s=0; s<mult_aln.nseqs;s++) {
                len=mult_aln.seqs[0].len;
                is=0;
                ie=len;
                        for(k=0;k<len;k++) {
                                c = mult_aln.seqs[s].data[k];
                                if (isalpha(c)) {
                                        is=k;
                                        break;
                                }
                        }
                        for(k=len-1;k>=0;k--) {
                                c = mult_aln.seqs[s].data[k];
                                if (isalpha(c)) {
                                        ie=k;
                                        break;
                                }
                        }
                for(p=is;p<=ie;p++) ntot[p]++;
        }



	seqvector = (sint *) ckalloc( (NUMRES+2) * sizeof(sint) );
	seqvector1 = (sint *) ckalloc( (NUMRES+2) * sizeof(sint) );

/* calculate the distance between each pair of residues in the space */
	resdist=(float **)ckalloc((NUMRES+1)*sizeof(float *));
	for(i=0;i<NUMRES;i++)
		resdist[i]=(float *)ckalloc((NUMRES+1)*sizeof(float));
	for(i=0;i<NUMRES;i++) {
		for (r=0;r<NUMRES; r++)
			seqvector[r]=matrix.score[r][i];
		resdist[i][i]=0.0;
		for(j=i+1;j<NUMRES;j++) {
			for (r=0;r<NUMRES; r++)
				seqvector1[r]=matrix.score[r][j];
			resdist[i][j]=0.0;
			for(r=0;r<NUMRES;r++) {
				diff=seqvector1[r]-seqvector[r];
				diff/=1000.0;
				resdist[i][j]+=diff*diff;
			}
			resdist[i][j]=sqrt((double)resdist[i][j]);
			resdist[j][i]=resdist[i][j];
		}
	}
	ckfree(seqvector);
	ckfree(seqvector1);
    	for(p=0; p<mult_aln.seqs[0].len; p++)
	{
		if(ntot[p]<2) continue;
		mean=0.0;
    		for(s=0,n=0; s<mult_aln.nseqs; s++)
		{
			if (p<mult_aln.seqs[s].len && isalpha(mult_aln.seqs[s].data[p])) {
    				for(s1=s+1; s1<mult_aln.nseqs; s1++) {
					if (p<mult_aln.seqs[s1].len && isalpha(mult_aln.seqs[s1].data[p])) {
						seqdist=resdist[mult_aln.seqs[s].data[p]-'a'][mult_aln.seqs[s1].data[p]-'a'];
						mean+=seqdist;
						n++;
					}
				}
			}
		}
		if(n>0) mean/=(float)n;

/* calculate mean of seq distances */
		ns[p]=0;
    		for(s=0; s<mult_aln.nseqs; s++)
			if(p<mult_aln.seqs[s].len && isalpha(mult_aln.seqs[s].data[p]))
				ns[p]++;

		colscore[p]=exp((double)(-mean/2.0))*100.0;
		colscore[p]=normalise_score(colscore[p],ns[p],ntot[p],nseqs);
		/*colscore[p]=normalise_score(colscore[p],ns,nseqs,nseqs);*/


	}
	for(i=0;i<NUMRES;i++)
		ckfree(resdist[i]);
	ckfree(resdist);

	half_win = score_win / 2;
	if (half_win == 0) half_win = 1;

        for (i=0;i<half_win+1 && i<mult_aln.seqs[0].len ;i++) {
                for(j=0;j<half_win;j++)
                	winscore[i]+=(float)colscore[i+j];
		winscore[i]=(float)winscore[i]/(float)score_win;
        	if(winscore[i]>100) winscore[i]=100;
        	if(ns[i]<ntot[i]/2.0) winscore[i]=0;
	}

        for (i=0; i+score_win<mult_aln.seqs[0].len; i++) {
                for(j=0;j<score_win;j++)
                        winscore[i+half_win]+=colscore[i+j];
                winscore[i+half_win]=(float)winscore[i+half_win]/(float)score_win;
                if(winscore[i+half_win]>100) winscore[i+half_win]=100;
                if(ns[i+half_win]<ntot[i+half_win]/2.0) winscore[i+half_win]=0;
        }

        for (i=mult_aln.seqs[0].len-score_win;i+half_win<mult_aln.seqs[0].len; i++) {
                for(j=0;j<half_win;j++)
                        winscore[i+half_win]+=colscore[i+j];
                winscore[i+half_win]=(float)winscore[i+half_win]/(float)score_win;
                if(winscore[i+half_win]>100) winscore[i+half_win]=100;
                if(ns[i+half_win]<ntot[i+half_win]/2.0) winscore[i+half_win]=0;
        }

	if(conserved>1.0) {
		for(i=0;i<mult_aln.seqs[0].len; i++) {
			if(ns[i]<nseqs*(conserved-1.0)) winscore[i]=0;
		}
	}
	else if(conserved>0) {
		for(i=0;i<mult_aln.seqs[0].len; i++) {
			if(ns[i]<ntot[i]*conserved) winscore[i]=0;
		}
	}

	ckfree(ns);
	ckfree(ntot);
	ckfree(colscore);
}

int get_blocks_for_subgroup(ALN mult_aln,double **tmat,GROUP group,BLOCK *blocks,sint window,float conserved)
{
	int i,nblocks,nseqs;
	sint *colscore;
	Boolean inregion;
	int first,last;
	int block_cutoff;

	nseqs=group.len;
	if(nseqs<10) block_cutoff=15;
	else block_cutoff=10;

        colscore = (sint *) ckalloc( (mult_aln.seqs[group.seqs[0]].len+1) * sizeof (sint) );
        calc_normd_for_subgroup(mult_aln,tmat,group,colscore,window,conserved);

        inregion=FALSE;
        nblocks=0;
        for(i=0;i<mult_aln.seqs[group.seqs[0]].len;i++) {
                if(colscore[i]>=block_cutoff) {
                        if(inregion==FALSE) {
                                first=i;
                                inregion=TRUE;
                        }
                }
                else {
                        if(inregion==TRUE) {
                                last=i-1;
				if(last>=first+6) {
					blocks[nblocks].first=first;
					blocks[nblocks].last=last;
					nblocks++;
				}
                                inregion=FALSE;
                        }
                }
        }
        if(inregion==TRUE) {
                last=mult_aln.seqs[group.seqs[0]].len;
		if(last>=first+6) {
			blocks[nblocks].first=first;
			blocks[nblocks].last=last;
			nblocks++;
		}
        }

for(i=0;i<nblocks;i++) fprintf(stdout,"%d %d\n",blocks[i].first,blocks[i].last);
	ckfree(colscore);

	return nblocks;

}

static void calc_normd_for_subgroup(ALN mult_aln,double **tmat,GROUP group,sint *winscore,sint score_win,float conserved)
{
	char c;
	sint n,i,j,k,s,s1,p,r,r1;
	sint nseqs;
	sint is,ie,len,iseq,jseq;
	sint *colscore,last_score,next_score;
	sint *ns;
	sint half_win;
	float mean,meanid;
	float t,q1,q3,ul;
	float *seqdist,diff;
	float *ntot,tmp;
	double dscore;
	sint maxres;
	sint *seqvector;
	sint *seqvector1;
        COMP_MATRIX matrix;
	Boolean *fragment;
	sint *matptr,*gapptr=NULL;
	extern sint   swgapdnamt[],gon250mt[];

	matrix.format=0;

	if(mult_aln.dnaflag)
		matptr=swgapdnamt;
	else 
		matptr=gon250mt;
	maxres = get_cl_matrix(mult_aln.dnaflag,matptr,gapptr, TRUE,100,&matrix);
	if (maxres == 0)
	{
		error("matrix not found for column score");
		return;
	}

/* check for fragments */
        fragment = (Boolean *) ckalloc( (mult_aln.nseqs+1) * sizeof (Boolean) );
	nseqs=group.len;

        colscore = (int *) ckalloc( (mult_aln.seqs[0].len+1) * sizeof (int) );
	meanid=n=0;
        for(iseq=0;iseq<group.len;iseq++) {
        	i=group.seqs[iseq];
        	for(jseq=iseq+1;jseq<group.len;jseq++) {
        		j=group.seqs[jseq];
                       	dscore = (1.0-tmat[i][j])*100;
			meanid+=dscore;
			n++;
                        if(dscore>60) {
                                tmp=(float)mult_aln.seqs[i].reslen/(float)mult_aln.seqs[j].reslen;
                                if(tmp<0.8) fragment[i]=TRUE;
                                else if(tmp>1.25) fragment[j]=TRUE;
                        }

               	}
        }
	meanid/=100.0*(float)n;
        for(iseq=0;iseq<group.len;iseq++) {
        	i=group.seqs[iseq];
	}


	s=group.seqs[0];
	ntot = (float *) ckalloc( (mult_aln.seqs[s].len+2) * sizeof(float) );
	ns = (sint *) ckalloc( (mult_aln.seqs[s].len+2) * sizeof(sint) );
        for(p=0; p<mult_aln.seqs[s].len; p++) {
                ntot[p]=0;
        }
/* make ntot the number of sequences at this position, excluding fragments */
        for(iseq=0;iseq<group.len;iseq++) {
        	s=group.seqs[iseq];
                len=mult_aln.seqs[0].len;
                is=0;
                ie=len;
                        for(k=0;k<len;k++) {
                                c = mult_aln.seqs[s].data[k];
                                if (isalpha(c)) {
                                        is=k;
                                        break;
                                }
                        }
                        for(k=len-1;k>=0;k--) {
                                c = mult_aln.seqs[s].data[k];
                                if (isalpha(c)) {
                                        ie=k;
                                        break;
                                }
                        }
                for(p=is;p<=ie;p++) ntot[p]++;
        }



	seqvector = (sint *) ckalloc( (NUMRES+2) * sizeof(sint) );
	seqvector1 = (sint *) ckalloc( (NUMRES+2) * sizeof(sint) );
	seqdist=(float *)ckalloc((mult_aln.nseqs+1)*sizeof(float));

    	for(p=0; p<mult_aln.seqs[group.seqs[0]].len; p++)
	{
		if(ntot[p]<1) continue;
		for(iseq=0;iseq<group.len;iseq++) seqdist[group.seqs[iseq]]=0.0;
        	for(iseq=0;iseq<group.len;iseq++) {
        		s=group.seqs[iseq];
			if (p<mult_aln.seqs[s].len && isalpha(mult_aln.seqs[s].data[p])) {
				for (r=0;r<NUMRES; r++)
					seqvector[r]=matrix.score[r][(int)mult_aln.seqs[s].data[p]-'a'];
        			for(jseq=iseq+1;jseq<group.len;jseq++) {
        				s1=group.seqs[jseq];
					if (p<mult_aln.seqs[s1].len && isalpha(mult_aln.seqs[s1].data[p])) {
						for (r=0;r<NUMRES; r++)
							seqvector1[r]=matrix.score[r][(int)mult_aln.seqs[s1].data[p]-'a'];
						seqdist[s]=0.0;
						for(r=0;r<NUMRES;r++) {
							diff=seqvector1[r]-seqvector[r];
							diff/=1000.0;
							seqdist[s]+=diff*diff;
						}
						seqdist[s]=sqrt((double)seqdist[s]);
					}
				}
			}
		}

/* calculate mean of seq distances */
		mean=0.0;
		ns[p]=0;
        	for(iseq=0;iseq<group.len;iseq++) {
        		s=group.seqs[iseq];
			if(p<mult_aln.seqs[s].len && isalpha(mult_aln.seqs[s].data[p]))
			{
				mean+=seqdist[s];
				ns[p]++;
			}
		}
		if(ns[p]>0) mean/=(float)ns[p];

		colscore[p]=exp((double)(-mean/2.0))*100.0;
		if(conserved>0) colscore[p]=normalise_score(colscore[p],ns[p],ntot[p],nseqs);
	}
	ckfree(fragment);
	ckfree(seqvector);
	ckfree(seqvector1);
	ckfree(seqdist);

	half_win = score_win / 2;
	if (half_win == 0) half_win = 1;

        for (i=0;i<half_win+1 && i<mult_aln.seqs[0].len ;i++) {
		for(j=0;j<score_win;j++)
                	winscore[i]+=(float)colscore[i+j];
        	winscore[i]=(float)winscore[i]/(float)score_win;;
        	if(winscore[i]>100) winscore[i]=100;
        	if(ns[i]<1) winscore[i]=0;
	}

        for (i=0; i+score_win<mult_aln.seqs[0].len; i++) {
                for(j=0;j<score_win;j++)
                        winscore[i+half_win]+=colscore[i+j];
                winscore[i+half_win]=(float)winscore[i+half_win]/(float)score_win;
                if(winscore[i+half_win]>100) winscore[i+half_win]=100;
                if(ns[i+half_win]<1) winscore[i+half_win]=0;
        }

        for (i=mult_aln.seqs[0].len-score_win;i+half_win<mult_aln.seqs[0].len; i++) {
                for(j=0;j<half_win;j++)
                        winscore[i+half_win]+=colscore[i+j];
                winscore[i+half_win]=(float)winscore[i+half_win]/(float)score_win;
                if(winscore[i+half_win]>100) winscore[i+half_win]=100;
                if(ns[i+half_win]<1) winscore[i+half_win]=0;
        }
	if(conserved>0) {
	for(i=0;i<mult_aln.seqs[0].len; i++) {
		if(ns[i]==0) continue;
		if(ns[i]<ntot[i]*conserved || ntot[i]<nseqs*conserved) winscore[i]=0;
	}
	}

	ckfree(ns);
	ckfree(ntot);
	ckfree(colscore);
}


/* 
   Routine for producing unrooted NJ trees from seperately aligned
   pairwise distances.  This produces the GUIDE DENDROGRAMS in
   PHYLIP format.
	- ofile		the file to write the tree to
	- seqs		the sequences
	- nseqs		the number of sequences
	- tmat		pairwise distance matrix
*/
void guide_tree(FILE *ofile,SEQ *seqs,sint nseqs,double **tmat,sint algo)
{
        sint i;
	float dist;

        if(nseqs==2) {
                dist=tmat[0][1]/2.0;
                fprintf(ofile,"(%s:%0.5f,%s:%0.5f);\n",
                        seqs[0].name,dist,seqs[1].name,dist);
        }
        else {
		quicknj(ofile,seqs,nseqs,tmat);
	}
        fclose(ofile);

}

void calc_blockprf1(ALN mult_aln,sint *seqweight,sint firstcol,sint lastcol,GROUP group,COMP_MATRIX matrix,sint *is,sint *ie,PROF *prf,float meanpcid)
{
	sint l,s,i;
	sint res,d;
	float f;
	float *weight,sum;
	Boolean *fragment;
	float **freq;

/* normalise the sequence weights for this group to sum to 100, so that
   the profiles for each group are directly comparable. */

	if(group.len==0) return;

        weight=(float *)ckalloc((group.len+1)*sizeof(float));
        fragment=(Boolean *)ckalloc((group.len+1)*sizeof(Boolean));
	sum=0;

	for(s=0;s<group.len;s++) {
		if(is[group.seqs[s]]>lastcol || ie[group.seqs[s]]<firstcol)
                	fragment[s]=TRUE;
	}

	for(s=0;s<group.len;s++) {
		if(!fragment[s]) sum+=seqweight[group.seqs[s]];
	}

	if(sum<=0) return;
	for(s=0;s<group.len;s++) 
                weight[s]=100*((float)seqweight[group.seqs[s]]/(float)sum);

        freq = (float **) ckalloc( (lastcol-firstcol+2) * sizeof (float *) );
        for(i=0; i<lastcol-firstcol+2; i++)
                freq[i] = (float *) ckalloc( (LENCOL+2) * sizeof(float) );

	for(l=firstcol;l<=lastcol;l++) {
		for(s=0;s<group.len;s++) {
			if(!fragment[s] && isalpha(mult_aln.seqs[group.seqs[s]].data[l])) {
				i=mult_aln.seqs[group.seqs[s]].data[l]-'a';
				freq[l-firstcol][i]+=weight[s];
			}
		}
	}

	for(l=firstcol;l<=lastcol;l++) {
		for (res=0; res<NUMRES; res++) {
			f=0;
			for (d=0; d<NUMRES; d++)
				f += (freq[l-firstcol][d] * matrix.score[d][res]);
			prf->data[l-firstcol][res]=(sint)(f/meanpcid);
		}
	}
	prf->len=lastcol-firstcol+1;
	prf->nseqs=group.len;

	ckfree(fragment);
	ckfree(weight);
        for(i=0; i<lastcol-firstcol+2; i++)
		ckfree(freq[i]);
	ckfree(freq);

}

void calc_blockprf2(ALN mult_aln,sint *seqweight,sint firstcol,sint lastcol,GROUP group,sint *is,sint *ie,PROF *prf)
{
        sint l,s,i;
        sint d;
        float *weight,sum,sum1;
        Boolean *fragment;

/* normalise the sequence weights for this group to sum to 100, so that
   the profiles for each group are directly comparable. */

        if(group.len==0) return;

        weight=(float *)ckalloc((group.len+1)*sizeof(float));
        fragment=(Boolean *)ckalloc((group.len+1)*sizeof(Boolean));
        sum=0;

        for(s=0;s<group.len;s++) {
                if(is[group.seqs[s]]>firstcol || ie[group.seqs[s]]<lastcol)
                        fragment[s]=TRUE;
        }

        for(s=0;s<group.len;s++) {
                if(!fragment[s]) {
                        sum+=seqweight[group.seqs[s]];
		}
        }

        if(sum<=0) return;
        for(s=0;s<group.len;s++)
                if(!fragment[s]) {
                weight[s]=100.0*((float)seqweight[group.seqs[s]]/(float)sum);
		}

        for(l=firstcol;l<=lastcol;l++) {
                for(d=0;d<NUMRES;d++) {
                        sum1=0;
                        for(s=0;s<group.len;s++) {
                                if(!fragment[s] && isalpha(mult_aln.seqs[group.seqs[s]].data[l])) {
                                        i=mult_aln.seqs[group.seqs[s]].data[l]-'a';
                                        if (d == i)  {
                                                sum1 += weight[s];
                                        }
                                }
                        }
                        prf->data[l-firstcol][d]=(10*sum1);
                }
        }

        prf->len=lastcol-firstcol+1;
	prf->nseqs=group.len;
        ckfree(fragment);
        ckfree(weight);

}

int get_groups(char *filename,ALNPTR mult_aln,sint *secgroup,sint *orggroup)
{
        FILE *fin;
        char line[MAXLINE+1];
        sint i,og,g,ngroups;
        sint nseqs;
        char tmp[MAXLINE+1];
        char sname[MAXLINE+1];
        Boolean found;

        nseqs=mult_aln->nseqs;
        if((fin=fopen(filename,"r"))==NULL) {
                error("Could not open file %s",filename);
                return 0;
        }

        ngroups=0;
        while(fgets(line,MAXLINE+1,fin)) {
                if(keyword(line,"SECGROUP")) {
                        sscanf(line,"%s %s %d %d\n",tmp,sname,&g,&og);
                        found=FALSE;
                        for(i=0;i<nseqs;i++) {
                                if(strcmp(sname,mult_aln->seqs[i].name)==0) {
                                        secgroup[i]=g;
                                        orggroup[i]=og;
                                        found=TRUE;
                                        break;
                                }
                        }
                        if(g>ngroups) ngroups=g;
                        if(found==FALSE) {
                                error("sequence %s is not in alignment",sname);
                                return 0;
                        }
                }
        }

        fclose(fin);
        return ngroups;

}

sint score_block_vs_block(sint s1,sint s2,PROF prf1,PROF prf2)
{
        sint i,s;
        sint offset1,offset2,len;
        sint score;

        if(s1>s2) {
                offset1=0;
                offset2=s1-s2;
        }
        else if(s2>s1) {
                offset2=0;
                offset1=s2-s1;
        }
        else {
                offset1=0;
                offset2=0;
        }
        len=MIN((prf1.len-offset1),(prf2.len-offset2));

        score=0;
        for(i=0;i<len;i++) {
                s=prfscore(prf1,prf2,i+offset1,i+offset2);
                score+=s;
        }
	if(len>0)
        	score=(float)score/(float)len;
	else score=0.0;

        return score;
}

sint prfscore(PROF prf1, PROF prf2,sint n,sint m)
{
        sint    ix;
        sint  score;

        score = 0;
        for (ix=0; ix<26; ix++)
                score += (prf1.data[n][ix] * prf2.data[m][ix]);
        return(score/1000);


}

sint score_sequence(char *seq,PROF prf,sint firstcol,sint lastcol)
{
        sint l,i;
        sint score;
	sint ngaps;

        score=ngaps=0;
        for(l=firstcol;l<=lastcol;l++) {
                if(isalpha(seq[l])) {
                        i=seq[l]-'a';
                        score+=prf.data[l-firstcol][i];
                }
        }

/* normalise score by length of block */
	l=lastcol-firstcol;
	if(l>0)
        	score/=l;
	else score=0.0;

        return score;
}


sint score_partial_sequence(char *seq,PROF prf,sint firstcol,sint lastcol,sint *first,sint *last)
{
        sint l,i,j;
	sint f;
        sint *wscore,score;
	sint window=20;
	Boolean in_block=FALSE;

	wscore=(sint *)ckalloc((lastcol-firstcol+2)*sizeof(sint));
        for(l=firstcol+window;l<=lastcol-window;l++) {
		wscore[l-firstcol]=0;
		for(j=l-window;j<l+window;j++) {
                	if(isalpha(seq[l])) {
                        	i=seq[l]-'a';
                        	wscore[l-firstcol]+=prf.data[l-firstcol][i];
                	}
		}
        }
	(*first)=(*last)=0;
        for(i=0;i<=lastcol-firstcol;i++) {
		if(wscore[i]>5000000) {
			if(in_block==FALSE) {
				f=i;
				in_block=TRUE;
			}
		}
		else {
			if(in_block==TRUE) {
				l=i;
				if((l-f) > (*last)-(*first)) {
					(*last)=l;
					(*first)=f;
				}
				in_block=FALSE;
			}
		}
	}
	ckfree(wscore);

	score=0;
	for(l=firstcol+(*first);l<=firstcol+(*last);l++) {
               	if(isalpha(seq[l])) {
                       	i=seq[l]-'a';
                       	score+=prf.data[l-firstcol][i];
               	}
	}

/* normalise score by length of block */
	l=(*last)-(*first);
	if(l>20)
        	score/=l;
	else score=0.0;

        return score;
}
