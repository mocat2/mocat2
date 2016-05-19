/* Phyle of filogenetic tree calculating functions for CLUSTAL W */

/* The nj_tree function has been modified to run in O(N3) time, instead
   of O(N4). The code imitates a NJ function by Bill Bruno in his
   nnneighbor.c code (which in turn includes code copyrighted by
   Joseph Felsenstein!).

	Julie Thompson 4.5.99
*/

/* version 3.52c. (c) Copyright 1993 by Joseph Felsenstein and Mary Kuhner.
   Written by Mary Kuhner, Jon Yamato, Joseph Felsenstein, Akiko Fuseki,
   Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

/* modified by Bill Bruno to replace some negative distances with
   zero and treat adjacent branches sensibly */
/* "non-negative neighbor" */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "score.h"
#include "dayhoff.h"    /* set correction for amino acid distances >= 75% */

#define MAX(a,b) ((a)>(b)?(a):(b))
/*
 *   Prototypes
 */
static Boolean transition(char base1, char base2);
static char * tree_gap_delete(SEQ *seqs,sint nseqs); 
static void compare_tree(PTREE *tree1, PTREE *tree2, sint *hits, sint n);
static sint two_way_split(FILE *ofile,SEQ *seqs,sint nseqs,PTREE *phy_tree, sint start_row, sint flag, sint *boot_totals,sint bootstrap);
static sint two_way_split_nexus(FILE *ofile,SEQ *seqs,sint nseqs,PTREE *phy_tree, sint start_row, sint flag, sint *boot_totals,sint bootstrap);
static void print_tree(FILE *ofile,sint nseqs,PTREE *phy_tree, sint *totals);
static Boolean is_ambiguity(char c);
sint calc_distance_matrix(FILE *ofile,ALN mult_aln,double **tmat,COMP_MATRIX matrix,TREE_OPT tree_opt);
static sint dna_distance_matrix(FILE *ofile,SEQ *seqs,sint nseqs,double **tmat,Boolean tossgaps,Boolean kimura,Boolean use_ambiguities,sint *boot_positions);
static sint prot_distance_matrix(FILE *ofile,SEQ *seqs, sint nseqs,double **tmat, Boolean tossgaps, Boolean kimura,sint *boot_positions);
static sint prot_distance1_matrix(FILE *tree,SEQ *seqs, sint nseqs,double **tmat, COMP_MATRIX matrix, sint cutoff,Boolean tossgaps, Boolean kimura,sint *boot_positions);



static sint debug=0;

/* 
   Calculate pairwise sequence distances using percent identity scores.
   Calculate a NJ tree from the pairwise distances (in the nseqs*nseqs array tmat).
   Write the tree to one of three file formats: clustal,phylip or distance matrix.
   (The 3 formats may be selected simultaneously).

   This is the routine for getting the REAL trees after alignment.
*/
void phylogenetic_tree(ALN mult_aln,char *phylip_name,char *clustal_name,char *dist_name,char *nexus_name,Boolean clustalw_menu,TREE_OPT tree_opt)
{	
	FILE *ofile;
	char path[FILENAMELEN+1];
	sint i,j;
	sint overspill = 0;
	sint total_dists;
	static  PTREE *standard_tree;
	static  PTREE *save_tree;
	char lin2[10];
	sint nseqs;
	sint maxres;
	COMP_MATRIX matrix;
	double **tmat;
	FILE 	*phylip_tree_file;
	FILE 	*clustal_tree_file;
	FILE 	*distances_tree_file;
	FILE 	*nexus_tree_file;

	nseqs=mult_aln.nseqs;
	if(nseqs<=0) {
		error("You must load an alignment first");
		return;
	}

	if(nseqs<2) {
		error("Alignment has only %d sequences",nseqs);
		return;
	}

if(tree_opt.output_tree_clustal) {
        if (clustal_name[0]!=EOS) {
                if((clustal_tree_file = open_explicit_file(
                clustal_name))==NULL) return;
        }
        else {
		if((clustal_tree_file = open_output_file(mult_aln.filename,
		"CLUSTAL tree output","nj",clustal_name)) == NULL) return;
        }
}

if(tree_opt.output_tree_phylip) {
        if (tree_opt.phylip_outname[0]!=EOS) {
                strcpy(phylip_name,tree_opt.phylip_outname);
        }
        if (phylip_name[0]!=EOS) {
                if((phylip_tree_file = open_explicit_file(
                phylip_name))==NULL) return;
        }
        else {
		if((phylip_tree_file = open_output_file(mult_aln.filename,
		"PHYLIP tree output","ph",phylip_name)) == NULL) return;
        }
}

if(tree_opt.output_tree_distances)
{
        if (tree_opt.dist_outname[0]!=EOS) {
                strcpy(dist_name,tree_opt.dist_outname);
        }
        if (dist_name[0]!=EOS) {
                if((distances_tree_file = open_explicit_file(
                dist_name))==NULL) return;
        }
        else {
		if((distances_tree_file = open_output_file(mult_aln.filename,
		"distance matrix output","dst",dist_name)) == NULL) return;
        }
}

if(tree_opt.output_tree_nexus)
{
        if (nexus_name[0]!=EOS) {
                if((nexus_tree_file = open_explicit_file(
                nexus_name))==NULL) return;
        }
        else {
                if((nexus_tree_file = open_output_file(mult_aln.filename,
                "NEXUS tree output","tre",nexus_name)) == NULL) return;
        }
}


        tmat = (double **) ckalloc( (nseqs+1) * sizeof (double *) );
        for(i=0;i<nseqs;i++)
                tmat[i] = (double *)ckalloc( (nseqs+1) * sizeof (double) );

	if(tree_opt.output_tree_clustal) 
		ofile=clustal_tree_file;     /* Turn on file output */
	else
		ofile=NULL;

        if(tree_opt.use_matrix) {
                /*maxres = get_cl_matrix(mult_aln.dnaflag, gon250mt, TRUE, 10, &matrix);*/
                maxres = get_user_matrix(tree_opt.matrixfile,TRUE, 1, &matrix);
                if(maxres<=0) {
                        error("No matrix in file %s\n",tree_opt.matrixfile);
                        return;
                }
        }
	calc_distance_matrix(ofile,mult_aln,tmat,matrix,tree_opt);


/* check if any distances overflowed the distance corrections */
	if ( overspill > 0 ) {
		total_dists = (nseqs*(nseqs-1))/2;
		fprintf(stdout,"\n");
		fprintf(stdout,"\n WARNING: %ld of the distances out of a total of %ld",
		(long)overspill,(long)total_dists);
		fprintf(stdout,"\n were out of range for the distance correction.");
		fprintf(stdout,"\n");
		fprintf(stdout,"\n SUGGESTIONS: 1) remove the most distant sequences");
		fprintf(stdout,"\n           or 2) use the PHYLIP package");
		fprintf(stdout,"\n           or 3) turn off the correction.");
        	fprintf(stdout,"\n Note: Use option 3 with caution! With this degree");
        	fprintf(stdout,"\n of divergence you will have great difficulty");
        	fprintf(stdout,"\n getting robust and reliable trees.");
		fprintf(stdout,"\n\n");
		if (clustalw_menu) 
			getstr("Press [RETURN] to continue",lin2);
	}

	if(tree_opt.output_tree_distances) {
                print_distance_matrix(distances_tree_file,mult_aln.seqs,nseqs,tmat);
	}
	if(tree_opt.output_tree_clustal || tree_opt.output_tree_phylip || tree_opt.output_tree_nexus) {


		standard_tree   = (PTREE *) ckalloc( (nseqs+1) * sizeof (PTREE) );
		for(i=0; i<nseqs; i++) 
			standard_tree[i].description  = (char *) ckalloc( (nseqs+1) * sizeof(char) );
		save_tree   = (PTREE *) ckalloc( (nseqs+1) * sizeof (PTREE) );
		for(i=0; i<nseqs; i++) 
			save_tree[i].description  = (char *) ckalloc( (nseqs+1) * sizeof(char) );

/* for clustal tree option, tree is written to output file during nj_tree(),
   for phylip tree option, nj_tree() calculates the tree which is written to file 
   by print_phylip_tree */ 
		if(tree_opt.output_tree_clustal) 
			ofile = clustal_tree_file;
		else 
			ofile = NULL;
		if(tree_opt.treealgo==NJ)
			nj_tree(ofile,nseqs,standard_tree,tmat);
		else if(tree_opt.treealgo==BIONJ)
			bionj(phylip_tree_file,mult_aln.seqs,nseqs,tmat);
        	for(i=0; i<nseqs; i++)
                	for(j=0; j<nseqs; j++)
                        	save_tree[i].description[j]  = standard_tree[i].description[j];

		if(tree_opt.output_tree_phylip && tree_opt.treealgo==NJ) 
			print_phylip_tree(phylip_tree_file,mult_aln.seqs,nseqs,standard_tree,tmat);
        	for(i=0; i<nseqs; i++)
                	for(j=0; j<nseqs; j++)
                        	standard_tree[i].description[j]  = save_tree[i].description[j];

        	if(tree_opt.output_tree_nexus)
                	print_nexus_tree(nexus_tree_file,mult_aln.seqs,nseqs,standard_tree,tmat);


		for (i=0;i<nseqs;i++)
			ckfree((void *)standard_tree[i].description);
		standard_tree=ckfree((void *)standard_tree);
		for (i=0;i<nseqs;i++)
			ckfree((void *)save_tree[i].description);
		save_tree=ckfree((void *)save_tree);
	}


if(tree_opt.output_tree_clustal) {
	fclose(clustal_tree_file);	
	info("Phylogenetic tree file created:   [%s]",clustal_name);
}

if(tree_opt.output_tree_phylip) {
	fclose(phylip_tree_file);	
	info("Phylogenetic tree file created:   [%s]",phylip_name);
}

if(tree_opt.output_tree_distances) {
	fclose(distances_tree_file);	
	info("Distance matrix  file  created:   [%s]",dist_name);
}

if(tree_opt.output_tree_nexus) {
	fclose(nexus_tree_file);	
	info("Nexus tree file  created:   [%s]",nexus_name);
}

       for(i=0;i<nseqs;i++)
                tmat[i] = ckfree(tmat[i]);
        tmat = ckfree(tmat);


}

/* 
   Calculate pairwise sequence distances using percent identity scores.
*/
sint calc_distance_matrix(FILE *ofile,ALN mult_aln,double **tmat,COMP_MATRIX matrix,TREE_OPT tree_opt)
{
	sint i,j,overspill;
	sint 	*boot_positions;

/* boot_positions is an array that is needed for compatability with the bootstrap_tree()
   function. Just set bootstrap[i]=i to indicate that the distances should be calculated
   over the entire length of the sequences */

	boot_positions = (sint *)ckalloc( (mult_aln.seqs[0].len+2) * sizeof (sint) );

	for(j=0; j<mult_aln.seqs[0].len; j++) 
		boot_positions[j] = j;		

	if(mult_aln.dnaflag)
		overspill = dna_distance_matrix(ofile,mult_aln.seqs,mult_aln.nseqs,tmat,tree_opt.tossgaps,tree_opt.kimura,tree_opt.use_ambiguities,boot_positions);
        else if(tree_opt.use_matrix)
                overspill = prot_distance1_matrix(ofile,mult_aln.seqs,mult_aln.nseqs,tmat,matrix,tree_opt.cutoff,tree_opt.tossgaps,tree_opt.kimura,boot_positions);
        else
                overspill = prot_distance_matrix(ofile,mult_aln.seqs,mult_aln.nseqs,tmat,tree_opt.tossgaps,tree_opt.kimura,boot_positions);


	boot_positions=ckfree((void *)boot_positions);

	return overspill;
}



static Boolean transition(char base1, char base2) /* TRUE if transition; else FALSE */
/* 

   assumes that the bases of DNA sequences have been translated as
   a,A = 0;   c,C = 1;   g,G = 2;   t,T,u,U = 3;  N = 4;  
   a,A = 0;   c,C = 2;   g,G = 6;   t,T,u,U =17;  

   A <--> G  and  T <--> C  are transitions;  all others are transversions.

*/
{
	if( ((base1 == 0) && (base2 == 6)) || ((base1 == 6) && (base2 == 0)) )
		return TRUE;                                     /* A <--> G */
	if( ((base1 ==17) && (base2 == 2)) || ((base1 == 2) && (base2 ==17)) )
		return TRUE;                                     /* T <--> C */
    return FALSE;
}


static char * tree_gap_delete(SEQ *seqs,sint nseqs)   /* flag all positions in alignment that have a gap */
{			  /* in ANY sequence */
	sint seqn;
	sint posn;
	char *tree_gaps;

	tree_gaps = (char *)ckalloc( (seqs[0].len+1) * sizeof (char) );
        
	for(posn=0; posn<seqs[0].len; posn++) {
		tree_gaps[posn] = 0;
     		for(seqn=0; seqn<nseqs; seqn++)  {
			if(!isalpha(seqs[seqn].data[posn])) {
			   tree_gaps[posn] = 1;
				break;
			}
		}
	}
	return tree_gaps;
}

/*
   Write the pairwise sequence distance matrix to a file.
*/
void print_distance_matrix(FILE *ofile,SEQ *seqs,sint nseqs,double **tmat)
{
	sint i,j,k;
	sint max_names;
	
	max_names=0;
	for(i=0;i<nseqs;i++) {
                if(strlen(seqs[i].name)>max_names)
                   max_names=strlen(seqs[i].name);
	}

	fprintf(ofile,"%6d",(pint)nseqs);
	for(i=0;i<nseqs;i++) {
		fprintf(ofile,"\n%-*s ",max_names,seqs[i].name);
		for(j=0;j<nseqs;j++) {
			fprintf(ofile,"%6.3f ",tmat[i][j]);
			if((j+1) % 8 == 0) {
				if(j!=nseqs-1) {
					fprintf(ofile,"\n"); 
					for(k=0;k<max_names;k++) fprintf(ofile," ");
				}
			}
		}
	}
}

/*
    Calculate a neighbour-joining tree from the pairwise distance matrix.
    The tree is stored in structure phy_tree 
*/
void nj_tree(FILE *tree, sint nseqs,PTREE *phy_tree, double **tmat)
{
	register int i;
	sint l[4],nude,k;
	sint nc,mini,minj,j,ii,jj;
	double fnseqs,fnseqs2=0,sumd;
	double diq,djq,dij,d2r,dr,dio,djo,da;
	double tmin,total,dmin;
	double bi,bj,b1,b2,b3,branch[4];
	sint typei,typej;             /* 0 = node; 1 = OTU */
	sint 	*tkill;
	double 	*av;

        /* IMPROVEMENT 1, STEP 0 : declare  variables */
        double *sum_cols, *sum_rows, *join;

        /* IMPROVEMENT 2, STEP 0 : declare  variables */
        sint loop_limit;
        typedef struct _ValidNodeID {
            sint n;
            struct _ValidNodeID *prev;
            struct _ValidNodeID *next;
        } ValidNodeID;
        ValidNodeID *tvalid, *lpi, *lpj, *lpii, *lpjj, *lp_prev, *lp_next;

        /*
         * correspondence of the loop counter variables.
         *   i .. lpi->n,       ii .. lpii->n
         *   j .. lpj->n,       jj .. lpjj->n
         */

	
	fnseqs = (double)nseqs;

/*********************** First initialisation ***************************/
	
	if(tree!=NULL)  {
		fprintf(tree,"\n\n\t\t\tNeighbor-joining Method\n");
		fprintf(tree,"\n Saitou, N. and Nei, M. (1987)");
		fprintf(tree," The Neighbor-joining Method:");
		fprintf(tree,"\n A New Method for Reconstructing Phylogenetic Trees.");
		fprintf(tree,"\n Mol. Biol. Evol., 4(4), 406-425\n");
		fprintf(tree,"\n\n This is an UNROOTED tree\n");
		fprintf(tree,"\n Numbers in parentheses are branch lengths\n\n");
	}	

	if (fnseqs == 2) {
		if (tree!=NULL) fprintf(tree,"Cycle   1     =  SEQ:   1 (%9.5f) joins  SEQ:   2 (%9.5f)",tmat[0][1],tmat[0][1]);
		return;
	}

	mini = minj = 0;

	tkill 		= (sint *) ckalloc( (nseqs+1) * sizeof (sint) );
	av   		= (double *) ckalloc( (nseqs+1) * sizeof (double)   );

        /* IMPROVEMENT 1, STEP 1 : Allocate memory */
        sum_cols        = (double *) ckalloc( (nseqs+1) * sizeof (double)   );
        sum_rows        = (double *) ckalloc( (nseqs+1) * sizeof (double)   );
        join            = (double *) ckalloc( (nseqs+1) * sizeof (double)   );

        /* IMPROVEMENT 2, STEP 1 : Allocate memory */
        tvalid  = (ValidNodeID *) ckalloc( (nseqs+1) * sizeof (ValidNodeID) );
        /* tvalid[0] is special entry in array. it points a header of valid entry list */
        tvalid[0].n = 0;
        tvalid[0].prev = NULL;
        tvalid[0].next = &tvalid[1];

        /* IMPROVEMENT 2, STEP 2 : Construct and initialize the entry chain list */
        for(i=1, loop_limit = nseqs,
                lpi=&tvalid[1], lp_prev=&tvalid[0], lp_next=&tvalid[2] ;
                i<=loop_limit ;
                ++i, ++lpi, ++lp_prev, ++lp_next) {
		tmat[i-1][i-1] = av[i] = 0.0;
		tkill[i] = 0;
                lpi->n = i;
                lpi->prev = lp_prev;
                lpi->next = lp_next;

                /* IMPROVEMENT 1, STEP 2 : Initialize arrays */
                sum_cols[i] = sum_rows[i] = join[i] = 0.0;
	}
        tvalid[loop_limit].next = NULL;

        /*
         * IMPROVEMENT 1, STEP 3 : Calculate the sum of score value that
         * is sequence[i] to others.
         */
        sumd = 0.0;
        for (lpj=tvalid[0].next ; lpj!=NULL ; lpj = lpj->next) {
                double tmp_sum = 0.0;
                j = lpj->n;
                /* calculate sum_rows[j] */
                for (lpi=tvalid[0].next ; (lpi->n) < (lpj->n) ; lpi = lpi->next) {
                        i = lpi->n;
                        tmp_sum += tmat[i-1][j-1];
                        /* tmat[j][i] = tmat[i][j]; */
                }
                sum_rows[j] = tmp_sum;

                tmp_sum = 0.0;
                /* Set lpi to that lpi->n is greater than j */
                if ((lpi != NULL) && (lpi->n == j)) {
                        lpi = lpi->next;
                }
                /* calculate sum_cols[j] */
                for( ; lpi!=NULL ; lpi = lpi->next) {
                        i = lpi->n;
                        tmp_sum += tmat[j-1][i-1];
                        /* tmat[i][j] = tmat[j][i]; */
                }
                sum_cols[j] = tmp_sum;
        }

/*********************** Enter The Main Cycle ***************************/
 	for(nc=1, loop_limit = nseqs-3; nc<=loop_limit; ++nc) {              	/**start main cycle**/
                sumd = 0.0;
                /* IMPROVEMENT 1, STEP 4 : use sum value */
                for(lpj=tvalid[0].next ; lpj!=NULL ; lpj = lpj->next) {
                        sumd += sum_cols[lpj->n];
                }

                /* IMPROVEMENT 3, STEP 0 : multiply tmin and 2*fnseqs2 */
                fnseqs2 = fnseqs - 2.0;         /* Set fnseqs2 at this point. */
                tmin = 99999.0 * 2.0 * fnseqs2;

/*.................compute SMATij values and find the smallest one ........*/

                mini = minj = 0;

                /* jj must starts at least 2 */
                if ((tvalid[0].next != NULL) && (tvalid[0].next->n == 1)) {
                        lpjj = tvalid[0].next->next;
                } else {
                        lpjj = tvalid[0].next;
                }

                for( ; lpjj != NULL; lpjj = lpjj->next) {
                        jj = lpjj->n;
                        for(lpii=tvalid[0].next ; (lpii->n) < (lpjj->n) ; lpii = lpii->next) {
                                ii = lpii->n;
                                diq = djq = 0.0;

                                /* IMPROVEMENT 1, STEP 4 : use sum value */
                                diq = sum_cols[ii] + sum_rows[ii];
                                djq = sum_cols[jj] + sum_rows[jj];
                                /*
                                 * always ii < jj in this point. Use upper
                                 * triangle of score matrix.
                                 */
                                dij = tmat[ii-1][jj-1];

                                /*
                                 * IMPROVEMENT 3, STEP 1 : fnseqs2 is
                                 * already calculated.
                                 */
                                /* fnseqs2 = fnseqs - 2.0 */

                                /* IMPROVEMENT 4 : transform the equation */
  /*-------------------------------------------------------------------*
   * OPTIMIZE of expression 'total = d2r + fnseqs2*dij + dr*2.0'       *
   * total = d2r + fnseq2*dij + 2.0*dr                                 *
   *       = d2r + fnseq2*dij + 2(sumd - dij - d2r)                    *
   *       = d2r + fnseq2*dij + 2*sumd - 2*dij - 2*d2r                 *
   *       =       fnseq2*dij + 2*sumd - 2*dij - 2*d2r + d2r           *
   *       = fnseq2*dij + 2*sumd - 2*dij - d2r                         *
   *       = fnseq2*dij + 2*sumd - 2*dij - (diq + djq - 2*dij)         *
   *       = fnseq2*dij + 2*sumd - 2*dij - diq - djq + 2*dij           *
   *       = fnseq2*dij + 2*sumd - 2*dij + 2*dij - diq - djq           *
   *       = fnseq2*dij + 2*sumd  - diq - djq                          *
   *-------------------------------------------------------------------*/
                                total = fnseqs2*dij + 2.0*sumd  - diq - djq;

                                /*
                                 * IMPROVEMENT 3, STEP 2 : abbrevlate
                                 * the division on comparison between
                                 * total and tmin.
                                 */
                                /* total = total / (2.0*fnseqs2); */

                                if(total < tmin) {
                                        tmin = total;
                                        mini = ii;
                                        minj = jj;
                                }
                        }
                }

                /* MEMO: always ii < jj in avobe loop, so mini < minj */

/*.................compute branch lengths and print the results ........*/


		dio = djo = 0.0;
                /* IMPROVEMENT 1, STEP 4 : use sum value */
                dio = sum_cols[mini] + sum_rows[mini];
                djo = sum_cols[minj] + sum_rows[minj];

                dmin = tmat[mini-1][minj-1];
                dio = (dio - dmin) / fnseqs2;
                djo = (djo - dmin) / fnseqs2;
                bi = (dmin + dio - djo) * 0.5;
                bj = dmin - bi;
                bi = bi - av[mini];
                bj = bj - av[minj];

                if( av[mini] > 0.0 )
                        typei = 0;
                else
                        typei = 1;
                if( av[minj] > 0.0 )
                        typej = 0;
                else
                        typej = 1;

		if(tree!=NULL) 
	 	    fprintf(tree,"\n Cycle%4d     = ",(pint)nc);

/*
   set negative branch lengths to zero.  Also set any tiny positive
   branch lengths to zero.
*/		if( fabs(bi) < 0.0001) bi = 0.0;
		if( fabs(bj) < 0.0001) bj = 0.0;

		if(tree!=NULL) {
		    if(typei == 0) 
			fprintf(tree,"Node:%4d (%9.5f) joins ",(pint)mini,bi);
		    else 
			fprintf(tree," SEQ:%4d (%9.5f) joins ",(pint)mini,bi);

		    if(typej == 0) 
			fprintf(tree,"Node:%4d (%9.5f)",(pint)minj,bj);
		    else 
			fprintf(tree," SEQ:%4d (%9.5f)",(pint)minj,bj);

		    fprintf(tree,"\n");
	    	}	


	    	phy_tree[nc-1].left_branch = bi;
	    	phy_tree[nc-1].right_branch = bj;

		for(i=0; i<nseqs; i++)
			phy_tree[nc-1].description[i] = 0;

	     	if(typei == 0) { 
			for(i=nc-2; i>=0; i--)
				if(phy_tree[i].description[mini-1] == 1) {
					for(j=0; j<nseqs; j++)  
					     if(phy_tree[i].description[j] == 1)
						    phy_tree[nc-1].description[j] = 1;
					break;
				}
		}
		else
			phy_tree[nc-1].description[mini-1] = 1;

		if(typej == 0) {
			for(i=nc-2; i>=0; i--) 
				if(phy_tree[i].description[minj-1] == 1) {
					for(j=0; j<nseqs; j++)  
					     if(phy_tree[i].description[j] == 1)
						    phy_tree[nc-1].description[j] = 1;
					break;
				}
		}
		else
			phy_tree[nc-1].description[minj-1] = 1;
			

/* 
   Here is where the -0.00005 branch lengths come from for 3 or more
   identical seqs.
*/
/*		if(dmin <= 0.0) dmin = 0.0001; */
                if(dmin <= 0.0) dmin = 0.000001;
		av[mini] = dmin * 0.5;

/*........................Re-initialisation................................*/

		fnseqs = fnseqs - 1.0;
		tkill[minj] = 1;

                /* IMPROVEMENT 2, STEP 3 : Remove tvalid[minj] from chain list. */
                /* [ Before ]
                 *  +---------+        +---------+        +---------+
                 *  |prev     |<-------|prev     |<-------|prev     |<---
                 *  |    n    |        | n(=minj)|        |    n    |
                 *  |     next|------->|     next|------->|     next|----
                 *  +---------+        +---------+        +---------+
                 *
                 * [ After ]
                 *  +---------+                           +---------+
                 *  |prev     |<--------------------------|prev     |<---
                 *  |    n    |                           |    n    |
                 *  |     next|-------------------------->|     next|----
                 *  +---------+                           +---------+
                 *                     +---------+
                 *              NULL---|prev     |
                 *                     | n(=minj)|
                 *                     |     next|---NULL
                 *                     +---------+
                 */
                (tvalid[minj].prev)->next = tvalid[minj].next;
                if (tvalid[minj].next != NULL) {
                        (tvalid[minj].next)->prev = tvalid[minj].prev;
                }
                tvalid[minj].prev = tvalid[minj].next = NULL;

                /* IMPROVEMENT 1, STEP 5 : re-calculate sum values. */
                for(lpj=tvalid[0].next ; lpj != NULL ; lpj = lpj->next) {
                        double tmp_di = 0.0;
                        double tmp_dj = 0.0;
                        j = lpj->n;

                        /*
                         * subtrace a score value related with 'minj' from
                         * sum arrays .
                         */
                        if (j < minj) {
                                tmp_dj = tmat[j-1][minj-1];
                                sum_cols[j] -= tmp_dj;
                        } else if (j > minj) {
                                tmp_dj = tmat[minj-1][j-1];
                                sum_rows[j] -= tmp_dj;
                        } /* nothing to do when j is equal to minj. */


                        /*
                         * subtrace a score value related with 'mini' from
                         * sum arrays .
                         */
                        if (j < mini) {
                                tmp_di = tmat[j-1][mini-1];
                                sum_cols[j] -= tmp_di;
                        } else if (j > mini) {
                                tmp_di = tmat[mini-1][j-1];
                                sum_rows[j] -= tmp_di;
                        } /* nothing to do when j is equal to mini. */

                        /*
                         * calculate a score value of the new inner node.
                         * then, store it temporary to join[] array.
                         */
                        join[j] = (tmp_dj + tmp_di) * 0.5;
                }

                /*
                 * 1)
                 * Set the score values (stored in join[]) into the matrix,
                 * row/column position is 'mini'.
                 * 2)
                 * Add a score value of the new inner node to sum arrays.
                 */
                for(lpj=tvalid[0].next ; lpj != NULL; lpj = lpj->next) {
                        j = lpj->n;
                        if (j < mini) {
                                tmat[j-1][mini-1] = join[j];
                                sum_cols[j] += join[j];
                        } else if (j > mini) {
                                tmat[mini-1][j-1] = join[j];
                                sum_rows[j] += join[j];
                        } /* nothing to do when j is equal to mini. */
                }

                /* Re-calculate sum_rows[mini],sum_cols[mini]. */
                sum_cols[mini] = sum_rows[mini] = 0.0;

                /* calculate sum_rows[mini] */
                da = 0.0;
                for(lpj=tvalid[0].next ; lpj->n < mini ; lpj = lpj->next) {
                        da += join[(lpj->n)];
                }
                sum_rows[mini] = da;

                /* skip if 'lpj->n' is equal to 'mini' */
                if ((lpj != NULL) && (lpj->n == mini)) {
                        lpj = lpj->next;
                }

                /* calculate sum_cols[mini] */
                da = 0.0;
                for( ; lpj != NULL; lpj = lpj->next) {
                        da += join[lpj->n];
                }
                sum_cols[mini] = da;

                /*
                 * Clean up sum_rows[minj], sum_cols[minj] and score matrix
                 * related with 'minj'.
                 */
                sum_cols[minj] = sum_rows[minj] = 0.0;

		for(j=1; j<=nseqs; ++j)
			tmat[minj-1][j-1] = tmat[j-1][minj-1] = join[j] = 0.0;


/****/	}						/**end main cycle**/

/******************************Last Cycle (3 Seqs. left)********************/

	nude = 1;

        for(lpi=tvalid[0].next; lpi != NULL; lpi = lpi->next) {
                l[nude] = lpi->n;
                ++nude;
        }


	b1 = (tmat[l[1]-1][l[2]-1] + tmat[l[1]-1][l[3]-1] - tmat[l[2]-1][l[3]-1]) * 0.5;
	b2 =  tmat[l[1]-1][l[2]-1] - b1;
	b3 =  tmat[l[1]-1][l[3]-1] - b1;
 
	branch[0] = b1 - av[l[1]];
	branch[1] = b2 - av[l[2]];
	branch[2] = b3 - av[l[3]];

/* Reset tiny negative and positive branch lengths to zero */
	if( fabs(branch[0]) < 0.0001) branch[0] = 0.0;
	if( fabs(branch[1]) < 0.0001) branch[1] = 0.0;
	if( fabs(branch[2]) < 0.0001) branch[2] = 0.0;

	phy_tree[nseqs-3].left_branch = branch[0];
	phy_tree[nseqs-2].left_branch = branch[1];
	phy_tree[nseqs-1].left_branch = branch[2];

	for(i=0; i<nseqs; i++)
		phy_tree[nseqs-3].description[i] = 0;

	if(tree!=NULL) 
		fprintf(tree,"\n Cycle%4d (Last cycle, trichotomy):\n",(pint)nc);

	for(i=1; i<=3; i++) {
	   if( av[l[i]] > 0.0) {
	      	if(tree!=NULL)
	      	    fprintf(tree,"\n\t\t Node:%4d (%9.5f) ",(pint)l[i],branch[i-1]);
		for(k=nseqs-4; k>=0; k--)
			if(phy_tree[k].description[l[i]-1] == 1) {
				for(j=0; j<nseqs; j++)
				 	if(phy_tree[k].description[j] == 1)
					    phy_tree[nseqs-3].description[j] = i;
				break;
			}
	   }
	   else  {
	      	if(tree!=NULL)
	   	    fprintf(tree,"\n\t\t  SEQ:%4d (%9.5f) ",(pint)l[i],branch[i-1]);
		phy_tree[nseqs-3].description[l[i]-1] = i;
	   }
	   if(i < 2) {
	      	if(tree!=NULL)
	            fprintf(tree,"joins");
	   }
	}

	if(tree!=NULL)
		fprintf(tree,"\n");

	tkill=ckfree((void *)tkill);
	av=ckfree((void *)av);

        /* IMPROVEMENT 1, STEP 6 : release memory area */
        ckfree(sum_cols);
        ckfree(sum_rows);
        ckfree(join);

        /* IMPROVEMENT 2, STEP 4 : release memory area */
        ckfree(tvalid);

}


/*
    Perform bootstrapping of a phylogenetic tree.
    First calculates the usual NJ tree using the full length of the sequences.
    Then, compares this 'standard' tree with a 'sample' tree calculated
    from randomly selected regions of the sequences.
    The tree is written to one of 2 file formats: clustal or phylip.
    (The 2 formats may be selected simultaneously).
*/

void bootstrap_tree(ALN mult_aln,char *phylip_name,char *clustal_name,char *nexus_name,Boolean clustalw_menu,TREE_OPT tree_opt, BSTREE_OPT bstree_opt)
{
	FILE *ofile;
	sint i,j;
	int ranno;
	char path[MAXLINE+1];
    char dummy[10];
	static PTREE *sample_tree;
	static PTREE *standard_tree;
	static  PTREE *save_tree;
	sint total_dists, overspill = 0, total_overspill = 0;
	sint nfails = 0;
	sint nseqs;
	sint maxres;
	COMP_MATRIX matrix;
	double **tmat;
	FILE 	*phylip_tree_file;
	FILE 	*clustal_tree_file;
	FILE 	*nexus_tree_file;
	sint 	*boot_positions;
	sint	*boot_totals;

	nseqs=mult_aln.nseqs;
	if(nseqs<=0) {
		error("You must load an alignment first");
		return;
	}

        if(nseqs<4) {
                error("Alignment has only %d sequences",nseqs);
                return;
        }

	if(!tree_opt.output_tree_clustal && !tree_opt.output_tree_phylip && !tree_opt.output_tree_nexus) {
		error("You must select either clustal or phylip or nexus tree output format");
		return;
	}
	get_path(mult_aln.filename, path);
	
	if (tree_opt.output_tree_clustal) {
        if (clustal_name[0]!=EOS) {
                if((clustal_tree_file = open_explicit_file(
                clustal_name))==NULL) return;
        }
        else {
		if((clustal_tree_file = open_output_file(path,
		"\nEnter name for bootstrap output file  ","njb",
		clustal_name)) == NULL) return;
        }
	}

	if (tree_opt.output_tree_phylip) {
        if (phylip_name[0]!=EOS) {
                if((phylip_tree_file = open_explicit_file(
                phylip_name))==NULL) return;
        }
	else {
		if((phylip_tree_file = open_output_file(path,
		"\nEnter name for bootstrap output file  ","phb",
		phylip_name)) == NULL) return;
	}
	}

	if (tree_opt.output_tree_nexus) {
        if (nexus_name[0]!=EOS) {
                if((nexus_tree_file = open_explicit_file(
                nexus_name))==NULL) return;
        }
	else {
		if((nexus_tree_file = open_output_file(path,
		"\nEnter name for bootstrap output file  ","treb",
		nexus_name)) == NULL) return;
	}
	}
	boot_totals    = (sint *)ckalloc( (nseqs+1) * sizeof (sint) );
	for(i=0;i<nseqs;i++)
		boot_totals[i]=0;
		
        tmat = (double **) ckalloc( (nseqs+1) * sizeof (double *) );
        for(i=0;i<nseqs;i++)
                tmat[i] = (double *)ckalloc( (nseqs+1) * sizeof (double) );
        for(i=0;i<nseqs;i++)
                for(j=0;j<nseqs;j++)
                        tmat[i][j]=0.0;

	if(tree_opt.output_tree_clustal) 
		ofile=clustal_tree_file;     /* Turn on file output */
	else
		ofile=NULL;

/* First, calculate the standard tree */

        if(tree_opt.use_matrix) {
                /*maxres = get_cl_matrix(mult_aln.dnaflag, gon250mt, TRUE, 10, &matrix);*/
                maxres = get_user_matrix(tree_opt.matrixfile,TRUE, 1, &matrix);
                if(maxres<=0) {
                        error("No matrix in file %s\n",tree_opt.matrixfile);
                        return;
                }
        }
	calc_distance_matrix(ofile,mult_aln,tmat,matrix,tree_opt);

/* check if any distances overflowed the distance corrections */
	if ( overspill > 0 ) {
		total_dists = (nseqs*(nseqs-1))/2;
		fprintf(stdout,"\n");
		fprintf(stdout,"\n WARNING: %d of the distances out of a total of %d",
		(pint)overspill,(pint)total_dists);
		fprintf(stdout,"\n were out of range for the distance correction.");
		fprintf(stdout,"\n This may not be fatal but you have been warned!");
		fprintf(stdout,"\n");
		fprintf(stdout,"\n SUGGESTIONS: 1) turn off the correction");
		fprintf(stdout,"\n           or 2) remove the most distant sequences");
		fprintf(stdout,"\n           or 3) use the PHYLIP package.");
		fprintf(stdout,"\n\n");
		if (clustalw_menu) 
			getstr("Press [RETURN] to continue",dummy);
	}


	standard_tree   = (PTREE *) ckalloc( (nseqs+1) * sizeof (PTREE) );
	for(i=0; i<nseqs; i++) 
		standard_tree[i].description  = (char *) ckalloc( (nseqs+1) * sizeof(char) );

/* for clustal tree option, tree is written to output file during nj_tree() */
	if(tree_opt.output_tree_clustal) 
		ofile = clustal_tree_file;
	else 
		ofile = NULL;
	nj_tree(ofile,nseqs,standard_tree,tmat);

	if (tree_opt.output_tree_clustal)
		fprintf(clustal_tree_file,"\n\n\t\t\tBootstrap Confidence Limits\n\n");

/* start the random number generator */

	if(clustalw_menu) 
   		bstree_opt.boot_ran_seed = 
getint("\n\nEnter seed no. for random number generator ",1,1000,bstree_opt.boot_ran_seed);

       	addrandinit((unsigned long) bstree_opt.boot_ran_seed);

	if (tree_opt.output_tree_clustal)
		fprintf(clustal_tree_file,"\n Random number generator seed = %7u\n",
		bstree_opt.boot_ran_seed);

	if(clustalw_menu) 
  		bstree_opt.boot_ntrials = 
getint("\n\nEnter number of bootstrap trials ",1,10000,bstree_opt.boot_ntrials);

	if (tree_opt.output_tree_clustal) {
  		fprintf(clustal_tree_file,"\n Number of bootstrap trials   = %7d\n",
			(pint)bstree_opt.boot_ntrials);

		fprintf(clustal_tree_file,
		"\n\n Diagrammatic representation of the above tree: \n");
		fprintf(clustal_tree_file,"\n Each row represents 1 tree cycle;");
		fprintf(clustal_tree_file," defining 2 groups.\n");
		fprintf(clustal_tree_file,"\n Each column is 1 sequence; ");
		fprintf(clustal_tree_file,"the stars in each line show 1 group; ");
		fprintf(clustal_tree_file,"\n the dots show the other\n");
		fprintf(clustal_tree_file,"\n Numbers show occurences in bootstrap samples.");
	}

	sample_tree   = (PTREE *) ckalloc( (nseqs+1) * sizeof (PTREE) );
	for(i=0; i<nseqs; i++) 
		sample_tree[i].description   = (char *) ckalloc( (nseqs+1) * sizeof(char) );

	if (clustalw_menu)
	fprintf(stdout,"\n\nEach dot represents 10 trials\n\n");
        total_overspill = 0;
	nfails = 0;
	boot_positions = (sint *)ckalloc( (mult_aln.seqs[0].len+2) * sizeof (sint) );
	for(i=0; i<bstree_opt.boot_ntrials; i++) {
		for(j=0; j<mult_aln.seqs[0].len; j++) { /* select alignment */
							    /* positions for */
			ranno = addrand( (unsigned long) mult_aln.seqs[0].len);
			boot_positions[j] = ranno; 	    /* bootstrap sample */
		}
		if(mult_aln.dnaflag)
			overspill = dna_distance_matrix(NULL,mult_aln.seqs,nseqs,tmat,tree_opt.tossgaps,tree_opt.kimura,tree_opt.use_ambiguities,boot_positions);
                else if(tree_opt.use_matrix)
                        overspill = prot_distance1_matrix(NULL,mult_aln.seqs,nseqs,tmat,matrix,tree_opt.cutoff,tree_opt.tossgaps,tree_opt.kimura,boot_positions);
                else
                        overspill = prot_distance_matrix(NULL,mult_aln.seqs,nseqs,tmat,tree_opt.tossgaps,tree_opt.kimura,boot_positions);


		if( overspill > 0) {
			total_overspill = total_overspill + overspill;
			nfails++;
		}			


/* calculate the sample tree without writing to the output file */
		nj_tree(NULL,nseqs,sample_tree,tmat);

		compare_tree(standard_tree, sample_tree, boot_totals, nseqs);
		if (clustalw_menu) {
			if(i % 10  == 0) fprintf(stdout,".");
			if(i % 100 == 0) fprintf(stdout,"\n");
		}
	}

	boot_positions=ckfree((void *)boot_positions);

/* check if any distances overflowed the distance corrections */
	if ( nfails > 0 ) {
		total_dists = (nseqs*(nseqs-1))/2;
		fprintf(stdout,"\n");
		fprintf(stdout,"\n WARNING: %ld of the distances out of a total of %ld times %ld",
		(long)total_overspill,(long)total_dists,(long)bstree_opt.boot_ntrials);
		fprintf(stdout,"\n were out of range for the distance correction.");
		fprintf(stdout,"\n This affected %d out of %d bootstrap trials.",
		(pint)nfails,(pint)bstree_opt.boot_ntrials);
		fprintf(stdout,"\n This may not be fatal but you have been warned!");
		fprintf(stdout,"\n");
		fprintf(stdout,"\n SUGGESTIONS: 1) turn off the correction");
		fprintf(stdout,"\n           or 2) remove the most distant sequences");
		fprintf(stdout,"\n           or 3) use the PHYLIP package.");
		fprintf(stdout,"\n\n");
		if (clustalw_menu) 
			getstr("Press [RETURN] to continue",dummy);
	}



	for (i=0;i<nseqs;i++)
		ckfree((void *)sample_tree[i].description);
	sample_tree=ckfree((void *)sample_tree);

	if (tree_opt.output_tree_clustal)
		print_tree(clustal_tree_file,nseqs,standard_tree,boot_totals);

	save_tree   = (PTREE *) ckalloc( (nseqs+1) * sizeof (PTREE) );
	for(i=0; i<nseqs; i++) 
		save_tree[i].description  = (char *) ckalloc( (nseqs+1) * sizeof(char) );

       	for(i=0; i<nseqs; i++)
               	for(j=0; j<nseqs; j++)
                       	save_tree[i].description[j]  = standard_tree[i].description[j];

	if(tree_opt.output_tree_phylip) {
		print_phylip_bootstrap_tree(phylip_tree_file,mult_aln.seqs,nseqs,standard_tree,tmat,
						boot_totals,
						 bstree_opt.bootstrap_format);
	}
        for(i=0; i<nseqs; i++)
               	for(j=0; j<nseqs; j++)
                       	standard_tree[i].description[j]  = save_tree[i].description[j];

        if(tree_opt.output_tree_nexus) {
		print_nexus_bootstrap_tree(nexus_tree_file,mult_aln.seqs,nseqs,standard_tree,tmat,
						boot_totals,
						 bstree_opt.bootstrap_format);
        }


	boot_totals=ckfree((void *)boot_totals);

	for (i=0;i<nseqs;i++)
		ckfree((void *)standard_tree[i].description);
	standard_tree=ckfree((void *)standard_tree);

        for(i=0;i<nseqs;i++)
                tmat[i] = ckfree(tmat[i]);
        tmat = ckfree(tmat);

	if (tree_opt.output_tree_clustal) {
		fclose(clustal_tree_file);
		info("Bootstrap output file completed       [%s]",clustal_name);
	}

	if (tree_opt.output_tree_phylip) {
		fclose(phylip_tree_file);
		info("Bootstrap output file completed       [%s]",phylip_name);
	}

	if (tree_opt.output_tree_nexus) {
		fclose(nexus_tree_file);
		info("Bootstrap output file completed       [%s]",nexus_name);
	}
}


static void compare_tree(PTREE *tree1, PTREE *tree2, sint *hits, sint n)
{	
	sint i,j,k;
	sint nhits1, nhits2;

	for(i=0; i<n-3; i++)  {
		for(j=0; j<n-3; j++)  {
			nhits1 = 0;
			nhits2 = 0;
			for(k=0; k<n; k++) {
				if(tree1[i].description[k] == tree2[j].description[k]) nhits1++;
				if(tree1[i].description[k] != tree2[j].description[k]) nhits2++;
			}
			if((nhits1 == n) || (nhits2 == n)) hits[i]++;
		}
	}
}

/* Write a tree to a file in nexus format */

void print_nexus_tree(FILE *tree,SEQ *seqs, sint nseqs, PTREE *phy_tree, double **tmat)
{
	sint i;
	sint old_row;
	
        fprintf(tree,"#NEXUS\n\n");

        fprintf(tree,"BEGIN TREES;\n\n");
        fprintf(tree,"\tTRANSLATE\n");
        for(i=0;i<nseqs-1;i++) {
                fprintf(tree,"\t\t%d    %s,\n",(pint)(i+1),seqs[i].name);
        }
                fprintf(tree,"\t\t%d    %s\n",(pint)nseqs,seqs[nseqs-1].name);
        fprintf(tree,"\t\t;\n");

        fprintf(tree,"\tUTREE PAUP_1= ");

	if(nseqs==2) {
		fprintf(tree,"(1:%7.5f,2:%7.5f);",tmat[0][1],tmat[0][1]);
		return;
	}

	fprintf(tree,"(");
 
	old_row=two_way_split_nexus(tree,seqs,nseqs,phy_tree, nseqs-3,1,NULL,0);
	fprintf(tree,":%7.5f",phy_tree[nseqs-3].left_branch);
	fprintf(tree,",");

	old_row=two_way_split_nexus(tree,seqs,nseqs,phy_tree, nseqs-3,2,NULL,0);
	fprintf(tree,":%7.5f",phy_tree[nseqs-2].left_branch);
	fprintf(tree,",");

	old_row=two_way_split_nexus(tree,seqs,nseqs,phy_tree, nseqs-3,3,NULL,0);
	fprintf(tree,":%7.5f",phy_tree[nseqs-1].left_branch);
	fprintf(tree,")");
	fprintf(tree,";");

	fprintf(tree,"\nENDBLOCK;\n");
}

/* Write a tree to a file in nexus format, with bootstrap information */

void print_nexus_bootstrap_tree(FILE *tree,SEQ *seqs, sint nseqs, PTREE *phy_tree, double **tmat, sint *boot_totals,sint bootstrap)
{
	sint i;
	sint old_row;
	
        fprintf(tree,"#NEXUS\n\n");
        fprintf(tree,"BEGIN TREES;\n\n");
        fprintf(tree,"\tTRANSLATE\n");
        for(i=0;i<nseqs-1;i++) {
                fprintf(tree,"\t\t%d    %s;\n",(pint)(i+1),seqs[i].name);
        }
        fprintf(tree,"\t\t%d    %s\n",(pint)nseqs,seqs[nseqs-1].name);
        fprintf(tree,"\t\t;\n");

        fprintf(tree,"\tUTREE PAUP_1= ");
	if(nseqs==2) {
		fprintf(tree,"(1:%7.5f,2:%7.5f);",tmat[0][1],tmat[0][1]);
		return;
	}

	fprintf(tree,"(");
 
	old_row=two_way_split_nexus(tree,seqs,nseqs,phy_tree, nseqs-3,1,boot_totals,bootstrap);
	fprintf(tree,":%7.5f",phy_tree[nseqs-3].left_branch);
	if ((bootstrap==BS_BRANCH_LABELS) && (old_row>=0) && (boot_totals[old_row]>0))
		fprintf(tree,"[%d]",(pint)boot_totals[old_row]);
	fprintf(tree,",");

	old_row=two_way_split_nexus(tree,seqs,nseqs,phy_tree, nseqs-3,2,boot_totals,bootstrap);
	fprintf(tree,":%7.5f",phy_tree[nseqs-2].left_branch);
	if ((bootstrap==BS_BRANCH_LABELS) && (old_row>=0) && (boot_totals[old_row]>0))
		fprintf(tree,"[%d]",(pint)boot_totals[old_row]);
	fprintf(tree,",");

	old_row=two_way_split_nexus(tree,seqs,nseqs,phy_tree, nseqs-3,3,boot_totals,bootstrap);
	fprintf(tree,":%7.5f",phy_tree[nseqs-1].left_branch);
	if ((bootstrap==BS_BRANCH_LABELS) && (old_row>=0) && (boot_totals[old_row]>0))
		fprintf(tree,"[%d]",(pint)boot_totals[old_row]);
	fprintf(tree,")");
        if (bootstrap==BS_NODE_LABELS) fprintf(tree,"TRICHOTOMY");
	fprintf(tree,";");
	fprintf(tree,"\nENDBLOCK;\n");
}


sint two_way_split_nexus
(FILE *ofile,SEQ *seqs, sint nseqs, PTREE *phy_tree,sint start_row, sint flag, sint *boot_totals,sint bootstrap)
{
	sint row, new_row = 0, old_row, col, test_col = 0;
	Boolean single_seq;

	if(start_row != nseqs-3) fprintf(ofile,"("); 

	for(col=0; col<nseqs; col++) {
		if(phy_tree[start_row].description[col] == flag) {
			test_col = col;
			break;
		}
	}

	single_seq = TRUE;
	for(row=start_row-1; row>=0; row--) 
		if(phy_tree[row].description[test_col] == 1) {
			single_seq = FALSE;
			new_row = row;
			break;
		}

	if(single_seq) {
		phy_tree[start_row].description[test_col] = 0;
		fprintf(ofile,"%d",test_col+1);
		if(start_row == nseqs-3) {
			return(-1);
		}

		fprintf(ofile,":%7.5f,",phy_tree[start_row].left_branch);
	}
	else {
		for(col=0; col<nseqs; col++) {
		    if((phy_tree[start_row].description[col]==1)&&
		       (phy_tree[new_row].description[col]==1))
				phy_tree[start_row].description[col] = 0;
		}
		old_row=two_way_split_nexus(ofile,seqs,nseqs,phy_tree, new_row, (sint)1, boot_totals,bootstrap);
		if(start_row == nseqs-3) {
			return(new_row);
		}

		fprintf(ofile,":%7.5f",phy_tree[start_row].left_branch);
		if ((bootstrap==BS_BRANCH_LABELS) && (boot_totals[old_row]>0))
			fprintf(ofile,"[%d]",(pint)boot_totals[old_row]);

		fprintf(ofile,",");
	}


	for(col=0; col<nseqs; col++) 
		if(phy_tree[start_row].description[col] == flag) {
			test_col = col;
			break;
		}
	
	single_seq = TRUE;
	new_row = 0;
	for(row=start_row-1; row>=0; row--) 
		if(phy_tree[row].description[test_col] == 1) {
			single_seq = FALSE;
			new_row = row;
			break;
		}

	if(single_seq) {
		phy_tree[start_row].description[test_col] = 0;
		fprintf(ofile,"%d",test_col+1);
		fprintf(ofile,":%7.5f)",phy_tree[start_row].right_branch);
	}
	else {
		for(col=0; col<nseqs; col++) {
		    if((phy_tree[start_row].description[col]==1)&&
		       (phy_tree[new_row].description[col]==1))
				phy_tree[start_row].description[col] = 0;
		}
		old_row=two_way_split_nexus(ofile,seqs,nseqs,phy_tree, new_row, (sint)1, boot_totals,bootstrap);
		fprintf(ofile,":%7.5f",phy_tree[start_row].right_branch);
		if ((bootstrap==BS_BRANCH_LABELS) && (boot_totals[old_row]>0))
			fprintf(ofile,"[%d]",(pint)boot_totals[old_row]);

		fprintf(ofile,")");
	}
	if ((bootstrap==BS_NODE_LABELS) && (boot_totals[start_row]>0))
			fprintf(ofile,"%d",(pint)boot_totals[start_row]);
	
	return(start_row);
}


/* Write a tree to a file in phylip format */

void print_phylip_tree(FILE *ofile,SEQ *seqs, sint nseqs, PTREE *phy_tree, double **tmat)
{
	sint old_row;
	
	if(nseqs==2) {
		fprintf(ofile,"(%s:%7.5f,%s:%7.5f);",seqs[0].name,tmat[0][1],seqs[1].name,tmat[0][1]);
		return;
	}

	fprintf(ofile,"(\n");
 
	old_row=two_way_split(ofile,seqs,nseqs,phy_tree, nseqs-3,1,NULL,0);
	fprintf(ofile,":%7.5f",phy_tree[nseqs-3].left_branch);
	fprintf(ofile,",\n");

	old_row=two_way_split(ofile,seqs,nseqs,phy_tree, nseqs-3,2,NULL,0);
	fprintf(ofile,":%7.5f",phy_tree[nseqs-2].left_branch);
	fprintf(ofile,",\n");

	old_row=two_way_split(ofile,seqs,nseqs,phy_tree, nseqs-3,3,NULL,0);
	fprintf(ofile,":%7.5f",phy_tree[nseqs-1].left_branch);
	fprintf(ofile,")");
	fprintf(ofile,";\n");
}

/* Write a tree to a file in phylip format, with bootstrap information */

void print_phylip_bootstrap_tree(FILE *ofile,SEQ *seqs, sint nseqs, PTREE *phy_tree, double **tmat, sint *boot_totals,sint bootstrap)
{
	sint old_row;
	
	if(nseqs==2) {
		fprintf(ofile,"(%s:%7.5f,%s:%7.5f);",seqs[0].name,tmat[0][1],seqs[1].name,tmat[0][1]);
		return;
	}

	fprintf(ofile,"(\n");
 
	old_row=two_way_split(ofile,seqs,nseqs,phy_tree, nseqs-3,1,boot_totals,bootstrap);
	fprintf(ofile,":%7.5f",phy_tree[nseqs-3].left_branch);
	if ((bootstrap==BS_BRANCH_LABELS) && (old_row>=0) && (boot_totals[old_row]>0))
		fprintf(ofile,"[%d]",(pint)boot_totals[old_row]);
	fprintf(ofile,",\n");

	old_row=two_way_split(ofile,seqs,nseqs,phy_tree, nseqs-3,2,boot_totals,bootstrap);
	fprintf(ofile,":%7.5f",phy_tree[nseqs-2].left_branch);
	if ((bootstrap==BS_BRANCH_LABELS) && (old_row>=0) && (boot_totals[old_row]>0))
		fprintf(ofile,"[%d]",(pint)boot_totals[old_row]);
	fprintf(ofile,",\n");

	old_row=two_way_split(ofile,seqs,nseqs,phy_tree, nseqs-3,3,boot_totals,bootstrap);
	fprintf(ofile,":%7.5f",phy_tree[nseqs-1].left_branch);
	if ((bootstrap==BS_BRANCH_LABELS) && (old_row>=0) && (boot_totals[old_row]>0))
		fprintf(ofile,"[%d]",(pint)boot_totals[old_row]);
	fprintf(ofile,")");
        if (bootstrap==BS_NODE_LABELS) fprintf(ofile,"TRICHOTOMY");
	fprintf(ofile,";\n");
}


sint two_way_split
(FILE *ofile,SEQ *seqs, sint nseqs, PTREE *phy_tree,sint start_row, sint flag, sint *boot_totals,sint bootstrap)
{
	sint row, new_row = 0, old_row, col, test_col = 0;
	Boolean single_seq;

	if(start_row != nseqs-3) fprintf(ofile,"(\n"); 

	for(col=0; col<nseqs; col++) {
		if(phy_tree[start_row].description[col] == flag) {
			test_col = col;
			break;
		}
	}

	single_seq = TRUE;
	for(row=start_row-1; row>=0; row--) 
		if(phy_tree[row].description[test_col] == 1) {
			single_seq = FALSE;
			new_row = row;
			break;
		}

	if(single_seq) {
		phy_tree[start_row].description[test_col] = 0;
		fprintf(ofile,"%s",seqs[test_col].name);
		if(start_row == nseqs-3) {
			return(-1);
		}

		fprintf(ofile,":%7.5f,\n",phy_tree[start_row].left_branch);
	}
	else {
		for(col=0; col<nseqs; col++) {
		    if((phy_tree[start_row].description[col]==1)&&
		       (phy_tree[new_row].description[col]==1))
				phy_tree[start_row].description[col] = 0;
		}
		old_row=two_way_split(ofile,seqs,nseqs,phy_tree, new_row, (sint)1, boot_totals,bootstrap);
		if(start_row == nseqs-3) {
			return(new_row);
		}

		fprintf(ofile,":%7.5f",phy_tree[start_row].left_branch);
		if ((bootstrap==BS_BRANCH_LABELS) && (boot_totals[old_row]>0))
			fprintf(ofile,"[%d]",(pint)boot_totals[old_row]);

		fprintf(ofile,",\n");
	}


	for(col=0; col<nseqs; col++) 
		if(phy_tree[start_row].description[col] == flag) {
			test_col = col;
			break;
		}
	
	single_seq = TRUE;
	new_row = 0;
	for(row=start_row-1; row>=0; row--) 
		if(phy_tree[row].description[test_col] == 1) {
			single_seq = FALSE;
			new_row = row;
			break;
		}

	if(single_seq) {
		phy_tree[start_row].description[test_col] = 0;
		fprintf(ofile,"%s",seqs[test_col].name);
		fprintf(ofile,":%7.5f)\n",phy_tree[start_row].right_branch);
	}
	else {
		for(col=0; col<nseqs; col++) {
		    if((phy_tree[start_row].description[col]==1)&&
		       (phy_tree[new_row].description[col]==1))
				phy_tree[start_row].description[col] = 0;
		}
		old_row=two_way_split(ofile,seqs,nseqs,phy_tree, new_row, (sint)1, boot_totals,bootstrap);
		fprintf(ofile,":%7.5f",phy_tree[start_row].right_branch);
		if ((bootstrap==BS_BRANCH_LABELS) && (boot_totals[old_row]>0))
			fprintf(ofile,"[%d]",(pint)boot_totals[old_row]);

		fprintf(ofile,")\n");
	}
	if ((bootstrap==BS_NODE_LABELS) && (boot_totals[start_row]>0))
			fprintf(ofile,"%d",(pint)boot_totals[start_row]);
	
	return(start_row);
}



static void print_tree(FILE *ofile,sint nseqs,PTREE *phy_tree,sint *totals)
{
	sint row,col;

	fprintf(ofile,"\n");

	for(row=0; row<nseqs-3; row++)  {
		fprintf(ofile," \n");
		for(col=0; col<nseqs; col++) { 
			if(phy_tree[row].description[col] == 0)
				fprintf(ofile,"*");
			else
				fprintf(ofile,".");
		}
		if(totals[row] > 0)
			fprintf(ofile,"%7d",(pint)totals[row]);
	}
	fprintf(ofile," \n");
	for(col=0; col<nseqs; col++) 
		fprintf(ofile,"%1d",(pint)phy_tree[nseqs-3].description[col]);
	fprintf(ofile,"\n");
}



static sint dna_distance_matrix(FILE *tree,SEQ *seqs,sint nseqs,double **tmat,Boolean tossgaps,Boolean kimura,Boolean use_ambiguities,sint *boot_positions)
{   
	sint m,n;
	sint j,i;
	char res1, res2;
    sint overspill = 0;
	double p,q,e,a,b,k;	
	char 	*tree_gaps;

	tree_gaps=tree_gap_delete(seqs,nseqs); /* flag positions with gaps (tree_gaps[i] = 1 ) */
	
	if(tree!=NULL) {
		fprintf(tree,"\n");
		fprintf(tree,"\n DIST   = percentage divergence (/100)");
		fprintf(tree,"\n p      = rate of transition (A <-> G; C <-> T)");
		fprintf(tree,"\n q      = rate of transversion");
		fprintf(tree,"\n Length = number of sites used in comparison");
		fprintf(tree,"\n");
	    if(tossgaps) {
		fprintf(tree,"\n All sites with gaps (in any sequence) deleted!");
		fprintf(tree,"\n");
	    }
	    if(kimura) {
		fprintf(tree,"\n Distances corrected by Kimura's 2 parameter model:");
		fprintf(tree,"\n\n Kimura, M. (1980)");
		fprintf(tree," A simple method for estimating evolutionary ");
		fprintf(tree,"rates of base");
		fprintf(tree,"\n substitutions through comparative studies of ");
		fprintf(tree,"nucleotide sequences.");
		fprintf(tree,"\n J. Mol. Evol., 16, 111-120.");
		fprintf(tree,"\n\n");
	    }
	}

	for(m=0;   m<nseqs;  m++)     /* for every pair of sequence */
	for(n=m+1; n<nseqs; n++) {
		p = q = e = 0.0;
		tmat[m][n] = tmat[n][m] = 0.0;
		for(i=0; i<seqs[0].len; i++) {
			j = boot_positions[i];
                    	if(!tossgaps || (tree_gaps[j] <= 0) ) {
				res1 = seqs[m].data[j];
				res2 = seqs[n].data[j];
				if( isalpha(res1) && isalpha(res2)) {
					if(!use_ambiguities)
					if(!is_ambiguity(res1) && !is_ambiguity(res2)) {
						e++;
                        			if(res1 != res2) {
							if(transition(res1,res2))
								p++;
							else
								q++;
						}
		        		}
				}
			}
		}


	/* Kimura's 2 parameter correction for multiple substitutions */

		if(!kimura) {
			if (e == 0) {
				fprintf(stdout,"\n WARNING: sequences %d and %d are non-overlapping\n",m,n);
				k = 0.0;
				p = 0.0;
				q = 0.0;
			}
			else {
				k = (p+q)/e;
				if(p > 0.0)
					p = p/e;
				else
					p = 0.0;
				if(q > 0.0)
					q = q/e;
				else
					q = 0.0;
			}
			tmat[m][n] = tmat[n][m] = k;
			if(tree!=NULL)                    /* if screen output */
				fprintf(tree,        
 	     "%4d vs.%4d:  DIST = %7.4f; p = %6.4f; q = %6.4f; length = %6.0f\n"
        	                 ,(pint)m,(pint)n,k,p,q,e);
		}
		else {
			if (e == 0) {
				fprintf(stdout,"\n WARNING: sequences %d and %d are non-overlapping\n",m,n);
				p = 0.0;
				q = 0.0;
			}
			else {
				if(p > 0.0)
					p = p/e;
				else
					p = 0.0;
				if(q > 0.0)
					q = q/e;
				else
					q = 0.0;
			}

			if( ((2.0*p)+q) == 1.0 )
				a = 0.0;
			else
				a = 1.0/(1.0-(2.0*p)-q);

			if( q == 0.5 )
				b = 0.0;
			else
				b = 1.0/(1.0-(2.0*q));

/* watch for values going off the scale for the correction. */
			if( (a<=0.0) || (b<=0.0) ) {
				overspill++;
				k = 3.5;  /* arbitrary high score */ 
			}
			else 
				k = 0.5*log(a) + 0.25*log(b);
			tmat[m][n] = tmat[n][m] = k;
			if(tree!=NULL)                      /* if screen output */
	   			fprintf(tree,
             "%4d vs.%4d:  DIST = %7.4f; p = %6.4f; q = %6.4f; length = %6.0f\n"
        	                ,(pint)m,(pint)n,k,p,q,e);

		}
	}
	tree_gaps=ckfree((void *)tree_gaps);
	return overspill;	/* return the number of off-scale values */
}

static sint prot_distance1_matrix(FILE *tree,SEQ *seqs, sint nseqs,double **tmat, COMP_MATRIX matrix, sint cutoff,Boolean tossgaps, Boolean kimura,sint *boot_positions)
{
        sint m,n;
        sint j,i;
        sint res1, res2;
    sint overspill = 0;
        double p,e,k, table_entry;
        char *tree_gaps;

        tree_gaps=tree_gap_delete(seqs,nseqs);  /* flag positions with gaps (tree_gaps[i] = 1 ) */

        if(tree!=NULL) {
                fprintf(tree,"\n");
                fprintf(tree,"\n DIST   = percentage divergence (/100)");
                fprintf(tree,"\n Length = number of sites used in comparison");
                fprintf(tree,"\n\n");
                if(tossgaps) {
                        fprintf(tree,"\n All sites with gaps (in any sequence) deleted");
                        fprintf(tree,"\n");
                }
                if(kimura) {
                        fprintf(tree,"\n Distances up tp 0.75 corrected by Kimura's empirical method:");
                        fprintf(tree,"\n\n Kimura, M. (1983)");
                        fprintf(tree," The Neutral Theory of Molecular Evolution.");
                        fprintf(tree,"\n Page 75. Cambridge University Press, Cambridge, England.");
                        fprintf(tree,"\n\n");
                }
        }

        for(m=0;   m<nseqs;  m++)     /* for every pair of sequence */
        for(n=m+1; n<nseqs; n++) {
                p = e = 0.0;
                tmat[m][n] = tmat[n][m] = 0.0;
                for(i=0; i<seqs[0].len;i++) {
                        j = boot_positions[i];
                        if(!(tossgaps && (tree_gaps[j] > 0))) {
                                res1 = seqs[m].data[j];
                                res2 = seqs[n].data[j];
                                if( isalpha(res1) && isalpha(res2)) {
                                        e++;
                                        if(res1!=res2) p++;
                                        if(matrix.score[res1-'a'][res2-'a']<=cutoff) p+=0.5;
                                }
                        }
                }

                if(e <= 0.0)
                        k = 0.0;
                else
                        k = p/e;

/* DES debug */
/* fprintf(stdout,"Seq1=%4d Seq2=%4d  k =%7.4f \n",(pint)m,(pint)n,k); */
/* DES debug */

                if(kimura) {
                        if(k < 0.75) { /* use Kimura's formula */
                                if(k > 0.0) k = - log(1.0 - k - (k * k/5.0) );
                        }
                        else {
                                if(k > 0.930) {
                                   overspill++;
                                   k = 10.0; /* arbitrarily set to 1000% */
                                }
                                else {
                                   table_entry = (k*1000.0) - 750.0;
                                   k = (double)dayhoff_pams[(int)table_entry];
                                   k = k/100.0;
                                }
                        }
                }

                tmat[m][n] = tmat[n][m] = k;
                    if(tree!=NULL)                    /* if screen output */
                        fprintf(tree,
                         "%4d vs.%4d  DIST = %6.4f;  length = %6.0f\n",
                         (pint)m+1,(pint)n+1,k,e);
        }
        tree_gaps=ckfree((void *)tree_gaps);
        return overspill;
}


/* calculates pairwise distance matrix - puts distances in array tmat */
static sint prot_distance_matrix(FILE *tree,SEQ *seqs, sint nseqs,double **tmat, Boolean tossgaps, Boolean kimura,sint *boot_positions)
{
	sint m,n;
	sint j,i;
	sint res1, res2;
    sint overspill = 0;
	double p,e,k, table_entry;	 
	char *tree_gaps;

	tree_gaps=tree_gap_delete(seqs,nseqs);  /* flag positions with gaps (tree_gaps[i] = 1 ) */
	
	if(tree!=NULL) {
		fprintf(tree,"\n");
		fprintf(tree,"\n DIST   = percentage divergence (/100)");
		fprintf(tree,"\n Length = number of sites used in comparison");
		fprintf(tree,"\n\n");
	        if(tossgaps) {
			fprintf(tree,"\n All sites with gaps (in any sequence) deleted");
			fprintf(tree,"\n");
		}
	    	if(kimura) {
			fprintf(tree,"\n Distances up tp 0.75 corrected by Kimura's empirical method:");
			fprintf(tree,"\n\n Kimura, M. (1983)");
			fprintf(tree," The Neutral Theory of Molecular Evolution.");
			fprintf(tree,"\n Page 75. Cambridge University Press, Cambridge, England.");
			fprintf(tree,"\n\n");
	    	}
	}

	for(m=0;   m<nseqs;  m++)     /* for every pair of sequence */
	for(n=m+1; n<nseqs; n++) {
		p = e = 0.0;
		tmat[m][n] = tmat[n][m] = 0.0;
		for(i=0; i<seqs[0].len;i++) {
			j = boot_positions[i];
	            	if(!(tossgaps && (tree_gaps[j] > 0))) {
				res1 = seqs[m].data[j];
				res2 = seqs[n].data[j];
				if( isalpha(res1) && isalpha(res2)) {
					e++;
                        		if(res1 != res2) p++;
		        	}
			}
		}
		if(e <= 0.0) 
			k = 1.0;
		else
			k = p/e;

/* DES debug */
/* fprintf(stdout,"Seq1=%4d Seq2=%4d  k =%7.4f \n",(pint)m,(pint)n,k); */
/* DES debug */

		if(kimura) {
			if(k < 0.75) { /* use Kimura's formula */
				if(k > 0.0) k = - log(1.0 - k - (k * k/5.0) );
			}
			else {
				if(k > 0.930) {
				   overspill++;
				   k = 10.0; /* arbitrarily set to 1000% */
				}
				else {
				   table_entry = (k*1000.0) - 750.0;
                                   k = (double)dayhoff_pams[(int)table_entry];
                                   k = k/100.0;
				}
			}
		}

		tmat[m][n] = tmat[n][m] = k;
		    if(tree!=NULL)                    /* if screen output */
			fprintf(tree,        
 	                 "%4d vs.%4d  DIST = %6.4f;  length = %6.0f\n",
 	                 (pint)m+1,(pint)n+1,k,e);
	}
	tree_gaps=ckfree((void *)tree_gaps);
	return overspill;
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
        static PTREE *standard_tree;
        sint i;
	float dist;

        if(nseqs==2) {
                dist=tmat[0][1]/2.0;
                fprintf(ofile,"(%s:%0.5f,%s:%0.5f);\n",
                        seqs[0].name,dist,seqs[1].name,dist);
        }
        else {
		if(algo==BIONJ) {
			bionj(ofile,seqs,nseqs,tmat);
		}
		else if (algo==NJ) {
        		standard_tree   = (PTREE *) ckalloc( (nseqs) * sizeof (PTREE) );
        		for(i=0; i<nseqs; i++)
                		standard_tree[i].description  = (char *) ckalloc( (nseqs+1) * sizeof(char));

        		nj_tree(NULL,nseqs,standard_tree,tmat);

        		print_phylip_tree(ofile,seqs,nseqs,standard_tree,tmat);

        		for (i=0;i<nseqs;i++)
                		ckfree((void *)standard_tree[i].description);
        		standard_tree=ckfree((void *)standard_tree);

		}
	}
        fclose(ofile);

}

static Boolean is_ambiguity(char c)
{
        int i;
	char codes[]="ACGTU";

	for(i=0;i<5;i++)
        	if(c==codes[i])
             	   return FALSE;

        return TRUE;
}

