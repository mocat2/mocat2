#include "mltaln.h"


void profilealignment2( int n0, int n2, char **aln0, char **aln2, int alloclen, char alg ) // n1 ha allgap
{
	int i, newlen;
	double *effarr0, *effarr2;
	float dumfl;
	double eff;
	effarr0 = AllocateDoubleVec( n0 );
	effarr2 = AllocateDoubleVec( n2 );

	commongappick( n0, aln0 );
	commongappick( n2, aln2 );

	eff = 1.0 / (double)n0; for( i=0; i<n0; i++ ) effarr0[i] = eff;
	eff = 1.0 / (double)n2; for( i=0; i<n2; i++ ) effarr2[i] = eff;

	newgapstr = "-";
	if( alg == 'M' )
		MSalignmm( n_dis_consweight_multi, aln0, aln2, effarr0, effarr2, n0, n2, alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 ); //outgap=1??
	else
		A__align( n_dis_consweight_multi, aln0, aln2, effarr0, effarr2, n0, n2, alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 ); //outgap=1??

	newlen = strlen( aln0[0] );

#if 0 // tabun hitsuyou
	for( j=0; j<newlen; j++ )
	{
//		fprintf( stderr, "j=%d\n", j );
		for( i=0; i<n0; i++ )
		{
			if( aln0[i][j] != '-' ) break;
		}
		if( i == n0 ) 
		{
			for( i=0; i<n1; i++ ) 
			{
				if( aln1[i][j] != '-' ) break;
			}
		}
		else i = -1;

		if( i == n1 ) 
		{
			for( i=0; i<n1; i++ ) aln1[i][j] = '=';
		}
	}
	fprintf( stderr, "in profilealignment,\n" );
	for( i=0; i<n0; i++ ) fprintf( stderr, "\n>aln0[%d] = \n%s\n", i, aln0[i] );
	for( i=0; i<n1; i++ ) fprintf( stderr, "\n>aln1[%d] = \n%s\n", i, aln1[i] );
	for( i=0; i<n2; i++ ) fprintf( stderr, "\n>aln2[%d] = \n%s\n", i, aln2[i] );
#endif

	free( effarr0 );
	free( effarr2 );
}

void profilealignment( int n0, int n1, int n2, char **aln0, char **aln1, char **aln2, int alloclen, char alg ) // n1 ha allgap
{
	int i, j, newlen;
	double *effarr0, *effarr2;
	float dumfl;
	double eff;
	effarr0 = AllocateDoubleVec( n0 );
	effarr2 = AllocateDoubleVec( n2 );

	commongappick( n0, aln0 );
	commongappick( n2, aln2 );

	eff = 1.0 / (double)n0; for( i=0; i<n0; i++ ) effarr0[i] = eff;
	eff = 1.0 / (double)n2; for( i=0; i<n2; i++ ) effarr2[i] = eff;

	newgapstr = "-";
	if( alg == 'M' )
		MSalignmm( n_dis_consweight_multi, aln0, aln2, effarr0, effarr2, n0, n2, alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 ); //outgap=1??
	else
		A__align( n_dis_consweight_multi, aln0, aln2, effarr0, effarr2, n0, n2, alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 ); //outgap=1??

	newlen = strlen( aln0[0] );

	for( i=0; i<newlen; i++ ) aln1[0][i] = '-';
	aln1[0][i] = 0;
	for( i=1; i<n1; i++ ) strcpy( aln1[i], aln1[0] );

	for( j=0; j<newlen; j++ )
	{
//		fprintf( stderr, "j=%d\n", j );
		for( i=0; i<n0; i++ )
		{
			if( aln0[i][j] != '-' ) break;
		}
		if( i == n0 ) 
		{
			for( i=0; i<n1; i++ ) 
			{
				if( aln1[i][j] != '-' ) break;
			}
		}
		else i = -1;

		if( i == n1 ) 
		{
			for( i=0; i<n1; i++ ) aln1[i][j] = '=';
		}
	}
#if 0
	fprintf( stderr, "in profilealignment,\n" );
	for( i=0; i<n0; i++ ) fprintf( stderr, "\n>aln0[%d] = \n%s\n", i, aln0[i] );
	for( i=0; i<n1; i++ ) fprintf( stderr, "\n>aln1[%d] = \n%s\n", i, aln1[i] );
	for( i=0; i<n2; i++ ) fprintf( stderr, "\n>aln2[%d] = \n%s\n", i, aln2[i] );
#endif

	free( effarr0 );
	free( effarr2 );
}

void eq2dash( char *s )
{
	while( *s )
	{
		if( *s == '=' ) *s = '-';
		s++;
	}
}

void findnewgaps( int n, int rep, char **seq, int *gaplen )
{
	int i, pos, len, len1;

	len = strlen( seq[0] );	
//	for( i=0; i<len; i++ ) gaplen[i] = 0; // calloc de shokika sareteirukara hontou ha iranai
	len1 = len + 1;
	for( i=0; i<len1; i++ ) gaplen[i] = 0; // reallo de shokika sareteirukara iru!
	
	pos = 0;
	for( i=0; i<len; i++ )
	{
		if( seq[rep][i] == '=' ) 
		{
//			fprintf( stderr, "Newgap! pos = %d\n", pos );
			gaplen[pos]++;
		}
		else
			pos++;
	}
}

void findcommongaps( int n, char **seq, int *gaplen )
{
	int i, j, pos, len, len1;
	len = strlen( seq[0] );	
	len1 = len+1;

//	fprintf( stderr, "seq[0] = %s\n", seq[0] );
	for( i=0; i<len1; i++ ) gaplen[i] = 0;
	
	pos = 0;
	for( i=0; i<len; i++ )
	{
		for( j=0; j<n; j++ )
			if( seq[j][i] != '-' ) break;

		if( j == n ) gaplen[pos]++;
		else
			pos++;
	}
#if 0
	for( i=0; i<pos; i++ )
	{
		fprintf( stderr, "vec[%d] = %d\n", i, gaplen[i] );
	}
#endif
}

void adjustgapmap( int newlen, int *gapmap, char *seq )
{
	int j;
	int pos;
	int newlen1 = newlen+1;
	int *tmpmap;


	tmpmap = AllocateIntVec( newlen+2 );
	j = 0;
	pos = 0;
	while( *seq )
	{
//		fprintf( stderr, "j=%d *seq = %c\n", j, *seq );
		if( *seq++ == '=' )
			tmpmap[j++] = 0;
		else
		{
			tmpmap[j++] = gapmap[pos++];
		}
	}
	tmpmap[j++] = gapmap[pos];

	for(j=0; j<newlen1; j++)
		gapmap[j] = tmpmap[j];

	free( tmpmap );
}

static void strncpy0( char *s1, char *s2, int n )
{
	while( n-- ) *s1++ = *s2++;
	*s1 = 0;
}

static int countnogaplen( int *gaplen, int *term )
{
	int v = 0;
	while( gaplen < term )
	{
		if( *gaplen++ == 0 ) v++;
		else break;
	}
	return( v );
}

void insertnewgaps( int njob, int *alreadyaligned, char **seq, int *ex1, int *ex2, int *gaplen, int *gapmap, int alloclen, char alg, char gapchar )
{
	int *mar;
	char *gaps;
	char *cptr;
	int i, j, k, len, rep, len0, lp, blocklen;
	char **mseq2, **mseq0, **mseq1;
	char **aseq, *newchar;
	int ngroup2, ngroup0, ngroup1;
	int *list0, *list1, *list2;
	int posin12, gapshift, newpos;
	int mlen1, mlen0, mlen2;
	int tmptmptmpmark = 0;

	mar = calloc( njob, sizeof( int ) );
	list0 = calloc( njob, sizeof( int ) );
	list1 = calloc( njob, sizeof( int ) );
	list2 = calloc( njob, sizeof( int ) );

	for( i=0; i<njob; i++ ) mar[i] = 0;
	for( i=0; i<njob; i++ ) 
	{
		if( alreadyaligned[i]==0 ) mar[i] = 3;
	}
	for( i=0; (k=ex1[i])>-1; i++ ) 
	{
		mar[k] = 1;
//		fprintf( stderr, "excluding %d\n", ex1[i] );
	}
	for( i=0; (k=ex2[i])>-1; i++ ) 
	{
		mar[k] = 2;
//		fprintf( stderr, "excluding %d\n", ex2[i] );
	}

	ngroup2 = ngroup1 = ngroup0 = 0;
	for( i=0; i<njob; i++ )
	{
		if( mar[i] == 2 ) 
		{
			list2[ngroup2] = i;
			ngroup2++;
		}
		if( mar[i] == 1 ) 
		{
			list1[ngroup1] = i;
			ngroup1++;
		}
		if( mar[i] == 0 ) 
		{
			list0[ngroup0] = i;
//			fprintf( stderr, "inserting new gaps to %d\n", i );
			ngroup0++;
		}
	}
	list0[ngroup0] = list1[ngroup1] = list2[ngroup2] = -1;
	if( ngroup0 == 0 )
	{
//		fprintf( stderr, "Nothing to do\n" );
		free( mar );
		free( list0 );
		free( list1 );
		free( list2 );
		return;
	}

	for( i=0; i<njob; i++ ) if( mar[i] == 0 ) break;
	rep = i;
	len = strlen( seq[rep] );
	len0 = len+1;

//
//	if( i == njob )
//	{
////		fprintf( stderr, "Nothing to do\n" );
//		free( mar );
//		return;
//	}

	mseq2 = AllocateCharMtx( ngroup2, alloclen );
	mseq1 = AllocateCharMtx( ngroup1, alloclen );
	mseq0 = AllocateCharMtx( ngroup0, alloclen );
	aseq = AllocateCharMtx( njob, alloclen );
	gaps = calloc( alloclen, sizeof( char ) );


	for( i=0; i<njob; i++ ) aseq[i][0] = 0;
	newpos = 0;
	posin12 = 0;
//	fprintf( stderr, "\ngaplen[] = \n" );
//	for(i=0; i<len0; i++ ) fprintf( stderr, "%d ", gaplen[i] );
//	fprintf( stderr, "\n" );

	for( j=0; j<len0; j++ )
	{
//		fprintf( stderr, "\nj=%d, gaplen[%d]=%d\n", j, j, gaplen[j] );
		if( gaplen[j] )
		{
//			fprintf( stderr, "j=%d GAP!\n", j );
			tmptmptmpmark = 1;
			for( i=0; i<ngroup0; i++ ) mseq0[i][0] = 0;
			for( i=0; i<ngroup1; i++ ) mseq1[i][0] = 0;
			for( i=0; i<ngroup2; i++ ) mseq2[i][0] = 0;
			mlen0 = mlen1 = mlen2 = 0;

			gapshift = gaplen[j];
			cptr = gaps;
			while( gapshift-- ) *cptr++ = gapchar;
			*cptr = 0;
			gapshift = gaplen[j];

			for( i=0; i<ngroup0; i++ ) strncpy0( mseq0[i]+mlen0, gaps, gapshift );
			for( i=0; i<ngroup1; i++ ) strncpy0( mseq1[i]+mlen1, seq[list1[i]]+posin12, gapshift );
			for( i=0; i<ngroup2; i++ ) strncpy0( mseq2[i]+mlen2, seq[list2[i]]+posin12, gapshift );
			posin12 += gapshift;
			mlen0 += gapshift;
			mlen1 += gapshift;
			mlen2 += gapshift;

			gapshift = gapmap[posin12];
//			fprintf( stderr, "gapmap[%d] kouho = %d\n", posin12, gapmap[posin12] );


			for( i=0; i<ngroup0; i++ ) strncpy0( mseq0[i]+mlen0, seq[list0[i]]+j, gapshift );
			for( i=0; i<ngroup1; i++ ) strncpy0( mseq1[i]+mlen1, seq[list1[i]]+posin12, gapshift );
			for( i=0; i<ngroup2; i++ ) strncpy0( mseq2[i]+mlen2, seq[list2[i]]+posin12, gapshift );
			mlen0 += gapshift;
			mlen1 += gapshift;
			mlen2 += gapshift;
#if 0
			for( i=0; i<ngroup0; i++ ) fprintf( stderr, "### mseq0[%d] = %s\n", i, mseq0[i] );
			for( i=0; i<ngroup1; i++ ) fprintf( stderr, "### mseq1[%d] = %s\n", i, mseq1[i] );
			for( i=0; i<ngroup2; i++ ) fprintf( stderr, "### mseq2[%d] = %s\n", i, mseq2[i] );
#endif

			if( gapshift ) 
				profilealignment( ngroup0, ngroup1, ngroup2, mseq0, mseq1, mseq2, alloclen, alg );

			j += gapshift;
			posin12 += gapshift;

			newpos = strlen( aseq[rep] ); // kufuu?
			for( i=0; i<ngroup0; i++ ) strcpy( aseq[list0[i]]+newpos, mseq0[i] );
			for( i=0; i<ngroup1; i++ ) strcpy( aseq[list1[i]]+newpos, mseq1[i] );
			for( i=0; i<ngroup2; i++ ) strcpy( aseq[list2[i]]+newpos, mseq2[i] );

//			fprintf( stderr, "gapshift = %d\n", gapshift );
		}
		blocklen = 1 + countnogaplen( gaplen+j+1, gaplen+len0 );
//		fprintf( stderr, "\nj=%d, blocklen=%d, len0=%d\n", j, blocklen, len0 );
//		blocklen = 1;
//		if( tmptmptmpmark ) exit( 1 );

		newpos = strlen( aseq[rep] );

#if 0
		for( i=0; i<ngroup0; i++ ) aseq[list0[i]][newpos] = seq[list0[i]][j];
		for( i=0; i<ngroup1; i++ ) aseq[list1[i]][newpos] = seq[list1[i]][posin12];
		for( i=0; i<ngroup2; i++ ) aseq[list2[i]][newpos] = seq[list2[i]][posin12];
		for( i=0; i<ngroup0; i++ ) aseq[list0[i]][newpos+1] = 0;
		for( i=0; i<ngroup1; i++ ) aseq[list1[i]][newpos+1] = 0;
		for( i=0; i<ngroup2; i++ ) aseq[list2[i]][newpos+1] = 0;
#else

		for( i=0; i<ngroup0; i++ )
		{
			lp = list0[i];
			newchar = aseq[lp] + newpos;
			strncpy0( newchar, seq[lp]+j, blocklen );
		}
		for( i=0; i<ngroup1; i++ )
		{
			lp = list1[i];
			newchar = aseq[lp] + newpos;
			strncpy0( newchar, seq[lp]+posin12, blocklen );
		}
		for( i=0; i<ngroup2; i++ )
		{
			lp = list2[i];
			newchar = aseq[lp] + newpos;
			strncpy0( newchar, seq[lp]+posin12, blocklen );
		}
//		fprintf( stderr, "### aseq[l0] = %s\n", aseq[list0[0]] );
//		fprintf( stderr, "### aseq[l1] = %s\n", aseq[list1[0]] );
//		fprintf( stderr, "### aseq[l2] = %s\n", aseq[list2[0]] );
//		exit( 1 );
#endif

//		fprintf( stderr, "j=%d -> %d\n", j, j+blocklen-1 );
		j += (blocklen-1);


		posin12 += (blocklen-1);


		posin12++;
	}
#if 0
	fprintf( stderr, "\n" );
	fprintf( stderr, " seq[l0] = %s\n", seq[list0[0]] );
	fprintf( stderr, " seq[l1] = %s\n", seq[list1[0]] );
	fprintf( stderr, " seq[l2] = %s\n", seq[list2[0]] );
	fprintf( stderr, "=====>\n" );
	fprintf( stderr, "aseq[l0] = %s\n", aseq[list0[0]] );
	fprintf( stderr, "aseq[l1] = %s\n", aseq[list1[0]] );
	fprintf( stderr, "aseq[l2] = %s\n", aseq[list2[0]] );
//if( tmptmptmpmark ) exit( 1 );
#endif

//	for( i=0; i<njob; i++ ) if( mar[i] != 3 ) strcpy( seq[i], aseq[i] );
	for( i=0; i<ngroup0; i++ ) strcpy( seq[list0[i]], aseq[list0[i]] );
	for( i=0; i<ngroup1; i++ ) strcpy( seq[list1[i]], aseq[list1[i]] );
	for( i=0; i<ngroup2; i++ ) strcpy( seq[list2[i]], aseq[list2[i]] );


	free( mar );
	free( gaps );
	free( list0 );
	free( list1 );
	free( list2 );
	FreeCharMtx( mseq2 );
	FreeCharMtx( mseq1 ); // ? added 2012/02/12
	FreeCharMtx( mseq0 );
	FreeCharMtx( aseq ); // ? added 2012/02/12
}


void restorecommongaps( int njob, char **seq, int *ex1, int *ex2, int *gaplen, int alloclen, char gapchar )
{
	int *mar;
	char *tmpseq;
	char *cptr;
	int *iptr;
	int *tmpgaplen;
	int i, j, k, len, rep, len1;

	mar = calloc( njob, sizeof( int ) );
	tmpseq = calloc( alloclen, sizeof( char ) );
	tmpgaplen = calloc( alloclen, sizeof( int ) );
//	tmpseq = calloc( alloclen+2, sizeof( char ) );
//	tmpgaplen = calloc( alloclen+2, sizeof( int ) );


	for( i=0; i<njob; i++ ) mar[i] = 0;
	for( i=0; (k=ex1[i])>-1; i++ ) 
	{
		mar[k] = 1;
//		fprintf( stderr, "excluding %d\n", ex1[i] );
	}
	for( i=0; (k=ex2[i])>-1; i++ ) 
	{
		mar[k] = 1;
//		fprintf( stderr, "excluding %d\n", ex2[i] );
	}

	for( i=0; i<njob; i++ )
		if( mar[i] ) break;

	if( i == njob )
	{
//		fprintf( stderr, "Nothing to do\n" );
		free( mar );
		free( tmpseq );
		free( tmpgaplen );
		return;
	}
	rep = i;
	len = strlen( seq[rep] );
	len1 = len+1;


	for( i=0; i<njob; i++ )
	{
		if( mar[i] == 0 ) continue;
		cptr = tmpseq;
		for( j=0; j<len1; j++ )
		{
			for( k=0; k<gaplen[j]; k++ )
				*(cptr++) = gapchar; // ???
			*(cptr++) = seq[i][j];
		}
		*cptr = 0;
		strcpy( seq[i], tmpseq );
	}

	iptr = tmpgaplen;
	for( j=0; j<len1; j++ )
	{
		*(iptr++) = gaplen[j];
		for( k=0; k<gaplen[j]; k++ )
			*(iptr++) = 0;
	}
	*iptr = -1;

	iptr = tmpgaplen;
	while( *iptr != -1 ) *gaplen++ = *iptr++;

	free( mar );
	free( tmpseq );
	free( tmpgaplen );
}
