#include "mltaln.h"
#include "dp.h"

#define DEBUG 0
#define XXXXXXX    0
#define USE_PENALTY_EX  1


#if 1
static void match_calc_mtx( double **mtx, float *match, char **s1, char **s2, int i1, int lgth2 ) 
{
	char *seq2 = s2[0];
	double *doubleptr = mtx[(int)s1[0][i1]];

	while( lgth2-- )
		*match++ = doubleptr[(int)*seq2++];
}
#else
static void match_calc( float *match, char **s1, char **s2, int i1, int lgth2 )
{
	int j;

	for( j=0; j<lgth2; j++ )
		match[j] = amino_dis[(*s1)[i1]][(*s2)[j]];
}
#endif

static float Atracking( float *lasthorizontalw, float *lastverticalw, 
						char **seq1, char **seq2, 
                        char **mseq1, char **mseq2, 
                        int **ijp,
						int tailgp )
{
	int i, j, l, iin, jin, ifi, jfi, lgth1, lgth2, k, limk;
//	char gap[] = "-";
	char *gap;
	gap = newgapstr;
	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );
	float wm;


#if 0
	for( i=0; i<lgth1; i++ ) 
	{
		fprintf( stderr, "lastverticalw[%d] = %f\n", i, lastverticalw[i] );
	}
#endif
 
    for( i=0; i<lgth1+1; i++ ) 
    {
        ijp[i][0] = i + 1;
    }
    for( j=0; j<lgth2+1; j++ ) 
    {
        ijp[0][j] = -( j + 1 );
    }

	if( tailgp == 1 )
		;
	else
	{
		wm = lastverticalw[0];
		for( i=0; i<lgth1; i++ )
		{
			if( lastverticalw[i] >= wm )
			{
				wm = lastverticalw[i];
				iin = i; jin = lgth2-1;
				ijp[lgth1][lgth2] = +( lgth1 - i );
			}
		}
		for( j=0; j<lgth2; j++ )
		{
			if( lasthorizontalw[j] >= wm )
			{
				wm = lasthorizontalw[j];
				iin = lgth1-1; jin = j;
				ijp[lgth1][lgth2] = -( lgth2 - j );
			}
		}
	}



	mseq1[0] += lgth1+lgth2;
	*mseq1[0] = 0;
	mseq2[0] += lgth1+lgth2;
	*mseq2[0] = 0;


	iin = lgth1; jin = lgth2;
	limk = lgth1+lgth2 + 1;
	for( k=0; k<limk; k++ ) 
	{
		if( ijp[iin][jin] < 0 ) 
		{
			ifi = iin-1; jfi = jin+ijp[iin][jin];
		}
		else if( ijp[iin][jin] > 0 )
		{
			ifi = iin-ijp[iin][jin]; jfi = jin-1;
		}
		else
		{
			ifi = iin-1; jfi = jin-1;
		}
		l = iin - ifi;
		while( --l ) 
		{
			*--mseq1[0] = seq1[0][ifi+l];
			*--mseq2[0] = *gap;
			k++;
		}
		l= jin - jfi;
		while( --l )
		{
			*--mseq1[0] = *gap;
			*--mseq2[0] = seq2[0][jfi+l];
			k++;
		}
		if( iin <= 0 || jin <= 0 ) break;
		*--mseq1[0] = seq1[0][ifi];
		*--mseq2[0] = seq2[0][jfi];
		k++;
		iin = ifi; jin = jfi;
	}
	return( 0.0 );
}


float G__align11( double **n_dynamicmtx, char **seq1, char **seq2, int alloclen, int headgp, int tailgp )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{
//	int k;
	register int i, j;
	int lasti;                      /* outgap == 0 -> lgth1, outgap == 1 -> lgth1+1 */
	int lgth1, lgth2;
	int resultlen;
	float wm;   /* int ?????? */
	float g;
	float *currentw, *previousw;
	float fpenalty = (float)penalty;
#if USE_PENALTY_EX
	float fpenalty_ex = (float)penalty_ex;
#endif
#if 1
	float *wtmp;
	int *ijppt;
	float *mjpt, *prept, *curpt;
	int *mpjpt;
#endif
	static TLS float mi = 0.0;
	static TLS float *m = NULL;
	static TLS int **ijp = NULL;
	static TLS int mpi = 0;
	static TLS int *mp = NULL;
	static TLS float *w1 = NULL;
	static TLS float *w2 = NULL;
	static TLS float *match = NULL;
	static TLS float *initverticalw = NULL;    /* kufuu sureba iranai */
	static TLS float *lastverticalw = NULL;    /* kufuu sureba iranai */
	static TLS char **mseq1 = NULL;
	static TLS char **mseq2 = NULL;
	static TLS char **mseq = NULL;
	static TLS int **intwork = NULL;
	static TLS float **floatwork = NULL;
	static TLS int orlgth1 = 0, orlgth2 = 0;
	static TLS double **amino_dynamicmtx = NULL; // ??


	if( seq1 == NULL )
	{
		if( orlgth1 > 0 && orlgth2 > 0 )
		{
			orlgth1 = 0;
			orlgth2 = 0;
			if( mseq1 ) free( mseq1 ); mseq1 = NULL;
			if( mseq2 ) free( mseq2 ); mseq2 = NULL;
			if( w1 ) FreeFloatVec( w1 ); w1 = NULL;
			if( w2 ) FreeFloatVec( w2 ); w2 = NULL;
			if( match ) FreeFloatVec( match ); match = NULL;
			if( initverticalw ) FreeFloatVec( initverticalw ); initverticalw = NULL;
			if( lastverticalw ) FreeFloatVec( lastverticalw ); lastverticalw = NULL;

			if( m ) FreeFloatVec( m ); m = NULL;
			if( mp ) FreeIntVec( mp ); mp = NULL;

			if( mseq ) FreeCharMtx( mseq ); mseq = NULL;



			if( floatwork ) FreeFloatMtx( floatwork ); floatwork = NULL;
			if( intwork ) FreeIntMtx( intwork ); intwork = NULL;

			if( amino_dynamicmtx ) FreeDoubleMtx( amino_dynamicmtx ); amino_dynamicmtx = NULL;
		}
		return( 0.0 );
	}


	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );

#if 0
	if( lgth1 <= 0 || lgth2 <= 0 )
	{
		fprintf( stderr, "WARNING (g11): lgth1=%d, lgth2=%d\n", lgth1, lgth2 );
	}
#endif

#if 1
	if( lgth1 == 0 && lgth2 == 0 )
		return( 0.0 );

	if( lgth1 == 0 )
	{
		seq1[0][lgth2] = 0;
		while( lgth2 ) seq1[0][--lgth2] = *newgapstr;
//		fprintf( stderr, "seq1[0] = %s\n", seq1[0] );
		return( 0.0 );
	}

	if( lgth2 == 0 )
	{
		seq2[0][lgth1] = 0;
		while( lgth1 ) seq2[0][--lgth1] = *newgapstr;
//		fprintf( stderr, "seq2[0] = %s\n", seq2[0] );
		return( 0.0 );
	}
#endif


	wm = 0.0;

	if( orlgth1 == 0 )
	{
		mseq1 = AllocateCharMtx( njob, 0 );
		mseq2 = AllocateCharMtx( njob, 0 );
	}



	if( lgth1 > orlgth1 || lgth2 > orlgth2 )
	{
		int ll1, ll2;

		if( orlgth1 > 0 && orlgth2 > 0 )
		{
			FreeFloatVec( w1 );
			FreeFloatVec( w2 );
			FreeFloatVec( match );
			FreeFloatVec( initverticalw );
			FreeFloatVec( lastverticalw );

			FreeFloatVec( m );
			FreeIntVec( mp );

			FreeCharMtx( mseq );



			FreeFloatMtx( floatwork );
			FreeIntMtx( intwork );
			FreeDoubleMtx( amino_dynamicmtx );
		}

		ll1 = MAX( (int)(1.3*lgth1), orlgth1 ) + 100;
		ll2 = MAX( (int)(1.3*lgth2), orlgth2 ) + 100;

#if DEBUG
		fprintf( stderr, "\ntrying to allocate (%d+%d)xn matrices ... ", ll1, ll2 );
#endif

		w1 = AllocateFloatVec( ll2+2 );
		w2 = AllocateFloatVec( ll2+2 );
		match = AllocateFloatVec( ll2+2 );

		initverticalw = AllocateFloatVec( ll1+2 );
		lastverticalw = AllocateFloatVec( ll1+2 );

		m = AllocateFloatVec( ll2+2 );
		mp = AllocateIntVec( ll2+2 );

		mseq = AllocateCharMtx( njob, ll1+ll2 );


		floatwork = AllocateFloatMtx( nalphabets, MAX( ll1, ll2 )+2 ); 
		intwork = AllocateIntMtx( nalphabets, MAX( ll1, ll2 )+2 ); 
		amino_dynamicmtx = AllocateDoubleMtx( 0x80, 0x80 );

#if DEBUG
		fprintf( stderr, "succeeded\n" );
#endif

		orlgth1 = ll1 - 100;
		orlgth2 = ll2 - 100;
	}
    for( i=0; i<nalphabets; i++) for( j=0; j<nalphabets; j++ )
		amino_dynamicmtx[(int)amino[i]][(int)amino[j]] = (double)n_dynamicmtx[i][j];


	mseq1[0] = mseq[0];
	mseq2[0] = mseq[1];


	if( orlgth1 > commonAlloc1 || orlgth2 > commonAlloc2 )
	{
		int ll1, ll2;

		if( commonAlloc1 && commonAlloc2 )
		{
			FreeIntMtx( commonIP );
		}

		ll1 = MAX( orlgth1, commonAlloc1 );
		ll2 = MAX( orlgth2, commonAlloc2 );

#if DEBUG
		fprintf( stderr, "\n\ntrying to allocate %dx%d matrices ... ", ll1+1, ll2+1 );
#endif

		commonIP = AllocateIntMtx( ll1+10, ll2+10 );

#if DEBUG
		fprintf( stderr, "succeeded\n\n" );
#endif

		commonAlloc1 = ll1;
		commonAlloc2 = ll2;
	}
	ijp = commonIP;


#if 0
	for( i=0; i<lgth1; i++ ) 
		fprintf( stderr, "ogcp1[%d]=%f\n", i, ogcp1[i] );
#endif

	currentw = w1;
	previousw = w2;


	match_calc_mtx( amino_dynamicmtx, initverticalw, seq2, seq1, 0, lgth1 );
	match_calc_mtx( amino_dynamicmtx, currentw, seq1, seq2, 0, lgth2 );

	if( headgp == 1 )
	{
		for( i=1; i<lgth1+1; i++ )
		{
			initverticalw[i] += fpenalty;
		}
		for( j=1; j<lgth2+1; j++ )
		{
			currentw[j] += fpenalty;
		}
	}

	for( j=1; j<lgth2+1; ++j ) 
	{
		m[j] = currentw[j-1]; mp[j] = 0;
	}

	if( lgth2 == 0 )
		lastverticalw[0] = 0.0;               // lgth2==0 no toki error
	else
		lastverticalw[0] = currentw[lgth2-1]; // lgth2==0 no toki error

	if( tailgp ) lasti = lgth1+1; else lasti = lgth1;

#if XXXXXXX
fprintf( stderr, "currentw = \n" );
for( i=0; i<lgth1+1; i++ )
{
	fprintf( stderr, "%5.2f ", currentw[i] );
}
fprintf( stderr, "\n" );
fprintf( stderr, "initverticalw = \n" );
for( i=0; i<lgth2+1; i++ )
{
	fprintf( stderr, "%5.2f ", initverticalw[i] );
}
fprintf( stderr, "\n" );
#endif

	for( i=1; i<lasti; i++ )
	{
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

		match_calc_mtx( amino_dynamicmtx, currentw, seq1, seq2, i, lgth2 );
#if XXXXXXX
fprintf( stderr, "\n" );
fprintf( stderr, "i=%d\n", i );
fprintf( stderr, "currentw = \n" );
for( j=0; j<lgth2; j++ )
{
	fprintf( stderr, "%5.2f ", currentw[j] );
}
fprintf( stderr, "\n" );
#endif
#if XXXXXXX
fprintf( stderr, "\n" );
fprintf( stderr, "i=%d\n", i );
fprintf( stderr, "currentw = \n" );
for( j=0; j<lgth2; j++ )
{
	fprintf( stderr, "%5.2f ", currentw[j] );
}
fprintf( stderr, "\n" );
#endif
		currentw[0] = initverticalw[i];

		mi = previousw[0]; mpi = 0;

		ijppt = ijp[i] + 1;
		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;
		mpjpt = mp + 1;
		for( j=1; j<lgth2+1; j++ )
		{
			wm = *prept;
			*ijppt = 0;

#if 0
			fprintf( stderr, "%5.0f->", wm );
#endif
#if 0
			fprintf( stderr, "%5.0f?", g );
#endif
			if( (g=mi+fpenalty) > wm )
			{
				wm = g;
				*ijppt = -( j - mpi );
			}
			if( (g=*prept) >= mi )
			{
				mi = g;
				mpi = j-1;
			}
#if USE_PENALTY_EX
			mi += fpenalty_ex;
#endif

#if 0 
			fprintf( stderr, "%5.0f?", g );
#endif
			if( (g=*mjpt + fpenalty) > wm )
			{
				wm = g;
				*ijppt = +( i - *mpjpt );
			}
			if( (g=*prept) >= *mjpt )
			{
				*mjpt = g;
				*mpjpt = i-1;
			}
#if USE_PENALTY_EX
			m[j] += fpenalty_ex;
#endif

#if 0
			fprintf( stderr, "%5.0f ", wm );
#endif
			*curpt++ += wm;
			ijppt++;
			mjpt++;
			prept++;
			mpjpt++;
		}
		lastverticalw[i] = currentw[lgth2-1]; // lgth2==0 no toki error
	}

	Atracking( currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, tailgp );

	resultlen = strlen( mseq1[0] );
	if( alloclen < resultlen || resultlen > N )
	{
		fprintf( stderr, "alloclen=%d, resultlen=%d, N=%d\n", alloclen, resultlen, N );
		ErrorExit( "LENGTH OVER!\n" );
	}


	strcpy( seq1[0], mseq1[0] );
	strcpy( seq2[0], mseq2[0] );
#if 0
	fprintf( stderr, "\n" );
	fprintf( stderr, ">\n%s\n", mseq1[0] );
	fprintf( stderr, ">\n%s\n", mseq2[0] );
	fprintf( stderr, "wm = %f\n", wm );
#endif

	return( wm );
}

float G__align11_noalign( double **n_dynamicmtx, int penal, int penal_ex, char **seq1, char **seq2, int alloclen )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{
//	int k;
	register int i, j;
	int lasti;                      /* outgap == 0 -> lgth1, outgap == 1 -> lgth1+1 */
	int lgth1, lgth2;
//	int resultlen;
	float wm;   /* int ?????? */
	float g;
	float *currentw, *previousw;
	float fpenalty = (float)penal;
#if USE_PENALTY_EX
	float fpenalty_ex = (float)penal_ex;
#endif
#if 1
	float *wtmp;
	float *mjpt, *prept, *curpt;
//	int *mpjpt;
#endif
	static TLS float mi, *m;
	static TLS float *w1, *w2;
	static TLS float *match;
	static TLS float *initverticalw;    /* kufuu sureba iranai */
	static TLS float *lastverticalw;    /* kufuu sureba iranai */
	static TLS int **intwork;
	static TLS float **floatwork;
	static TLS int orlgth1 = 0, orlgth2 = 0;
	static TLS double **amino_dynamicmtx;

	if( seq1 == NULL )
	{
		if( orlgth1 > 0 && orlgth2 > 0 )
		{
			orlgth1 = 0;
			orlgth2 = 0;
			FreeFloatVec( w1 );
			FreeFloatVec( w2 );
			FreeFloatVec( match );
			FreeFloatVec( initverticalw );
			FreeFloatVec( lastverticalw );
			free( m );


			FreeFloatMtx( floatwork );
			FreeIntMtx( intwork );
			FreeDoubleMtx( amino_dynamicmtx );
		}
		return( 0.0 );
	}


	wm = 0.0;



	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );



#if 0
	if( lgth1 <= 0 || lgth2 <= 0 )
	{
		fprintf( stderr, "WARNING (g11): lgth1=%d, lgth2=%d\n", lgth1, lgth2 );
	}
#endif

	if( lgth1 > orlgth1 || lgth2 > orlgth2 )
	{
		int ll1, ll2;

		if( orlgth1 > 0 && orlgth2 > 0 )
		{
			FreeFloatVec( w1 );
			FreeFloatVec( w2 );
			FreeFloatVec( match );
			FreeFloatVec( initverticalw );
			FreeFloatVec( lastverticalw );

			FreeFloatVec( m );




			FreeFloatMtx( floatwork );
			FreeIntMtx( intwork );

			FreeDoubleMtx( amino_dynamicmtx );
		}

		ll1 = MAX( (int)(1.3*lgth1), orlgth1 ) + 100;
		ll2 = MAX( (int)(1.3*lgth2), orlgth2 ) + 100;

#if DEBUG
		fprintf( stderr, "\ntrying to allocate (%d+%d)xn matrices ... ", ll1, ll2 );
#endif

		w1 = AllocateFloatVec( ll2+2 );
		w2 = AllocateFloatVec( ll2+2 );
		match = AllocateFloatVec( ll2+2 );

		initverticalw = AllocateFloatVec( ll1+2 );
		lastverticalw = AllocateFloatVec( ll1+2 );

		m = AllocateFloatVec( ll2+2 );



		floatwork = AllocateFloatMtx( nalphabets, MAX( ll1, ll2 )+2 ); 
		intwork = AllocateIntMtx( nalphabets, MAX( ll1, ll2 )+2 ); 


		amino_dynamicmtx = AllocateDoubleMtx( 0x80, 0x80 );
#if DEBUG
		fprintf( stderr, "succeeded\n" );
#endif

		orlgth1 = ll1 - 100;
		orlgth2 = ll2 - 100;
	}


    for( i=0; i<nalphabets; i++) for( j=0; j<nalphabets; j++ )
		amino_dynamicmtx[(int)amino[i]][(int)amino[j]] = (double)n_dynamicmtx[i][j];




#if 0
	for( i=0; i<lgth1; i++ ) 
		fprintf( stderr, "ogcp1[%d]=%f\n", i, ogcp1[i] );
#endif

	currentw = w1;
	previousw = w2;


	match_calc_mtx( amino_dynamicmtx, initverticalw, seq2, seq1, 0, lgth1 );


	match_calc_mtx( amino_dynamicmtx, currentw, seq1, seq2, 0, lgth2 );

	if( 1 ) // tsuneni outgap-1
	{
		for( i=1; i<lgth1+1; i++ )
		{
			initverticalw[i] += fpenalty;
		}
		for( j=1; j<lgth2+1; j++ )
		{
			currentw[j] += fpenalty;
		}
	}

	for( j=1; j<lgth2+1; ++j ) 
	{
		m[j] = currentw[j-1];
	}

	if( lgth2 == 0 )
		lastverticalw[0] = 0.0;               // lgth2==0 no toki error
	else
		lastverticalw[0] = currentw[lgth2-1]; // lgth2==0 no toki error

	if( 1 ) lasti = lgth1+1; else lasti = lgth1; // tsuneni outgap-1

#if XXXXXXX
fprintf( stderr, "currentw = \n" );
for( i=0; i<lgth1+1; i++ )
{
	fprintf( stderr, "%5.2f ", currentw[i] );
}
fprintf( stderr, "\n" );
fprintf( stderr, "initverticalw = \n" );
for( i=0; i<lgth2+1; i++ )
{
	fprintf( stderr, "%5.2f ", initverticalw[i] );
}
fprintf( stderr, "\n" );
#endif

	for( i=1; i<lasti; i++ )
	{
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

		match_calc_mtx( amino_dynamicmtx, currentw, seq1, seq2, i, lgth2 );
#if XXXXXXX
fprintf( stderr, "\n" );
fprintf( stderr, "i=%d\n", i );
fprintf( stderr, "currentw = \n" );
for( j=0; j<lgth2; j++ )
{
	fprintf( stderr, "%5.2f ", currentw[j] );
}
fprintf( stderr, "\n" );
#endif
#if XXXXXXX
fprintf( stderr, "\n" );
fprintf( stderr, "i=%d\n", i );
fprintf( stderr, "currentw = \n" );
for( j=0; j<lgth2; j++ )
{
	fprintf( stderr, "%5.2f ", currentw[j] );
}
fprintf( stderr, "\n" );
#endif
		currentw[0] = initverticalw[i];

		mi = previousw[0];

		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;
		for( j=1; j<lgth2+1; j++ )
		{
			wm = *prept;

#if 0
			fprintf( stderr, "%5.0f->", wm );
#endif
#if 0
			fprintf( stderr, "%5.0f?", g );
#endif
			if( (g=mi+fpenalty) > wm )
			{
				wm = g;
			}
			if( (g=*prept) >= mi )
			{
				mi = g;
			}
#if USE_PENALTY_EX
			mi += fpenalty_ex;
#endif

#if 0 
			fprintf( stderr, "%5.0f?", g );
#endif
			if( (g=*mjpt + fpenalty) > wm )
			{
				wm = g;
			}
			if( (g=*prept) >= *mjpt )
			{
				*mjpt = g;
			}
#if USE_PENALTY_EX
			m[j] += fpenalty_ex;
#endif

#if 0
			fprintf( stderr, "%5.0f ", wm );
#endif
			*curpt++ += wm;
			mjpt++;
			prept++;
		}
		lastverticalw[i] = currentw[lgth2-1]; // lgth2==0 no toki error
	}

#if 0
	fprintf( stderr, "\n" );
	fprintf( stderr, ">\n%s\n", mseq1[0] );
	fprintf( stderr, ">\n%s\n", mseq2[0] );
	fprintf( stderr, "wm = %f\n", wm );
#endif

	return( wm );
}
