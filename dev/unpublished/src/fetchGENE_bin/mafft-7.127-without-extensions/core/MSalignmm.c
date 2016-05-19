#include "mltaln.h"
#include "dp.h"

#define MEMSAVE 1

#define DEBUG 0
#define USE_PENALTY_EX  0
#define STOREWM 0

#define DPTANNI 100


static TLS int reccycle = 0;


static void match_calc_add( double **scoringmtx, float *match, float **cpmx1, float **cpmx2, int i1, int lgth2, float **floatwork, int **intwork, int initialize )
{
	int j, k, l;
//	float scarr[26];
	float **cpmxpd = floatwork;
	int **cpmxpdn = intwork;
	int count = 0;
	float *matchpt;
	float **cpmxpdpt;
	int **cpmxpdnpt;
	int cpkd;
	float *scarr;
	scarr = calloc( nalphabets, sizeof( float ) );
	if( initialize )
	{
		for( j=0; j<lgth2; j++ )
		{
			count = 0;
			for( l=0; l<nalphabets; l++ )
			{
				if( cpmx2[j][l] )
				{
					cpmxpd[j][count] = cpmx2[j][l];
					cpmxpdn[j][count] = l;
					count++;
				}
			}
			cpmxpdn[j][count] = -1;
		}
	}

	for( l=0; l<nalphabets; l++ )
	{
		scarr[l] = 0.0;
		for( k=0; k<nalphabets; k++ )
		{
//			scarr[l] += n_dis[k][l] * cpmx1[i1][k];
			scarr[l] += scoringmtx[k][l] * cpmx1[i1][k];
		}
	}
#if 0 /* これを使うときはfloatworkのアロケートを逆にする */
	{
		float *fpt, **fptpt, *fpt2;
		int *ipt, **iptpt;
		fpt2 = match;
		iptpt = cpmxpdn;
		fptpt = cpmxpd;
		while( lgth2-- )
		{
			*fpt2 = 0.0;
			ipt=*iptpt,fpt=*fptpt;
			while( *ipt > -1 )
				*fpt2 += scarr[*ipt++] * *fpt++;
			fpt2++,iptpt++,fptpt++;
		} 
	}
	for( j=0; j<lgth2; j++ )
	{
		match[j] = 0.0;
		for( k=0; cpmxpdn[j][k]>-1; k++ )
			match[j] += scarr[cpmxpdn[j][k]] * cpmxpd[j][k];
	} 
#else
	matchpt = match;
	cpmxpdnpt = cpmxpdn;
	cpmxpdpt = cpmxpd;
	while( lgth2-- )
	{
//		*matchpt = 0.0; // add dakara
		for( k=0; (cpkd=(*cpmxpdnpt)[k])>-1; k++ )
			*matchpt += scarr[cpkd] * (*cpmxpdpt)[k];
		matchpt++;
		cpmxpdnpt++;
		cpmxpdpt++;
	}
#endif
	free( scarr );
}

static void match_calc( double **n_dynamicmtx, float *match, float **cpmx1, float **cpmx2, int i1, int lgth2, float **floatwork, int **intwork, int initialize )
{
	int j, k, l;
//	float scarr[26];
	float **cpmxpd = floatwork;
	int **cpmxpdn = intwork;
	int count = 0;
	float *matchpt;
	float **cpmxpdpt;
	int **cpmxpdnpt;
	int cpkd;
	float *scarr;
	scarr = calloc( nalphabets, sizeof( float ) );
	if( initialize )
	{
		for( j=0; j<lgth2; j++ )
		{
			count = 0;
			for( l=0; l<nalphabets; l++ )
			{
				if( cpmx2[j][l] )
				{
					cpmxpd[j][count] = cpmx2[j][l];
					cpmxpdn[j][count] = l;
					count++;
				}
			}
			cpmxpdn[j][count] = -1;
		}
	}

	for( l=0; l<nalphabets; l++ )
	{
		scarr[l] = 0.0;
		for( k=0; k<nalphabets; k++ )
		{
//			scarr[l] += n_dis[k][l] * cpmx1[i1][k];
			scarr[l] += n_dynamicmtx[k][l] * cpmx1[i1][k];
		}
	}
#if 0 /* これを使うときはfloatworkのアロケートを逆にする */
	{
		float *fpt, **fptpt, *fpt2;
		int *ipt, **iptpt;
		fpt2 = match;
		iptpt = cpmxpdn;
		fptpt = cpmxpd;
		while( lgth2-- )
		{
			*fpt2 = 0.0;
			ipt=*iptpt,fpt=*fptpt;
			while( *ipt > -1 )
				*fpt2 += scarr[*ipt++] * *fpt++;
			fpt2++,iptpt++,fptpt++;
		} 
	}
	for( j=0; j<lgth2; j++ )
	{
		match[j] = 0.0;
		for( k=0; cpmxpdn[j][k]>-1; k++ )
			match[j] += scarr[cpmxpdn[j][k]] * cpmxpd[j][k];
	} 
#else
	matchpt = match;
	cpmxpdnpt = cpmxpdn;
	cpmxpdpt = cpmxpd;
	while( lgth2-- )
	{
		*matchpt = 0.0;
		for( k=0; (cpkd=(*cpmxpdnpt)[k])>-1; k++ )
			*matchpt += scarr[cpkd] * (*cpmxpdpt)[k];
		matchpt++;
		cpmxpdnpt++;
		cpmxpdpt++;
	}
#endif
	free( scarr );
}

static float Atracking( float *lasthorizontalw, float *lastverticalw,
						char **seq1, char **seq2, 
                        char **mseq1, char **mseq2, 
                        int **ijp, int icyc, int jcyc,
						int ist, int ien, int jst, int jen, 
						int fulllen1, int fulllen2, int tailgp )
{
	int i, j, l, iin, jin, ifi, jfi, lgth1, lgth2, k, klim;
	char *gaptable1, *gt1bk;
	char *gaptable2, *gt2bk;
	float wm;
	lgth1 = ien-ist+1;
	lgth2 = jen-jst+1;

	gt1bk = AllocateCharVec( lgth1+lgth2+1 );
	gt2bk = AllocateCharVec( lgth1+lgth2+1 );

#if 0
	for( i=0; i<lgth1; i++ ) 
	{
		fprintf( stderr, "lastverticalw[%d] = %f\n", i, lastverticalw[i] );
	}
#endif


//	fprintf( stderr, "in Atracking, tailgp=%d, ien=%d, fulllen1=%d, jen=%d, fulllen2=%d\n", tailgp, ien, fulllen1, jen, fulllen2 );

	if( tailgp == 1 )
		;
	else if( ien == fulllen1-1 || jen == fulllen2-1 )
	{
//		fprintf( stderr, "searching lasthorizontalw and lasthorizontalw\n" );
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
#if 0
	else if( jen == fulllen2-1 )
	{
		fprintf( stderr, "searching lastverticalw\n" );
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
	}
	else if( ien == fulllen1-1 )
	{
		fprintf( stderr, "searching lasthorizontalw\n" );
		wm = lasthorizontalw[0];
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
#endif
 
    for( i=0; i<lgth1+1; i++ ) 
    {
        ijp[i][0] = i + 1;
    }
    for( j=0; j<lgth2+1; j++ ) 
    {
        ijp[0][j] = -( j + 1 );
    }


	gaptable1 = gt1bk + lgth1+lgth2;
	*gaptable1 = 0;
	gaptable2 = gt2bk + lgth1+lgth2;
	*gaptable2 = 0;

//	if( lgth2 == 1 ) fprintf( stderr, "in Atracking, mseq1 = %s, mseq2 = %s\n", mseq1[0], mseq2[0] );

	iin = lgth1; jin = lgth2;
	klim = lgth1+lgth2;
	for( k=0; k<=klim; k++ ) 
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
			*--gaptable1 = 'o';
			*--gaptable2 = '-';
			k++;
		}
		l= jin - jfi;
		while( --l )
		{
			*--gaptable1 = '-';
			*--gaptable2 = 'o';
			k++;
		}
		if( iin <= 0 || jin <= 0 ) break;
		*--gaptable1 = 'o';
		*--gaptable2 = 'o';
		k++;
		iin = ifi; jin = jfi;

	}

	for( i=0; i<icyc; i++ ) gapireru( mseq1[i], seq1[i]+ist, gaptable1 );
	for( j=0; j<jcyc; j++ ) gapireru( mseq2[j], seq2[j]+jst, gaptable2 );

	free( gt1bk );
	free( gt2bk );

//	fprintf( stderr, "in Atracking (owari), mseq1 = %s\n", mseq1[0] );
//	fprintf( stderr, "in Atracking (owari), mseq2 = %s\n", mseq2[0] );
	return( 0.0 );
}

static float MSalignmm_tanni( double **n_dynamicmtx, int icyc, int jcyc, double *eff1, double *eff2, char **seq1, char **seq2, float **cpmx1, float **cpmx2, int ist, int ien, int jst, int jen, int alloclen, int fulllen1, int fulllen2, char **mseq1, char **mseq2, float **gapinfo, int headgp, int tailgp, float headgapfreq1_g, float headgapfreq2_g )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{
//	int k;
	register int i, j;
	int ll1, ll2;
	int lasti, lastj;
	float wm = 0.0;   /* int ?????? */
	float g;
	float *currentw, *previousw;
#if USE_PENALTY_EX
	float fpenalty_ex = (float)penalty_ex;
#endif
#if 1
	float *wtmp;
	int *ijppt;
	float *mjpt, *prept, *curpt;
	int *mpjpt;
#endif
	float mi, *m;
	int **ijp;
	int mpi, *mp;
	float *w1, *w2;
	float *initverticalw;    /* kufuu sureba iranai */
	float *lastverticalw;    /* kufuu sureba iranai */
	int **intwork;
	float **floatwork;
	int **intmtx;
	int lgth1, lgth2;
	float *ogcp1;
	float *fgcp1;
	float *ogcp2;
	float *fgcp2;
	float *gapfreq1f;
	float *gapfreq2f;
//	char **aseq1;
//	char **aseq2;
//	char **aseq1bk, **aseq2bk;
	float headgapfreq1;
	float headgapfreq2;


	ogcp1 = gapinfo[0] + ist;
	fgcp1 = gapinfo[1] + ist;
	ogcp2 = gapinfo[2] + jst;
	fgcp2 = gapinfo[3] + jst;
	gapfreq1f = gapinfo[4] + ist;
	gapfreq2f = gapinfo[5] + jst;

	if( ist > 0 ) headgapfreq1 = gapfreq1f[-1];
	else headgapfreq1 = headgapfreq1_g;
	if( jst > 0 ) headgapfreq2 = gapfreq2f[-1];
	else headgapfreq2 = headgapfreq2_g;

#if STOREWM
	char ttt1[10000], ttt2[10000];
#endif


	lgth1 = ien-ist+1;
	lgth2 = jen-jst+1;

#if STOREWM
	strncpy( ttt1, seq1[0]+ist, lgth1 ); ttt1[lgth1] = 0;
	strncpy( ttt2, seq2[0]+jst, lgth2 ); ttt2[lgth2] = 0;

	fprintf( stderr, "in _tanni ist,ien = %d,%d, lgth1=%d\n", ist, ien, lgth1 );
	fprintf( stderr, "in _tanni jst,jen = %d,%d, lgth2=%d\n", jst, jen, lgth2 );
	fprintf( stderr, "ttt1 = %s\n", ttt1 );
	fprintf( stderr, "ttt2 = %s\n", ttt2 );
#endif

#if 0
	fprintf( stderr, "in _tanni ist,ien = %d,%d, fulllen1=%d\n", ist, ien, fulllen1 );
	fprintf( stderr, "in _tanni jst,jen = %d,%d, fulllen2=%d\n", jst, jen, fulllen2 );
	fprintf( stderr, "in _tanni seq1[0] = %-*.*s\n", ien-ist+1, ien-ist+1, seq1[0]+ist );
	fprintf( stderr, "in _tanni seq2[0] = %-*.*s\n", jen-jst+1, jen-jst+1, seq2[0]+jst );
#endif


	ll1 = ( (int)(lgth1) ) + 100;
	ll2 = ( (int)(lgth2) ) + 100;

//	aseq1 = AllocateCharMtx( icyc, 0 );
//	aseq2 = AllocateCharMtx( jcyc, 0 );
//	aseq1bk = AllocateCharMtx( icyc, lgth1+lgth2+100 );
//	aseq2bk = AllocateCharMtx( jcyc, lgth1+lgth2+100 );
//	for( i=0; i<icyc; i++ ) aseq1[i] = aseq1bk[i];
//	for( i=0; i<jcyc; i++ ) aseq2[i] = aseq2bk[i];

	w1 = AllocateFloatVec( ll2+2 );
	w2 = AllocateFloatVec( ll2+2 );

	initverticalw = AllocateFloatVec( ll1+2 );
	lastverticalw = AllocateFloatVec( ll1+2 );

	m = AllocateFloatVec( ll2+2 );
	mp = AllocateIntVec( ll2+2 );

	floatwork = AllocateFloatMtx( MAX( ll1, ll2 )+2, nalphabets+1 ); 
	intwork = AllocateIntMtx( MAX( ll1, ll2 )+2, nalphabets+1 ); 


	intmtx = AllocateIntMtx( ll1+1, ll2+1 );

	ijp = intmtx;

	currentw = w1;
	previousw = w2;

	match_calc( n_dynamicmtx, initverticalw, cpmx2+jst, cpmx1+ist, 0, lgth1, floatwork, intwork, 1 );

	match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, 0, lgth2, floatwork, intwork, 1 );

	if( headgp || ist != 0 )
	{
		for( i=1; i<lgth1+1; i++ )
		{
			initverticalw[i] += ( ogcp1[0] * headgapfreq2 + fgcp1[i-1] * gapfreq2f[0] ) ;
		}
	}
	if( headgp || jst != 0 )
	{
		for( j=1; j<lgth2+1; j++ )
		{
			currentw[j] += ( ogcp2[0] * headgapfreq1 + fgcp2[j-1] * gapfreq1f[0] ) ;
		}
	}

	for( j=1; j<lgth2+1; ++j ) 
	{
		m[j] = currentw[j-1] + ogcp1[1] * gapfreq2f[j-1]; mp[j] = 0;;
	}

	lastverticalw[0] = currentw[lgth2-1];



	if( tailgp || jen != fulllen2-1 ) lasti = lgth1+1; else lasti = lgth1;
//	if( 1 ) lasti = lgth1+1; else lasti = lgth1;
	for( i=1; i<lasti; i++ )
	{
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

		match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, i, lgth2, floatwork, intwork, 0 );
		currentw[0] = initverticalw[i];

		mi = previousw[0] + ogcp2[1] * gapfreq1f[i-1];
		mpi = 0;

		ijppt = ijp[i] + 1;
		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;
		mpjpt = mp + 1;
//		if( tailgp && jen != fulllen2-1 ) lastj = lgth2+1; else lastj = lgth2;
		lastj = lgth2+1; 
		for( j=1; j<lastj; j++ )
		{
			wm = *prept;
			*ijppt = 0;

#if 0
			fprintf( stderr, "%5.0f->", wm );
#endif
			g = mi + fgcp2[j-1]  * gapfreq1f[i];
#if 0
			fprintf( stderr, "%5.0f?", g );
#endif
			if( g > wm )
			{
				wm = g;
				*ijppt = -( j - mpi );
			}
			g = *prept + ogcp2[j] * gapfreq1f[i-1];
			if( g >= mi )
			{
				mi = g;
				mpi = j-1;
			}
#if USE_PENALTY_EX
			mi += fpenalty_ex;
#endif

			g = *mjpt + fgcp1[i-1] * gapfreq2f[j];
#if 0 
			fprintf( stderr, "%5.0f?", g );
#endif
			if( g > wm )
			{
				wm = g;
				*ijppt = +( i - *mpjpt );
			}
			g = *prept + ogcp1[i] * gapfreq2f[j-1];
			if( g >= *mjpt )
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
			*curpt += wm;


			ijppt++;
			mjpt++;
			prept++;
			mpjpt++;
			curpt++;
		}
		lastverticalw[i] = currentw[lgth2-1];
	}

//	fprintf( stderr, "wm = %f\n", wm );

	Atracking( currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, icyc, jcyc, ist, ien, jst, jen, fulllen1, fulllen2, tailgp );
#if 0
	fprintf( stderr, "res in _tanni mseq1[0] = %s\n", mseq1[0] );
	fprintf( stderr, "res in _tanni mseq2[0] = %s\n", mseq2[0] );
#endif

//	for( i=0; i<icyc; i++ ) strcpy( mseq1[i], aseq1[i] );
//	for( i=0; i<jcyc; i++ ) strcpy( mseq2[i], aseq2[i] );

//	fprintf( stderr, "in _tanni, aseq1 = %s\n", aseq1[0] );
//	fprintf( stderr, "in _tanni, mseq1 = %s\n", mseq1[0] );

	FreeFloatVec( w1 );
	FreeFloatVec( w2 );
	FreeFloatVec( initverticalw );
	FreeFloatVec( lastverticalw );

	FreeFloatVec( m );
	FreeIntVec( mp );


	FreeFloatMtx( floatwork );
	FreeIntMtx( intwork );

	FreeIntMtx( intmtx );


//	FreeCharMtx( aseq1bk );
//	FreeCharMtx( aseq2bk );

//	free( aseq1 );
//	free( aseq2 );

	return( wm );

}

static void freearrays_rec1(
	float *w1, float *w2, float *initverticalw, float *lastverticalw,
	float *midw, float *midm, float *midn,
	int *jumpbacki, int *jumpbackj, int *jumpforwi, int *jumpforwj, int *jumpdummi, int *jumpdummj,
	float *m, int *mp,
	float **floatwork, int **intwork
#if STOREWM
	, float **WMMTX, float **WMMTX2
#endif
)
{

	FreeFloatVec( w1 );
	FreeFloatVec( w2 );
	FreeFloatVec( initverticalw );
	FreeFloatVec( lastverticalw );
	FreeFloatVec( midw );
	FreeFloatVec( midm );
	FreeFloatVec( midn );

	FreeIntVec( jumpbacki );
	FreeIntVec( jumpbackj );
	FreeIntVec( jumpforwi );
	FreeIntVec( jumpforwj );
	FreeIntVec( jumpdummi );
	FreeIntVec( jumpdummj );

	FreeFloatVec( m );
	FreeIntVec( mp );

	FreeFloatMtx( floatwork );
	FreeIntMtx( intwork );

#if STOREWM
	FreeFloatMtx( WMMTX );
	FreeFloatMtx( WMMTX2 );
#endif
}

static void freearrays_rec2( char *gaps, char **aseq1, char **aseq2 )
{
	free( gaps );
#if MEMSAVE
	free( aseq1 );
	free( aseq2 );
#else
	FreeCharMtx( aseq1 );
	FreeCharMtx( aseq2 );
#endif
}

static float MSalignmm_rec( double **n_dynamicmtx, int icyc, int jcyc, double *eff1, double *eff2, char **seq1, char **seq2, float **cpmx1, float **cpmx2, int ist, int ien, int jst, int jen, int alloclen, int fulllen1, int fulllen2, char **mseq1, char **mseq2, int depth, float **gapinfo, int *chudanpt, int chudanref, int *chudanres, int headgp, int tailgp, float headgapfreq1_g, float headgapfreq2_g )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{
//	int k;
	int alnlen;
	float value = 0.0;
	register int i, j;
	char **aseq1, **aseq2;
	int ll1, ll2, l, len;
	int lasti, lastj, imid;
	int jmid = 0; // by D.Mathog, a guess
	float wm = 0.0;   /* int ?????? */
	float g;
	float *currentw, *previousw;
#if USE_PENALTY_EX
	float fpenalty_ex = (float)penalty_ex;
#endif
	float *wtmp;
//	short *ijppt;
	int *mpjpt;
//	short **ijp;
	int *mp;
	int mpi;
	float *mjpt, *prept, *curpt;
	float mi;
	float *m;
	float *w1, *w2;
//	float *match;
	float *initverticalw;    /* kufuu sureba iranai */
	float *lastverticalw;    /* kufuu sureba iranai */
	int **intwork;
	float **floatwork;
//	short **shortmtx;
#if STOREWM
	float **WMMTX;
	float **WMMTX2;
#endif
	float *midw;
	float *midm;
	float *midn;
	int lgth1, lgth2;
	float maxwm;
	int *jumpforwi;
	int *jumpforwj;
	int *jumpbacki;
	int *jumpbackj;
	int *jumpdummi; //muda
	int *jumpdummj = NULL; // by D.Mathog, a guess
	int jumpi, jumpj = 0; // by D.Mathog, a guess
	char *gaps;
	int ijpi, ijpj;
	float *ogcp1;
	float *fgcp1;
	float *ogcp2;
	float *fgcp2;
	float firstm;
	int firstmp;
#if STOREWM
	static TLS char ttt1[50000];
	static TLS char ttt2[50000];
#endif
	float *gapfreq1f;
	float *gapfreq2f;
	float headgapfreq1;
	float headgapfreq2;

#if 0
	int nglen1, nglen2;
	nglen1 = seqlen( seq1[0] );
	nglen2 = seqlen( seq2[0] );
#endif

//	fprintf( stderr, "fulllen1 = %d, fulllen2 = %d, headgp = %d, tailgp = %d\n", fulllen1, fulllen2, headgp, tailgp );

	ogcp1 = gapinfo[0] + ist;
	fgcp1 = gapinfo[1] + ist;
	ogcp2 = gapinfo[2] + jst;
	fgcp2 = gapinfo[3] + jst;
	gapfreq1f = gapinfo[4] + ist;
	gapfreq2f = gapinfo[5] + jst;

	if( ist > 0 ) headgapfreq1 = gapfreq1f[-1];
	else headgapfreq1 = headgapfreq1_g;
	if( jst > 0 ) headgapfreq2 = gapfreq2f[-1];
	else headgapfreq2 = headgapfreq2_g;

	depth++;
	reccycle++;

	lgth1 = ien-ist+1;
	lgth2 = jen-jst+1;

//	if( lgth1 < 5 )
//		fprintf( stderr, "\nWARNING: lgth1 = %d\n", lgth1 );
//	if( lgth2 < 5 )
//		fprintf( stderr, "\nWARNING: lgth2 = %d\n", lgth2 );
//


#if STOREWM
	fprintf( stderr, "==== MSalign (depth=%d, reccycle=%d), ist=%d, ien=%d, jst=%d, jen=%d\n", depth, reccycle, ist, ien, jst, jen );
	strncpy( ttt1, seq1[0]+ist, lgth1 );
	strncpy( ttt2, seq2[0]+jst, lgth2 );
	ttt1[lgth1] = 0;
	ttt2[lgth2] = 0;
	fprintf( stderr, "seq1 = %s\n", ttt1 );
	fprintf( stderr, "seq2 = %s\n", ttt2 );
#endif
	if( lgth2 <= 0 ) // lgth1 <= 0 ha?
	{
//		fprintf( stderr, "\n\n==== jimei\n\n" );
//		exit( 1 );
		for( i=0; i<icyc; i++ ) 
		{
			strncpy( mseq1[i], seq1[i]+ist, lgth1 );
			mseq1[i][lgth1] = 0;
		}
		for( i=0; i<jcyc; i++ ) 
		{
			mseq2[i][0] = 0;
			for( j=0; j<lgth1; j++ )
//				strcat( mseq2[i], "-" );
				strcat( mseq2[i], newgapstr );
		}

//		fprintf( stderr, "==== mseq1[0] (%d) = %s\n", depth, mseq1[0] );
//		fprintf( stderr, "==== mseq2[0] (%d) = %s\n", depth, mseq2[0] );

		return( 0.0 );
	}

#if MEMSAVE
	aseq1 = AllocateCharMtx( icyc, 0 );
	aseq2 = AllocateCharMtx( jcyc, 0 );
	for( i=0; i<icyc; i++ ) aseq1[i] = mseq1[i];
	for( i=0; i<jcyc; i++ ) aseq2[i] = mseq2[i];
#else
	aseq1 = AllocateCharMtx( icyc, lgth1+lgth2+100 );
	aseq2 = AllocateCharMtx( jcyc, lgth1+lgth2+100 );
#endif

//	fprintf( stderr, "####(s) seq1[0] (%d) = \n%-*.*s\n (a%d-a%d)\n", depth, ien-ist+1, ien-ist+1, seq1[0]+ist, ist, ien );
//	fprintf( stderr, "####(s) seq2[0] (%d) = \n%-*.*s\n (b%d-b%d)\n", depth, jen-jst+1, jen-jst+1, seq2[0]+jst, jst, jen );

//  if( lgth1 < DPTANNI && lgth2 < DPTANNI ) // & dato lgth ==1 no kanousei ga arunode yokunai 
//    if( lgth1 < DPTANNI ) // kore mo lgth2 ga mijikasugiru kanousei ari
    if( lgth1 < DPTANNI || lgth2 < DPTANNI ) // zettai ni anzen ka?
	{
//		fprintf( stderr, "==== Going to _tanni\n" );

		value = MSalignmm_tanni( n_dynamicmtx, icyc, jcyc, eff1, eff2, seq1, seq2, cpmx1, cpmx2, ist, ien, jst, jen, alloclen, fulllen1, fulllen2, aseq1, aseq2, gapinfo, headgp, tailgp, headgapfreq1_g, headgapfreq2_g );	


#if MEMSAVE
		free( aseq1 );
		free( aseq2 );
#else
		for( i=0; i<icyc; i++ ) strcpy( mseq1[i], aseq1[i] );
		for( i=0; i<jcyc; i++ ) strcpy( mseq2[i], aseq2[i] );

		FreeCharMtx( aseq1 );
		FreeCharMtx( aseq2 );
#endif

//		fprintf( stderr, "value = %f\n", value );

		return( value );
	}
//	fprintf( stderr, "Trying to divide the mtx\n" );

	ll1 = ( (int)(lgth1) ) + 100;
	ll2 = ( (int)(lgth2) ) + 100;

//	fprintf( stderr, "ll1,ll2=%d,%d\n", ll1, ll2 );

	w1 = AllocateFloatVec( ll2+2 );
	w2 = AllocateFloatVec( ll2+2 );
//	match = AllocateFloatVec( ll2+2 );
	midw = AllocateFloatVec( ll2+2 );
	midn = AllocateFloatVec( ll2+2 );
	midm = AllocateFloatVec( ll2+2 );
	jumpbacki = AllocateIntVec( ll2+2 );
	jumpbackj = AllocateIntVec( ll2+2 );
	jumpforwi = AllocateIntVec( ll2+2 );
	jumpforwj = AllocateIntVec( ll2+2 );
	jumpdummi = AllocateIntVec( ll2+2 );
	jumpdummj = AllocateIntVec( ll2+2 );

	initverticalw = AllocateFloatVec( ll1+2 );
	lastverticalw = AllocateFloatVec( ll1+2 );

	m = AllocateFloatVec( ll2+2 );
	mp = AllocateIntVec( ll2+2 );
	gaps = AllocateCharVec( MAX( ll1, ll2 ) + 2 );

	floatwork = AllocateFloatMtx( MAX( ll1, ll2 )+2, nalphabets ); 
	intwork = AllocateIntMtx( MAX( ll1, ll2 )+2, nalphabets ); 

#if DEBUG
	fprintf( stderr, "succeeded\n" );
#endif

#if STOREWM
	WMMTX = AllocateFloatMtx( ll1, ll2 );
	WMMTX2 = AllocateFloatMtx( ll1, ll2 );
#endif
#if 0
	shortmtx = AllocateShortMtx( ll1, ll2 );

#if DEBUG
	fprintf( stderr, "succeeded\n\n" );
#endif

	ijp = shortmtx;
#endif

	currentw = w1;
	previousw = w2;

	match_calc( n_dynamicmtx, initverticalw, cpmx2+jst, cpmx1+ist, 0, lgth1, floatwork, intwork, 1 );

	match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, 0, lgth2, floatwork, intwork, 1 );

	for( i=1; i<lgth1+1; i++ )
	{
//		initverticalw[i] += ( ogcp1[0] + fgcp1[i-1] ) ;
		initverticalw[i] += ( ogcp1[0] * headgapfreq2 + fgcp1[i-1] * gapfreq2f[0] ) ;
	}
	for( j=1; j<lgth2+1; j++ )
	{
		currentw[j] += ( ogcp2[0] * headgapfreq1 + fgcp2[j-1] * gapfreq1f[0]) ;
	}

#if STOREWM
	WMMTX[0][0] = initverticalw[0];
	for( i=1; i<lgth1+1; i++ )
	{
		WMMTX[i][0] = initverticalw[i];
	}
	for( j=1; j<lgth2+1; j++ )
	{
		WMMTX[0][j] = currentw[j];
	}
#endif


	for( j=1; j<lgth2+1; ++j ) 
	{
		m[j] = currentw[j-1] + ogcp1[1] * gapfreq2f[j-1];
//		m[j] = currentw[j-1];
		mp[j] = 0;
	}

	lastverticalw[0] = currentw[lgth2-1];

	imid = lgth1 * 0.5;

	jumpi = 0; // atode kawaru.
	lasti = lgth1+1;
#if STOREWM
	for( i=1; i<lasti; i++ )
#else
	for( i=1; i<=imid; i++ )
#endif
	{
#ifdef enablemultithread
//		fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
		if( chudanpt && *chudanpt != chudanref ) 
		{
//			fprintf( stderr, "\n\n## CHUUDAN!!! zenhan\n" );
			*chudanres = 1;
			freearrays_rec1
			(
				w1, w2, initverticalw, lastverticalw, midw, midm, midn,
				jumpbacki, jumpbackj, jumpforwi, jumpforwj, jumpdummi, jumpdummj,
				m, mp,
				floatwork, intwork
#if STOREWM
				, WMMTX, WMMTX2
#endif
			);
			freearrays_rec2( gaps, aseq1, aseq2 );
			return( -1.0 );
		}
#endif
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

		match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, i, lgth2, floatwork, intwork, 0 );
		currentw[0] = initverticalw[i];

		m[0] = ogcp1[i];
#if STOREM
		WMMTX2[i][0] = m[0];
#endif
		if( i == imid ) midm[0] = m[0];

		mi = previousw[0] + ogcp2[1] * gapfreq1f[i-1]; 
//		mi = previousw[0];
		mpi = 0;


//		ijppt = ijp[i] + 1;
		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;
		mpjpt = mp + 1;


		lastj = lgth2+1;
		for( j=1; j<lastj; j++ )
		{

			wm = *prept;

#if 0
			fprintf( stderr, "%5.0f->", wm );
#endif
			g = mi + fgcp2[j-1] * gapfreq1f[i];
//			g = mi + fpenalty;
#if 0
			fprintf( stderr, "%5.0f?", g );
#endif
			if( g > wm )
			{
				wm = g;
//				*ijppt = -( j - mpi );
			}
			g = *prept + ogcp2[j] * gapfreq1f[i-1];
//			g = *prept;
			if( g >= mi )
			{
				mi = g;
				mpi = j-1;
			}
#if USE_PENALTY_EX
			mi += fpenalty_ex;
#endif

			g = *mjpt + fgcp1[i-1] * gapfreq2f[j];
//			g = *mjpt + fpenalty;
#if 0 
			fprintf( stderr, "%5.0f?", g );
#endif
			if( g > wm )
			{
				wm = g;
//				*ijppt = +( i - *mpjpt );
			}


			g = *prept + ogcp1[i] * gapfreq2f[j-1];
//			g = *prept;
			if( g >= *mjpt )
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
			*curpt += wm;

#if STOREWM
			WMMTX[i][j] = *curpt;
			WMMTX2[i][j] = *mjpt;
#endif

			if( i == imid ) //muda
			{	
				jumpbackj[j] = *mpjpt; // muda atode matomeru
				jumpbacki[j] = mpi; // muda atode matomeru
//				fprintf( stderr, "jumpbackj[%d] in forward dp is %d\n", j, *mpjpt );
//				fprintf( stderr, "jumpbacki[%d] in forward dp is %d\n", j, mpi );
				midw[j] = *curpt;
				midm[j] = *mjpt;
				midn[j] = mi;
			}

//			fprintf( stderr, "m[%d] = %f\n", j, m[j] );
			mjpt++;
			prept++;
			mpjpt++;
			curpt++;

		}
		lastverticalw[i] = currentw[lgth2-1];

#if STOREWM
		WMMTX2[i][lgth2] = m[lgth2-1];
#endif

#if 0  // ue
		if( i == imid )
		{
			for( j=0; j<lgth2; j++ ) midw[j] = currentw[j];
			for( j=0; j<lgth2; j++ ) midm[j] = m[j];
		}
#endif
	}
//	for( j=0; j<lgth2; j++ ) midw[j] = WMMTX[imid][j];
//	for( j=0; j<lgth2; j++ ) midm[j] = WMMTX2[imid][j];

#if 0
    for( i=0; i<lgth1; i++ )
    {
        for( j=0; j<lgth2; j++ )
        {
            fprintf( stderr, "% 10.2f ", WMMTX[i][j] );
        }
        fprintf( stderr, "\n" );
    }
	fprintf( stderr, "\n" );
	fprintf( stderr, "WMMTX2 = \n" );
    for( i=0; i<lgth1; i++ )
    {
        for( j=0; j<lgth2; j++ )
        {
            fprintf( stderr, "% 10.2f ", WMMTX2[i][j] );
        }
        fprintf( stderr, "\n" );
    }
	fprintf( stderr, "\n" );
#endif

// gyakudp

	match_calc( n_dynamicmtx, initverticalw, cpmx2+jst, cpmx1+ist, lgth2-1, lgth1, floatwork, intwork, 1 );
	match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, lgth1-1, lgth2, floatwork, intwork, 1 );

	for( i=0; i<lgth1-1; i++ )
	{
//		initverticalw[i] += ( fgcp1[lgth1-1] + ogcp1[i+1] );
		initverticalw[i] += ( fgcp1[lgth1-1] * gapfreq2f[lgth2] + ogcp1[i+1] * gapfreq2f[lgth2-1] );
	}
	for( j=0; j<lgth2-1; j++ )
	{
//		currentw[j] += ( fgcp2[lgth2-1] + ogcp2[j+1] );
		currentw[j] += ( fgcp2[lgth2-1] * gapfreq1f[lgth1] + ogcp2[j+1] * gapfreq1f[lgth1-1] );
	}

#if STOREWM
	for( i=0; i<lgth1-1; i++ )
	{
		WMMTX[i][lgth2-1] += ( fgcp1[lgth1-1] + ogcp1[i+1] );
		fprintf( stderr, "fgcp1[lgth1-1] + ogcp1[i+1] = %f\n", fgcp1[lgth1-1] + ogcp1[i+1] );
	}
	for( j=0; j<lgth2-1; j++ )
	{
		WMMTX[lgth1-1][j] += ( fgcp2[lgth2-1] + ogcp2[j+1] );
		fprintf( stderr, "fgcp2[lgth2-1] + ogcp2[j+1] = %f\n", fgcp2[lgth2-1] + ogcp2[j+1] );
	}
#endif






#if 0
	for( j=lgth2-1; j>0; --j )
	{
		m[j-1] = currentw[j] + fgcp2[lgth2-2];
//		m[j-1] = currentw[j];
		mp[j] = lgth1-1;
	}
#else
	for( j=lgth2-1; j>-1; --j )
	{
		m[j] = currentw[j+1] + fgcp1[lgth1-2] * gapfreq2f[j+1];
//		m[j-1] = currentw[j];
		mp[j] = lgth1-1;
	}
#endif

//	for( j=0; j<lgth2; j++ ) m[j] = 0.0;
	// m[lgth2-1] ha irunoka?


//	for( i=lgth1-2; i>=imid; i-- )
	firstm = -9999999.9;
//	firstmp = lgth1-1;
	firstmp = lgth1;
	for( i=lgth1-2; i>-1; i-- )
	{
#ifdef enablemultithread
//		fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
		if( chudanpt && *chudanpt != chudanref ) 
		{
//			fprintf( stderr, "\n\n## CHUUDAN!!! kouhan\n" );
			*chudanres = 1;
			freearrays_rec1
			(
				w1, w2, initverticalw, lastverticalw, midw, midm, midn,
				jumpbacki, jumpbackj, jumpforwi, jumpforwj, jumpdummi, jumpdummj,
				m, mp,
				floatwork, intwork
#if STOREWM
				, WMMTX, WMMTX2
#endif
			);
			freearrays_rec2( gaps, aseq1, aseq2 );
			return( -1.0 );
		}
#endif
		wtmp = previousw;
		previousw = currentw;
		currentw = wtmp;
		previousw[lgth2-1] = initverticalw[i+1];
//		match_calc( currentw, seq1, seq2, i, lgth2 );
		match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, i, lgth2, floatwork, intwork, 0 );

		currentw[lgth2-1] = initverticalw[i];

//		m[lgth2] = fgcp1[i];
//		WMMTX2[i][lgth2] += m[lgth2];
//		fprintf( stderr, "m[] = %f\n", m[lgth2] );

		mi = previousw[lgth2-1] + fgcp2[lgth2-2] * gapfreq1f[i+1];
//		mi = previousw[lgth2-1];
		mpi = lgth2 - 1;

		mjpt = m + lgth2 - 2;
		prept = previousw + lgth2 - 1;
		curpt = currentw + lgth2 - 2;
		mpjpt = mp + lgth2 - 2;


		for( j=lgth2-2; j>-1; j-- )
		{
			wm = *prept;
			ijpi = i+1;
			ijpj = j+1;

			g = mi + ogcp2[j+1] * gapfreq1f[i];
//			g = mi + fpenalty;
			if( g > wm )
			{
				wm = g;
				ijpj = mpi;
				ijpi = i+1;
			}

			g = *prept + fgcp2[j] * gapfreq1f[i+1];
//			g = *prept;
			if( g >= mi )
			{
//				fprintf( stderr, "i,j=%d,%d - renewed! mpi = %d\n", i, j, j+1 );
				mi = g;
				mpi = j + 1;
			}

#if USE_PENALTY_EX
			mi += fpenalty_ex;
#endif

//			fprintf( stderr, "i,j=%d,%d *mpjpt = %d\n", i, j, *mpjpt );
			g = *mjpt + ogcp1[i+1] * gapfreq2f[j];
//			g = *mjpt + fpenalty;
			if( g > wm )
			{
				wm = g;
				ijpi = *mpjpt;
				ijpj = j+1;
			}

//			if( i == imid )fprintf( stderr, "i,j=%d,%d \n", i, j );
			g = *prept + fgcp1[i] * gapfreq2f[j+1];
//			g = *prept;
			if( g >= *mjpt )
			{
				*mjpt = g;
				*mpjpt = i + 1;
			}

#if USE_PENALTY_EX
			m[j] += fpenalty_ex;
#endif

			if( i == jumpi || i == imid - 1 )
			{
				jumpforwi[j] = ijpi; //muda
				jumpforwj[j] = ijpj; //muda
//				fprintf( stderr, "jumpfori[%d] = %d\n", j, ijpi );
//				fprintf( stderr, "jumpforj[%d] = %d\n", j, ijpj );
			}
			if( i == imid ) // muda
			{
				midw[j] += wm;
//				midm[j+1] += *mjpt + fpenalty; //??????
				midm[j+1] += *mjpt; //??????
			}
			if( i == imid - 1 )
			{
//				midn[j] += mi + fpenalty;  //????
				midn[j] += mi;  //????
			}
#if STOREWM
			WMMTX[i][j] += wm;
//			WMMTX2[i][j+1] += *mjpt + fpenalty;
			WMMTX2[i][j+1] += *mjpt;
#endif
			*curpt += wm;

			mjpt--;
			prept--;
			mpjpt--;
			curpt--;
		}
//		fprintf( stderr, "adding *mjpt (=%f) to WMMTX2[%d][%d]\n", *mjpt, i, j+1 );
		g = *prept + fgcp1[i];
		if( firstm < g ) 
		{
			firstm = g;
			firstmp = i + 1;
		}
#if STOREWM
		WMMTX2[i][j+1] += firstm;
#endif
		if( i == imid ) midm[j+1] += firstm;


		if( i == imid - 1 )	
		{
			maxwm = midw[1];
			jmid = 0;
//			if( depth == 1 ) fprintf( stderr, "maxwm!! = %f\n", maxwm );
			for( j=2; j<lgth2-1; j++ )
			{
				wm = midw[j];
				if( wm > maxwm )
				{
					jmid = j;
					maxwm = wm;
				}
//				if( depth == 1 ) fprintf( stderr, "maxwm!! = %f\n", maxwm );
			}
			for( j=0; j<lgth2+1; j++ )
			{
				wm = midm[j];
				if( wm > maxwm )
				{
					jmid = j;
					maxwm = wm;
				}
//				if( depth == 1 ) fprintf( stderr, "maxwm!! = %f\n", maxwm );
			}

//			if( depth == 1 ) fprintf( stderr, "maxwm!! = %f\n", maxwm );


//			fprintf( stderr, "### imid=%d, jmid=%d\n", imid, jmid );
			wm = midw[jmid];
			jumpi = imid-1;
			jumpj = jmid-1;
			if( jmid > 0 && midn[jmid-1] > wm ) //060413
			{
				jumpi = imid-1;
				jumpj = jumpbacki[jmid];
				wm = midn[jmid-1];
//				fprintf( stderr, "rejump (n)\n" );
			}
			if( midm[jmid] > wm )
			{
				jumpi = jumpbackj[jmid];
				jumpj = jmid-1;
				wm = midm[jmid];
//				fprintf( stderr, "rejump (m) jumpi=%d\n", jumpi );
			}


//			fprintf( stderr, "--> imid=%d, jmid=%d\n", imid, jmid );
//			fprintf( stderr, "--> jumpi=%d, jumpj=%d\n", jumpi, jumpj );
#if STOREWM
			fprintf( stderr, "imid = %d\n", imid );
			fprintf( stderr, "midn = \n" );
			for( j=0; j<lgth2; j++ )
			{
				fprintf( stderr, "% 7.1f ", midn[j] );
			}
			fprintf( stderr, "\n" );
			fprintf( stderr, "midw = \n" );
			for( j=0; j<lgth2; j++ )
			{
				fprintf( stderr, "% 7.1f ", midw[j] );
			}
			fprintf( stderr, "\n" );
			fprintf( stderr, "midm = \n" );
			for( j=0; j<lgth2; j++ )
			{
				fprintf( stderr, "% 7.1f ", midm[j] );
			}
			fprintf( stderr, "\n" );
#endif
//			fprintf( stderr, "maxwm = %f\n", maxwm );
		}
		if( i == jumpi ) //saki?
		{
//			fprintf( stderr, "#### FIRST i=%d, jumpi<imid=%d<%d, ist=%d, ien=%d, firstmp=%d, lgth1=%d\n", i, jumpi, imid, ist, ien, firstmp, lgth1 );
//			fprintf( stderr, "#### mark 1\n" );
//			fprintf( stderr, "imid, jumpi = %d,%d\n", imid, jumpi );
//			fprintf( stderr, "jmid, jumpj = %d,%d\n", jmid, jumpj );

			if( jmid == 0 )
			{
//				fprintf( stderr, "#### CHUI2!\n" );
				jumpj = 0; jmid = 1;
#if 0 // v6.823 made
				jumpi = firstmp-1;
				imid = firstmp;
#endif
#if 0
				jumpi = 0;
				imid = 1;
#else
//				if( 1 || firstmp > 100 ) // naze 100
				if( imid < firstmp-1 ) // naze 100
				{
					jumpi = firstmp;
					imid = firstmp+1;
				}
#if 0
				else
				{
					jumpi = 0;
					imid = 1;
				}
#endif
#endif
			}

#if 0
			else if( jmid == lgth2 ) 
			{
				fprintf( stderr, "CHUI1!\n" );
				jumpi=0; jumpj=0;
				imid=jumpforwi[0]; jmid=lgth2-1;
			}
#else // 060414
			else if( jmid >= lgth2 ) 
			{
//				fprintf( stderr, "CHUI1!\n" );
				jumpi=imid-1; jmid=lgth2;
				jumpj = lgth2-1;
			}
#endif
			else
			{
//				fprintf( stderr, "#### CHUI3!\n" );
				imid = jumpforwi[jumpj];
				jmid = jumpforwj[jumpj];
				if( imid == jumpi ) jumpi = imid-1;
			}
#if 0
			fprintf( stderr, "jumpi -> %d\n", jumpi );
			fprintf( stderr, "jumpj -> %d\n", jumpj );
			fprintf( stderr, "imid -> %d\n", imid );
			fprintf( stderr, "jmid -> %d\n", jmid );
#endif
//			fprintf( stderr, "#### FINAL i=%d, jumpi<imid=%d<%d, ist=%d, ien=%d, firstmp=%d\n", i, jumpi, imid, ist, ien, firstmp );

#if STOREWM
			break;
#else
			break;
#endif
		}
	}

//	fprintf( stderr, "imid = %d, but jumpi = %d\n", imid, jumpi );
//	fprintf( stderr, "jmid = %d, but jumpj = %d\n", jmid, jumpj );

//	for( j=0; j<lgth2; j++ ) midw[j] += currentw[j];
//	for( j=0; j<lgth2; j++ ) midm[j] += m[j+1];
//	for( j=0; j<lgth2; j++ ) midw[j] += WMMTX[imid][j];
//	for( j=0; j<lgth2; j++ ) midw[j] += WMMTX[imid][j];


#if STOREWM
	fprintf( stderr, "WMMTX = \n" );
    for( i=0; i<lgth1; i++ )
    {
        fprintf( stderr, "%d ", i );
        for( j=0; j<lgth2; j++ )
        {
            fprintf( stderr, "% 7.2f ", WMMTX[i][j] );
        }
        fprintf( stderr, "\n" );
    }
//	fprintf( stderr, "WMMTX2 = (p = %f)\n", fpenalty );
    for( i=0; i<lgth1; i++ )
    {
        fprintf( stderr, "%d ", i );
        for( j=0; j<lgth2+1; j++ )
        {
            fprintf( stderr, "% 7.2f ", WMMTX2[i][j] );
        }
        fprintf( stderr, "\n" );
    }

	fprintf( stderr, "jumpbacki = \n" );
	for( j=0; j<lgth2; j++ )
	{
		fprintf( stderr, "% 7d ", jumpbacki[j] );
	}
	fprintf( stderr, "\n" );
	fprintf( stderr, "jumpbackj = \n" );
	for( j=0; j<lgth2; j++ )
	{
		fprintf( stderr, "% 7d ", jumpbackj[j] );
	}
	fprintf( stderr, "\n" );
	fprintf( stderr, "jumpforwi = \n" );
	for( j=0; j<lgth2; j++ )
	{
		fprintf( stderr, "% 7d ", jumpforwi[j] );
	}
	fprintf( stderr, "\n" );
	fprintf( stderr, "jumpforwj = \n" );
	for( j=0; j<lgth2; j++ )
	{
		fprintf( stderr, "% 7d ", jumpforwj[j] );
	}
	fprintf( stderr, "\n" );


#endif


//	Atracking( currentw, lastverticalw, seq1, seq2, mseq1, mseq2, cpmx1, cpmx2, ijp, icyc, jcyc );

#if 0 // irukamo
	resultlen = strlen( mseq1[0] );
	if( alloclen < resultlen || resultlen > N )
	{
		fprintf( stderr, "alloclen=%d, resultlen=%d, N=%d\n", alloclen, resultlen, N );
		ErrorExit( "LENGTH OVER!\n" );
	}
#endif



#if 0
	fprintf( stderr, "jumpi = %d, imid = %d\n", jumpi, imid );
	fprintf( stderr, "jumpj = %d, jmid = %d\n", jumpj, jmid );

	fprintf( stderr, "imid = %d\n", imid );
	fprintf( stderr, "jmid = %d\n", jmid );
#endif

	freearrays_rec1
	(
		w1, w2, initverticalw, lastverticalw, midw, midm, midn,
		jumpbacki, jumpbackj, jumpforwi, jumpforwj, jumpdummi, jumpdummj,
		m, mp,
		floatwork, intwork
#if STOREWM
		, WMMTX, WMMTX2
#endif
	);


//	fprintf( stderr, "==== calling myself (first)\n" );

	value = MSalignmm_rec( n_dynamicmtx, icyc, jcyc, eff1, eff2, seq1, seq2, cpmx1, cpmx2, ist, ist+jumpi, jst, jst+jumpj, alloclen, fulllen1, fulllen2, aseq1, aseq2, depth, gapinfo, NULL, 0, NULL, headgp, tailgp, headgapfreq1_g, headgapfreq2_g ); // chudan mada	
#if 0
		fprintf( stderr, "aseq1[0] = %s\n", aseq1[0] );
		fprintf( stderr, "aseq2[0] = %s\n", aseq2[0] );
#endif
#if MEMSAVE
#else
	for( i=0; i<icyc; i++ ) strcpy( mseq1[i], aseq1[i] );
	for( i=0; i<jcyc; i++ ) strcpy( mseq2[i], aseq2[i] );
#endif

//	fprintf( stderr, "====(f) aseq1[0] (%d) = %s (%d-%d)\n", depth, aseq1[0], ist, ien );
//	fprintf( stderr, "====(f) aseq2[0] (%d) = %s (%d-%d)\n", depth, aseq2[0], jst, jen );

	len = strlen( mseq1[0] );
//	fprintf( stderr, "len = %d\n", len );
	l = jmid - jumpj - 1;
//	fprintf( stderr, "l=%d\n", l );
	if( l > 0 )
	{
//		for( i=0; i<l; i++ ) gaps[i] = '-'; gaps[i] = 0;
		for( i=0; i<l; i++ ) gaps[i] = *newgapstr; gaps[i] = 0;
		for( i=0; i<icyc; i++ ) 
		{
			strcat( mseq1[i], gaps );
			mseq1[i][len+l] = 0;
		}
		for( j=0; j<jcyc; j++ )
		{
			strncat( mseq2[j], seq2[j]+jst+jumpj+1, l );
			mseq2[j][len+l] = 0;
		}
//		fprintf( stderr, "penalizing (2) .. %f(%d), %f(%d)\n", ogcp2[jumpj+1], jumpj+1, fgcp2[jmid-1], jmid-1 );
		value +=  ( ogcp2[jumpj+1] + fgcp2[jmid-1] );
//		value += fpenalty;
	}
	len = strlen( mseq1[0] );
	l = imid - jumpi - 1;
//	fprintf( stderr, "l=%d\n", l );
	if( l > 0 )
	{
//		for( i=0; i<l; i++ ) gaps[i] = '-'; gaps[i] = 0;
		for( i=0; i<l; i++ ) gaps[i] = *newgapstr; gaps[i] = 0;
		for( i=0; i<icyc; i++ )
		{
			strncat( mseq1[i], seq1[i]+ist+jumpi+1, l );
			mseq1[i][len+l] = 0;
		}
		for( j=0; j<jcyc; j++ ) 
		{
			strcat( mseq2[j], gaps );
			mseq2[j][len+l] = 0;
		}

//		for( i=0; i<lgth1; i++ ) fprintf( stderr, "ogcp1[%d] = %f\n", i, ogcp1[i] );
//		for( i=0; i<lgth1; i++ ) fprintf( stderr, "fgcp1[%d] = %f\n", i, fgcp1[i] );


//		fprintf( stderr, "penalizing (1) .. ogcp1[%d] = %f, fgcp1[%d] = %f\n", jumpi+1, ogcp1[jumpi+1], imid-1, fgcp1[imid-1] );
		value += ( ogcp1[jumpi+1] + fgcp1[imid-1] );
//		value += fpenalty;
	}
#if 0
	for( i=0; i<icyc; i++ ) fprintf( stderr, "after gapfill mseq1[%d]=%s\n", i, mseq1[i] );
	for( i=0; i<jcyc; i++ ) fprintf( stderr, "after gapfill mseq2[%d]=%s\n", i, mseq2[i] );
#endif

//	fprintf( stderr, "==== calling myself (second)\n" );

#if MEMSAVE
	alnlen = strlen( aseq1[0] );
	for( i=0; i<icyc; i++ ) aseq1[i] += alnlen;
	for( i=0; i<jcyc; i++ ) aseq2[i] += alnlen;
#endif

	value += MSalignmm_rec( n_dynamicmtx, icyc, jcyc, eff1, eff2, seq1, seq2, cpmx1, cpmx2, ist+imid, ien, jst+jmid, jen, alloclen, fulllen1, fulllen2, aseq1, aseq2, depth, gapinfo, NULL, 0, NULL, headgp, tailgp, headgapfreq1_g, headgapfreq2_g ); // chudan mada
#if 0
		fprintf( stderr, "aseq1[0] = %s\n", aseq1[0] );
		fprintf( stderr, "aseq2[0] = %s\n", aseq2[0] );
#endif



#if DEBUG
	if( value - maxwm > 1 || maxwm - value > 1 )
	{
		fprintf( stderr, "WARNING value  = %f, but maxwm = %f\n", value, maxwm );
		for( i=0; i<icyc; i++ )
		{
			fprintf( stderr, ">1-%d\n%s\n", i, mseq1[i] );
			fprintf( stderr, "%s\n", aseq1[i] );
		}
		for( i=0; i<jcyc; i++ )
		{
			fprintf( stderr, ">2-%d\n%s\n", i, mseq2[i] );
			fprintf( stderr, "%s\n", aseq2[i] );
		}

//		exit( 1 );
	}
	else
	{
		fprintf( stderr, "value = %.0f, maxwm = %.0f -> ok\n", value, maxwm );
	}
#endif

#if MEMSAVE
#else
	for( i=0; i<icyc; i++ ) strcat( mseq1[i], aseq1[i] );
	for( i=0; i<jcyc; i++ ) strcat( mseq2[i], aseq2[i] );
#endif

//	fprintf( stderr, "====(s) aseq1[0] (%d) = \n%s\n (a%d-a%d)\n", depth, mseq1[0], ist, ien );
//	fprintf( stderr, "====(s) aseq2[0] (%d) = \n%s\n (b%d-b%d)\n", depth, mseq2[0], jst, jen );

	freearrays_rec2( gaps, aseq1, aseq2 );

#if 0
	if( seqlen( seq1[0] ) != nglen1 )
	{
		fprintf( stderr, "bug! hairetsu ga kowareta! (nglen1) seqlen(seq1[0])=%d but nglen1=%d\n", seqlen( seq1[0] ), nglen1 );
		fprintf( stderr, "seq1[0] = %s\n", seq1[0] );
		exit( 1 );
	}
	else
		fprintf( stderr, "nglen1 is ok in _rec\n" );
	if( seqlen( seq2[0] ) != nglen2 )
	{
		fprintf( stderr, "bug! hairetsu ga kowareta! (nglen2) seqlen(seq2[0])=%d but nglen2=%d\n", seqlen( seq2[0] ), nglen2 );
		exit( 1 );
	}
	else
		fprintf( stderr, "nglen2 is ok in _rec\n" );
#endif

	return( value );
}

static void freearrays( 
	float *ogcp1, 
	float *ogcp2, 
	float *fgcp1, 
	float *fgcp2, 
	float **cpmx1,
	float **cpmx2, 
	float *gapfreq1f,
	float *gapfreq2f,
	float **gapinfo, 
	char **mseq1, 
	char **mseq2 
)
{
	FreeFloatVec( ogcp1 );
	FreeFloatVec( ogcp2 );
	FreeFloatVec( fgcp1 );
	FreeFloatVec( fgcp2 );
	FreeFloatMtx( cpmx1 );
	FreeFloatMtx( cpmx2 );
	FreeFloatVec( gapfreq1f );
	FreeFloatVec( gapfreq2f );
	free( (void *)gapinfo );

	FreeCharMtx( mseq1 );
	FreeCharMtx( mseq2 );
}

float MSalignmm( double **n_dynamicmtx, char **seq1, char **seq2, double *eff1, double *eff2, int icyc, int jcyc, int alloclen, char *sgap1, char *sgap2, char *egap1, char *egap2, int *chudanpt, int chudanref, int *chudanres, int headgp, int tailgp )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{
//	int k;
	int i, j;
	int ll1, ll2;
	int lgth1, lgth2;
	float wm = 0.0;   /* int ?????? */
	char **mseq1;
	char **mseq2;
	float *ogcp1;
	float *ogcp2;
	float *fgcp1;
	float *fgcp2;
	float **cpmx1;
	float **cpmx2;
	float **gapinfo;
	float fpenalty = (float)penalty;
	float *gapfreq1f;
	float *gapfreq2f;
	int nglen1, nglen2;
	float headgapfreq1;
	float headgapfreq2;

#if 0
	fprintf( stderr, "eff in SA+++align\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "eff1[%d] = %f\n", i, eff1[i] );
#endif

	nglen1 = seqlen( seq1[0] );
	nglen2 = seqlen( seq2[0] );

#if 0
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "seq1[%d] at root = %s\n", i, seq1[i] );
	for( j=0; j<jcyc; j++ ) fprintf( stderr, "seq2[%d] at root = %s\n", j, seq2[j] );
	fprintf( stderr, "\n" );
#endif

	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );

	ll1 = ( (int)(lgth1) ) + 100;
	ll2 = ( (int)(lgth2) ) + 100;

	mseq1 = AllocateCharMtx( icyc, ll1+ll2 );
	mseq2 = AllocateCharMtx( jcyc, ll1+ll2 );

	gapinfo = AllocateFloatMtx( 6, 0 );
	ogcp1 = AllocateFloatVec( ll1+2 );
	ogcp2 = AllocateFloatVec( ll2+2 );
	fgcp1 = AllocateFloatVec( ll1+2 );
	fgcp2 = AllocateFloatVec( ll2+2 );


	cpmx1 = AllocateFloatMtx( ll1+2, nalphabets+1 );
	cpmx2 = AllocateFloatMtx( ll2+2, nalphabets+1 );

	gapfreq1f = AllocateFloatVec( ll1+2 ); // must be filled with 0.0
	gapfreq2f = AllocateFloatVec( ll2+2 ); // must be filled with 0.0

	for( i=0; i<icyc; i++ ) 
	{
		if( strlen( seq1[i] ) != lgth1 )
		{
			fprintf( stderr, "i = %d / %d\n", i, icyc );
			fprintf( stderr, "bug! hairetsu ga kowareta!\n" );
			exit( 1 );
		}
	}
	for( j=0; j<jcyc; j++ )
	{
		if( strlen( seq2[j] ) != lgth2 )
		{
			fprintf( stderr, "j = %d / %d\n", j, jcyc );
			fprintf( stderr, "bug! hairetsu ga kowareta!\n" );
			exit( 1 );
		}
	}

	MScpmx_calc_new( seq1, cpmx1, eff1, lgth1, icyc );
	MScpmx_calc_new( seq2, cpmx2, eff2, lgth2, jcyc );


#if 1

	if( sgap1 )
	{
		new_OpeningGapCount( ogcp1, icyc, seq1, eff1, lgth1, sgap1 );
		new_OpeningGapCount( ogcp2, jcyc, seq2, eff2, lgth2, sgap2 );
		new_FinalGapCount( fgcp1, icyc, seq1, eff1, lgth1, egap2 );
		new_FinalGapCount( fgcp2, jcyc, seq2, eff2, lgth2, egap2 );
		outgapcount( &headgapfreq1, icyc, sgap1, eff1 );
		outgapcount( &headgapfreq2, jcyc, sgap2, eff2 );
		outgapcount( gapfreq1f+lgth1, icyc, egap1, eff1 );
		outgapcount( gapfreq2f+lgth2, jcyc, egap2, eff2 );
	}
	else
	{
		st_OpeningGapCount( ogcp1, icyc, seq1, eff1, lgth1 );
		st_OpeningGapCount( ogcp2, jcyc, seq2, eff2, lgth2 );
		st_FinalGapCount( fgcp1, icyc, seq1, eff1, lgth1 );
		st_FinalGapCount( fgcp2, jcyc, seq2, eff2, lgth2 );
		headgapfreq1 = 0.0;
		headgapfreq2 = 0.0;
		gapfreq1f[lgth1] = 0.0;
		gapfreq2f[lgth2] = 0.0;
	}

	if( legacygapcost == 0 )
	{
		gapcountf( gapfreq1f, seq1, icyc, eff1, lgth1 );
		gapcountf( gapfreq2f, seq2, jcyc, eff2, lgth2 );
		for( i=0; i<lgth1+1; i++ ) gapfreq1f[i] = 1.0 - gapfreq1f[i];
		for( i=0; i<lgth2+1; i++ ) gapfreq2f[i] = 1.0 - gapfreq2f[i];
		headgapfreq1 = 1.0 - headgapfreq1;
		headgapfreq2 = 1.0 - headgapfreq2;
	}
	else
	{
		for( i=0; i<lgth1+1; i++ ) gapfreq1f[i] = 1.0;
		for( i=0; i<lgth2+1; i++ ) gapfreq2f[i] = 1.0;
		headgapfreq1 = 1.0;
		headgapfreq2 = 1.0;
	}

#if 1
	for( i=0; i<lgth1; i++ ) 
	{
		ogcp1[i] = 0.5 * ( 1.0 - ogcp1[i] ) * fpenalty * ( gapfreq1f[i] );
		fgcp1[i] = 0.5 * ( 1.0 - fgcp1[i] ) * fpenalty * ( gapfreq1f[i] );
//		fprintf( stderr, "fgcp1[%d] = %f\n", i, fgcp1[i] );
	}
	for( i=0; i<lgth2; i++ ) 
	{
		ogcp2[i] = 0.5 * ( 1.0 - ogcp2[i] ) * fpenalty * ( gapfreq2f[i] );
		fgcp2[i] = 0.5 * ( 1.0 - fgcp2[i] ) * fpenalty * ( gapfreq2f[i] );
//		fprintf( stderr, "fgcp2[%d] = %f\n", i, fgcp2[i] );
	}
#else
	for( i=0; i<lgth1; i++ ) 
	{
		ogcp1[i] = 0.5 * fpenalty;
		fgcp1[i] = 0.5 * fpenalty;
	}
	for( i=0; i<lgth2; i++ ) 
	{
		ogcp2[i] = 0.5 * fpenalty;
		fgcp2[i] = 0.5 * fpenalty;
	}
#endif

	gapinfo[0] = ogcp1;
	gapinfo[1] = fgcp1;
	gapinfo[2] = ogcp2;
	gapinfo[3] = fgcp2;
	gapinfo[4] = gapfreq1f;
	gapinfo[5] = gapfreq2f;
#endif

#if 0
	fprintf( stdout, "in MSalignmm.c\n" );
	for( i=0; i<icyc; i++ )
	{
		fprintf( stdout, ">%d of GROUP1\n", i );
		fprintf( stdout, "%s\n", seq1[i] );
	}
	for( i=0; i<jcyc; i++ )
	{
		fprintf( stdout, ">%d of GROUP2\n", i );
		fprintf( stdout, "%s\n", seq2[i] );
	}
	fflush( stdout );
#endif

	wm = MSalignmm_rec( n_dynamicmtx, icyc, jcyc, eff1, eff2, seq1, seq2, cpmx1, cpmx2, 0, lgth1-1, 0, lgth2-1, alloclen, lgth1, lgth2, mseq1, mseq2, 0, gapinfo, chudanpt, chudanref, chudanres, headgp, tailgp, headgapfreq1, headgapfreq2 );
#ifdef enablemultithread
	if( chudanres && *chudanres ) 
	{
//		fprintf( stderr, "\n\n## CHUUDAN!!! relay\n" );
		*chudanres = 1;
		freearrays( ogcp1, ogcp2, fgcp1, fgcp2, cpmx1, cpmx2, gapfreq1f, gapfreq2f, gapinfo, mseq1, mseq2 );
		return( -1.0 );
	}
#endif

#if 0
		fprintf( stderr, "\n" );
		fprintf( stderr, " seq1[0] = %s\n", seq1[0] );
		fprintf( stderr, " seq2[0] = %s\n", seq2[0] );
		fprintf( stderr, "mseq1[0] = %s\n", mseq1[0] );
		fprintf( stderr, "mseq2[0] = %s\n", mseq2[0] );
		fprintf( stderr, "\n" );
#endif

//	fprintf( stderr, "wm = %f\n", wm );


	for( i=0; i<icyc; i++ ) strcpy( seq1[i], mseq1[i] );
	for( i=0; i<jcyc; i++ ) strcpy( seq2[i], mseq2[i] );

	if( seqlen( seq1[0] ) != nglen1 )
	{
		fprintf( stderr, "bug! hairetsu ga kowareta! (nglen1) seqlen(seq1[0])=%d but nglen1=%d\n", seqlen( seq1[0] ), nglen1 );
		fprintf( stderr, "seq1[0] = %s\n", seq1[0] );
		exit( 1 );
	}
	if( seqlen( seq2[0] ) != nglen2 )
	{
		fprintf( stderr, "bug! hairetsu ga kowareta! (nglen2) seqlen(seq2[0])=%d but nglen2=%d\n", seqlen( seq2[0] ), nglen2 );
		exit( 1 );
	}


	freearrays( ogcp1, ogcp2, fgcp1, fgcp2, cpmx1, cpmx2, gapfreq1f, gapfreq2f, gapinfo, mseq1, mseq2 );

	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );
	for( i=0; i<icyc; i++ ) 
	{
		if( strlen( seq1[i] ) != lgth1 )
		{
			fprintf( stderr, "i = %d / %d\n", i, icyc );
			fprintf( stderr, "hairetsu ga kowareta (end of MSalignmm) !\n" );
			exit( 1 );
		}
	}
	for( j=0; j<jcyc; j++ )
	{
		if( strlen( seq2[j] ) != lgth2 )
		{
			fprintf( stderr, "j = %d / %d\n", j, jcyc );
			fprintf( stderr, "hairetsu ga kowareta (end of MSalignmm) !\n" );
			exit( 1 );
		}
	}

#if 0
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, " seq1[i] = %s\n", seq1[i] );
	for( j=0; j<jcyc; j++ ) fprintf( stderr, " seq2[j] = %s\n", seq2[j] );
	fprintf( stderr, "\n" );
#endif

	return( wm );
}

// -------------------------------
//    variousdist
// -------------------------------

static void fillzero( float *s, int l )
{
	while( l-- ) *s++ = 0.0; 
}

static float MSalignmm_tanni_variousdist( double ***matrices, int icyc, int jcyc, char **seq1, char **seq2, float ***cpmx1s, float ***cpmx2s, int ist, int ien, int jst, int jen, int alloclen, int fulllen1, int fulllen2, char **mseq1, char **mseq2, float **gapinfo, int headgp, int tailgp, float headgapfreq1_g, float headgapfreq2_g )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{
//	int k;
	register int i, j, c;
	int ll1, ll2;
	int lasti, lastj;
	float wm = 0.0;   /* int ?????? */
	float g;
	float *currentw, *previousw;
#if USE_PENALTY_EX
	float fpenalty_ex = (float)penalty_ex;
#endif
#if 1
	float *wtmp;
	int *ijppt;
	float *mjpt, *prept, *curpt;
	int *mpjpt;
#endif
	float mi, *m;
	int **ijp;
	int mpi, *mp;
	float *w1, *w2;
	float *initverticalw;    /* kufuu sureba iranai */
	float *lastverticalw;    /* kufuu sureba iranai */
	int ***intwork;
	float ***floatwork;
	int **intmtx;
	int lgth1, lgth2;
	float *ogcp1;
	float *fgcp1;
	float *ogcp2;
	float *fgcp2;
	float *gapfreq1f;
	float *gapfreq2f;
//	char **aseq1;
//	char **aseq2;
//	char **aseq1bk, **aseq2bk;
	float headgapfreq1;
	float headgapfreq2;


	ogcp1 = gapinfo[0] + ist;
	fgcp1 = gapinfo[1] + ist;
	ogcp2 = gapinfo[2] + jst;
	fgcp2 = gapinfo[3] + jst;
	gapfreq1f = gapinfo[4] + ist;
	gapfreq2f = gapinfo[5] + jst;

	if( ist > 0 ) headgapfreq1 = gapfreq1f[-1];
	else headgapfreq1 = headgapfreq1_g;
	if( jst > 0 ) headgapfreq2 = gapfreq2f[-1];
	else headgapfreq2 = headgapfreq2_g;

#if STOREWM
	char ttt1[10000], ttt2[10000];
#endif


	lgth1 = ien-ist+1;
	lgth2 = jen-jst+1;

#if STOREWM
	strncpy( ttt1, seq1[0]+ist, lgth1 ); ttt1[lgth1] = 0;
	strncpy( ttt2, seq2[0]+jst, lgth2 ); ttt2[lgth2] = 0;

	fprintf( stderr, "in _tanni ist,ien = %d,%d, lgth1=%d\n", ist, ien, lgth1 );
	fprintf( stderr, "in _tanni jst,jen = %d,%d, lgth2=%d\n", jst, jen, lgth2 );
	fprintf( stderr, "ttt1 = %s\n", ttt1 );
	fprintf( stderr, "ttt2 = %s\n", ttt2 );
#endif

#if 0
	fprintf( stderr, "in _tanni ist,ien = %d,%d, fulllen1=%d\n", ist, ien, fulllen1 );
	fprintf( stderr, "in _tanni jst,jen = %d,%d, fulllen2=%d\n", jst, jen, fulllen2 );
	fprintf( stderr, "in _tanni seq1[0] = %-*.*s\n", ien-ist+1, ien-ist+1, seq1[0]+ist );
	fprintf( stderr, "in _tanni seq2[0] = %-*.*s\n", jen-jst+1, jen-jst+1, seq2[0]+jst );
#endif


	ll1 = ( (int)(lgth1) ) + 100;
	ll2 = ( (int)(lgth2) ) + 100;

//	aseq1 = AllocateCharMtx( icyc, 0 );
//	aseq2 = AllocateCharMtx( jcyc, 0 );
//	aseq1bk = AllocateCharMtx( icyc, lgth1+lgth2+100 );
//	aseq2bk = AllocateCharMtx( jcyc, lgth1+lgth2+100 );
//	for( i=0; i<icyc; i++ ) aseq1[i] = aseq1bk[i];
//	for( i=0; i<jcyc; i++ ) aseq2[i] = aseq2bk[i];

	w1 = AllocateFloatVec( ll2+2 );
	w2 = AllocateFloatVec( ll2+2 );

	initverticalw = AllocateFloatVec( ll1+2 );
	lastverticalw = AllocateFloatVec( ll1+2 );

	m = AllocateFloatVec( ll2+2 );
	mp = AllocateIntVec( ll2+2 );

	floatwork = AllocateFloatCub( maxdistclass, MAX( ll1, ll2 )+2, nalphabets+1 ); 
	intwork = AllocateIntCub( maxdistclass, MAX( ll1, ll2 )+2, nalphabets+1 ); 


	intmtx = AllocateIntMtx( ll1+1, ll2+1 );

	ijp = intmtx;

	currentw = w1;
	previousw = w2;

#if 0
	match_calc( n_dynamicmtx, initverticalw, cpmx2+jst, cpmx1+ist, 0, lgth1, floatwork, intwork, 1 );
	match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, 0, lgth2, floatwork, intwork, 1 );
#else
	fillzero( initverticalw, lgth1 );
	for( c=0; c<maxdistclass; c++ )
		match_calc_add( matrices[c], initverticalw, cpmx2s[c]+jst, cpmx1s[c]+ist, 0, lgth1, floatwork[c], intwork[c], 1 );

	fillzero( currentw, lgth2 );
	for( c=0; c<maxdistclass; c++ )
		match_calc_add( matrices[c], currentw, cpmx1s[c]+ist, cpmx2s[c]+jst, 0, lgth2, floatwork[c], intwork[c], 1 );
#endif

	if( headgp || ist != 0 )
	{
		for( i=1; i<lgth1+1; i++ )
		{
			initverticalw[i] += ( ogcp1[0] * headgapfreq2 + fgcp1[i-1] * gapfreq2f[0] ) ;
		}
	}
	if( headgp || jst != 0 )
	{
		for( j=1; j<lgth2+1; j++ )
		{
			currentw[j] += ( ogcp2[0] * headgapfreq1 + fgcp2[j-1] * gapfreq1f[0] ) ;
		}
	}

	for( j=1; j<lgth2+1; ++j ) 
	{
		m[j] = currentw[j-1] + ogcp1[1] * gapfreq2f[j-1]; mp[j] = 0;;
	}

	lastverticalw[0] = currentw[lgth2-1];



	if( tailgp || jen != fulllen2-1 ) lasti = lgth1+1; else lasti = lgth1;
//	if( 1 ) lasti = lgth1+1; else lasti = lgth1;
	for( i=1; i<lasti; i++ )
	{
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

#if 0
		match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, i, lgth2, floatwork, intwork, 0 );
#else
		fillzero( currentw, lgth2 );
		for( c=0; c<maxdistclass; c++ )
			match_calc_add( matrices[c], currentw, cpmx1s[c]+ist, cpmx2s[c]+jst, i, lgth2, floatwork[c], intwork[c], 0 );
#endif
		currentw[0] = initverticalw[i];

		mi = previousw[0] + ogcp2[1] * gapfreq1f[i-1];
		mpi = 0;

		ijppt = ijp[i] + 1;
		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;
		mpjpt = mp + 1;
//		if( tailgp && jen != fulllen2-1 ) lastj = lgth2+1; else lastj = lgth2;
		lastj = lgth2+1; 
		for( j=1; j<lastj; j++ )
		{
			wm = *prept;
			*ijppt = 0;

#if 0
			fprintf( stderr, "%5.0f->", wm );
#endif
			g = mi + fgcp2[j-1]  * gapfreq1f[i];
#if 0
			fprintf( stderr, "%5.0f?", g );
#endif
			if( g > wm )
			{
				wm = g;
				*ijppt = -( j - mpi );
			}
			g = *prept + ogcp2[j] * gapfreq1f[i-1];
			if( g >= mi )
			{
				mi = g;
				mpi = j-1;
			}
#if USE_PENALTY_EX
			mi += fpenalty_ex;
#endif

			g = *mjpt + fgcp1[i-1] * gapfreq2f[j];
#if 0 
			fprintf( stderr, "%5.0f?", g );
#endif
			if( g > wm )
			{
				wm = g;
				*ijppt = +( i - *mpjpt );
			}
			g = *prept + ogcp1[i] * gapfreq2f[j-1];
			if( g >= *mjpt )
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
			*curpt += wm;


			ijppt++;
			mjpt++;
			prept++;
			mpjpt++;
			curpt++;
		}
		lastverticalw[i] = currentw[lgth2-1];
	}

//	fprintf( stderr, "wm = %f\n", wm );

	Atracking( currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, icyc, jcyc, ist, ien, jst, jen, fulllen1, fulllen2, tailgp );
#if 0
	fprintf( stderr, "res in _tanni mseq1[0] = %s\n", mseq1[0] );
	fprintf( stderr, "res in _tanni mseq2[0] = %s\n", mseq2[0] );
#endif

//	for( i=0; i<icyc; i++ ) strcpy( mseq1[i], aseq1[i] );
//	for( i=0; i<jcyc; i++ ) strcpy( mseq2[i], aseq2[i] );

//	fprintf( stderr, "in _tanni, aseq1 = %s\n", aseq1[0] );
//	fprintf( stderr, "in _tanni, mseq1 = %s\n", mseq1[0] );

	FreeFloatVec( w1 );
	FreeFloatVec( w2 );
	FreeFloatVec( initverticalw );
	FreeFloatVec( lastverticalw );

	FreeFloatVec( m );
	FreeIntVec( mp );


	FreeFloatCub( floatwork );
	FreeIntCub( intwork );

	FreeIntMtx( intmtx );


//	FreeCharMtx( aseq1bk );
//	FreeCharMtx( aseq2bk );

//	free( aseq1 );
//	free( aseq2 );

	return( wm );

}

static void freearrays_rec1_variousdist(
	float *w1, float *w2, float *initverticalw, float *lastverticalw,
	float *midw, float *midm, float *midn,
	int *jumpbacki, int *jumpbackj, int *jumpforwi, int *jumpforwj, int *jumpdummi, int *jumpdummj,
	float *m, int *mp,
	float ***floatwork, int ***intwork
#if STOREWM
	, float **WMMTX, float **WMMTX2
#endif
)
{

	FreeFloatVec( w1 );
	FreeFloatVec( w2 );
	FreeFloatVec( initverticalw );
	FreeFloatVec( lastverticalw );
	FreeFloatVec( midw );
	FreeFloatVec( midm );
	FreeFloatVec( midn );

	FreeIntVec( jumpbacki );
	FreeIntVec( jumpbackj );
	FreeIntVec( jumpforwi );
	FreeIntVec( jumpforwj );
	FreeIntVec( jumpdummi );
	FreeIntVec( jumpdummj );

	FreeFloatVec( m );
	FreeIntVec( mp );

	FreeFloatCub( floatwork );
	FreeIntCub( intwork );

#if STOREWM
	FreeFloatMtx( WMMTX );
	FreeFloatMtx( WMMTX2 );
#endif
}

static float MSalignmm_rec_variousdist( double ***matrices, int icyc, int jcyc, char **seq1, char **seq2, float ***cpmx1s, float ***cpmx2s, int ist, int ien, int jst, int jen, int alloclen, int fulllen1, int fulllen2, char **mseq1, char **mseq2, int depth, float **gapinfo, int *chudanpt, int chudanref, int *chudanres, int headgp, int tailgp, float headgapfreq1_g, float headgapfreq2_g )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{
//	int k;
	int alnlen;
	float value = 0.0;
	register int i, j, c;
	char **aseq1, **aseq2;
	int ll1, ll2, l, len;
	int lasti, lastj, imid;
	int jmid = 0; // by D.Mathog, a guess
	float wm = 0.0;   /* int ?????? */
	float g;
	float *currentw, *previousw;
#if USE_PENALTY_EX
	float fpenalty_ex = (float)penalty_ex;
#endif
	float *wtmp;
//	short *ijppt;
	int *mpjpt;
//	short **ijp;
	int *mp;
	int mpi;
	float *mjpt, *prept, *curpt;
	float mi;
	float *m;
	float *w1, *w2;
//	float *match;
	float *initverticalw;    /* kufuu sureba iranai */
	float *lastverticalw;    /* kufuu sureba iranai */
	int ***intwork;
	float ***floatwork;
//	short **shortmtx;
#if STOREWM
	float **WMMTX;
	float **WMMTX2;
#endif
	float *midw;
	float *midm;
	float *midn;
	int lgth1, lgth2;
	float maxwm;
	int *jumpforwi;
	int *jumpforwj;
	int *jumpbacki;
	int *jumpbackj;
	int *jumpdummi; //muda
	int *jumpdummj = NULL; // by D.Mathog, a guess
	int jumpi, jumpj = 0; // by D.Mathog, a guess
	char *gaps;
	int ijpi, ijpj;
	float *ogcp1;
	float *fgcp1;
	float *ogcp2;
	float *fgcp2;
	float firstm;
	int firstmp;
#if STOREWM
	static TLS char ttt1[50000];
	static TLS char ttt2[50000];
#endif
	float *gapfreq1f;
	float *gapfreq2f;
	float headgapfreq1;
	float headgapfreq2;

#if 0
	int nglen1, nglen2;
	nglen1 = seqlen( seq1[0] );
	nglen2 = seqlen( seq2[0] );
#endif

//	fprintf( stderr, "fulllen1 = %d, fulllen2 = %d, headgp = %d, tailgp = %d\n", fulllen1, fulllen2, headgp, tailgp );

	ogcp1 = gapinfo[0] + ist;
	fgcp1 = gapinfo[1] + ist;
	ogcp2 = gapinfo[2] + jst;
	fgcp2 = gapinfo[3] + jst;
	gapfreq1f = gapinfo[4] + ist;
	gapfreq2f = gapinfo[5] + jst;

	if( ist > 0 ) headgapfreq1 = gapfreq1f[-1];
	else headgapfreq1 = headgapfreq1_g;
	if( jst > 0 ) headgapfreq2 = gapfreq2f[-1];
	else headgapfreq2 = headgapfreq2_g;

	depth++;
	reccycle++;

	lgth1 = ien-ist+1;
	lgth2 = jen-jst+1;

//	if( lgth1 < 5 )
//		fprintf( stderr, "\nWARNING: lgth1 = %d\n", lgth1 );
//	if( lgth2 < 5 )
//		fprintf( stderr, "\nWARNING: lgth2 = %d\n", lgth2 );
//


#if STOREWM
	fprintf( stderr, "==== MSalign (depth=%d, reccycle=%d), ist=%d, ien=%d, jst=%d, jen=%d\n", depth, reccycle, ist, ien, jst, jen );
	strncpy( ttt1, seq1[0]+ist, lgth1 );
	strncpy( ttt2, seq2[0]+jst, lgth2 );
	ttt1[lgth1] = 0;
	ttt2[lgth2] = 0;
	fprintf( stderr, "seq1 = %s\n", ttt1 );
	fprintf( stderr, "seq2 = %s\n", ttt2 );
#endif
	if( lgth2 <= 0 ) // lgth1 <= 0 ha?
	{
//		fprintf( stderr, "\n\n==== jimei\n\n" );
//		exit( 1 );
		for( i=0; i<icyc; i++ ) 
		{
			strncpy( mseq1[i], seq1[i]+ist, lgth1 );
			mseq1[i][lgth1] = 0;
		}
		for( i=0; i<jcyc; i++ ) 
		{
			mseq2[i][0] = 0;
			for( j=0; j<lgth1; j++ )
//				strcat( mseq2[i], "-" );
				strcat( mseq2[i], newgapstr );
		}

//		fprintf( stderr, "==== mseq1[0] (%d) = %s\n", depth, mseq1[0] );
//		fprintf( stderr, "==== mseq2[0] (%d) = %s\n", depth, mseq2[0] );

		return( 0.0 );
	}

#if MEMSAVE
	aseq1 = AllocateCharMtx( icyc, 0 );
	aseq2 = AllocateCharMtx( jcyc, 0 );
	for( i=0; i<icyc; i++ ) aseq1[i] = mseq1[i];
	for( i=0; i<jcyc; i++ ) aseq2[i] = mseq2[i];
#else
	aseq1 = AllocateCharMtx( icyc, lgth1+lgth2+100 );
	aseq2 = AllocateCharMtx( jcyc, lgth1+lgth2+100 );
#endif

//	fprintf( stderr, "####(s) seq1[0] (%d) = \n%-*.*s\n (a%d-a%d)\n", depth, ien-ist+1, ien-ist+1, seq1[0]+ist, ist, ien );
//	fprintf( stderr, "####(s) seq2[0] (%d) = \n%-*.*s\n (b%d-b%d)\n", depth, jen-jst+1, jen-jst+1, seq2[0]+jst, jst, jen );

//  if( lgth1 < DPTANNI && lgth2 < DPTANNI ) // & dato lgth ==1 no kanousei ga arunode yokunai 
//    if( lgth1 < DPTANNI ) // kore mo lgth2 ga mijikasugiru kanousei ari
    if( lgth1 < DPTANNI || lgth2 < DPTANNI ) // zettai ni anzen ka?
	{
//		fprintf( stderr, "==== Going to _tanni\n" );

		value = MSalignmm_tanni_variousdist( matrices, icyc, jcyc, seq1, seq2, cpmx1s, cpmx2s, ist, ien, jst, jen, alloclen, fulllen1, fulllen2, aseq1, aseq2, gapinfo, headgp, tailgp, headgapfreq1_g, headgapfreq2_g );	


#if MEMSAVE
		free( aseq1 );
		free( aseq2 );
#else
		for( i=0; i<icyc; i++ ) strcpy( mseq1[i], aseq1[i] );
		for( i=0; i<jcyc; i++ ) strcpy( mseq2[i], aseq2[i] );

		FreeCharMtx( aseq1 );
		FreeCharMtx( aseq2 );
#endif

//		fprintf( stderr, "value = %f\n", value );

		return( value );
	}
//	fprintf( stderr, "Trying to divide the mtx\n" );

	ll1 = ( (int)(lgth1) ) + 100;
	ll2 = ( (int)(lgth2) ) + 100;

//	fprintf( stderr, "ll1,ll2=%d,%d\n", ll1, ll2 );

	w1 = AllocateFloatVec( ll2+2 );
	w2 = AllocateFloatVec( ll2+2 );
//	match = AllocateFloatVec( ll2+2 );
	midw = AllocateFloatVec( ll2+2 );
	midn = AllocateFloatVec( ll2+2 );
	midm = AllocateFloatVec( ll2+2 );
	jumpbacki = AllocateIntVec( ll2+2 );
	jumpbackj = AllocateIntVec( ll2+2 );
	jumpforwi = AllocateIntVec( ll2+2 );
	jumpforwj = AllocateIntVec( ll2+2 );
	jumpdummi = AllocateIntVec( ll2+2 );
	jumpdummj = AllocateIntVec( ll2+2 );

	initverticalw = AllocateFloatVec( ll1+2 );
	lastverticalw = AllocateFloatVec( ll1+2 );

	m = AllocateFloatVec( ll2+2 );
	mp = AllocateIntVec( ll2+2 );
	gaps = AllocateCharVec( MAX( ll1, ll2 ) + 2 );

	floatwork = AllocateFloatCub( maxdistclass, MAX( ll1, ll2 )+2, nalphabets ); 
	intwork = AllocateIntCub( maxdistclass, MAX( ll1, ll2 )+2, nalphabets ); 

#if DEBUG
	fprintf( stderr, "succeeded\n" );
#endif

#if STOREWM
	WMMTX = AllocateFloatMtx( ll1, ll2 );
	WMMTX2 = AllocateFloatMtx( ll1, ll2 );
#endif
#if 0
	shortmtx = AllocateShortMtx( ll1, ll2 );

#if DEBUG
	fprintf( stderr, "succeeded\n\n" );
#endif

	ijp = shortmtx;
#endif

	currentw = w1;
	previousw = w2;

#if 0
	match_calc( n_dynamicmtx, initverticalw, cpmx2+jst, cpmx1+ist, 0, lgth1, floatwork, intwork, 1 );
	match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, 0, lgth2, floatwork, intwork, 1 );
#else
	fillzero( initverticalw, lgth1 );
	for( c=0; c<maxdistclass; c++ )
		match_calc_add( matrices[c], initverticalw, cpmx2s[c]+jst, cpmx1s[c]+ist, 0, lgth1, floatwork[c], intwork[c], 1 );

	fillzero( currentw, lgth2 );
	for( c=0; c<maxdistclass; c++ )
		match_calc_add( matrices[c], currentw, cpmx1s[c]+ist, cpmx2s[c]+jst, 0, lgth2, floatwork[c], intwork[c], 1 );
#endif


	for( i=1; i<lgth1+1; i++ )
	{
//		initverticalw[i] += ( ogcp1[0] + fgcp1[i-1] ) ;
		initverticalw[i] += ( ogcp1[0] * headgapfreq2 + fgcp1[i-1] * gapfreq2f[0] ) ;
	}
	for( j=1; j<lgth2+1; j++ )
	{
		currentw[j] += ( ogcp2[0] * headgapfreq1 + fgcp2[j-1] * gapfreq1f[0]) ;
	}

#if STOREWM
	WMMTX[0][0] = initverticalw[0];
	for( i=1; i<lgth1+1; i++ )
	{
		WMMTX[i][0] = initverticalw[i];
	}
	for( j=1; j<lgth2+1; j++ )
	{
		WMMTX[0][j] = currentw[j];
	}
#endif


	for( j=1; j<lgth2+1; ++j ) 
	{
		m[j] = currentw[j-1] + ogcp1[1] * gapfreq2f[j-1];
//		m[j] = currentw[j-1];
		mp[j] = 0;
	}

	lastverticalw[0] = currentw[lgth2-1];

	imid = lgth1 * 0.5;

	jumpi = 0; // atode kawaru.
	lasti = lgth1+1;
#if STOREWM
	for( i=1; i<lasti; i++ )
#else
	for( i=1; i<=imid; i++ )
#endif
	{
#ifdef enablemultithread
//		fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
		if( chudanpt && *chudanpt != chudanref ) 
		{
//			fprintf( stderr, "\n\n## CHUUDAN!!! zenhan\n" );
			*chudanres = 1;
			freearrays_rec1_variousdist
			(
				w1, w2, initverticalw, lastverticalw, midw, midm, midn,
				jumpbacki, jumpbackj, jumpforwi, jumpforwj, jumpdummi, jumpdummj,
				m, mp,
				floatwork, intwork
#if STOREWM
				, WMMTX, WMMTX2
#endif
			);
			freearrays_rec2( gaps, aseq1, aseq2 );
			return( -1.0 );
		}
#endif
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

#if 0
		match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, i, lgth2, floatwork, intwork, 0 );
#else
		fillzero( currentw, lgth2 );
		for( c=0; c<maxdistclass; c++ )
			match_calc_add( matrices[c], currentw, cpmx1s[c]+ist, cpmx2s[c]+jst, i, lgth2, floatwork[c], intwork[c], 0 );
#endif
		currentw[0] = initverticalw[i];

		m[0] = ogcp1[i];
#if STOREM
		WMMTX2[i][0] = m[0];
#endif
		if( i == imid ) midm[0] = m[0];

		mi = previousw[0] + ogcp2[1] * gapfreq1f[i-1]; 
//		mi = previousw[0];
		mpi = 0;


//		ijppt = ijp[i] + 1;
		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;
		mpjpt = mp + 1;


		lastj = lgth2+1;
		for( j=1; j<lastj; j++ )
		{

			wm = *prept;

#if 0
			fprintf( stderr, "%5.0f->", wm );
#endif
			g = mi + fgcp2[j-1] * gapfreq1f[i];
//			g = mi + fpenalty;
#if 0
			fprintf( stderr, "%5.0f?", g );
#endif
			if( g > wm )
			{
				wm = g;
//				*ijppt = -( j - mpi );
			}
			g = *prept + ogcp2[j] * gapfreq1f[i-1];
//			g = *prept;
			if( g >= mi )
			{
				mi = g;
				mpi = j-1;
			}
#if USE_PENALTY_EX
			mi += fpenalty_ex;
#endif

			g = *mjpt + fgcp1[i-1] * gapfreq2f[j];
//			g = *mjpt + fpenalty;
#if 0 
			fprintf( stderr, "%5.0f?", g );
#endif
			if( g > wm )
			{
				wm = g;
//				*ijppt = +( i - *mpjpt );
			}


			g = *prept + ogcp1[i] * gapfreq2f[j-1];
//			g = *prept;
			if( g >= *mjpt )
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
			*curpt += wm;

#if STOREWM
			WMMTX[i][j] = *curpt;
			WMMTX2[i][j] = *mjpt;
#endif

			if( i == imid ) //muda
			{	
				jumpbackj[j] = *mpjpt; // muda atode matomeru
				jumpbacki[j] = mpi; // muda atode matomeru
//				fprintf( stderr, "jumpbackj[%d] in forward dp is %d\n", j, *mpjpt );
//				fprintf( stderr, "jumpbacki[%d] in forward dp is %d\n", j, mpi );
				midw[j] = *curpt;
				midm[j] = *mjpt;
				midn[j] = mi;
			}

//			fprintf( stderr, "m[%d] = %f\n", j, m[j] );
			mjpt++;
			prept++;
			mpjpt++;
			curpt++;

		}
		lastverticalw[i] = currentw[lgth2-1];

#if STOREWM
		WMMTX2[i][lgth2] = m[lgth2-1];
#endif

#if 0  // ue
		if( i == imid )
		{
			for( j=0; j<lgth2; j++ ) midw[j] = currentw[j];
			for( j=0; j<lgth2; j++ ) midm[j] = m[j];
		}
#endif
	}
//	for( j=0; j<lgth2; j++ ) midw[j] = WMMTX[imid][j];
//	for( j=0; j<lgth2; j++ ) midm[j] = WMMTX2[imid][j];

#if 0
    for( i=0; i<lgth1; i++ )
    {
        for( j=0; j<lgth2; j++ )
        {
            fprintf( stderr, "% 10.2f ", WMMTX[i][j] );
        }
        fprintf( stderr, "\n" );
    }
	fprintf( stderr, "\n" );
	fprintf( stderr, "WMMTX2 = \n" );
    for( i=0; i<lgth1; i++ )
    {
        for( j=0; j<lgth2; j++ )
        {
            fprintf( stderr, "% 10.2f ", WMMTX2[i][j] );
        }
        fprintf( stderr, "\n" );
    }
	fprintf( stderr, "\n" );
#endif

// gyakudp

#if 0
	match_calc( n_dynamicmtx, initverticalw, cpmx2+jst, cpmx1+ist, lgth2-1, lgth1, floatwork, intwork, 1 );
	match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, lgth1-1, lgth2, floatwork, intwork, 1 );
#else
	fillzero( initverticalw, lgth1 );
	for( c=0; c<maxdistclass; c++ )
		match_calc_add( matrices[c], initverticalw, cpmx2s[c]+jst, cpmx1s[c]+ist, lgth2-1, lgth1, floatwork[c], intwork[c], 1 );

	fillzero( currentw, lgth2 );
	for( c=0; c<maxdistclass; c++ )
		match_calc_add( matrices[c], currentw, cpmx1s[c]+ist, cpmx2s[c]+jst, lgth1-1, lgth2, floatwork[c], intwork[c], 1 );
#endif

	for( i=0; i<lgth1-1; i++ )
	{
//		initverticalw[i] += ( fgcp1[lgth1-1] + ogcp1[i+1] );
		initverticalw[i] += ( fgcp1[lgth1-1] * gapfreq2f[lgth2] + ogcp1[i+1] * gapfreq2f[lgth2-1] );
	}
	for( j=0; j<lgth2-1; j++ )
	{
//		currentw[j] += ( fgcp2[lgth2-1] + ogcp2[j+1] );
		currentw[j] += ( fgcp2[lgth2-1] * gapfreq1f[lgth1] + ogcp2[j+1] * gapfreq1f[lgth1-1] );
	}

#if STOREWM
	for( i=0; i<lgth1-1; i++ )
	{
		WMMTX[i][lgth2-1] += ( fgcp1[lgth1-1] + ogcp1[i+1] );
		fprintf( stderr, "fgcp1[lgth1-1] + ogcp1[i+1] = %f\n", fgcp1[lgth1-1] + ogcp1[i+1] );
	}
	for( j=0; j<lgth2-1; j++ )
	{
		WMMTX[lgth1-1][j] += ( fgcp2[lgth2-1] + ogcp2[j+1] );
		fprintf( stderr, "fgcp2[lgth2-1] + ogcp2[j+1] = %f\n", fgcp2[lgth2-1] + ogcp2[j+1] );
	}
#endif






#if 0
	for( j=lgth2-1; j>0; --j )
	{
		m[j-1] = currentw[j] + fgcp2[lgth2-2];
//		m[j-1] = currentw[j];
		mp[j] = lgth1-1;
	}
#else
	for( j=lgth2-1; j>-1; --j )
	{
		m[j] = currentw[j+1] + fgcp1[lgth1-2] * gapfreq2f[j+1];
//		m[j-1] = currentw[j];
		mp[j] = lgth1-1;
	}
#endif

//	for( j=0; j<lgth2; j++ ) m[j] = 0.0;
	// m[lgth2-1] ha irunoka?


//	for( i=lgth1-2; i>=imid; i-- )
	firstm = -9999999.9;
//	firstmp = lgth1-1;
	firstmp = lgth1;
	for( i=lgth1-2; i>-1; i-- )
	{
#ifdef enablemultithread
//		fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
		if( chudanpt && *chudanpt != chudanref ) 
		{
//			fprintf( stderr, "\n\n## CHUUDAN!!! kouhan\n" );
			*chudanres = 1;
			freearrays_rec1_variousdist
			(
				w1, w2, initverticalw, lastverticalw, midw, midm, midn,
				jumpbacki, jumpbackj, jumpforwi, jumpforwj, jumpdummi, jumpdummj,
				m, mp,
				floatwork, intwork
#if STOREWM
				, WMMTX, WMMTX2
#endif
			);
			freearrays_rec2( gaps, aseq1, aseq2 );
			return( -1.0 );
		}
#endif
		wtmp = previousw;
		previousw = currentw;
		currentw = wtmp;
		previousw[lgth2-1] = initverticalw[i+1];
#if 0
		match_calc( n_dynamicmtx, currentw, cpmx1+ist, cpmx2+jst, i, lgth2, floatwork, intwork, 0 );
#else
		fillzero( currentw, lgth2 );
		for( c=0; c<maxdistclass; c++ )
			match_calc_add( matrices[c], currentw, cpmx1s[c]+ist, cpmx2s[c]+jst, i, lgth2, floatwork[c], intwork[c], 0 );
#endif

		currentw[lgth2-1] = initverticalw[i];

//		m[lgth2] = fgcp1[i];
//		WMMTX2[i][lgth2] += m[lgth2];
//		fprintf( stderr, "m[] = %f\n", m[lgth2] );

		mi = previousw[lgth2-1] + fgcp2[lgth2-2] * gapfreq1f[i+1];
//		mi = previousw[lgth2-1];
		mpi = lgth2 - 1;

		mjpt = m + lgth2 - 2;
		prept = previousw + lgth2 - 1;
		curpt = currentw + lgth2 - 2;
		mpjpt = mp + lgth2 - 2;


		for( j=lgth2-2; j>-1; j-- )
		{
			wm = *prept;
			ijpi = i+1;
			ijpj = j+1;

			g = mi + ogcp2[j+1] * gapfreq1f[i];
//			g = mi + fpenalty;
			if( g > wm )
			{
				wm = g;
				ijpj = mpi;
				ijpi = i+1;
			}

			g = *prept + fgcp2[j] * gapfreq1f[i+1];
//			g = *prept;
			if( g >= mi )
			{
//				fprintf( stderr, "i,j=%d,%d - renewed! mpi = %d\n", i, j, j+1 );
				mi = g;
				mpi = j + 1;
			}

#if USE_PENALTY_EX
			mi += fpenalty_ex;
#endif

//			fprintf( stderr, "i,j=%d,%d *mpjpt = %d\n", i, j, *mpjpt );
			g = *mjpt + ogcp1[i+1] * gapfreq2f[j];
//			g = *mjpt + fpenalty;
			if( g > wm )
			{
				wm = g;
				ijpi = *mpjpt;
				ijpj = j+1;
			}

//			if( i == imid )fprintf( stderr, "i,j=%d,%d \n", i, j );
			g = *prept + fgcp1[i] * gapfreq2f[j+1];
//			g = *prept;
			if( g >= *mjpt )
			{
				*mjpt = g;
				*mpjpt = i + 1;
			}

#if USE_PENALTY_EX
			m[j] += fpenalty_ex;
#endif

			if( i == jumpi || i == imid - 1 )
			{
				jumpforwi[j] = ijpi; //muda
				jumpforwj[j] = ijpj; //muda
//				fprintf( stderr, "jumpfori[%d] = %d\n", j, ijpi );
//				fprintf( stderr, "jumpforj[%d] = %d\n", j, ijpj );
			}
			if( i == imid ) // muda
			{
				midw[j] += wm;
//				midm[j+1] += *mjpt + fpenalty; //??????
				midm[j+1] += *mjpt; //??????
			}
			if( i == imid - 1 )
			{
//				midn[j] += mi + fpenalty;  //????
				midn[j] += mi;  //????
			}
#if STOREWM
			WMMTX[i][j] += wm;
//			WMMTX2[i][j+1] += *mjpt + fpenalty;
			WMMTX2[i][j+1] += *mjpt;
#endif
			*curpt += wm;

			mjpt--;
			prept--;
			mpjpt--;
			curpt--;
		}
//		fprintf( stderr, "adding *mjpt (=%f) to WMMTX2[%d][%d]\n", *mjpt, i, j+1 );
		g = *prept + fgcp1[i];
		if( firstm < g ) 
		{
			firstm = g;
			firstmp = i + 1;
		}
#if STOREWM
		WMMTX2[i][j+1] += firstm;
#endif
		if( i == imid ) midm[j+1] += firstm;


		if( i == imid - 1 )	
		{
			maxwm = midw[1];
			jmid = 0;
//			if( depth == 1 ) fprintf( stderr, "maxwm!! = %f\n", maxwm );
			for( j=2; j<lgth2-1; j++ )
			{
				wm = midw[j];
				if( wm > maxwm )
				{
					jmid = j;
					maxwm = wm;
				}
//				if( depth == 1 ) fprintf( stderr, "maxwm!! = %f\n", maxwm );
			}
			for( j=0; j<lgth2+1; j++ )
			{
				wm = midm[j];
				if( wm > maxwm )
				{
					jmid = j;
					maxwm = wm;
				}
//				if( depth == 1 ) fprintf( stderr, "maxwm!! = %f\n", maxwm );
			}

//			if( depth == 1 ) fprintf( stderr, "maxwm!! = %f\n", maxwm );


//			fprintf( stderr, "### imid=%d, jmid=%d\n", imid, jmid );
			wm = midw[jmid];
			jumpi = imid-1;
			jumpj = jmid-1;
			if( jmid > 0 && midn[jmid-1] > wm ) //060413
			{
				jumpi = imid-1;
				jumpj = jumpbacki[jmid];
				wm = midn[jmid-1];
//				fprintf( stderr, "rejump (n)\n" );
			}
			if( midm[jmid] > wm )
			{
				jumpi = jumpbackj[jmid];
				jumpj = jmid-1;
				wm = midm[jmid];
//				fprintf( stderr, "rejump (m) jumpi=%d\n", jumpi );
			}


//			fprintf( stderr, "--> imid=%d, jmid=%d\n", imid, jmid );
//			fprintf( stderr, "--> jumpi=%d, jumpj=%d\n", jumpi, jumpj );
#if STOREWM
			fprintf( stderr, "imid = %d\n", imid );
			fprintf( stderr, "midn = \n" );
			for( j=0; j<lgth2; j++ )
			{
				fprintf( stderr, "% 7.1f ", midn[j] );
			}
			fprintf( stderr, "\n" );
			fprintf( stderr, "midw = \n" );
			for( j=0; j<lgth2; j++ )
			{
				fprintf( stderr, "% 7.1f ", midw[j] );
			}
			fprintf( stderr, "\n" );
			fprintf( stderr, "midm = \n" );
			for( j=0; j<lgth2; j++ )
			{
				fprintf( stderr, "% 7.1f ", midm[j] );
			}
			fprintf( stderr, "\n" );
#endif
//			fprintf( stderr, "maxwm = %f\n", maxwm );
		}
		if( i == jumpi ) //saki?
		{
//			fprintf( stderr, "#### FIRST i=%d, jumpi<imid=%d<%d, ist=%d, ien=%d, firstmp=%d, lgth1=%d\n", i, jumpi, imid, ist, ien, firstmp, lgth1 );
//			fprintf( stderr, "#### mark 1\n" );
//			fprintf( stderr, "imid, jumpi = %d,%d\n", imid, jumpi );
//			fprintf( stderr, "jmid, jumpj = %d,%d\n", jmid, jumpj );

			if( jmid == 0 )
			{
//				fprintf( stderr, "#### CHUI2!\n" );
				jumpj = 0; jmid = 1;
#if 0 // v6.823 made
				jumpi = firstmp-1;
				imid = firstmp;
#endif
#if 0
				jumpi = 0;
				imid = 1;
#else
//				if( 1 || firstmp > 100 ) // naze 100
				if( imid < firstmp-1 ) // naze 100
				{
					jumpi = firstmp;
					imid = firstmp+1;
				}
#if 0
				else
				{
					jumpi = 0;
					imid = 1;
				}
#endif
#endif
			}

#if 0
			else if( jmid == lgth2 ) 
			{
				fprintf( stderr, "CHUI1!\n" );
				jumpi=0; jumpj=0;
				imid=jumpforwi[0]; jmid=lgth2-1;
			}
#else // 060414
			else if( jmid >= lgth2 ) 
			{
//				fprintf( stderr, "CHUI1!\n" );
				jumpi=imid-1; jmid=lgth2;
				jumpj = lgth2-1;
			}
#endif
			else
			{
//				fprintf( stderr, "#### CHUI3!\n" );
				imid = jumpforwi[jumpj];
				jmid = jumpforwj[jumpj];
				if( imid == jumpi ) jumpi = imid-1;
			}
#if 0
			fprintf( stderr, "jumpi -> %d\n", jumpi );
			fprintf( stderr, "jumpj -> %d\n", jumpj );
			fprintf( stderr, "imid -> %d\n", imid );
			fprintf( stderr, "jmid -> %d\n", jmid );
#endif
//			fprintf( stderr, "#### FINAL i=%d, jumpi<imid=%d<%d, ist=%d, ien=%d, firstmp=%d\n", i, jumpi, imid, ist, ien, firstmp );

#if STOREWM
			break;
#else
			break;
#endif
		}
	}

//	fprintf( stderr, "imid = %d, but jumpi = %d\n", imid, jumpi );
//	fprintf( stderr, "jmid = %d, but jumpj = %d\n", jmid, jumpj );

//	for( j=0; j<lgth2; j++ ) midw[j] += currentw[j];
//	for( j=0; j<lgth2; j++ ) midm[j] += m[j+1];
//	for( j=0; j<lgth2; j++ ) midw[j] += WMMTX[imid][j];
//	for( j=0; j<lgth2; j++ ) midw[j] += WMMTX[imid][j];


#if STOREWM
	fprintf( stderr, "WMMTX = \n" );
    for( i=0; i<lgth1; i++ )
    {
        fprintf( stderr, "%d ", i );
        for( j=0; j<lgth2; j++ )
        {
            fprintf( stderr, "% 7.2f ", WMMTX[i][j] );
        }
        fprintf( stderr, "\n" );
    }
//	fprintf( stderr, "WMMTX2 = (p = %f)\n", fpenalty );
    for( i=0; i<lgth1; i++ )
    {
        fprintf( stderr, "%d ", i );
        for( j=0; j<lgth2+1; j++ )
        {
            fprintf( stderr, "% 7.2f ", WMMTX2[i][j] );
        }
        fprintf( stderr, "\n" );
    }

	fprintf( stderr, "jumpbacki = \n" );
	for( j=0; j<lgth2; j++ )
	{
		fprintf( stderr, "% 7d ", jumpbacki[j] );
	}
	fprintf( stderr, "\n" );
	fprintf( stderr, "jumpbackj = \n" );
	for( j=0; j<lgth2; j++ )
	{
		fprintf( stderr, "% 7d ", jumpbackj[j] );
	}
	fprintf( stderr, "\n" );
	fprintf( stderr, "jumpforwi = \n" );
	for( j=0; j<lgth2; j++ )
	{
		fprintf( stderr, "% 7d ", jumpforwi[j] );
	}
	fprintf( stderr, "\n" );
	fprintf( stderr, "jumpforwj = \n" );
	for( j=0; j<lgth2; j++ )
	{
		fprintf( stderr, "% 7d ", jumpforwj[j] );
	}
	fprintf( stderr, "\n" );


#endif


//	Atracking( currentw, lastverticalw, seq1, seq2, mseq1, mseq2, cpmx1, cpmx2, ijp, icyc, jcyc );

#if 0 // irukamo
	resultlen = strlen( mseq1[0] );
	if( alloclen < resultlen || resultlen > N )
	{
		fprintf( stderr, "alloclen=%d, resultlen=%d, N=%d\n", alloclen, resultlen, N );
		ErrorExit( "LENGTH OVER!\n" );
	}
#endif



#if 0
	fprintf( stderr, "jumpi = %d, imid = %d\n", jumpi, imid );
	fprintf( stderr, "jumpj = %d, jmid = %d\n", jumpj, jmid );

	fprintf( stderr, "imid = %d\n", imid );
	fprintf( stderr, "jmid = %d\n", jmid );
#endif

	freearrays_rec1_variousdist
	(
		w1, w2, initverticalw, lastverticalw, midw, midm, midn,
		jumpbacki, jumpbackj, jumpforwi, jumpforwj, jumpdummi, jumpdummj,
		m, mp,
		floatwork, intwork
#if STOREWM
		, WMMTX, WMMTX2
#endif
	);


//	fprintf( stderr, "==== calling myself (first)\n" );

	value = MSalignmm_rec_variousdist( matrices, icyc, jcyc, seq1, seq2, cpmx1s, cpmx2s, ist, ist+jumpi, jst, jst+jumpj, alloclen, fulllen1, fulllen2, aseq1, aseq2, depth, gapinfo, NULL, 0, NULL, headgp, tailgp, headgapfreq1_g, headgapfreq2_g ); // chudan mada	
#if 0
		fprintf( stderr, "aseq1[0] = %s\n", aseq1[0] );
		fprintf( stderr, "aseq2[0] = %s\n", aseq2[0] );
#endif
#if MEMSAVE
#else
	for( i=0; i<icyc; i++ ) strcpy( mseq1[i], aseq1[i] );
	for( i=0; i<jcyc; i++ ) strcpy( mseq2[i], aseq2[i] );
#endif

//	fprintf( stderr, "====(f) aseq1[0] (%d) = %s (%d-%d)\n", depth, aseq1[0], ist, ien );
//	fprintf( stderr, "====(f) aseq2[0] (%d) = %s (%d-%d)\n", depth, aseq2[0], jst, jen );

	len = strlen( mseq1[0] );
//	fprintf( stderr, "len = %d\n", len );
	l = jmid - jumpj - 1;
//	fprintf( stderr, "l=%d\n", l );
	if( l > 0 )
	{
//		for( i=0; i<l; i++ ) gaps[i] = '-'; gaps[i] = 0;
		for( i=0; i<l; i++ ) gaps[i] = *newgapstr; gaps[i] = 0;
		for( i=0; i<icyc; i++ ) 
		{
			strcat( mseq1[i], gaps );
			mseq1[i][len+l] = 0;
		}
		for( j=0; j<jcyc; j++ )
		{
			strncat( mseq2[j], seq2[j]+jst+jumpj+1, l );
			mseq2[j][len+l] = 0;
		}
//		fprintf( stderr, "penalizing (2) .. %f(%d), %f(%d)\n", ogcp2[jumpj+1], jumpj+1, fgcp2[jmid-1], jmid-1 );
		value +=  ( ogcp2[jumpj+1] + fgcp2[jmid-1] );
//		value += fpenalty;
	}
	len = strlen( mseq1[0] );
	l = imid - jumpi - 1;
//	fprintf( stderr, "l=%d\n", l );
	if( l > 0 )
	{
//		for( i=0; i<l; i++ ) gaps[i] = '-'; gaps[i] = 0;
		for( i=0; i<l; i++ ) gaps[i] = *newgapstr; gaps[i] = 0;
		for( i=0; i<icyc; i++ )
		{
			strncat( mseq1[i], seq1[i]+ist+jumpi+1, l );
			mseq1[i][len+l] = 0;
		}
		for( j=0; j<jcyc; j++ ) 
		{
			strcat( mseq2[j], gaps );
			mseq2[j][len+l] = 0;
		}

//		for( i=0; i<lgth1; i++ ) fprintf( stderr, "ogcp1[%d] = %f\n", i, ogcp1[i] );
//		for( i=0; i<lgth1; i++ ) fprintf( stderr, "fgcp1[%d] = %f\n", i, fgcp1[i] );


//		fprintf( stderr, "penalizing (1) .. ogcp1[%d] = %f, fgcp1[%d] = %f\n", jumpi+1, ogcp1[jumpi+1], imid-1, fgcp1[imid-1] );
		value += ( ogcp1[jumpi+1] + fgcp1[imid-1] );
//		value += fpenalty;
	}
#if 0
	for( i=0; i<icyc; i++ ) fprintf( stderr, "after gapfill mseq1[%d]=%s\n", i, mseq1[i] );
	for( i=0; i<jcyc; i++ ) fprintf( stderr, "after gapfill mseq2[%d]=%s\n", i, mseq2[i] );
#endif

//	fprintf( stderr, "==== calling myself (second)\n" );

#if MEMSAVE
	alnlen = strlen( aseq1[0] );
	for( i=0; i<icyc; i++ ) aseq1[i] += alnlen;
	for( i=0; i<jcyc; i++ ) aseq2[i] += alnlen;
#endif

	value += MSalignmm_rec_variousdist( matrices, icyc, jcyc, seq1, seq2, cpmx1s, cpmx2s, ist+imid, ien, jst+jmid, jen, alloclen, fulllen1, fulllen2, aseq1, aseq2, depth, gapinfo, NULL, 0, NULL, headgp, tailgp, headgapfreq1_g, headgapfreq2_g ); // chudan mada
#if 0
		fprintf( stderr, "aseq1[0] = %s\n", aseq1[0] );
		fprintf( stderr, "aseq2[0] = %s\n", aseq2[0] );
#endif



#if DEBUG
	if( value - maxwm > 1 || maxwm - value > 1 )
	{
		fprintf( stderr, "WARNING value  = %f, but maxwm = %f\n", value, maxwm );
		for( i=0; i<icyc; i++ )
		{
			fprintf( stderr, ">1-%d\n%s\n", i, mseq1[i] );
			fprintf( stderr, "%s\n", aseq1[i] );
		}
		for( i=0; i<jcyc; i++ )
		{
			fprintf( stderr, ">2-%d\n%s\n", i, mseq2[i] );
			fprintf( stderr, "%s\n", aseq2[i] );
		}

//		exit( 1 );
	}
	else
	{
		fprintf( stderr, "value = %.0f, maxwm = %.0f -> ok\n", value, maxwm );
	}
#endif

#if MEMSAVE
#else
	for( i=0; i<icyc; i++ ) strcat( mseq1[i], aseq1[i] );
	for( i=0; i<jcyc; i++ ) strcat( mseq2[i], aseq2[i] );
#endif

//	fprintf( stderr, "====(s) aseq1[0] (%d) = \n%s\n (a%d-a%d)\n", depth, mseq1[0], ist, ien );
//	fprintf( stderr, "====(s) aseq2[0] (%d) = \n%s\n (b%d-b%d)\n", depth, mseq2[0], jst, jen );

	freearrays_rec2( gaps, aseq1, aseq2 );

#if 0
	if( seqlen( seq1[0] ) != nglen1 )
	{
		fprintf( stderr, "bug! hairetsu ga kowareta! (nglen1) seqlen(seq1[0])=%d but nglen1=%d\n", seqlen( seq1[0] ), nglen1 );
		fprintf( stderr, "seq1[0] = %s\n", seq1[0] );
		exit( 1 );
	}
	else
		fprintf( stderr, "nglen1 is ok in _rec\n" );
	if( seqlen( seq2[0] ) != nglen2 )
	{
		fprintf( stderr, "bug! hairetsu ga kowareta! (nglen2) seqlen(seq2[0])=%d but nglen2=%d\n", seqlen( seq2[0] ), nglen2 );
		exit( 1 );
	}
	else
		fprintf( stderr, "nglen2 is ok in _rec\n" );
#endif

	return( value );
}

static void freearrays_variousdist( 
	float *ogcp1, 
	float *ogcp2, 
	float *fgcp1, 
	float *fgcp2, 
	float ***cpmx1s,
	float ***cpmx2s, 
	float *gapfreq1f,
	float *gapfreq2f,
	float **gapinfo, 
	char **mseq1, 
	char **mseq2 
)
{
	FreeFloatVec( ogcp1 );
	FreeFloatVec( ogcp2 );
	FreeFloatVec( fgcp1 );
	FreeFloatVec( fgcp2 );
	FreeFloatCub( cpmx1s );
	FreeFloatCub( cpmx2s );
	FreeFloatVec( gapfreq1f );
	FreeFloatVec( gapfreq2f );
	free( (void *)gapinfo );

	FreeCharMtx( mseq1 );
	FreeCharMtx( mseq2 );
}


float MSalignmm_variousdist( double **pairoffset, double ***matrices, double **dummy_mtx, char **seq1, char **seq2, double *eff1, double *eff2, double **eff1s, double **eff2s, int icyc, int jcyc, int alloclen, char *sgap1, char *sgap2, char *egap1, char *egap2, int *chudanpt, int chudanref, int *chudanres, int headgp, int tailgp )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{
//	fprintf( stderr, "IN MSalignmm_variousdist\n" );
//	int k;
	int i, j, c;
	int ll1, ll2;
	int lgth1, lgth2;
	float wm = 0.0;   /* int ?????? */
	char **mseq1;
	char **mseq2;
	float *ogcp1;
	float *ogcp2;
	float *fgcp1;
	float *fgcp2;
	float ***cpmx1s;
	float ***cpmx2s;
	float **gapinfo;
	float fpenalty = (float)penalty;
	float *gapfreq1f;
	float *gapfreq2f;
	int nglen1, nglen2;
	float headgapfreq1;
	float headgapfreq2;

#if 0
	fprintf( stderr, "eff in SA+++align\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "eff1[%d] = %f\n", i, eff1[i] );
#endif

	nglen1 = seqlen( seq1[0] );
	nglen2 = seqlen( seq2[0] );

#if 0
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "seq1[%d] at root = %s\n", i, seq1[i] );
	for( j=0; j<jcyc; j++ ) fprintf( stderr, "seq2[%d] at root = %s\n", j, seq2[j] );
	fprintf( stderr, "\n" );
#endif

	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );

	ll1 = ( (int)(lgth1) ) + 100;
	ll2 = ( (int)(lgth2) ) + 100;

	mseq1 = AllocateCharMtx( icyc, ll1+ll2 );
	mseq2 = AllocateCharMtx( jcyc, ll1+ll2 );

	gapinfo = AllocateFloatMtx( 6, 0 );
	ogcp1 = AllocateFloatVec( ll1+2 );
	ogcp2 = AllocateFloatVec( ll2+2 );
	fgcp1 = AllocateFloatVec( ll1+2 );
	fgcp2 = AllocateFloatVec( ll2+2 );


	cpmx1s = AllocateFloatCub( maxdistclass, ll1+2, nalphabets+1 );
	cpmx2s = AllocateFloatCub( maxdistclass, ll2+2, nalphabets+1 );

	gapfreq1f = AllocateFloatVec( ll1+2 ); // must be filled with 0.0
	gapfreq2f = AllocateFloatVec( ll2+2 ); // must be filled with 0.0

	for( i=0; i<icyc; i++ ) 
	{
		if( strlen( seq1[i] ) != lgth1 )
		{
			fprintf( stderr, "i = %d / %d\n", i, icyc );
			fprintf( stderr, "bug! hairetsu ga kowareta!\n" );
			exit( 1 );
		}
	}
	for( j=0; j<jcyc; j++ )
	{
		if( strlen( seq2[j] ) != lgth2 )
		{
			fprintf( stderr, "j = %d / %d\n", j, jcyc );
			fprintf( stderr, "bug! hairetsu ga kowareta!\n" );
			exit( 1 );
		}
	}

	for( c=0; c<maxdistclass; c++ )
	{
		MScpmx_calc_new( seq1, cpmx1s[c], eff1s[c], lgth1, icyc );
		MScpmx_calc_new( seq2, cpmx2s[c], eff2s[c], lgth2, jcyc );
	}


#if 1

	if( sgap1 )
	{
		new_OpeningGapCount( ogcp1, icyc, seq1, eff1, lgth1, sgap1 );
		new_OpeningGapCount( ogcp2, jcyc, seq2, eff2, lgth2, sgap2 );
		new_FinalGapCount( fgcp1, icyc, seq1, eff1, lgth1, egap2 );
		new_FinalGapCount( fgcp2, jcyc, seq2, eff2, lgth2, egap2 );
		outgapcount( &headgapfreq1, icyc, sgap1, eff1 );
		outgapcount( &headgapfreq2, jcyc, sgap2, eff2 );
		outgapcount( gapfreq1f+lgth1, icyc, egap1, eff1 );
		outgapcount( gapfreq2f+lgth2, jcyc, egap2, eff2 );
	}
	else
	{
		st_OpeningGapCount( ogcp1, icyc, seq1, eff1, lgth1 );
		st_OpeningGapCount( ogcp2, jcyc, seq2, eff2, lgth2 );
		st_FinalGapCount( fgcp1, icyc, seq1, eff1, lgth1 );
		st_FinalGapCount( fgcp2, jcyc, seq2, eff2, lgth2 );
		headgapfreq1 = 0.0;
		headgapfreq2 = 0.0;
		gapfreq1f[lgth1] = 0.0;
		gapfreq2f[lgth2] = 0.0;
	}

	if( legacygapcost == 0 )
	{
		gapcountf( gapfreq1f, seq1, icyc, eff1, lgth1 );
		gapcountf( gapfreq2f, seq2, jcyc, eff2, lgth2 );
		for( i=0; i<lgth1+1; i++ ) gapfreq1f[i] = 1.0 - gapfreq1f[i];
		for( i=0; i<lgth2+1; i++ ) gapfreq2f[i] = 1.0 - gapfreq2f[i];
		headgapfreq1 = 1.0 - headgapfreq1;
		headgapfreq2 = 1.0 - headgapfreq2;
	}
	else
	{
		for( i=0; i<lgth1+1; i++ ) gapfreq1f[i] = 1.0;
		for( i=0; i<lgth2+1; i++ ) gapfreq2f[i] = 1.0;
		headgapfreq1 = 1.0;
		headgapfreq2 = 1.0;
	}

#if 1
	for( i=0; i<lgth1; i++ ) 
	{
		ogcp1[i] = 0.5 * ( 1.0 - ogcp1[i] ) * fpenalty * ( gapfreq1f[i] );
		fgcp1[i] = 0.5 * ( 1.0 - fgcp1[i] ) * fpenalty * ( gapfreq1f[i] );
//		fprintf( stderr, "fgcp1[%d] = %f\n", i, fgcp1[i] );
	}
	for( i=0; i<lgth2; i++ ) 
	{
		ogcp2[i] = 0.5 * ( 1.0 - ogcp2[i] ) * fpenalty * ( gapfreq2f[i] );
		fgcp2[i] = 0.5 * ( 1.0 - fgcp2[i] ) * fpenalty * ( gapfreq2f[i] );
//		fprintf( stderr, "fgcp2[%d] = %f\n", i, fgcp2[i] );
	}
#else
	for( i=0; i<lgth1; i++ ) 
	{
		ogcp1[i] = 0.5 * fpenalty;
		fgcp1[i] = 0.5 * fpenalty;
	}
	for( i=0; i<lgth2; i++ ) 
	{
		ogcp2[i] = 0.5 * fpenalty;
		fgcp2[i] = 0.5 * fpenalty;
	}
#endif

	gapinfo[0] = ogcp1;
	gapinfo[1] = fgcp1;
	gapinfo[2] = ogcp2;
	gapinfo[3] = fgcp2;
	gapinfo[4] = gapfreq1f;
	gapinfo[5] = gapfreq2f;
#endif

#if 0
	fprintf( stdout, "in MSalignmm.c\n" );
	for( i=0; i<icyc; i++ )
	{
		fprintf( stdout, ">%d of GROUP1\n", i );
		fprintf( stdout, "%s\n", seq1[i] );
	}
	for( i=0; i<jcyc; i++ )
	{
		fprintf( stdout, ">%d of GROUP2\n", i );
		fprintf( stdout, "%s\n", seq2[i] );
	}
	fflush( stdout );
#endif

	wm = MSalignmm_rec_variousdist( matrices, icyc, jcyc, seq1, seq2, cpmx1s, cpmx2s, 0, lgth1-1, 0, lgth2-1, alloclen, lgth1, lgth2, mseq1, mseq2, 0, gapinfo, chudanpt, chudanref, chudanres, headgp, tailgp, headgapfreq1, headgapfreq2 );
#ifdef enablemultithread
	if( chudanres && *chudanres ) 
	{
//		fprintf( stderr, "\n\n## CHUUDAN!!! relay\n" );
		*chudanres = 1;
		freearrays_variousdist( ogcp1, ogcp2, fgcp1, fgcp2, cpmx1s, cpmx2s, gapfreq1f, gapfreq2f, gapinfo, mseq1, mseq2 );
		return( -1.0 );
	}
#endif

#if 0
		fprintf( stderr, "\n" );
		fprintf( stderr, " seq1[0] = %s\n", seq1[0] );
		fprintf( stderr, " seq2[0] = %s\n", seq2[0] );
		fprintf( stderr, "mseq1[0] = %s\n", mseq1[0] );
		fprintf( stderr, "mseq2[0] = %s\n", mseq2[0] );
		fprintf( stderr, "\n" );
#endif

//	fprintf( stderr, "wm = %f\n", wm );


	for( i=0; i<icyc; i++ ) strcpy( seq1[i], mseq1[i] );
	for( i=0; i<jcyc; i++ ) strcpy( seq2[i], mseq2[i] );

	if( seqlen( seq1[0] ) != nglen1 )
	{
		fprintf( stderr, "bug! hairetsu ga kowareta! (nglen1) seqlen(seq1[0])=%d but nglen1=%d\n", seqlen( seq1[0] ), nglen1 );
		fprintf( stderr, "seq1[0] = %s\n", seq1[0] );
		exit( 1 );
	}
	if( seqlen( seq2[0] ) != nglen2 )
	{
		fprintf( stderr, "bug! hairetsu ga kowareta! (nglen2) seqlen(seq2[0])=%d but nglen2=%d\n", seqlen( seq2[0] ), nglen2 );
		exit( 1 );
	}


	freearrays_variousdist( ogcp1, ogcp2, fgcp1, fgcp2, cpmx1s, cpmx2s, gapfreq1f, gapfreq2f, gapinfo, mseq1, mseq2 );

	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );
	for( i=0; i<icyc; i++ ) 
	{
		if( strlen( seq1[i] ) != lgth1 )
		{
			fprintf( stderr, "i = %d / %d\n", i, icyc );
			fprintf( stderr, "hairetsu ga kowareta (end of MSalignmm) !\n" );
			exit( 1 );
		}
	}
	for( j=0; j<jcyc; j++ )
	{
		if( strlen( seq2[j] ) != lgth2 )
		{
			fprintf( stderr, "j = %d / %d\n", j, jcyc );
			fprintf( stderr, "hairetsu ga kowareta (end of MSalignmm) !\n" );
			exit( 1 );
		}
	}

#if 0
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, " seq1[i] = %s\n", seq1[i] );
	for( j=0; j<jcyc; j++ ) fprintf( stderr, " seq2[j] = %s\n", seq2[j] );
	fprintf( stderr, "\n" );
#endif

	return( wm );
}
