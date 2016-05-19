
/* 
	tree-dependent iteration   
    algorithm A+ when group-to-group, C when group-to-singleSeqence 
	                 OR
    algorithm A+
*/

#include "mltaln.h"


#define DEBUG 0
#define RECORD 0
#define MCD 0

extern char **seq_g;
extern char **res_g;

static int nwa;


#ifdef enablemultithread
typedef struct _threadarg
{
	int thread_no;
	int *jobposintpt;
	int *ndonept;
	int *ntrypt;
	int *collectingpt;
	int njob;
	int nbranch;
	int maxiter;
	int nkozo;
	int *subgenerationpt;
	float *basegainpt;
	float *gainlist;
	float *tscorelist;
	int *generationofinput;
	char *kozoarivec;	
	char **mastercopy;
	char ***candidates;
	int *generationofmastercopypt;
	int *branchtable;
	RNApair ***singlerna;
	LocalHom **localhomtable;
	int alloclen;
	Node *stopol;
	int ***topol;
//	double **len;
	float **tscorehistory_detail;
	int *finishpt;
	int **skipthisbranch;
	double **distmtx;
	pthread_mutex_t *mutex;
	pthread_cond_t *collection_end;
	pthread_cond_t *collection_start;
} threadarg_t;
#endif

#if 1
static void shuffle( int *arr, int n )
{
	int i;
	int x;
	int b;

	for( i=1; i<n; i++ )
	{
		x = rand() % (i+1);
		if( x != i )
		{
			b = arr[i];
			arr[i] = arr[x];
			arr[x] = b;
		}
	}
}
#endif


static void makescoringmatrices( double ***matrices, double **originalmtx )
{
	int c;
	float rep;
	for( c=0; c<maxdistclass; c++ )
	{
		rep = (double) 2 * c / ndistclass; // rep:0..2
//		fprintf( stderr, "rep = %f\n", rep );
		makedynamicmtx( matrices[c], originalmtx, rep * 0.5 ); // upgma ni awaseru node, 0..1
//		fprintf( stderr, "c=%d, score for %c-%c = %f\n", c, 'W', 'W', matrices[c][amino_n['W']][amino_n['W']] );
	}
}

static void classifypairs( int n1, double **eff1s, double *eff1, int n2, double **eff2s, double *eff2, double **smalldistmtx, int **matnum, int maxdistclass )
{
	int i, j, c;
	for( c=0; c<maxdistclass; c++ ) 
	{
		for( i=0; i<n1; i++ ) eff1s[c][i] = 0.0;
		for( j=0; j<n2; j++ ) eff2s[c][j] = 0.0;
	}
	
//	fprintf( stderr, "\n" );
	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
	{
		c = (int)( smalldistmtx[i][j] / 2.0 * ndistclass ); // dist:0..2
//		if( c >= ndistclass ) c = ndistclass-1;
		if( c >= maxdistclass ) c = maxdistclass-1;
//		fprintf( stderr, "pair %d-%d (%f), dist=%f -> c=%d\n", i, j, eff1[i] * eff2[j], smalldistmtx[i][j], c );
		eff1s[c][i] = 1.0;
		eff2s[c][j] = 1.0;
		matnum[i][j] = c;
	}
	for( c=0; c<maxdistclass; c++ ) for( i=0; i<n1; i++ ) eff1s[c][i] *= eff1[i];
	for( c=0; c<maxdistclass; c++ ) for( i=0; i<n2; i++ ) eff2s[c][i] *= eff2[i];
#if 0
	double totaleff;
	totaleff = 0.0; for( c=0; c<maxdistclass; c++ ) for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ ) totaleff += eff1s[c][i] * eff2s[c][j];
	fprintf( stderr, "totaleff1s-2s  = %f\n", totaleff );
	totaleff = 0.0; for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ ) totaleff += eff1[i] * eff2[j];
	fprintf( stderr, "totaleff1-2  = %f\n", totaleff );

	totaleff = 0.0; for( c=0; c<maxdistclass; c++ ) for( i=0; i<n1; i++ ) totaleff += eff1s[c][i]; 
	fprintf( stderr, "totaleff1s = %f\n", totaleff );
	totaleff = 0.0; for( c=0; c<maxdistclass; c++ ) for( i=0; i<n2; i++ ) totaleff += eff2s[c][i]; 
	fprintf( stderr, "totaleff2s = %f\n", totaleff );
	totaleff = 0.0; for( i=0; i<n1; i++ ) totaleff += eff1[i]; 
	fprintf( stderr, "totaleff1 = %f\n", totaleff );
	totaleff = 0.0; for( i=0; i<n2; i++ ) totaleff += eff2[i]; 
	fprintf( stderr, "totaleff2 = %f\n", totaleff );
	{
//		for( i=0; i<n1; i++ ) fprintf( stderr, "eff1s[%d][%d] = %f\n", c, i, eff1s[c][i] );
//		for( i=0; i<n2; i++ ) fprintf( stderr, "eff2s[%d][%d] = %f\n", c, i, eff2s[c][i] );
//		fprintf( stderr, "\n" );
	}
#endif
}

static void Writeoption2( FILE *fp, int cycle, double cut )
{
	fprintf( fp, "%dth cycle\n", cycle );
    fprintf( fp, "marginal score to search : current score * (100-%d) / 100\n", (int)cut );
}

static void Writeoptions( FILE *fp )
{
	fprintf( fp, "Tree-dependent-iteration\n" );
    if( scoremtx == 1 )
        fprintf( fp, "Blosum %d\n", nblosum );
    else if( scoremtx == -1 )
        fprintf( fp, "DNA\n" );
    else if( scoremtx == 2 )
        fprintf( fp, "Miyata-Yasunaga\n" );
	else
		fprintf( fp, "JTT %dPAM\n", pamN );

	if( scoremtx == 0 || scoremtx == 1 )
    	fprintf( fp, "Gap Penalty = %+5.3f, %5.2f, %+5.3f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );
	else
		fprintf( fp, "Gap Penalty = %+5.3f\n", (double)penalty/1000 );
		

    if( scmtd == 3 )
        fprintf( fp, "score of rnd or sco\n" );
    else if( scmtd == 4 )
        fprintf( fp, "score = sigma( score for a pair of homologous amino acids ) / ( number of amino acids pairs )\n" );
    else if( scmtd == 5 )
        fprintf( fp, "score : SP\n" );
    if( mix )
        fprintf( fp, "?\n" );
    else
    { 
        if( weight == 2 )
            fprintf( fp, "weighted rationale-1,  geta2 = %f\n", geta2 );
        else if( weight == 3 )
            fprintf( fp, "weighted like ClustalW," );
        else if( weight == 4 )
            fprintf( fp, "weighted rationale-2,  geta2 = %f\n", geta2 );
        else
            fprintf( fp, "unweighted\n" );
    }
    if( weight && utree )
        fprintf( fp, "using tree defined by the file hat2.\n" );
    if( weight && !utree )
        fprintf( fp, "using temporary tree.\n" );

	if( treemethod == 'n' )
		fprintf( fp, "Tree is calculated with Neighbor-Joining Method.\n" );
	else if( treemethod == 'q' )
		fprintf( fp, "Tree is calculated with Minimum linkage.\n" );
	else if( treemethod == 'X' )
		fprintf( fp, "Tree is calculated with simplified UPG Method and UPG Method.\n" );
	else if( treemethod == 'E' )
		fprintf( fp, "Tree is calculated with UPG Method.\n" );
	else
		fprintf( fp, "Tree is calculated with unknown Method.\n" );
		
	if( alg == 'C' ) 
		fprintf( fp, "Algorithm A+ / C\n" );
	else if( alg == 'A' ) 
		fprintf( fp, "Algorithm A+ \n" );
	else if( alg == 'a' ) 
		fprintf( fp, "Algorithm A \n" );
	else 
		fprintf( fp, "Algorithm ? \n" );

	if( use_fft )
	{
		if( scoremtx == -1 )
		{
			fprintf( fp, "Basis : 4 nucleotides\n" );
		}
		else
		{
			if( fftscore )
				fprintf( fp, "Basis : Polarity and Volume\n" );
			else
				fprintf( fp, "Basis : 20 amino acids\n" );
		}
		fprintf( fp, "Threshold   of anchors = %d%%\n", fftThreshold );
		fprintf( fp, "window size of anchors = %dsites\n", fftWinSize );
	}
}

#ifdef enablemultithread

static void freelocalarrays( 
	float *tscorehistory,
	RNApair ***grouprna1, RNApair ***grouprna2,
	RNApair *rnapairboth,
	char *indication1, char *indication2,
	double *effarr, double *effarrforlocalhom, double *effarr1, double *effarr2,
	char **mseq1, char **mseq2,
	char **localcopy,
	int *gapmap1, int *gapmap2,
	double *effarr1_kozo, double *effarr2_kozo, double *effarr_kozo,
	int **memlist,
	char *pairbuf,
	LocalHom *** localhomshrink,
	double **smalldistmtx,
#if MCD
	double **offsetmtx,
#endif
	double ***scoringmatrices,
	double **eff1s, double **eff2s,
	int **whichmtx
)
{
//	fprintf( stderr, "Skipping freelocalarrays\n" );
//	return;
	int i;
	if( commonIP ) FreeIntMtx( commonIP );
	commonIP = NULL;
	Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
	Falign_localhom( NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, NULL, NULL,NULL, 0, NULL );
	if( rnakozo && rnaprediction == 'm' )
	{
		free( grouprna1 ); // nakamimo?
		free( grouprna2 ); // nakamimo?
	}

	free( tscorehistory );
	free( indication1 );
	free( indication2 );
	free( effarr );
	free( effarrforlocalhom );
	free( effarr1 );
	free( effarr2 );
	free( mseq1 );
	free( mseq2 );
	FreeCharMtx( localcopy );
	free( gapmap1 );
	free( gapmap2 );

	free( effarr1_kozo );
	free( effarr2_kozo );
	free( effarr_kozo );

	FreeIntMtx( memlist );
	free( pairbuf );

	if( smalldistmtx ) FreeDoubleMtx( smalldistmtx );
#if MCD
	if( offsetmtx ) FreeDoubleMtx( offsetmtx );
#endif
	if( scoringmatrices ) FreeDoubleCub( scoringmatrices );
	if( eff1s ) FreeDoubleMtx( eff1s );
	if( eff2s ) FreeDoubleMtx( eff2s );
	if( whichmtx ) FreeIntMtx( whichmtx );

	if( rnakozo ) free( rnapairboth );

	if( constraint )
	{
		for( i=0; i<njob; i++)
		{
			free( localhomshrink[i] ); // nakamimo??
		}
		free( localhomshrink );
	}
}


static void *athread( void *arg )
{

	threadarg_t *targ = (threadarg_t *)arg;
	int thread_no = targ->thread_no;
	int njob = targ->njob;
	int nbranch = targ->nbranch;
	int maxiter = targ->maxiter;
	int *ndonept = targ->ndonept;
	int *ntrypt = targ->ntrypt;
	int *collectingpt = targ->collectingpt;
	int *jobposintpt = targ->jobposintpt;
	int nkozo = targ->nkozo;
	float *gainlist = targ->gainlist;
	float *tscorelist = targ->tscorelist;
	int *generationofinput = targ->generationofinput;
	int *subgenerationpt = targ->subgenerationpt;
	float *basegainpt = targ->basegainpt;
	char *kozoarivec = targ->kozoarivec;	
	char **mastercopy = targ->mastercopy;
	char ***candidates = targ->candidates;
	int *generationofmastercopypt = targ->generationofmastercopypt;
	int *branchtable = targ->branchtable;
	RNApair ***singlerna = targ->singlerna;
	LocalHom **localhomtable = targ->localhomtable;
	int alloclen = targ->alloclen;
	Node * stopol = targ->stopol;
	int ***topol = targ->topol;
//	double **len = targ->len;
	float **tscorehistory_detail = targ->tscorehistory_detail;
	int *finishpt = targ->finishpt;
	int **skipthisbranch = targ->skipthisbranch;
	double **distmtx = targ->distmtx;

	int i, k, l, ii, j;
	float gain;
	int iterate;
	int **memlist;
	char *pairbuf;
	int locnjob;
	int s1, s2;
	int clus1, clus2;
	char **localcopy;
	char **mseq1, **mseq2;
	double *effarr, *effarr_kozo; // re-calc 
	double *effarr1, *effarr2, *effarr1_kozo, *effarr2_kozo;
	char *indication1, *indication2;
	int length;
	RNApair ***grouprna1, ***grouprna2;
	RNApair *rnapairboth;
	LocalHom ***localhomshrink;
	int *gapmap1, *gapmap2;
	float tscore, mscore, oimpmatch, impmatch;
	int identity;
	double tmpdouble;
	float naivescore0 = 0, naivescore1;
	double *effarrforlocalhom;
	float *tscorehistory;
	int intdum;
#if 0
	int oscillating;
	int lin, ldf;
#endif
	float maxgain;
	int bestthread;
	int branchpos;
	int subgenerationatfirst;
	double unweightedspscore;
	int myjob;
	int converged2 = 0;
	int chudanres;
	double **smalldistmtx;
#if MCD
	double **offsetmtx = NULL; // check senyou
#endif
	double ***scoringmatrices;
	double **eff1s, **eff2s; 
	int **whichmtx;


	locnjob = njob;

	if( utree == 0 )
	{
		fprintf( stderr, "Dynamic tree is not supported in the multithread version.\n" );
		exit( 1 );
	}
	if( score_check == 2 )
	{
		fprintf( stderr, "Score_check 2 is not supported in the multithread version.\n" );
		exit( 1 );
	}

	if( weight == 2 )
	{
		fprintf( stderr, "Weight 2 is not supported in the multithread version.\n" );
		exit( 1 );
	}
	if( cooling &&  cut > 0.0 )
	{
		fprintf( stderr, "Cooling is not supported in the multithread version.\n" );
		exit( 1 );
	}

	tscorehistory = calloc( maxiter, sizeof( float ) );

	if( rnakozo && rnaprediction == 'm' )
	{
		grouprna1 = (RNApair ***)calloc( njob, sizeof( RNApair ** ) );
		grouprna2 = (RNApair ***)calloc( njob, sizeof( RNApair ** ) );
	}
	else
	{
		grouprna1 = grouprna2 = NULL;
	}

	indication1 = AllocateCharVec( 150 );
	indication2 = AllocateCharVec( 150 );
	effarr = AllocateDoubleVec( locnjob );
	effarrforlocalhom = AllocateDoubleVec( locnjob );
	effarr1 = AllocateDoubleVec( locnjob );
	effarr2 = AllocateDoubleVec( locnjob );
	mseq1 = AllocateCharMtx( locnjob, 0 );
	mseq2 = AllocateCharMtx( locnjob, 0 );
	localcopy = AllocateCharMtx( locnjob, alloclen );
	gapmap1 = AllocateIntVec( alloclen );
	gapmap2 = AllocateIntVec( alloclen );
	if( specificityconsideration != 0 ) 
	{
		smalldistmtx = AllocateDoubleMtx( locnjob, locnjob ); // ookii?
#if MCD
		offsetmtx = AllocateDoubleMtx( locnjob, locnjob ); // ookii?
#endif
		scoringmatrices = AllocateDoubleCub( maxdistclass, nalphabets, nalphabets );
		makescoringmatrices( scoringmatrices, n_dis_consweight_multi );
		eff1s = AllocateDoubleMtx( maxdistclass, locnjob );
		eff2s = AllocateDoubleMtx( maxdistclass, locnjob );
		whichmtx = AllocateIntMtx( locnjob, locnjob );
	}
	else
	{
		smalldistmtx = NULL;
#if MCD
		offsetmtx = NULL;
#endif
		scoringmatrices = NULL;
		eff1s = eff2s = NULL;
		whichmtx = NULL;
	}

	effarr1_kozo = AllocateDoubleVec( locnjob ); // tsuneni allocate suru.
	effarr2_kozo = AllocateDoubleVec( locnjob ); // tsuneni allocate suru.
	effarr_kozo = AllocateDoubleVec( locnjob ); 
	for( i=0; i<locnjob; i++ )
		effarr_kozo[i] = effarr1_kozo[i] = effarr2_kozo[i] = 0.0;

	memlist = AllocateIntMtx( 2, locnjob );
	pairbuf = AllocateCharVec( locnjob );


	if( rnakozo ) rnapairboth = (RNApair *)calloc( alloclen, sizeof( RNApair ) );

	if( constraint )
	{
		localhomshrink = (LocalHom ***)calloc( njob, sizeof( LocalHom ** ) );
		for( i=0; i<njob; i++)
		{
			localhomshrink[i] = (LocalHom **)calloc( njob, sizeof( LocalHom * ) );
		}
	}


	if( thread_no == 0 )
	{
		*ntrypt = 0;
		srand( randomseed );
		*finishpt = 0;
		for( iterate=0; iterate<maxiter; iterate++ ) 
		{
			pthread_mutex_lock( targ->mutex );

			if( *collectingpt == 1 )
			{
				*collectingpt = 0;
				*generationofmastercopypt = iterate;
				*subgenerationpt = 0;
				*basegainpt = 0.0;
				*ndonept = 0;
				*jobposintpt = 0;
				for( i=0; i<nwa; i++ ) gainlist[i] = 0;
				for( i=0; i<nwa; i++ ) tscorelist[i] = 0.0;
				for( i=0; i<nbranch; i++ ) generationofinput[i] = -1;
				if( parallelizationstrategy != BESTFIRST && randomseed != 0 ) shuffle( branchtable, nbranch );
				pthread_cond_broadcast( targ->collection_end );
				pthread_mutex_unlock( targ->mutex );
			}
			else
			{
				pthread_cond_broadcast( targ->collection_end );
				pthread_mutex_unlock( targ->mutex );
				freelocalarrays
				( 
					tscorehistory,
					grouprna1, grouprna2,
					rnapairboth,
					indication1, indication2,
					effarr, effarrforlocalhom, effarr1, effarr2,
					mseq1, mseq2,
					localcopy,
					gapmap1, gapmap2,
					effarr1_kozo, effarr2_kozo, effarr_kozo,
					memlist, pairbuf,
					localhomshrink,
					smalldistmtx,
#if MCD
					offsetmtx,
#endif
					scoringmatrices,
					eff1s, eff2s,
					whichmtx
				);
//				return( NULL );
				pthread_exit( NULL );
			}

			pthread_mutex_lock( targ->mutex );
			while( *ndonept < nbranch )
				pthread_cond_wait( targ->collection_start, targ->mutex );
			pthread_mutex_unlock( targ->mutex );
//			fprintf( stderr, "Thread 0 got a signal, *collectionpt = %d\n", *collectingpt );

/* 
			Hoka no thread ga keisan
*/

			pthread_mutex_lock( targ->mutex );
			*collectingpt = 1; // chofuku

#if 0
			for( i=0; i<nbranch; i++ )
			{
				if( generationofinput[i] != iterate )
				{
					fprintf( stderr, "Error! generationofinput[%d] = %d, but iterate=%d\n", i, generationofinput[i], iterate );
					exit( 1 );

				}
			}
#endif
	
			maxgain = gainlist[1];
			bestthread = 1;
			for( i=2; i<nwa; i++ )
			{
				if( gainlist[i] > maxgain )
				{
					maxgain = gainlist[i];
					bestthread = i;
				}
			}
	
			if( maxgain > 0.0 )
			{
//				fprintf( stderr, "\nGain = %f\n", maxgain );
//				fprintf( stderr, "best gain = %f by thread %d\n", gainlist[bestthread], bestthread );
//				fprintf( stderr, "tscorelist[best] = %f by thread %d\n", tscorelist[bestthread], bestthread );
				if( parallelizationstrategy == BESTFIRST )
				{
					for( i=0; i<locnjob; i++ ) strcpy( mastercopy[i], candidates[bestthread][i] );
					if( scoreout )
					{
						unweightedspscore = plainscore( locnjob, mastercopy );
						fprintf( stderr, "\nSCORE %d = %.0f, ", iterate * nbranch, unweightedspscore );
						fprintf( stderr, "SCORE / residue = %f", unweightedspscore / ( locnjob * strlen( mastercopy[0] ) ) );
						if( weight || constraint ) fprintf( stderr, " (differs from the objective score)" );
						fprintf( stderr, "\n" );
					}
				}
#if 1
//				fprintf( stderr, "gain(%d, by %d) = %f\n", iterate, bestthread, maxgain );
				for( i=iterate-1; i>0; i-- )
				{
//					if( iterate-i < 15 ) fprintf( stderr, "hist[%d] = %f\n", i, tscorehistory[i] );
					if( tscorehistory[i] == tscorelist[bestthread] )
					{
						fprintf( stderr, "\nOscillating? %f == %f\n", tscorehistory[i], tscorelist[bestthread] );
						*collectingpt = -1;
						break;
					}
				}
				tscorehistory[iterate] = tscorelist[bestthread];
#endif
			}
			else
			{
				fprintf( stderr, "\nConverged.\n" );
				*collectingpt = -1;
//				pthread_cond_broadcast( targ->collection_end );
//				pthread_mutex_unlock( targ->mutex );
//				freelocalarrays();
//				return( NULL );
//				pthread_exit( NULL );
			}

#if 1
			if( *finishpt )
			{
				fprintf( stderr, "\nConverged2.\n" );
				*collectingpt = -1;
			}
#endif

			pthread_mutex_unlock( targ->mutex );
		}
		pthread_mutex_lock( targ->mutex );
		fprintf( stderr, "\nReached %d\n", maxiter );
		*collectingpt = -1;
		pthread_cond_broadcast( targ->collection_end );
		pthread_mutex_unlock( targ->mutex );
		freelocalarrays
		( 
			tscorehistory,
			grouprna1, grouprna2,
			rnapairboth,
			indication1, indication2,
			effarr, effarrforlocalhom, effarr1, effarr2,
			mseq1, mseq2,
			localcopy,
			gapmap1, gapmap2,
			effarr1_kozo, effarr2_kozo, effarr_kozo,
			memlist, pairbuf,
			localhomshrink,
			smalldistmtx,
#if MCD
			offsetmtx,
#endif
			scoringmatrices,
			eff1s, eff2s,
			whichmtx
		);
		return( NULL );
		pthread_exit( NULL );
	}
	else
	{
		while( 1 )
		{
#if 0
			if( iterate % 2 == 0 ) 
			{
				lin = 0; ldf = +1;
			}
			else
			{
				lin = locnjob - 2; ldf = -1;
			}	
			for( l=lin; l < locnjob-1 && l >= 0 ; l+=ldf )
				for( k=0; k<2; k++ ) 
#endif

			pthread_mutex_lock( targ->mutex );
			while( *collectingpt > 0 )
				pthread_cond_wait( targ->collection_end, targ->mutex );
			if( *collectingpt == -1 )
			{
				pthread_mutex_unlock( targ->mutex );
				freelocalarrays
				( 
					tscorehistory,
					grouprna1, grouprna2,
					rnapairboth,
					indication1, indication2,
					effarr, effarrforlocalhom, effarr1, effarr2,
					mseq1, mseq2,
					localcopy,
					gapmap1, gapmap2,
					effarr1_kozo, effarr2_kozo, effarr_kozo,
					memlist, pairbuf,
					localhomshrink,
					smalldistmtx,
#if MCD
					offsetmtx,
#endif
					scoringmatrices,
					eff1s, eff2s,
					whichmtx
				);
				return( NULL );
				pthread_exit( NULL );
			}
//			pthread_mutex_unlock( targ->mutex );


//			pthread_mutex_lock( targ->mutex );
			if( *jobposintpt == nbranch )
			{
				if( *collectingpt != -1 ) *collectingpt = 1; // chofuku
				pthread_mutex_unlock( targ->mutex );
				continue;
			}
//			fprintf( stderr, "JOB jobposintpt=%d\n", *jobposintpt );
			myjob = branchtable[*jobposintpt];
			l = myjob / 2;
			if( l == locnjob-2 ) k = 1;
			else k = myjob - l * 2;
//			fprintf( stderr, "JOB l=%d, k=%d\n", l, k );


			branchpos = myjob;
			(*jobposintpt)++;
			iterate = *generationofmastercopypt;
			(*ntrypt)++;
			pthread_mutex_unlock( targ->mutex );

//			fprintf( stderr, "\n IRANAI IRANAI *jobposintpt=%d, nbranch = %d\n", *jobposintpt, nbranch );

//			fprintf( stderr, "branchpos = %d (thread %d)\n", branchpos, thread_no );

//			fprintf( stderr, "iterate=%d, l=%d, k=%d (thread %d)\n", iterate, l, k, thread_no );

#if 0
			fprintf( stderr, "STEP %03d-%03d-%d (Thread %d)    ", iterate+1, l+1, k, thread_no );
			fprintf( stderr, "STEP %03d-%03d-%d (thread %d) %s       ", iterate+1, l+1, k, thread_no, use_fft?"\n":"\n" );
#endif

//			for( i=0; i<2; i++ ) for( j=0; j<locnjob; j++ ) pair[i][j] = 0;
			OneClusterAndTheOther_fast( locnjob, memlist[0], memlist[1], &s1, &s2, pairbuf, topol, l, k, smalldistmtx, distmtx );


#if 0
			fprintf( stderr, "STEP%d-%d\n", l, k );
			for( i=0; i<locnjob; i++ ) 
			{
				for( j=0; j<locnjob; j++ ) 
				{
					fprintf( stderr, "%#3d", pair[i][j] );
				}
				fprintf( stderr, "\n" );
			}
#endif
			if( !weight )
			{
				for( i=0; i<locnjob; i++ ) effarr[i] = 1.0;
				if( nkozo )
				{
					for( i=0; i<locnjob; i++ ) 
					{
						if( kozoarivec[i] )
							effarr_kozo[i] = 1.0;
						else
							effarr_kozo[i] = 0.0;
					}
				}
			}
			else if( weight == 4 )
			{
				weightFromABranch( locnjob, effarr, stopol, topol, l, k );
				if( nkozo ) // hitomadu single weight.
				{
					for( i=0; i<locnjob; i++ ) 
					{
						if( kozoarivec[i] ) effarr_kozo[i] = effarr[i]; 
						else effarr_kozo[i] = 0.0;
					}
				}
			}
			else
			{
				fprintf( stderr, "weight error!\n" );
				exit( 1 );
			}

			yarinaoshi:

			pthread_mutex_lock( targ->mutex );
			for( i=0; i<locnjob; i++ ) strcpy( localcopy[i], mastercopy[i] );
			subgenerationatfirst = *subgenerationpt;
			pthread_mutex_unlock( targ->mutex );
			length = strlen( localcopy[0] );

			if( nkozo )
			{
//				double tmptmptmp;
//				tmptmptmp = 0.0;
//				clus1 = conjuctionfortbfast_kozo( &tmptmptmp, pair[0], s1, localcopy, mseq1, effarr1, effarr, effarr1_kozo, effarr_kozo, indication1 );
				clus1 = fastconjuction_noname_kozo( memlist[0], localcopy, mseq1, effarr1, effarr, effarr1_kozo, effarr_kozo, indication1 );
				for( i=0; i<clus1; i++ ) effarr1_kozo[i] *= 1.0; // 0.5 ga sairyo ?
//				tmptmptmp = 0.0;
//				clus2 = conjuctionfortbfast_kozo( &tmptmptmp, pair[1], s2, localcopy, mseq2, effarr2, effarr, effarr2_kozo, effarr_kozo, indication2 );
				clus2 = fastconjuction_noname_kozo( memlist[1], localcopy, mseq2, effarr2, effarr, effarr2_kozo, effarr_kozo, indication2 );
				for( i=0; i<clus2; i++ ) effarr2_kozo[i] *= 1.0; // 0.5 ga sairyo ?

#if 0
				fprintf( stderr, "\ngroup1 = %s\n", indication1 );
				for( i=0; i<clus1; i++ ) fprintf( stderr, "effarr1_kozo[%d], effarr1[]  = %f, %f\n", i, effarr1_kozo[i], effarr1[i] );
				fprintf( stderr, "\ngroup2 = %s\n", indication2 );
				for( i=0; i<clus2; i++ ) fprintf( stderr, "effarr2_kozo[%d], effarr2[]  = %f, %f\n", i, effarr2_kozo[i], effarr2[i] );
#endif
			}
			else
			{
//				clus1 = conjuctionfortbfast( pair[0], s1, localcopy, mseq1, effarr1, effarr, indication1 );
//				clus2 = conjuctionfortbfast( pair[1], s2, localcopy, mseq2, effarr2, effarr, indication2 );
				clus1 = fastconjuction_noname( memlist[0], localcopy, mseq1, effarr1, effarr, indication1 );
				clus2 = fastconjuction_noname( memlist[1], localcopy, mseq2, effarr2, effarr, indication2 );
			}

	        if( rnakozo && rnaprediction == 'm' )
	        {       
//				makegrouprnait( grouprna1, singlerna, pair[0], s1 );
//				makegrouprnait( grouprna2, singlerna, pair[1], s2 );
				makegrouprna( grouprna1, singlerna, memlist[0] );
				makegrouprna( grouprna2, singlerna, memlist[1] );
	        }       

			if( smalldistmtx )
			{
				classifypairs( clus1, eff1s, effarr1, clus2, eff2s, effarr2, smalldistmtx, whichmtx, maxdistclass );
#if MCD
				for( i=0; i<clus1; i++ )for( j=0; j<clus2; j++ )
				{
//					offsetmtx[i][j] = dist2offset( smalldistmtx[i][j] );
//					int digit = (int)( smalldistmtx[i][j] / 2.0 * (double)ndistclass );
//					fprintf( stderr, "digit = %d\n", digit );
//					offsetmtx[i][j] = dist2offset( (double)digit * 2.0 / ndistclass );

					offsetmtx[i][j] = dist2offset( smalldistmtx[i][j] );
//					fprintf( stderr, "offsetmtx[%d][%d] = %f -> %f\n", i, j, smalldistmtx[i][j], offsetmtx[i][j] );
				}
#endif
			}

			if( score_check == 2 )
			{
				fprintf( stderr, "Score_check 2 is not supported in the multithread version.\n" );
				exit( 1 );
			}
			else if( score_check )
			{
				if( constraint )
				{
					if( RNAscoremtx == 'r' )
						intergroup_score_gapnomi( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble ); // gappick mae denaito dame
					else
						intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble ); // gappick mae denaito dame

//					shrinklocalhom( pair, int s1, int s2, localhomtable, localhomshrink );
//					msshrinklocalhom( pair[0], pair[1], s1, s2, localhomtable, localhomshrink );
					msshrinklocalhom_fast( memlist[0], memlist[1], localhomtable, localhomshrink );
					oimpmatch = 0.0;
					if( use_fft )
					{
						if( alg == 'Q' )
						{
							part_imp_match_init_strictQ( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
							if( rnakozo ) part_imp_rnaQ( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, gapmap1, gapmap2, NULL );
							for(  i=length-1; i>=0; i-- )
							{
								oimpmatch += part_imp_match_out_scQ( i, i );
//								fprintf( stderr, "#### i=%d, initial impmatch = %f seq1 = %c, seq2 = %c\n", i, oimpmatch, mseq1[0][i], mseq2[0][i] );
							}
						}
						else
						{
							part_imp_match_init_strict( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );
							if( rnakozo ) part_imp_rna( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, gapmap1, gapmap2, NULL );
							for(  i=length-1; i>=0; i-- )
							{
								oimpmatch += part_imp_match_out_sc( i, i );
//								fprintf( stderr, "#### i=%d, initial impmatch = %f seq1 = %c, seq2 = %c\n", i, oimpmatch, mseq1[0][i], mseq2[0][i] );
							}
						}
//						fprintf( stderr, "otmpmatch = %f\n", oimpmatch );
					}
					else
					{
						if( alg == 'Q' )
						{
							imp_match_init_strictQ( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
							if( rnakozo ) imp_rnaQ( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, gapmap1, gapmap2, NULL );

							for(  i=length-1; i>=0; i-- )
							{
								oimpmatch += imp_match_out_scQ( i, i );
//								fprintf( stderr, "#### i=%d, initial impmatch = %f\n", i, oimpmatch );
							}
						}
						else
						{
							imp_match_init_strict( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );

							fprintf( stderr, "not supported\n" );
							exit( 1 );

							for(  i=length-1; i>=0; i-- )
							{
								oimpmatch += imp_match_out_sc( i, i );
//								fprintf( stderr, "#### i=%d, initial impmatch = %f seq1 = %c, seq2 = %c\n", i, oimpmatch, mseq1[0][i], mseq2[0][i] );
							}
						}
//						fprintf( stderr, "otmpmatch = %f\n", oimpmatch );
					}
//					fprintf( stderr, "#### initial impmatch = %f\n", oimpmatch );
				}
				else
				{
					if( RNAscoremtx == 'r' )
						intergroup_score_gapnomi( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble ); // gappick mae denaito dame
					else
					{
						if( smalldistmtx )
#if 1
							intergroup_score_multimtx( whichmtx, scoringmatrices, mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble );
#else
							intergroup_score_dynmtx( smalldistmtx, amino_dis, mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble ); // gappick mae denaito dame
#endif
						else
							intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble ); // gappick mae denaito dame
					}
					oimpmatch = 0.0;
				}


//				fprintf( stderr, "#### tmpdouble = %f\n", tmpdouble );
				mscore = (double)oimpmatch + tmpdouble;
			}
			else
			{
				fprintf( stderr, "score_check = %d\n", score_check );
				fprintf( stderr, "Not supported.  Please add --threadit 0 to disable the multithreading in the iterative refinement calculation.\n" );
				exit( 1 );
			}


//			if( rnakozo ) foldalignedrna( clus1, clus2, mseq1, mseq2, effarr1, effarr2, rnapairboth );

//			if( !use_fft && !rnakozo )
			if( !use_fft )
			{
				commongappick_record( clus1, mseq1, gapmap1 );
				commongappick_record( clus2, mseq2, gapmap2 );
			}

#if 0
			fprintf( stderr, "##### mscore = %f\n", mscore );
#endif

#if DEBUG
			if( !devide )
			{
	       		fprintf( trap_g, "\n%d-%d-%d\n", iterate+1, l+1, k );
      		    	fprintf( trap_g, "group1 = %s\n", indication1 );
      		    	fprintf( trap_g, "group2 = %s\n", indication2 );
				fflush( trap_g );
			}

#endif
#if 0
			printf( "STEP %d-%d-%d\n", iterate, l, k );
			for( i=0; i<clus2; i++ ) printf( "%f ", effarr2[i] );
			printf( "\n" );
#endif
			if( !skipthisbranch[l][k] )
			{

				if( constraint == 2 )
				{
					if( use_fft )
					{
//						if( alg == 'Q' )
//							part_imp_match_init_strictQ( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
//						else
//							part_imp_match_init_strict( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );
						chudanres = 0;
						Falign_localhom( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, localhomshrink, &impmatch, gapmap1, gapmap2, subgenerationpt, subgenerationatfirst, &chudanres );
//						fprintf( stderr, "##### impmatch = %f\n", impmatch );
						if( chudanres && parallelizationstrategy == BAATARI2 )
						{
//							fprintf( stderr, "#### yarinaoshi!!! INS-i\n" );
							goto yarinaoshi;
						}
					}
					else
					{
						fprintf( stderr, "Not supported\n" );
						exit( 1 );
					}
				}
				else if( use_fft )
				{
					float totalwm;
					chudanres = 0;
//					totalwm = Falign( smalldistmtx, NULL, n_dis_consweight_multi, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, alloclen, &intdum, subgenerationpt, subgenerationatfirst, &chudanres );
#if MCD
					totalwm = Falign( offsetmtx, scoringmatrices, n_dis_consweight_multi, mseq1, mseq2, effarr1, effarr2, eff1s, eff2s, clus1, clus2, alloclen, &intdum, subgenerationpt, subgenerationatfirst, &chudanres );
#else
					totalwm = Falign( NULL, scoringmatrices, n_dis_consweight_multi, mseq1, mseq2, effarr1, effarr2, eff1s, eff2s, clus1, clus2, alloclen, &intdum, subgenerationpt, subgenerationatfirst, &chudanres );
#endif
					if( chudanres && parallelizationstrategy == BAATARI2 )
					{
//						fprintf( stderr, "#### yarinaoshi!!! FFT-NS-i\n" );
						goto yarinaoshi;
					}

				}
				else
				{
					fprintf( stderr, "\n\nUnexpected error.  Please contact kazutaka.katoh@aist.go.jp\n\n\n" );
					exit( 1 );
				}
//				fprintf( stderr, "## impmatch = %f\n", impmatch );

#if 1
				if( parallelizationstrategy == BAATARI2 && *subgenerationpt != subgenerationatfirst )
				{
//					fprintf( stderr, "\nYarinaoshi2!! (Thread %d)\n", thread_no );
					goto yarinaoshi;
				}
#endif
							
				identity = !strcmp( localcopy[s1], mastercopy[s1] );
				identity *= !strcmp( localcopy[s2], mastercopy[s2] );
				fprintf( stderr, "%03d-%04d-%d (thread %4d) identical     \r", iterate+1, *ndonept, k, thread_no );

			}
			else
			{
				identity = 1;
				fprintf( stderr, "%03d-%04d-%d (thread %4d) skip     \r", iterate+1, *ndonept, k, thread_no );
			}


/* Bug?  : idnetitcal but score change when scoreing mtx != JTT  */

			length = strlen( mseq1[0] );

			if( identity )
			{
				tscore = mscore;
			}
			else
			{
				if( score_check )
				{
					if( constraint == 2 )
					{
#if 1
						if( RNAscoremtx == 'r' )
							intergroup_score_gapnomi( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble );
						else
							intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble );
#else
						intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble );
#endif

						tscore = impmatch + tmpdouble;

//						fprintf( stderr, "tmpdouble=%f, impmatch = %f -> %f, tscore = %f\n", tmpdouble, oimpmatch, impmatch, tscore );
					}
					else
					{
						if( smalldistmtx )
#if 1
							intergroup_score_multimtx( whichmtx, scoringmatrices, mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble );
#else
							intergroup_score_dynmtx( smalldistmtx, amino_dis, mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble );
#endif
						else
							intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble );
						tscore = tmpdouble;
					}
//					fprintf( stderr, "#######ii=%d, iterate=%d score = %f -> %f \n", ii, iterate , mscore, tscore );
	#if 0
					for( i=0; i<1; i++ )
						fprintf( stderr, "%s\n", mseq1[i] );
					fprintf( stderr, "+++++++\n" );
					for( i=0; i<1; i++ )
						fprintf( stderr, "%s\n", mseq2[i] );
	#endif

				}
				else
				{
					tscore = mscore + 1.0;
//					tscore = 0.0;
//					fprintf( stderr, "in line 705, tscore=%f\n", tscore );
//					for( i=0; i<length; i++ )
//						tscore = tscore + (double)mseq1[0][i];
//					mscore = tscore - 1.0;
				}

				if( isnan( mscore ) )
				{
					fprintf( stderr, "\n\nmscore became NaN\n" );
					exit( 1 );
				}
				if( isnan( tscore ) )
				{
					fprintf( stderr, "\n\ntscore became NaN\n" );
					exit( 1 );
				}



//				fprintf( stderr, "@@@@@ mscore,tscore = %f,%f\n", mscore, tscore );

#if 1
				if( parallelizationstrategy == BAATARI1 && *subgenerationpt != subgenerationatfirst )
				{
//					fprintf( stderr, "\nYarinaoshi1!! (Thread %d)\n", thread_no );
					goto yarinaoshi;
				}
#endif
				gain = tscore - ( mscore - cut/100.0*mscore );
				if( gain > 0 )
				{
					if( parallelizationstrategy == BESTFIRST )
					{
						if( gain > gainlist[thread_no] ) 
						{
							gainlist[thread_no] = gain;
							for( i=0; i<locnjob; i++ ) strcpy( candidates[thread_no][i], localcopy[i] );
							tscorelist[thread_no] = tscore;
//						if( iterate == 0 ) fprintf( stderr, "hist %d-%d-%d, gain=%f (Thread %d)\n", iterate, l, k, gain, thread_no );
						}
					}
					else
					{
						pthread_mutex_lock( targ->mutex );
						for( i=0; i<locnjob; i++ ) strcpy( mastercopy[i], localcopy[i] );
						*subgenerationpt += 1;
						gainlist[thread_no] = *basegainpt + gain;
						*basegainpt += gain;

						if( scoreout )
						{
							unweightedspscore = plainscore( locnjob, localcopy );
							fprintf( stderr, "\nSCORE %d = %.0f, ", *ntrypt, unweightedspscore );
							fprintf( stderr, "SCORE / residue = %f", unweightedspscore / ( locnjob * strlen( mastercopy[0] ) ) );
							if( weight || constraint ) fprintf( stderr, " (differs from the objective score)" );
							fprintf( stderr, "\n" );
						}

						pthread_mutex_unlock( targ->mutex );
						tscorelist[thread_no] = tscore;
					}
#if 0
					fprintf( stderr, "tscore =  %f   mscore = %f  accepted.\n", tscore, mscore );
					fprintf( stderr, "\nbetter! gain = %f (thread %d)\r", gain, thread_no );
#else
					fprintf( stderr, "%03d-%04d-%d (thread %4d) better     \r", iterate+1, *ndonept, k, thread_no );
#endif

				}
				else 
				{
#if 0
					fprintf( stderr, "tscore =  %f   mscore = %f  rejected.\r", tscore, mscore );
					fprintf( stderr, "worse! gain = %f", gain );
#else
					fprintf( stderr, "%03d-%04d-%d (thread %4d) worse      \r", iterate+1, *ndonept, k, thread_no );
#endif
					tscore = mscore;
				}
			}
			converged2 = 0;
			for( ii=iterate-2; ii>=0; ii-=1 ) 
			{
//				fprintf( stderr, "Checking tscorehistory %f ?= %f\n", tscore, tscorehistory_detail[ii][branchpos] );
				if( tscore == tscorehistory_detail[ii][branchpos] )
				{
					converged2 = 1;
					break;
				}
			}
			if( parallelizationstrategy != BESTFIRST && converged2 ) 
			{
//				fprintf( stderr, "\nFINISH!\n" );
				pthread_mutex_lock( targ->mutex );
				*finishpt = 1;
				pthread_mutex_unlock( targ->mutex );
			}

			tscorehistory_detail[iterate][branchpos] = tscore;
			fprintf( stderr, "\r" );

			pthread_mutex_lock( targ->mutex );
			(*ndonept)++;
//			fprintf( stderr, "*ndonept = %d, nbranch = %d (thread %d) iterate=%d\n", *ndonept, nbranch, thread_no, iterate );
			generationofinput[branchpos] = iterate;
			if( *ndonept == nbranch )
			{
				if( *collectingpt != -1 ) *collectingpt = 1; // chofuku
//				fprintf( stderr, "Thread %d sends a signal, *ndonept = %d\n", thread_no, *ndonept );
				pthread_cond_signal( targ->collection_start );
			}
			pthread_mutex_unlock( targ->mutex );
		}            /* while( 1 ) */

	}                  /* for( iterate ) */
//	return( NULL );
}
#endif


int TreeDependentIteration( int locnjob, char **name, int nlen[M], 
                             char **aseq, char **bseq, int ***topol, double **len, 
							 double **distmtx,
							 int **skipthisbranch,
                             int alloclen, LocalHom **localhomtable, 
							 RNApair ***singlerna,
							 int nkozo, char *kozoarivec )
{
	int i, j, k, l, iterate, ii, iu, ju;
	int lin, ldf, length; 
	int clus1, clus2;
	int s1, s2;
	static double **imanoten;
	static Node *stopol;
	static double *effarrforlocalhom = NULL;
	static double *effarr = NULL;
	static double *effarr1 = NULL;
	static double *effarr2 = NULL;
	static double *effarr_kozo = NULL;
	static double *effarr1_kozo = NULL;
	static double *effarr2_kozo = NULL;
	static double **mtx = NULL;
	static int **node = NULL;
	static int *branchnode = NULL;
	static double **branchWeight = NULL;
	static char **mseq1, **mseq2;
	static float ***history;
	FILE *trap;
	double tscore, mscore;
	int identity;
	int converged;
	int oscillating;
	float naivescore0 = 0.0; // by D.Mathog, a guess
	float naivescore1;
#if 0
	char pair[njob][njob];
#else
	static int **memlist;
	static char *pairbuf;
#endif
#if DEBUG + RECORD
	double score_for_check0, score_for_check1;
	static double **effmtx = NULL;
	extern double score_calc0();
#endif
	static char *indication1, *indication2;
	static LocalHom ***localhomshrink = NULL;
	float impmatch = 0.0, oimpmatch = 0.0;
	static int *gapmap1;
	static int *gapmap2;
	double tmpdouble;
	int intdum;
	static RNApair *rnapairboth;
	RNApair ***grouprna1, ***grouprna2;
	double unweightedspscore;
	static double **smalldistmtx;
#if MCD
	static double **offsetmtx = NULL; // check senyou
#endif
	static double ***scoringmatrices;
	static double **eff1s, **eff2s; 
	static int **whichmtx;

	if( rnakozo && rnaprediction == 'm' )
	{
		grouprna1 = (RNApair ***)calloc( njob, sizeof( RNApair ** ) );
		grouprna2 = (RNApair ***)calloc( njob, sizeof( RNApair ** ) );
	}
	else
	{
		grouprna1 = grouprna2 = NULL;
	}

	Writeoptions( trap_g );
	fflush( trap_g );

	if( effarr == NULL ) /* locnjob == njob ni kagiru */
	{
		indication1 = AllocateCharVec( 150 );
		indication2 = AllocateCharVec( 150 );
		effarr = AllocateDoubleVec( locnjob );
		effarrforlocalhom = AllocateDoubleVec( locnjob );
		effarr1 = AllocateDoubleVec( locnjob );
		effarr2 = AllocateDoubleVec( locnjob );
		mseq1 = AllocateCharMtx( locnjob, 0 );
		mseq2 = AllocateCharMtx( locnjob, 0 );
		mtx = AllocateDoubleMtx( locnjob, locnjob );
		node = AllocateIntMtx( locnjob, locnjob );
		branchnode = AllocateIntVec( locnjob );
		branchWeight = AllocateDoubleMtx( locnjob, 2 );
		history = AllocateFloatCub( niter, locnjob, 2 );
		stopol = (Node *)calloc( locnjob * 2, sizeof( Node ) );
		gapmap1 = AllocateIntVec( alloclen );
		gapmap2 = AllocateIntVec( alloclen );
		if( score_check == 2 ) imanoten = AllocateDoubleMtx( njob, njob );
		if( specificityconsideration != 0 ) 
		{
			smalldistmtx = AllocateDoubleMtx( locnjob, locnjob ); // ookii?
#if MCD
			offsetmtx = AllocateDoubleMtx( locnjob, locnjob ); // ookii?
#endif
			scoringmatrices = AllocateDoubleCub( maxdistclass, nalphabets, nalphabets );
			makescoringmatrices( scoringmatrices, n_dis_consweight_multi );
			eff1s = AllocateDoubleMtx( maxdistclass, locnjob );
			eff2s = AllocateDoubleMtx( maxdistclass, locnjob );
			whichmtx = AllocateIntMtx( locnjob, locnjob );
		}
		else
		{
			smalldistmtx = NULL;
#if MCD
			offsetmtx = NULL;
#endif
			scoringmatrices = NULL;
			eff1s = eff2s = NULL;
			whichmtx = NULL;
		}

		effarr1_kozo = AllocateDoubleVec( locnjob ); // tsuneni allocate suru.
		effarr2_kozo = AllocateDoubleVec( locnjob ); // tsuneni allocate suru.
		effarr_kozo = AllocateDoubleVec( locnjob ); 
		for( i=0; i<locnjob; i++ )
			effarr_kozo[i] = effarr1_kozo[i] = effarr2_kozo[i] = 0.0;

#if 0
#else
		pairbuf = AllocateCharVec( locnjob );
		memlist = AllocateIntMtx( 2, locnjob );
		if( rnakozo ) rnapairboth = (RNApair *)calloc( alloclen, sizeof( RNApair ) );

		if( constraint )
		{
			localhomshrink = (LocalHom ***)calloc( njob, sizeof( LocalHom ** ) );
			for( i=0; i<njob; i++)
			{
				localhomshrink[i] = (LocalHom **)calloc( njob, sizeof( LocalHom * ) );
			}
		}
#endif
	}
#if DEBUG + RECORD
	if( !effmtx ) effmtx = AllocateDoubleMtx( locnjob, locnjob );
	for( i=0; i<locnjob; i++ ) for( j=0; j<locnjob; j++ ) effmtx[i][j] = 1.0;
#endif

	for( i=0; i<locnjob; i++ ) strcpy( bseq[i], aseq[i] );

	writePre( locnjob, name, nlen, aseq, 0 );

	if( utree )
	{
		if( constraint )
		{
			counteff_simple( locnjob, topol, len, effarrforlocalhom );
			calcimportance( locnjob, effarrforlocalhom, aseq, localhomtable );
		}

		if( weight == 2 ) 
		{
			countnode_int( locnjob, topol, node );
			if( nkozo )
			{
				fprintf( stderr, "Not supported, weight=%d nkozo=%d.\n", weight, nkozo );
			}
		}
		else if( weight == 4 )
		{
			treeCnv( stopol, locnjob, topol, len, branchWeight );
			calcBranchWeight( branchWeight, locnjob, stopol, topol, len );
		}
	}

#ifdef enablemultithread
	if( nthread > 0 )
	{
		threadarg_t *targ;
		pthread_t *handle;
		pthread_mutex_t mutex;
		pthread_cond_t collection_end;
		pthread_cond_t collection_start;
		int jobposint;
		int generationofmastercopy;
		int subgeneration;
		float basegain;
		int *generationofinput;
		float *gainlist;
		float *tscorelist;
		int ndone;
		int ntry;
		int collecting;
		int nbranch;
		int maxiter;
		char ***candidates;
		int *branchtable;
		float **tscorehistory_detail;
		int finish;

		nwa = nthread + 1;
		nbranch = (njob-1) * 2 - 1;
		maxiter = niter;

		targ = calloc( nwa, sizeof( threadarg_t ) );
		handle = calloc( nwa, sizeof( pthread_t ) );
		pthread_mutex_init( &mutex, NULL );
		pthread_cond_init( &collection_end, NULL );
		pthread_cond_init( &collection_start, NULL );

		gainlist = calloc( nwa, sizeof( float ) );
		tscorelist = calloc( nwa, sizeof( float ) );
		branchtable = calloc( nbranch, sizeof( int ) );
		generationofinput = calloc( nbranch, sizeof( int ) );
		if( parallelizationstrategy == BESTFIRST )
			candidates = AllocateCharCub( nwa, locnjob, alloclen );
		for( i=0; i<nbranch; i++ ) branchtable[i] = i;
		tscorehistory_detail = AllocateFloatMtx( maxiter, nbranch );

		collecting = 1;

		for( i=0; i<nwa; i++ )
		{
			targ[i].thread_no = i;
			targ[i].njob = njob;
			targ[i].nbranch = nbranch;
			targ[i].maxiter = maxiter;
			targ[i].ndonept = &ndone;
			targ[i].ntrypt = &ntry;
			targ[i].collectingpt = &collecting;
			targ[i].jobposintpt = &jobposint;
			targ[i].gainlist = gainlist;
			targ[i].tscorelist = tscorelist;
			targ[i].nkozo = nkozo;
			targ[i].kozoarivec = kozoarivec;
			targ[i].mastercopy = bseq;
			targ[i].candidates = candidates;
			targ[i].subgenerationpt = &subgeneration;
			targ[i].basegainpt = &basegain;
			targ[i].generationofmastercopypt = &generationofmastercopy;
			targ[i].generationofinput = generationofinput;
			targ[i].branchtable = branchtable;
			targ[i].singlerna = singlerna;
			targ[i].localhomtable = localhomtable;
			targ[i].alloclen = alloclen;
			targ[i].stopol = stopol;
			targ[i].topol = topol;
			targ[i].skipthisbranch = skipthisbranch;
			targ[i].distmtx = distmtx;
//			targ[i].len = len;
			targ[i].mutex = &mutex;
			targ[i].collection_end = &collection_end;
			targ[i].collection_start = &collection_start;
			targ[i].tscorehistory_detail = tscorehistory_detail;
			targ[i].finishpt = &finish;
	
			pthread_create( handle+i, NULL, athread, (void *)(targ+i) );
		}

		for( i=0; i<nwa; i++ )
		{
			pthread_join( handle[i], NULL );
		}

		pthread_mutex_destroy( &mutex );
		pthread_cond_destroy( &collection_end );
		pthread_cond_destroy( &collection_start );

		free( targ );
		free( handle );
		free( gainlist );
		free( tscorelist );
		free( branchtable );
		free( generationofinput );
		if( parallelizationstrategy == BESTFIRST )
			FreeCharCub( candidates );
		FreeFloatMtx( tscorehistory_detail );
	}
	else
#endif
	{
#if 0
		int *branchtable;
		int jobpos;
		int myjob;

		int nbranch;
		nbranch = (njob-1) * 2 - 1;

		branchtable = calloc( nbranch, sizeof( int ) );
		for( i=0; i<nbranch; i++ ) branchtable[i] = i;

		srand( randomseed );
#endif

		if( parallelizationstrategy == BESTFIRST )
		{
			fprintf( stderr, "Not implemented.  Try --thread 1 --bestfirst\n" );
			exit( 1 );
		}
		converged = 0;
		if( cooling ) cut *= 2.0;
		for( iterate = 0; iterate<niter; iterate++ ) 
		{
			if( cooling ) cut *= 0.5; /* ... */

#if 0
			if( randomseed != 0 ) shuffle( branchtable, nbranch );
#endif
	
			fprintf( trap_g, "\n" );
			Writeoption2( trap_g, iterate, cut );
			fprintf( trap_g, "\n" );

	
			if( utree == 0 )
			{
				if( nkozo )
				{
					fprintf( stderr, "The combination of dynamic tree and kozo is not supported.\n" );
					exit( 1 );
				}
				if( devide )
				{
					static char *buff1 = NULL;
					static char *buff2 = NULL;
					if( !buff1 )
					{
						buff1 = AllocateCharVec( alloclen );
						buff2 = AllocateCharVec( alloclen );
					}	
	
					for( i=0; i<locnjob-1; i++ ) for( j=i+1; j<locnjob; j++ ) 	
					{
						buff1[0] = buff2[0] = 0;
						strcat( buff1, res_g[i] ); strcat( buff2, res_g[j] );
						strcat( buff1,  bseq[i] ); strcat( buff2,  bseq[j] );
						strcat( buff1, seq_g[i] ); strcat( buff2, seq_g[j] );
	
						mtx[i][j] = (double)substitution_hosei( buff1, buff2 );
					}
				}
				else
				{
					for( i=0; i<locnjob-1; i++ ) for( j=i+1; j<locnjob; j++ ) 	
						mtx[i][j] = (double)substitution_hosei( bseq[i], bseq[j] );
				}
	
				if     ( treemethod == 'n' )
					nj( locnjob, mtx, topol, len );
				else if( treemethod == 's' )
					spg( locnjob, mtx, topol, len );
				else if( treemethod == 'X' )
					supg( locnjob, mtx, topol, len );
				else if( treemethod == 'p' )
					upg2( locnjob, mtx, topol, len );
				/* veryfastsupgは、今のところ使えません。*/
				/* 順番の問題があるので                  */
	
				if( weight == 2 )
					countnode_int( locnjob, topol, node );
				else if( weight == 4 )
				{
					treeCnv( stopol, locnjob, topol, len, branchWeight );
					calcBranchWeight( branchWeight, locnjob, stopol, topol, len );
				}
				trap = fopen( "hat2", "w" );
				if( !trap ) ErrorExit( "Cannot open hat2." );
				WriteHat2_pointer( trap, locnjob, name, mtx );
				fclose( trap );
				if( constraint )
				{
					counteff_simple( locnjob, topol, len, effarrforlocalhom );
					calcimportance( locnjob, effarrforlocalhom, aseq, localhomtable );
				}
			}
	
			if( iterate % 2 == 0 ) 
			{
				lin = 0; ldf = +1;
			}
			else
			{
				lin = locnjob - 2; ldf = -1;
			}	
	
			if( score_check == 2 )
			{
				effarr1[0] = 1.0;
				effarr2[0] = 1.0;
				length = strlen( bseq[0] );
				for( i=0; i<locnjob-1; i++ )
					for( j=i+1; j<locnjob; j++ )
						intergroup_score( bseq+i, bseq+j, effarr1, effarr2, 1, 1, length, imanoten[i]+j );
			}
	
#if 1
			for( l=lin; l < locnjob-1 && l >= 0 ; l+=ldf )
			{


				for( k=0; k<2; k++ ) 
				{
					if( l == locnjob-2 ) k = 1;
#else

			for( jobpos=0; jobpos<nbranch; jobpos++)
			{
				{
					myjob = branchtable[jobpos];
					l = myjob / 2;
					if( l == locnjob-2 ) k = 1;
					else k = myjob - l * 2;
#endif
	#if 0 // IRANAI!!!!
					fprintf( stderr, "\nSTEP %d-%d\n", l, k );
					for( i=0; topol[l][k][i]!=-1; i++ ) fprintf( stderr, " %d ", topol[l][k][i]+1 );
					fprintf( stderr, "\n" );
					fprintf( stderr, "SKIP %d\n", skipthisbranch[l][k] );
	#endif
	#if 1
					fprintf( stderr, "STEP %03d-%03d-%d ", iterate+1, l+1, k );
					fflush( stderr );
	#else
					fprintf( stderr, "STEP %03d-%03d-%d %s", iterate+1, l+1, k, use_fft?"\n":"\n" );
	#endif
					if( skipthisbranch[l][k] ) 
					{
						fprintf( stderr, " skip.      \r" );
						continue;
					}

//					for( i=0; i<2; i++ ) for( j=0; j<locnjob; j++ ) pair[i][j] = 0;
					OneClusterAndTheOther_fast( locnjob, memlist[0], memlist[1], &s1, &s2, pairbuf, topol, l, k, smalldistmtx, distmtx );
	#if 0 // IRANAI!!!!
					fprintf( stderr, "STEP%d-%d\n", l, k );
					for( i=0; i<2; i++ ) 
					{
						for( j=0; j<locnjob; j++ ) 
						{
							fprintf( stderr, "%#3d", pair[i][j] );
						}
						fprintf( stderr, "\n" );
					}
	#endif
					if( !weight )
					{
						for( i=0; i<locnjob; i++ ) effarr[i] = 1.0;
						if( nkozo )
						{
							for( i=0; i<locnjob; i++ ) 
							{
								if( kozoarivec[i] )
									effarr_kozo[i] = 1.0;
								else
									effarr_kozo[i] = 0.0;
							}
						}
					}
					else if( weight == 2 ) 
					{
						nodeFromABranch( locnjob, branchnode, node, topol, len, l, k );
						node_eff( locnjob, effarr, branchnode );
					}
					else if( weight == 4 )
					{
						weightFromABranch( locnjob, effarr, stopol, topol, l, k );
	#if 0
						if( nkozo )
						{
							assignstrweight( locnjob, effarr_kozo, stopol, topol, l, k, kozoarivec, effarr );
						}
	
	#else
						if( nkozo ) // hitomadu single weight.
							for( i=0; i<locnjob; i++ ) 
							{
								if( kozoarivec[i] ) effarr_kozo[i] = effarr[i]; 
								else effarr_kozo[i] = 0.0;
							}
	#endif
	#if 0
						fprintf( stderr, "\n" );
						fprintf( stderr, "effarr_kozo = \n" );
						for( i=0; i<locnjob; i++ ) fprintf( stderr, "%5.3f ", effarr_kozo[i] );
						fprintf( stderr, "\n" );
						fprintf( stderr, "effarr = \n" );
						for( i=0; i<locnjob; i++ ) fprintf( stderr, "%5.3f ", effarr[i] );
						fprintf( stderr, "\n\n" );
	#endif
					}
	
					for( i=0; i<locnjob; i++ ) strcpy( aseq[i], bseq[i] );
					length = strlen( aseq[0] );
	
					if( nkozo )
					{
//						double tmptmptmp;
//						tmptmptmp = 0.0;
//						clus1 = conjuctionfortbfast_kozo( &tmptmptmp, pair[0], s1, aseq, mseq1, effarr1, effarr, effarr1_kozo, effarr_kozo, indication1 );
						clus1 = fastconjuction_noname_kozo( memlist[0], aseq, mseq1, effarr1, effarr, effarr1_kozo, effarr_kozo, indication1 );
						for( i=0; i<clus1; i++ ) effarr1_kozo[i] *= 1.0; // 0.5 ga sairyo ?
//						tmptmptmp = 0.0;
//						clus2 = conjuctionfortbfast_kozo( &tmptmptmp, pair[1], s2, aseq, mseq2, effarr2, effarr, effarr2_kozo, effarr_kozo, indication2 );
						clus2 = fastconjuction_noname_kozo( memlist[1], aseq, mseq2, effarr2, effarr, effarr2_kozo, effarr_kozo, indication2 );
						for( i=0; i<clus2; i++ ) effarr2_kozo[i] *= 1.0; // 0.5 ga sairyo ?
	
	#if 0
						fprintf( stderr, "\ngroup1 = %s\n", indication1 );
						for( i=0; i<clus1; i++ ) fprintf( stderr, "effarr1_kozo[%d], effarr1[]  = %f, %f\n", i, effarr1_kozo[i], effarr1[i] );
						fprintf( stderr, "\ngroup2 = %s\n", indication2 );
						for( i=0; i<clus2; i++ ) fprintf( stderr, "effarr2_kozo[%d], effarr2[]  = %f, %f\n", i, effarr2_kozo[i], effarr2[i] );
	#endif
	
					}
					else
					{
//						clus1 = conjuctionfortbfast( pair[0], s1, aseq, mseq1, effarr1, effarr, indication1 );
//						clus2 = conjuctionfortbfast( pair[1], s2, aseq, mseq2, effarr2, effarr, indication2 );
						clus1 = fastconjuction_noname( memlist[0], aseq, mseq1, effarr1, effarr, indication1 );
						clus2 = fastconjuction_noname( memlist[1], aseq, mseq2, effarr2, effarr, indication2 );
					}
	
	
			        if( rnakozo && rnaprediction == 'm' )
			        {       
//						makegrouprnait( grouprna1, singlerna, pair[0], s1 );
//						makegrouprnait( grouprna2, singlerna, pair[1], s2 );
						makegrouprna( grouprna1, singlerna, memlist[0] );
						makegrouprna( grouprna2, singlerna, memlist[1] );
			        }       
	
					if( smalldistmtx )
					{
						classifypairs( clus1, eff1s, effarr1, clus2, eff2s, effarr2, smalldistmtx, whichmtx, maxdistclass );
#if MCD
						for( i=0; i<clus1; i++ )for( j=0; j<clus2; j++ )
						{
							offsetmtx[i][j] = dist2offset( smalldistmtx[i][j] );
//							fprintf( stderr, "offsetmtx[%d][%d] = %f -> %f\n", i, j, smalldistmtx[i][j], offsetmtx[i][j] );
//							smalldistmtx[i][j] = smalldistmtx[i][j] * 1.0 - specificityconsideration;
//							if( smalldistmtx[i][j] > 0 ) smalldistmtx[i][j] = 0.0;
//							fprintf( stderr, "offset[%d][%d] = %f\n", i, j, smalldistmtx[i][j] );
						}
#endif
					}

					if( score_check == 2 )
					{
						fprintf( stderr, "Not supported\n" );
						exit( 1 );
						if( constraint )
						{
//							msshrinklocalhom( pair[0], pair[1], s1, s2, localhomtable, localhomshrink );
							msshrinklocalhom_fast( memlist[0], memlist[1], localhomtable, localhomshrink );
							oimpmatch = 0.0;
							if( use_fft )
							{
								if( alg == 'Q' )
								{
									part_imp_match_init_strictQ( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
									if( rnakozo ) part_imp_rnaQ( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, gapmap1, gapmap2, NULL );
									for(  i=length-1; i>=0; i-- ) oimpmatch += part_imp_match_out_scQ( i, i );
								}
								else
								{
									part_imp_match_init_strict( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );
									if( rnakozo ) part_imp_rna( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, gapmap1, gapmap2, NULL );
									for(  i=length-1; i>=0; i-- ) oimpmatch += part_imp_match_out_sc( i, i );
								}
							}
							else
							{
								if( alg == 'Q' )
								{
									imp_match_init_strictQ( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
									if( rnakozo ) imp_rnaQ( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, gapmap1, gapmap2, NULL );
									for(  i=length-1; i>=0; i-- ) oimpmatch += imp_match_out_scQ( i, i );
								}
								else
								{
									imp_match_init_strict( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );
									fprintf( stderr, "not supported\n" );
									exit( 1 );
								}
							}
	//						fprintf( stderr, "### oimpmatch = %f\n", oimpmatch );
						}
						else
						{
							oimpmatch = 0.0;
						}
#if 0
						tmpdouble = 0.0;
						iu=0; 
						for( i=s1; i<locnjob; i++ ) 
						{
							if( !pair[0][i] ) continue;
							ju=0;
							for( j=s2; j<locnjob; j++ )
							{
								if( !pair[1][j] ) continue;
	//							fprintf( stderr, "i = %d, j = %d, effarr1=%f, effarr2=%f\n", i, j, effarr1[iu], effarr2[ju] );
								tmpdouble += effarr1[iu] * effarr2[ju] * imanoten[MIN(i,j)][MAX(i,j)];
								ju++;
							}
							iu++;
						}
#else // not yet checked
						fprintf( stderr, "##### NOT YET CHECKED!!!!\n" );
						exit( 1 );
						tmpdouble = 0.0;
						iu=0; 
						for( i=0; (s1=memlist[0][i])!=-1; i++ ) 
						{
							ju=0;
							for( j=0; (s2=memlist[1][j])!=-1; j++ ) 
							{
	//							fprintf( stderr, "i = %d, j = %d, effarr1=%f, effarr2=%f\n", i, j, effarr1[iu], effarr2[ju] );
								tmpdouble += effarr1[iu] * effarr2[ju] * imanoten[MIN(s1,s2)][MAX(s1,s2)];
								ju++;
							}
							iu++;
						}
#endif
						mscore = oimpmatch + tmpdouble;
					}
					else if( score_check )
					{
						if( constraint )
						{
							if( RNAscoremtx == 'r' )
								intergroup_score_gapnomi( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble ); // gappick mae denaito dame
							else
								intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble ); // gappick mae denaito dame
	
//							shrinklocalhom( pair, s1, s2, localhomtable, localhomshrink );
							msshrinklocalhom_fast( memlist[0], memlist[1], localhomtable, localhomshrink );
		//					weightimportance4( clus1, clus2,  effarr1, effarr2, localhomshrink ); // >>>
							oimpmatch = 0.0;
							if( use_fft )
							{
								if( alg == 'Q' )
								{
									part_imp_match_init_strictQ( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
									if( rnakozo ) part_imp_rnaQ( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, gapmap1, gapmap2, NULL );
									for(  i=length-1; i>=0; i-- )
									{
										oimpmatch += part_imp_match_out_scQ( i, i );
	//									fprintf( stderr, "#### i=%d, initial impmatch = %f seq1 = %c, seq2 = %c\n", i, oimpmatch, mseq1[0][i], mseq2[0][i] );
									}
								}
								else
								{
									part_imp_match_init_strict( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );
									if( rnakozo ) part_imp_rna( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, gapmap1, gapmap2, NULL );
									for(  i=length-1; i>=0; i-- )
									{
										oimpmatch += part_imp_match_out_sc( i, i );
		//								fprintf( stderr, "#### i=%d, initial impmatch = %f seq1 = %c, seq2 = %c\n", i, oimpmatch, mseq1[0][i], mseq2[0][i] );
									}
								}
		//						fprintf( stderr, "otmpmatch = %f\n", oimpmatch );
							}
							else
							{
								if( alg == 'Q' )
								{
									imp_match_init_strictQ( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
									if( rnakozo ) imp_rnaQ( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, gapmap1, gapmap2, NULL );
	
									for(  i=length-1; i>=0; i-- )
									{
										oimpmatch += imp_match_out_scQ( i, i );
	//									fprintf( stderr, "#### i=%d, initial impmatch = %f\n", i, oimpmatch );
									}
								}
								else
								{
									imp_match_init_strict( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );
	
									fprintf( stderr, "not supported\n" );
									exit( 1 );
	
									for(  i=length-1; i>=0; i-- )
									{
										oimpmatch += imp_match_out_sc( i, i );
		//								fprintf( stderr, "#### i=%d, initial impmatch = %f seq1 = %c, seq2 = %c\n", i, oimpmatch, mseq1[0][i], mseq2[0][i] );
									}
								}
		//						fprintf( stderr, "otmpmatch = %f\n", oimpmatch );
							}
		//					fprintf( stderr, "#### initial impmatch = %f\n", oimpmatch );
						}
						else
						{
							if( RNAscoremtx == 'r' )
								intergroup_score_gapnomi( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble ); // gappick mae denaito dame
							else
							{
								if( smalldistmtx )
#if 1
									intergroup_score_multimtx( whichmtx, scoringmatrices, mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble );
#else
									intergroup_score_dynmtx( offsetmtx, amino_dis, mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble ); // n_dis ha machigai
#endif
								else
									intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble ); // gappick mae denaito dame
							}
							oimpmatch = 0.0;
						}
	
	
	//					fprintf( stderr, "#### tmpdouble = %f\n", tmpdouble );
						mscore = (double)oimpmatch + tmpdouble;
					}
					else
					{
	//					fprintf( stderr, "score_check = %d\n", score_check );
	/* atode kousokuka */
						intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble );
						mscore = tmpdouble;
	/* atode kousokuka */
	
						if( constraint )
						{
							oimpmatch = 0.0;
//							shrinklocalhom( pair, s1, s2, localhomtable, localhomshrink );
							msshrinklocalhom_fast( memlist[0], memlist[1], localhomtable, localhomshrink );
							if( use_fft )
							{
								if( alg == 'Q' )
								{
									part_imp_match_init_strictQ( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
									if( rnakozo ) part_imp_rnaQ( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, gapmap1, gapmap2, NULL );
								}
								else
								{
									part_imp_match_init_strict( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );
									if( rnakozo ) part_imp_rna( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, gapmap1, gapmap2, NULL );
								}
							}
							else
							{
								if( alg == 'Q' )
								{
									imp_match_init_strictQ( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
									if( rnakozo ) imp_rnaQ( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, gapmap1, gapmap2, NULL );
								}
								else
								{
									imp_match_init_strict( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );
									fprintf( stderr, "Not supported\n" );
									exit( 1 );
								}
							}
						}
					}
	
	//				oimpmatch = 0.0;
					if( constraint )
					{
	#if 0 // iranai
						if( alg == 'Q' )
						{
							imp_match_init_strictQ( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
							for(  i=length-1; i>=0; i-- )
							{
								oimpmatch += imp_match_out_scQ( i, i );
	//							fprintf( stderr, "#### i=%d, initial impmatch = %f seq1 = %c, seq2 = %c\n", i, oimpmatch, mseq1[0][i], mseq2[0][i] );
							}
						}
						else
						{
							imp_match_init_strict( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );
							for(  i=length-1; i>=0; i-- )
							{
								oimpmatch += imp_match_out_sc( i, i );
	//							fprintf( stderr, "#### i=%d, initial impmatch = %f seq1 = %c, seq2 = %c\n", i, oimpmatch, mseq1[0][i], mseq2[0][i] );
							}
						}
	#endif
					}
	#if 0
					if( alg == 'H' )
						naivescore0 = naivepairscore( clus1, clus2, mseq1, mseq2, effarr1, effarr2, penalty ) + oimpmatch;
					else if( alg == 'Q' )
						naivescore0 = naiveQpairscore( clus1, clus2, mseq1, mseq2, effarr1, effarr2, penalty ) + oimpmatch;
					else if( alg == 'R' )
						naivescore0 = naiveRpairscore( clus1, clus2, mseq1, mseq2, effarr1, effarr2, penalty ) + oimpmatch;
	#endif
	
	//				if( rnakozo ) foldalignedrna( clus1, clus2, mseq1, mseq2, effarr1, effarr2, rnapairboth );
	
	//				if( !use_fft && !rnakozo )
					if( !use_fft )
					{
						commongappick_record( clus1, mseq1, gapmap1 );
						commongappick_record( clus2, mseq2, gapmap2 );
					}
	
	#if 0
					fprintf( stderr, "##### mscore = %f\n", mscore );
	#endif
		
	#if DEBUG
					if( !devide )
					{
			       		fprintf( trap_g, "\nSTEP%d-%d-%d\n", iterate+1, l+1, k );
	       		    	fprintf( trap_g, "group1 = %s\n", indication1 );
	   	   		    	fprintf( trap_g, "group2 = %s\n", indication2 );
						fflush( trap_g );
					}
		
	#endif
	#if 0
					printf( "STEP %d-%d-%d\n", iterate, l, k );
					for( i=0; i<clus2; i++ ) printf( "%f ", effarr2[i] );
					printf( "\n" );
	#endif
					if( constraint == 2 )
					{
						if( use_fft )
						{
	//						if( alg == 'Q' )
	//							part_imp_match_init_strictQ( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
	//						else
	//							part_imp_match_init_strict( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );
							Falign_localhom( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, localhomshrink, &impmatch, gapmap1, gapmap2, NULL, 0, NULL );
	//						fprintf( stderr, "##### impmatch = %f\n", impmatch );
						}
						else
						{
							fprintf( stderr, "\n\nUnexpected error.  Please contact kazutaka.katoh@aist.go.jp\n\n\n" );
							exit( 1 );
						}
					}
					else if( use_fft )
					{
						float totalwm;
#if MCD
						totalwm = Falign( offsetmtx, scoringmatrices, n_dis_consweight_multi, mseq1, mseq2, effarr1, effarr2, eff1s, eff2s, clus1, clus2, alloclen, &intdum, NULL, 0, NULL );
#else
						totalwm = Falign( NULL, scoringmatrices, n_dis_consweight_multi, mseq1, mseq2, effarr1, effarr2, eff1s, eff2s, clus1, clus2, alloclen, &intdum, NULL, 0, NULL );
#endif
	
	//					fprintf( stderr, "totalwm = %f\n", totalwm );
					}
					else
					{
						fprintf( stderr, "\n\nUnexpected error.  Please contact kazutaka.katoh@aist.go.jp\n\n\n" );
						exit( 1 );
					}
	//				fprintf( stderr, "## impmatch = %f\n", impmatch );
								
						if( checkC )
						{
							extern double DSPscore();
							extern double SSPscore();
							static double cur;
							static double pre;
		
	/*
							pre = SSPscore( locnjob, bseq );
							cur = SSPscore( locnjob, aseq );
	*/
							pre = DSPscore( locnjob, bseq );
							cur = DSPscore( locnjob, aseq );
		
							fprintf( stderr, "Previous Sscore = %f\n", pre );
							fprintf( stderr, "Currnet  Sscore = %f\n\n", cur );
						}
						
	//				fprintf( stderr, "## impmatch = %f\n", impmatch );
					identity = !strcmp( aseq[s1], bseq[s1] );
					identity *= !strcmp( aseq[s2], bseq[s2] );
	
	
	/* Bug?  : idnetitcal but score change when scoreing mtx != JTT  */
	
					length = strlen( mseq1[0] );
		
					if( identity )
					{
						tscore = mscore;
						if( !devide ) fprintf( trap_g, "tscore =  %f   identical.\n", tscore );
						fprintf( stderr, " identical.   " );
						converged++;
					}
					else
					{
						if( score_check )
						{
							if( constraint == 2 )
							{
	#if 1
								if( RNAscoremtx == 'r' )
									intergroup_score_gapnomi( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble );
								else
									intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble );
	#else
								intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble );
	#endif
	
								tscore = impmatch + tmpdouble;
	
	//							fprintf( stderr, "tmpdouble=%f, impmatch = %f -> %f, tscore = %f\n", tmpdouble, oimpmatch, impmatch, tscore );
							}
							else
							{
								if( smalldistmtx )
#if 1
									intergroup_score_multimtx( whichmtx, scoringmatrices, mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble );
#else
									intergroup_score_dynmtx( offsetmtx, amino_dis, mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble ); // n_dis ha machigai
#endif
								else
									intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble );
								tscore = tmpdouble;
							}
	//						fprintf( stderr, "#######ii=%d, iterate=%d score = %f -> %f \n", ii, iterate , mscore, tscore );
		#if 0
							for( i=0; i<1; i++ )
								fprintf( stderr, "%s\n", mseq1[i] );
							fprintf( stderr, "+++++++\n" );
							for( i=0; i<1; i++ )
								fprintf( stderr, "%s\n", mseq2[i] );
		#endif
		
						}
						else
						{
							tscore = mscore + 1.0;
	//						tscore = 0.0;
	//						fprintf( stderr, "in line 705, tscore=%f\n", tscore );
	//						for( i=0; i<length; i++ )
	//							tscore = tscore + (double)mseq1[0][i];
	//						mscore = tscore - 1.0;
						}
	
						if( isnan( mscore ) )
						{
							fprintf( stderr, "\n\nmscore became NaN\n" );
							exit( 1 );
						}
						if( isnan( tscore ) )
						{
							fprintf( stderr, "\n\ntscore became NaN\n" );
							exit( 1 );
						}
	
	
	
	//					fprintf( stderr, "@@@@@ mscore,tscore = %f,%f\n", mscore, tscore );
	
						if( tscore > mscore - cut/100.0*mscore ) 
						{
							writePre( locnjob, name, nlen, aseq, 0 );
							for( i=0; i<locnjob; i++ ) strcpy( bseq[i], aseq[i] );
							if( score_check == 2 )
							{
								effarr1[0] = 1.0;
								effarr2[0] = 1.0;
								for( i=0; i<locnjob-1; i++ )
									for( j=i+1; j<locnjob; j++ )
										intergroup_score( bseq+i, bseq+j, effarr1, effarr2, 1, 1, length, imanoten[i]+j );
							}
		
	#if 0
							fprintf( stderr, "tscore =  %f   mscore = %f  accepted.\n", tscore, mscore );
	#endif
							fprintf( stderr, " accepted." );
							converged = 0;
		
						}
						else 
						{
	#if 0
							fprintf( stderr, "tscore =  %f   mscore = %f  rejected.\n", tscore, mscore );
	#endif
							fprintf( stderr, " rejected." );
							tscore = mscore;
							converged++;
						}
					}
					fprintf( stderr, "\r" );
	
	
					history[iterate][l][k] = (float)tscore;
	
	//				fprintf( stderr, "tscore = %f\n", tscore );
		
					if( converged >= locnjob * 2 )
					{
						fprintf( trap_g, "Converged.\n\n" );
						fprintf( stderr, "\nConverged.\n\n" );
						if( scoreout )
						{
							unweightedspscore = plainscore( njob, bseq );
							fprintf( stderr, "\nSCORE %d = %.0f, ", iterate * ( (njob-1)*2-1 ), unweightedspscore );
							fprintf( stderr, "SCORE / residue = %f", unweightedspscore / ( njob * strlen( bseq[0] ) ) );
							if( weight || constraint ) fprintf( stderr, " (differs from the objective score)" );
							fprintf( stderr, "\n\n" );
						}
						if( grouprna1 ) free( grouprna1 );
						if( grouprna2 ) free( grouprna2 );
						return( 0 );
					}
					if( iterate >= 1 )
					{
		/*   oscillation check    */
						oscillating = 0;
						for( ii=iterate-2; ii>=0; ii-=2 ) 
						{
							if( (float)tscore == history[ii][l][k] )
							{
								oscillating = 1;
								break;
							}
						}
						if( ( oscillating && !cooling ) || ( oscillating && cut < 0.001 && cooling ) )
						{
							fprintf( trap_g, "Oscillating.\n" );
							fprintf( stderr, "\nOscillating.\n\n" );
							if( scoreout )
							{
								unweightedspscore = plainscore( njob, bseq );
								fprintf( stderr, "\nSCORE %d = %.0f, ", iterate * ( (njob-1)*2-1 ), unweightedspscore );
								fprintf( stderr, "SCORE / residue = %f", unweightedspscore / ( njob * strlen( bseq[0] ) ) );
								if( weight || constraint ) fprintf( stderr, " (differs from the objective score)" );
								fprintf( stderr, "\n\n" );
							}
	#if 1 /* hujuubun */
							if( grouprna1 ) free( grouprna1 );
							if( grouprna2 ) free( grouprna2 );
							return( -1 );
	#endif
						}
					}      /* if( iterate ) */
				}          /* for( k ) */
			}              /* for( l ) */
			if( scoreout )
			{
				unweightedspscore = plainscore( njob, bseq );
				fprintf( stderr, "\nSCORE %d = %.0f, ", iterate * ( (njob-1)*2-1 ), unweightedspscore );
				fprintf( stderr, "SCORE / residue = %f", unweightedspscore / ( njob * strlen( bseq[0] ) ) );
				if( weight || constraint ) fprintf( stderr, " (differs from the objective score)" );
				fprintf( stderr, "\n\n" );
			}
		}                  /* for( iterate ) */
	}
	if( grouprna1 ) free( grouprna1 );
	if( grouprna2 ) free( grouprna2 );
	return( 2 );
}                  	   /* int Tree... */
