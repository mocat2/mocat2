#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0

static int nadd;
static int treein;
static int topin;
static int treeout;
static int distout;
static int noalign;
static int multidist;
static int subalignment;
static int subalignmentoffset;

#ifdef enablemultithread
typedef struct _jobtable
{
    int i;  
    int j;  
} Jobtable;

typedef struct _distancematrixthread_arg
{
	int njob;
	int thread_no;
	float *selfscore;
	float **iscore;
	char **seq;
	Jobtable *jobpospt;
	pthread_mutex_t *mutex;
} distancematrixthread_arg_t;

typedef struct _treebasethread_arg
{
	int thread_no;
	int *nrunpt;
	int njob;
	int *nlen;
	int *jobpospt;
	int ***topol;
	Treedep *dep;
	char **aseq;
	double *effarr;
	int *alloclenpt;
	LocalHom **localhomtable;
	RNApair ***singlerna;
	double *effarr_kozo;
	int *fftlog;
	char *mergeoralign;
	pthread_mutex_t *mutex;
	pthread_cond_t *treecond;
} treebasethread_arg_t;
#endif

void arguments( int argc, char *argv[] )
{
    int c;

	nthread = 1;
	outnumber = 0;
	scoreout = 0;
	treein = 0;
	topin = 0;
	rnaprediction = 'm';
	rnakozo = 0;
	nevermemsave = 0;
	inputfile = NULL;
	addfile = NULL;
	addprofile = 1;
	fftkeika = 0;
	constraint = 0;
	nblosum = 62;
	fmodel = 0;
	calledByXced = 0;
	devide = 0;
	use_fft = 0; // chuui
	force_fft = 0;
	fftscore = 1;
	fftRepeatStop = 0;
	fftNoAnchStop = 0;
    weight = 3;
    utree = 1;
	tbutree = 1;
    refine = 0;
    check = 1;
    cut = 0.0;
    disp = 0;
    outgap = 1;
    alg = 'A';
    mix = 0;
	tbitr = 0;
	scmtd = 5;
	tbweight = 0;
	tbrweight = 3;
	checkC = 0;
	treemethod = 'X';
	contin = 0;
	scoremtx = 1;
	kobetsubunkatsu = 0;
	dorp = NOTSPECIFIED;
	ppenalty = NOTSPECIFIED;
	ppenalty_ex = NOTSPECIFIED;
	poffset = NOTSPECIFIED;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	geta2 = GETA2;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;
	RNAppenalty = NOTSPECIFIED;
	RNAppenalty_ex = NOTSPECIFIED;
	RNApthr = NOTSPECIFIED;
	TMorJTT = JTT;
	consweight_multi = 1.0;
	consweight_rna = 0.0;
	multidist = 0;
	subalignment = 0;
	subalignmentoffset = 0;
	legacygapcost = 0;
	specificityconsideration = 0.0;

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( ( c = *++argv[0] ) )
		{
            switch( c )
            {
				case 'i':
					inputfile = *++argv;
					fprintf( stderr, "inputfile = %s\n", inputfile );
					--argc;
                    goto nextoption;
				case 'I':
					nadd = myatoi( *++argv );
					fprintf( stderr, "nadd = %d\n", nadd );
					--argc;
					goto nextoption;
				case 'e':
					RNApthr = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'o':
					RNAppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'f':
					ppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( atof( *++argv ) * 1000 - 0.5 );
					fprintf( stderr, "ppenalty_ex = %d\n", ppenalty_ex );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "poffset = %d\n", poffset );
					--argc;
					goto nextoption;
				case 'k':
					kimuraR = myatoi( *++argv );
					fprintf( stderr, "kappa = %d\n", kimuraR );
					--argc;
					goto nextoption;
				case 'b':
					nblosum = myatoi( *++argv );
					scoremtx = 1;
					fprintf( stderr, "blosum %d / kimura 200\n", nblosum );
					--argc;
					goto nextoption;
				case 'j':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = JTT;
					fprintf( stderr, "jtt/kimura %d\n", pamN );
					--argc;
					goto nextoption;
				case 'm':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = TM;
					fprintf( stderr, "tm %d\n", pamN );
					--argc;
					goto nextoption;
				case 'l':
					fastathreshold = atof( *++argv );
					constraint = 2;
					--argc;
					goto nextoption;
				case 'r':
					consweight_rna = atof( *++argv );
					rnakozo = 1;
					--argc;
					goto nextoption;
				case 'c':
					consweight_multi = atof( *++argv );
					--argc;
					goto nextoption;
				case 'C':
					nthread = myatoi( *++argv );
					fprintf( stderr, "nthread = %d\n", nthread );
					--argc; 
					goto nextoption;
				case 'H':
					subalignment = 1;
					subalignmentoffset = myatoi( *++argv );
					--argc;
					goto nextoption;
				case 's':
					specificityconsideration = (double)myatof( *++argv );
//					fprintf( stderr, "specificityconsideration = %f\n", specificityconsideration );
					--argc; 
					goto nextoption;
				case 'R':
					rnaprediction = 'r';
#if 1
				case 'a':
					fmodel = 1;
					break;
#endif
				case 'K':
					addprofile = 0;
					break;
				case 'y':
					distout = 1;
					break;
				case 't':
					treeout = 1;
					break;
				case 'T':
					noalign = 1;
					break;
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
					break;
				case 'L':
					legacygapcost = 1;
					break;
#if 1
				case 'O':
					outgap = 0;
					break;
#else
				case 'O':
					fftNoAnchStop = 1;
					break;
#endif
				case 'S':
					scoreout = 1;
					break;
#if 0
				case 'e':
					fftscore = 0;
					break;
				case 'r':
					fmodel = -1;
					break;
				case 'R':
					fftRepeatStop = 1;
					break;
				case 's':
					treemethod = 's';
					break;
#endif
				case 'X':
					treemethod = 'X';
					break;
				case 'E':
					treemethod = 'E';
					break;
				case 'q':
					treemethod = 'q';
					break;
				case 'n' :
					outnumber = 1;
					break;
#if 0
				case 'a':
					alg = 'a';
					break;
				case 'H':
					alg = 'H';
					break;
#endif
				case 'Q':
					alg = 'Q';
					break;
				case 'A':
					alg = 'A';
					break;
				case 'M':
					alg = 'M';
					break;
				case 'N':
					nevermemsave = 1;
					break;
				case 'B':
					break;
				case 'F':
					use_fft = 1;
					break;
				case 'G':
					force_fft = 1;
					use_fft = 1;
					break;
				case 'U':
					treein = 1;
					break;
				case 'V':
					topin = 1;
					break;
				case 'u':
					tbrweight = 0;
					weight = 0;
					break;
				case 'v':
					tbrweight = 3;
					break;
				case 'd':
					multidist = 1;
					break;
#if 0
				case 'd':
					disp = 1;
					break;
#endif
/* Modified 01/08/27, default: user tree */
				case 'J':
					tbutree = 0;
					break;
/* modification end. */
				case 'z':
					fftThreshold = myatoi( *++argv );
					--argc; 
					goto nextoption;
				case 'w':
					fftWinSize = myatoi( *++argv );
					--argc;
					goto nextoption;
				case 'Z':
					checkC = 1;
					break;
                default:
                    fprintf( stderr, "illegal option %c\n", c );
                    argc = 0;
                    break;
            }
		}
		nextoption:
			;
	}
    if( argc == 1 )
    {
        cut = atof( (*argv) );
        argc--;
    }
    if( argc != 0 ) 
    {
        fprintf( stderr, "options: Check source file !\n" );
        exit( 1 );
    }
	if( tbitr == 1 && outgap == 0 )
	{
		fprintf( stderr, "conflicting options : o, m or u\n" );
		exit( 1 );
	}
	if( alg == 'C' && outgap == 0 )
	{
		fprintf( stderr, "conflicting options : C, o\n" );
		exit( 1 );
	}
}

#if 0
static void *distancematrixthread2( void *arg )
{
	distancematrixthread_arg_t *targ = (distancematrixthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	float *selfscore = targ->selfscore;
	float **iscore = targ->iscore;
	char **seq = targ->seq;
	Jobtable *jobpospt = targ->jobpospt;

	float ssi, ssj, bunbo;
	int i, j;

	while( 1 )
	{
		pthread_mutex_lock( targ->mutex );
		i = jobpospt->i;
		i++;
		if( i == njob-1 )
		{
			pthread_mutex_unlock( targ->mutex );
			return( NULL );
		}
		jobpospt->i = i;
		pthread_mutex_unlock( targ->mutex );

		ssi = selfscore[i];
		if( i % 10 == 0 ) fprintf( stderr, "\r% 5d / %d (thread %4d)", i, njob, thread_no );
		for( j=i+1; j<njob; j++ )
		{
			ssj = selfscore[j];
			bunbo = MIN( ssi, ssj );
			if( bunbo == 0.0 )
				iscore[i][j-i] = 1.0;
			else
				iscore[i][j-i] = 1.0 - naivepairscore11( seq[i], seq[j], penalty ) / bunbo;
		}
	}
}
#endif

#ifdef enablemultithread
static void *distancematrixthread( void *arg )
{
	distancematrixthread_arg_t *targ = (distancematrixthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	float *selfscore = targ->selfscore;
	float **iscore = targ->iscore;
	char **seq = targ->seq;
	Jobtable *jobpospt = targ->jobpospt;

	float ssi, ssj, bunbo;
	int i, j;

	while( 1 )
	{
		pthread_mutex_lock( targ->mutex );
		j = jobpospt->j;
		i = jobpospt->i;
		j++;
		if( j == njob )
		{
			i++;
			j = i + 1;
			if( i == njob-1 )
			{
				pthread_mutex_unlock( targ->mutex );
				return( NULL );
			}
		}
		jobpospt->j = j;
		jobpospt->i = i;
		pthread_mutex_unlock( targ->mutex );


		if( j==i+1 && i % 10 == 0 ) fprintf( stderr, "\r% 5d / %d (thread %4d)", i, njob, thread_no );
		ssi = selfscore[i];
		ssj = selfscore[j];
		bunbo = MIN( ssi, ssj );
		if( bunbo == 0.0 )
			iscore[i][j-i] = 2.0; // 2013/Oct/17
		else
			iscore[i][j-i] = ( 1.0 - naivepairscore11( seq[i], seq[j], penalty ) / bunbo ) * 2.0; // 2013/Oct/17
	}
}

static void *treebasethread( void *arg )
{
	treebasethread_arg_t *targ = (treebasethread_arg_t *)arg;
	int *nrunpt = targ->nrunpt;
	int thread_no = targ->thread_no;
	int njob = targ->njob;
	int *nlen = targ->nlen;
	int *jobpospt = targ->jobpospt;
	int ***topol = targ->topol;
	Treedep *dep = targ->dep;
	char **aseq = targ->aseq;
	double *effarr = targ->effarr;
	int *alloclen = targ->alloclenpt;
	LocalHom **localhomtable = targ->localhomtable;
	RNApair ***singlerna = targ->singlerna;
	double *effarr_kozo = targ->effarr_kozo;
	int *fftlog = targ->fftlog;
	char *mergeoralign = targ->mergeoralign;

	char **mseq1, **mseq2;
	char **localcopy;
	int i, j, l;
	int len1, len2;
	int clus1, clus2;
	float pscore;
	char *indication1, *indication2;
	double *effarr1 = NULL;
	double *effarr2 = NULL;
	double *effarr1_kozo = NULL;
	double *effarr2_kozo = NULL;
	LocalHom ***localhomshrink = NULL;
	int m1, m2;
	float dumfl = 0.0;
	int ffttry;
	RNApair ***grouprna1 = NULL, ***grouprna2 = NULL;
	double **dynamicmtx;


	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );
	localcopy = calloc( njob, sizeof( char * ) );
	dynamicmtx = AllocateDoubleMtx( nalphabets, nalphabets );

	if( rnakozo && rnaprediction == 'm' )
	{
		grouprna1 = (RNApair ***)calloc( njob, sizeof( RNApair ** ) );
		grouprna2 = (RNApair ***)calloc( njob, sizeof( RNApair ** ) );
	}
	else
	{
		grouprna1 = grouprna2 = NULL;
	}

	effarr1 = AllocateDoubleVec( njob );
	effarr2 = AllocateDoubleVec( njob );
	indication1 = AllocateCharVec( 150 );
	indication2 = AllocateCharVec( 150 );
#if 0
#else
	if( constraint )
	{
		localhomshrink = (LocalHom ***)calloc( njob, sizeof( LocalHom ** ) );
		for( i=0; i<njob; i++ )
			localhomshrink[i] = (LocalHom **)calloc( njob, sizeof( LocalHom *) );
	}
#endif
	effarr1_kozo = AllocateDoubleVec( njob ); //tsuneni allocate sareru.
	effarr2_kozo = AllocateDoubleVec( njob ); //tsuneni allocate sareru.
	for( i=0; i<njob; i++ ) effarr1_kozo[i] = 0.0;
	for( i=0; i<njob; i++ ) effarr2_kozo[i] = 0.0;


#if 0
#endif

#if 0
	for( i=0; i<njob; i++ )
		fprintf( stderr, "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif

#if 0 //-> main thread
	if( constraint )
		calcimportance( njob, effarr, aseq, localhomtable );
#endif


//	writePre( njob, name, nlen, aseq, 0 );


//	for( l=0; l<njob-1; l++ )
	while( 1 )
	{

		pthread_mutex_lock( targ->mutex );
		l = *jobpospt;
		if( l == njob-1 )
		{
			pthread_mutex_unlock( targ->mutex );
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
			Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
			A__align( NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
			free( mseq1 );
			free( mseq2 );
			free( localcopy );
			free( effarr1 );
			free( effarr2 );
			free( effarr1_kozo );
			free( effarr2_kozo );
			free( indication1 );
			free( indication2 );
			FreeDoubleMtx( dynamicmtx );
			if( rnakozo && rnaprediction == 'm' )
			{
				if( grouprna1 ) free( grouprna1 ); // nakami ha?
				if( grouprna2 ) free( grouprna2 ); // nakami ha?
				grouprna1 = grouprna2 = NULL;
			}
			if( constraint )
			{
				if( localhomshrink ) // nen no tame
				{
					for( i=0; i<njob; i++ )
					{
						free( localhomshrink[i] );
						localhomshrink[i] = NULL;
					}
					free( localhomshrink );
					localhomshrink = NULL;
				}
			}
			return( NULL );
		}
		*jobpospt = l+1;

		if( dep[l].child0 != -1 )
		{
			while( dep[dep[l].child0].done == 0 )
				pthread_cond_wait( targ->treecond, targ->mutex );
		}
		if( dep[l].child1 != -1 )
		{
			while( dep[dep[l].child1].done == 0 )
				pthread_cond_wait( targ->treecond, targ->mutex );
		}
//		while( *nrunpt >= nthread )
//			pthread_cond_wait( targ->treecond, targ->mutex );
		(*nrunpt)++;

//		pthread_mutex_unlock( targ->mutex );

		if( mergeoralign[l] == 'n' )
		{
//			fprintf( stderr, "SKIP!\n" );
			dep[l].done = 1;
			(*nrunpt)--;
			pthread_cond_broadcast( targ->treecond );
			free( topol[l][0] );
			free( topol[l][1] );
			free( topol[l] );
			pthread_mutex_unlock( targ->mutex );
			continue;
		}



		m1 = topol[l][0][0];
		m2 = topol[l][1][0];

//		fprintf( stderr, "\ndistfromtip = %f\n", dep[l].distfromtip );
//		makedynamicmtx( dynamicmtx, n_dis_consweight_multi, dep[l].distfromtip - 0.5 );
		makedynamicmtx( dynamicmtx, n_dis_consweight_multi, dep[l].distfromtip );

//		pthread_mutex_lock( targ->mutex );



        len1 = strlen( aseq[m1] );
        len2 = strlen( aseq[m2] );
        if( *alloclen <= len1 + len2 )
        {
			fprintf( stderr, "\nReallocating (by thread %d) ..", thread_no );
			*alloclen = ( len1 + len2 ) + 1000;
			ReallocateCharMtx( aseq, njob, *alloclen + 10  );
			fprintf( stderr, "done. *alloclen = %d\n", *alloclen );
		}

		for( i=0; (j=topol[l][0][i])!=-1; i++ )
		{
			localcopy[j] = calloc( *alloclen, sizeof( char ) );
			strcpy( localcopy[j], aseq[j] );
		}
		for( i=0; (j=topol[l][1][i])!=-1; i++ )
		{
			localcopy[j] = calloc( *alloclen, sizeof( char ) );
			strcpy( localcopy[j], aseq[j] );
		}

		pthread_mutex_unlock( targ->mutex );



		if( effarr_kozo )
		{
			clus1 = fastconjuction_noname_kozo( topol[l][0], localcopy, mseq1, effarr1, effarr, effarr1_kozo, effarr_kozo, indication1 );
			clus2 = fastconjuction_noname_kozo( topol[l][1], localcopy, mseq2, effarr2, effarr, effarr2_kozo, effarr_kozo, indication2 );
		}
		else
		{
			clus1 = fastconjuction_noname( topol[l][0], localcopy, mseq1, effarr1, effarr, indication1 );
			clus2 = fastconjuction_noname( topol[l][1], localcopy, mseq2, effarr2, effarr, indication2 );
		}



#if 1
		fprintf( stderr, "\rSTEP % 5d /%d (thread %4d) ", l+1, njob-1, thread_no );
#else
		fprintf( stderr, "STEP %d /%d (thread %d) \n", l+1, njob-1, thread_no );
		fprintf( stderr, "group1 = %.66s", indication1 );
		if( strlen( indication1 ) > 66 ) fprintf( stderr, "..." );
		fprintf( stderr, ", child1 = %d\n", dep[l].child0 );
		fprintf( stderr, "group2 = %.66s", indication2 );
		if( strlen( indication2 ) > 66 ) fprintf( stderr, "..." );
		fprintf( stderr, ", child2 = %d\n", dep[l].child1 );

		fprintf( stderr, "Group1's lengths = " );
		for( i=0; i<clus1; i++ ) fprintf( stderr, "%d ", strlen( mseq1[i] ) );
		fprintf( stderr, "\n" );
		fprintf( stderr, "Group2's lengths = " );
		for( i=0; i<clus2; i++ ) fprintf( stderr, "%d ", strlen( mseq2[i] ) );
		fprintf( stderr, "\n" );
		
#endif



//		for( i=0; i<clus1; i++ ) fprintf( stderr, "## STEP%d-eff for mseq1-%d %f\n", l+1, i, effarr1[i] );

		if( constraint )
		{
			fastshrinklocalhom( topol[l][0], topol[l][1], localhomtable, localhomshrink );
//			msfastshrinklocalhom( topol[l][0], topol[l][1], localhomtable, localhomshrink );
//			fprintf( stderr, "localhomshrink =\n" );
//			outlocalhompt( localhomshrink, clus1, clus2 );
//			weightimportance4( clus1, clus2, effarr1, effarr2, localhomshrink );
//			fprintf( stderr, "after weight =\n" );
//			outlocalhompt( localhomshrink, clus1, clus2 );
		}
		if( rnakozo && rnaprediction == 'm' )
		{
			makegrouprna( grouprna1, singlerna, topol[l][0] );
			makegrouprna( grouprna2, singlerna, topol[l][1] );
		}


/*
		fprintf( stderr, "before align all\n" );
		display( localcopy, njob );
		fprintf( stderr, "\n" );
		fprintf( stderr, "before align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		fprintf( stderr, "\n" );
		fprintf( stderr, "before align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		fprintf( stderr, "\n" );
*/

		if( !nevermemsave && ( constraint != 2  && alg != 'M'  && ( len1 > 30000 || len2 > 30000 ) ) )
		{
			fprintf( stderr, "\nlen1=%d, len2=%d, Switching to the memsave mode.\n", len1, len2 );
			alg = 'M';
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			commonAlloc1 = 0;
			commonAlloc2 = 0;
		}
		

//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000 );
		else						   ffttry = 0;
//		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 5000 && clus2 < 5000 ); // v6.708
//		fprintf( stderr, "f=%d, len1/fftlog[m1]=%f, clus1=%d, len2/fftlog[m2]=%f, clus2=%d\n", ffttry, (float)len1/fftlog[m1], clus1, (float)len2/fftlog[m2], clus2 );
//		fprintf( stderr, "f=%d, clus1=%d, fftlog[m1]=%d, clus2=%d, fftlog[m2]=%d\n", ffttry, clus1, fftlog[m1], clus2, fftlog[m2] );
		if( constraint == 2 )
		{
			if( alg == 'M' )
			{
				fprintf( stderr, "\n\nMemory saving mode is not supported.\n\n" );
				exit( 1 );
			}
			fprintf( stderr, "c" );
			if( alg == 'A' )
			{
				imp_match_init_strict( NULL, clus1, clus2, strlen( mseq1[0] ), strlen( mseq2[0] ), mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );
				if( rnakozo ) imp_rna( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, NULL, NULL, NULL );
				pscore = A__align( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, localhomshrink, &dumfl, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
			}
			else if( alg == 'H' )
			{
				imp_match_init_strictH( NULL, clus1, clus2, strlen( mseq1[0] ), strlen( mseq2[0] ), mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
				pscore = H__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, localhomshrink, &dumfl, NULL, NULL, NULL, NULL );
			}
			else if( alg == 'Q' )
			{
				imp_match_init_strictQ( NULL, clus1, clus2, strlen( mseq1[0] ), strlen( mseq2[0] ), mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
				if( rnakozo ) imp_rnaQ( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, NULL, NULL, NULL );
				pscore = Q__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, localhomshrink, &dumfl, NULL, NULL, NULL, NULL );
			}
			else if( alg == 'R' )
			{
				imp_match_init_strictR( NULL, clus1, clus2, strlen( mseq1[0] ), strlen( mseq2[0] ), mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
				pscore = R__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, localhomshrink, &dumfl, NULL, NULL, NULL, NULL );
			}
		}
		else if( force_fft || ( use_fft && ffttry ) )
		{
			fprintf( stderr, "f" );
			if( alg == 'M' )
			{
				fprintf( stderr, "m" );
				pscore = Falign_udpari_long( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			}
			else
				pscore = Falign( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1, NULL, 0, NULL );
		}
		else
		{
			fprintf( stderr, "d" );
			fftlog[m1] = 0;
			switch( alg )
			{
				case( 'a' ):
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen );
					break;
				case( 'M' ):
					fprintf( stderr, "m" );
					pscore = MSalignmm( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					break;
				case( 'A' ):
					pscore = A__align( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					break;
				case( 'Q' ):
					pscore = Q__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL );
					break;
				case( 'R' ):
					pscore = R__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL );
					break;
				case( 'H' ):
					pscore = H__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL );
					break;
				default:
					ErrorExit( "ERROR IN SOURCE FILE" );
			}
		}

		nlen[m1] = 0.5 * ( nlen[m1] + nlen[m2] );

#if SCOREOUT
		fprintf( stderr, "score = %10.2f\n", pscore );
#endif

/*
		fprintf( stderr, "after align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		fprintf( stderr, "\n" );
		fprintf( stderr, "after align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		fprintf( stderr, "\n" );
*/

//		writePre( njob, name, nlen, localcopy, 0 );

		if( disp ) display( localcopy, njob );

		pthread_mutex_lock( targ->mutex );
		dep[l].done = 1;
		(*nrunpt)--;
		pthread_cond_broadcast( targ->treecond );

//		pthread_mutex_unlock( targ->mutex );
//		pthread_mutex_lock( targ->mutex );

		for( i=0; (j=topol[l][0][i])!=-1; i++ )
			strcpy( aseq[j], localcopy[j] );
		for( i=0; (j=topol[l][1][i])!=-1; i++ )
			strcpy( aseq[j], localcopy[j] );
		pthread_mutex_unlock( targ->mutex );

		for( i=0; (j=topol[l][0][i])!=-1; i++ )
			free( localcopy[j] );
		for( i=0; (j=topol[l][1][i])!=-1; i++ )
			free( localcopy[j] );
		free( topol[l][0] );
		free( topol[l][1] );
		free( topol[l] );

	}
}
#endif

void treebase( int *nlen, char **aseq, int nadd, char *mergeoralign, char **mseq1, char **mseq2, int ***topol, Treedep *dep, double *effarr, int *alloclen, LocalHom **localhomtable, RNApair ***singlerna, double *effarr_kozo )
{
	int i, l, m;
	int len1nocommongap, len2nocommongap;
	int len1, len2;
	int clus1, clus2;
	float pscore, tscore;
	static char *indication1, *indication2;
	static double *effarr1 = NULL;
	static double *effarr2 = NULL;
	static double *effarr1_kozo = NULL;
	static double *effarr2_kozo = NULL;
	static LocalHom ***localhomshrink = NULL;
	static int *fftlog;
	int m1, m2;
	static int *gaplen;
	static int *gapmap;
	static int *alreadyaligned;
	float dumfl = 0.0;
	int ffttry;
	RNApair ***grouprna1 = NULL, ***grouprna2 = NULL;
	double **dynamicmtx;

	if( rnakozo && rnaprediction == 'm' )
	{
		grouprna1 = (RNApair ***)calloc( njob, sizeof( RNApair ** ) );
		grouprna2 = (RNApair ***)calloc( njob, sizeof( RNApair ** ) );
	}
	else
	{
		grouprna1 = grouprna2 = NULL;
	}

	if( effarr1 == NULL ) 
	{
		fftlog = AllocateIntVec( njob );
		effarr1 = AllocateDoubleVec( njob );
		effarr2 = AllocateDoubleVec( njob );
		indication1 = AllocateCharVec( 150 );
		indication2 = AllocateCharVec( 150 );
		gaplen = AllocateIntVec( *alloclen+10 );
		gapmap = AllocateIntVec( *alloclen+10 );
		alreadyaligned = AllocateIntVec( njob );
		dynamicmtx = AllocateDoubleMtx( nalphabets, nalphabets );
#if 0
#else
		if( constraint )
		{
			localhomshrink = (LocalHom ***)calloc( njob, sizeof( LocalHom ** ) );
			for( i=0; i<njob; i++ )
				localhomshrink[i] = (LocalHom **)calloc( njob, sizeof( LocalHom *) );
		}
#endif
		effarr1_kozo = AllocateDoubleVec( njob ); //tsuneni allocate sareru.
		effarr2_kozo = AllocateDoubleVec( njob ); //tsuneni allocate sareru.
		for( i=0; i<njob; i++ ) effarr1_kozo[i] = 0.0;
		for( i=0; i<njob; i++ ) effarr2_kozo[i] = 0.0;
	}

	for( i=0; i<njob-nadd; i++ ) alreadyaligned[i] = 1;
	for( i=njob-nadd; i<njob; i++ ) alreadyaligned[i] = 0;

	for( l=0; l<njob; l++ ) fftlog[l] = 1;

#if 0
	fprintf( stderr, "##### fftwinsize = %d, fftthreshold = %d\n", fftWinSize, fftThreshold );
#endif

#if 0
	for( i=0; i<njob; i++ )
		fprintf( stderr, "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif


	if( constraint )
		calcimportance( njob, effarr, aseq, localhomtable );
//		dontcalcimportance( njob, effarr, aseq, localhomtable ); // CHUUIII!!!!!


//	writePre( njob, name, nlen, aseq, 0 );

	tscore = 0.0;
	for( l=0; l<njob-1; l++ )
	{
//		fprintf( stderr, "\ndistfromtip = %f\n", dep[l].distfromtip );
		makedynamicmtx( dynamicmtx, n_dis_consweight_multi, dep[l].distfromtip );
//		makedynamicmtx( dynamicmtx, n_dis_consweight_multi, ( dep[l].distfromtip - 0.2 ) * 3 );
		if( mergeoralign[l] == 'n' )
		{
//			fprintf( stderr, "SKIP!\n" );
			free( topol[l][0] );
			free( topol[l][1] );
			free( topol[l] );
			continue;
		}

		m1 = topol[l][0][0];
		m2 = topol[l][1][0];
        len1 = strlen( aseq[m1] );
        len2 = strlen( aseq[m2] );
        if( *alloclen < len1 + len2 )
        {
			fprintf( stderr, "\nReallocating.." );
			*alloclen = ( len1 + len2 ) + 1000;
			ReallocateCharMtx( aseq, njob, *alloclen + 10  );
			gaplen = realloc( gaplen, ( *alloclen + 10 ) * sizeof( int ) );
			if( gaplen == NULL )
			{
				fprintf( stderr, "Cannot realloc gaplen\n" );
				exit( 1 );
			}
			gapmap = realloc( gapmap, ( *alloclen + 10 ) * sizeof( int ) );
			if( gapmap == NULL )
			{
				fprintf( stderr, "Cannot realloc gapmap\n" );
				exit( 1 );
			}
			fprintf( stderr, "done. *alloclen = %d\n", *alloclen );
		}

		if( effarr_kozo )
		{
			clus1 = fastconjuction_noname_kozo( topol[l][0], aseq, mseq1, effarr1, effarr, effarr1_kozo, effarr_kozo, indication1 );
			clus2 = fastconjuction_noname_kozo( topol[l][1], aseq, mseq2, effarr2, effarr, effarr2_kozo, effarr_kozo, indication2 );
		}
		else
		{
			clus1 = fastconjuction_noname( topol[l][0], aseq, mseq1, effarr1, effarr, indication1 );
			clus2 = fastconjuction_noname( topol[l][1], aseq, mseq2, effarr2, effarr, indication2 );
		}

		if( mergeoralign[l] == '1' || mergeoralign[l] == '2' ) // only in serial version
		{
			newgapstr = "=";
		}
		else
			newgapstr = "-";


		len1nocommongap = len1;
		len2nocommongap = len2;
		if( mergeoralign[l] == '1' ) // nai
		{
			findcommongaps( clus2, mseq2, gapmap );
			commongappick( clus2, mseq2 );
			len2nocommongap = strlen( mseq2[0] );
		}
		else if( mergeoralign[l] == '2' )
		{
			findcommongaps( clus1, mseq1, gapmap );
			commongappick( clus1, mseq1 );
			len1nocommongap = strlen( mseq1[0] );
		}
		

		fprintf( trap_g, "\nSTEP-%d\n", l );
		fprintf( trap_g, "group1 = %s\n", indication1 );
		fprintf( trap_g, "group2 = %s\n", indication2 );

#if 1
		fprintf( stderr, "\rSTEP % 5d /%d ", l+1, njob-1 );
		fflush( stderr );
#else
		fprintf( stdout, "STEP %d /%d\n", l+1, njob-1 );
		fprintf( stderr, "STEP %d /%d\n", l+1, njob-1 );
		fprintf( stderr, "group1 = %.66s", indication1 );
		if( strlen( indication1 ) > 66 ) fprintf( stderr, "..." );
		fprintf( stderr, "\n" );
		fprintf( stderr, "group2 = %.66s", indication2 );
		if( strlen( indication2 ) > 66 ) fprintf( stderr, "..." );
		fprintf( stderr, "\n" );
#endif



//		for( i=0; i<clus1; i++ ) fprintf( stderr, "## STEP%d-eff for mseq1-%d %f\n", l+1, i, effarr1[i] );

		if( constraint )
		{
			fastshrinklocalhom( topol[l][0], topol[l][1], localhomtable, localhomshrink );
//			msfastshrinklocalhom( topol[l][0], topol[l][1], localhomtable, localhomshrink );
//			fprintf( stderr, "localhomshrink =\n" );
//			outlocalhompt( localhomshrink, clus1, clus2 );
//			weightimportance4( clus1, clus2, effarr1, effarr2, localhomshrink );
//			fprintf( stderr, "after weight =\n" );
//			outlocalhompt( localhomshrink, clus1, clus2 );
		}
		if( rnakozo && rnaprediction == 'm' )
		{
			makegrouprna( grouprna1, singlerna, topol[l][0] );
			makegrouprna( grouprna2, singlerna, topol[l][1] );
		}


/*
		fprintf( stderr, "before align all\n" );
		display( aseq, njob );
		fprintf( stderr, "\n" );
		fprintf( stderr, "before align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		fprintf( stderr, "\n" );
		fprintf( stderr, "before align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		fprintf( stderr, "\n" );
*/

		if( !nevermemsave && ( constraint != 2  && alg != 'M'  && ( len1 > 30000 || len2 > 30000 ) ) )
		{
			fprintf( stderr, "\nlen1=%d, len2=%d, Switching to the memsave mode.\n", len1, len2 );
			alg = 'M';
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			commonAlloc1 = 0;
			commonAlloc2 = 0;
		}
		

//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000 );
		else						   ffttry = 0;
//		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 5000 && clus2 < 5000 ); // v6.708
//		fprintf( stderr, "f=%d, len1/fftlog[m1]=%f, clus1=%d, len2/fftlog[m2]=%f, clus2=%d\n", ffttry, (float)len1/fftlog[m1], clus1, (float)len2/fftlog[m2], clus2 );
//		fprintf( stderr, "f=%d, clus1=%d, fftlog[m1]=%d, clus2=%d, fftlog[m2]=%d\n", ffttry, clus1, fftlog[m1], clus2, fftlog[m2] );
		if( constraint == 2 )
		{
			if( alg == 'M' )
			{
				fprintf( stderr, "\n\nMemory saving mode is not supported.\n\n" );
				exit( 1 );
			}
			fprintf( stderr, "c" );
			if( alg == 'A' )
			{
				imp_match_init_strict( NULL, clus1, clus2, strlen( mseq1[0] ), strlen( mseq2[0] ), mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );
				if( rnakozo ) imp_rna( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, NULL, NULL, NULL );
				pscore = A__align( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, localhomshrink, &dumfl, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
			}
			else if( alg == 'H' )
			{
				imp_match_init_strictH( NULL, clus1, clus2, strlen( mseq1[0] ), strlen( mseq2[0] ), mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
				pscore = H__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, localhomshrink, &dumfl, NULL, NULL, NULL, NULL );
			}
			else if( alg == 'Q' )
			{
				imp_match_init_strictQ( NULL, clus1, clus2, strlen( mseq1[0] ), strlen( mseq2[0] ), mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
				if( rnakozo ) imp_rnaQ( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, NULL, NULL, NULL );
				pscore = Q__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, localhomshrink, &dumfl, NULL, NULL, NULL, NULL );
			}
			else if( alg == 'R' )
			{
				imp_match_init_strictR( NULL, clus1, clus2, strlen( mseq1[0] ), strlen( mseq2[0] ), mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
				pscore = R__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, localhomshrink, &dumfl, NULL, NULL, NULL, NULL );
			}
		}
		else if( force_fft || ( use_fft && ffttry ) )
		{
			fprintf( stderr, "f" );
			if( alg == 'M' )
			{
				fprintf( stderr, "m" );
				pscore = Falign_udpari_long( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			}
			else
				pscore = Falign( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1, NULL, 0, NULL );
		}
		else
		{
			fprintf( stderr, "d" );
			fftlog[m1] = 0;
			switch( alg )
			{
				case( 'a' ):
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen );
					break;
				case( 'M' ):
					fprintf( stderr, "m" );
					pscore = MSalignmm( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					break;
				case( 'A' ):
					pscore = A__align( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					break;
				case( 'Q' ):
					pscore = Q__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL );
					break;
				case( 'R' ):
					pscore = R__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL );
					break;
				case( 'H' ):
					pscore = H__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL );
					break;
				default:
					ErrorExit( "ERROR IN SOURCE FILE" );
			}
		}

		nlen[m1] = 0.5 * ( nlen[m1] + nlen[m2] );

#if SCOREOUT
		fprintf( stderr, "score = %10.2f\n", pscore );
#endif
		tscore += pscore;
/*
		fprintf( stderr, "after align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		fprintf( stderr, "\n" );
		fprintf( stderr, "after align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		fprintf( stderr, "\n" );
*/

//		writePre( njob, name, nlen, aseq, 0 );

		if( disp ) display( aseq, njob );

		if( mergeoralign[l] == '1' ) // jissainiha nai. atarashii hairetsu ha saigo dakara.
		{
			adjustgapmap( strlen( mseq2[0] )-len2nocommongap+len2, gapmap, mseq2[0] );
			restorecommongaps( njob, aseq, topol[l][0], topol[l][1], gapmap, *alloclen, '-' );
			findnewgaps( clus2, 0, mseq2, gaplen );
			insertnewgaps( njob, alreadyaligned, aseq, topol[l][1], topol[l][0], gaplen, gapmap, *alloclen, alg, '-' );
			for( i=0; i<njob; i++ ) eq2dash( aseq[i] );
			for( i=0; (m=topol[l][0][i])>-1; i++ ) alreadyaligned[m] = 1;
		}
		if( mergeoralign[l] == '2' )
		{
//			fprintf( stderr, ">mseq1[0] = \n%s\n", mseq1[0] );
//			fprintf( stderr, ">mseq2[0] = \n%s\n", mseq2[0] );
			adjustgapmap( strlen( mseq1[0] )-len1nocommongap+len1, gapmap, mseq1[0] );
			restorecommongaps( njob, aseq, topol[l][0], topol[l][1], gapmap, *alloclen, '-' );
			findnewgaps( clus1, 0, mseq1, gaplen );
			insertnewgaps( njob, alreadyaligned, aseq, topol[l][0], topol[l][1], gaplen, gapmap, *alloclen, alg, '-' );
			for( i=0; i<njob; i++ ) eq2dash( aseq[i] );
			for( i=0; (m=topol[l][1][i])>-1; i++ ) alreadyaligned[m] = 1;
		}

		free( topol[l][0] );
		free( topol[l][1] );
		free( topol[l] );
	}
#if SCOREOUT
	fprintf( stderr, "totalscore = %10.2f\n\n", tscore );
#endif
	if( rnakozo && rnaprediction == 'm' )
	{
		if( grouprna1 ) free( grouprna1 ); // nakami ha?
		if( grouprna2 ) free( grouprna2 ); // nakami ha?
		grouprna1 = grouprna2 = NULL;
	}
	if( constraint )
	{
		if( localhomshrink ) // nen no tame
		{
			for( i=0; i<njob; i++ )
			{
				free( localhomshrink[i] );
				localhomshrink[i] = NULL;
			}
			free( localhomshrink );
			localhomshrink = NULL;
		}
	}
}

static void WriteOptions( FILE *fp )
{

	if( dorp == 'd' ) fprintf( fp, "DNA\n" );
	else
	{
		if     ( scoremtx ==  0 ) fprintf( fp, "JTT %dPAM\n", pamN );
		else if( scoremtx ==  1 ) fprintf( fp, "BLOSUM %d\n", nblosum );
		else if( scoremtx ==  2 ) fprintf( fp, "M-Y\n" );
	}
    fprintf( stderr, "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );
    if( use_fft ) fprintf( fp, "FFT on\n" );

	fprintf( fp, "tree-base method\n" );
	if( tbrweight == 0 ) fprintf( fp, "unweighted\n" );
	else if( tbrweight == 3 ) fprintf( fp, "clustalw-like weighting\n" );
	if( tbitr || tbweight ) 
	{
		fprintf( fp, "iterate at each step\n" );
		if( tbitr && tbrweight == 0 ) fprintf( fp, "  unweighted\n" ); 
		if( tbitr && tbrweight == 3 ) fprintf( fp, "  reversely weighted\n" ); 
		if( tbweight ) fprintf( fp, "  weighted\n" ); 
		fprintf( fp, "\n" );
	}

   	 fprintf( fp, "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );

	if( alg == 'a' )
		fprintf( fp, "Algorithm A\n" );
	else if( alg == 'A' ) 
		fprintf( fp, "Algorithm A+\n" );
	else if( alg == 'C' ) 
		fprintf( fp, "Apgorithm A+/C\n" );
	else
		fprintf( fp, "Unknown algorithm\n" );

	if( treemethod == 'X' )
		fprintf( fp, "Tree = UPGMA (mix).\n" );
	else if( treemethod == 'E' )
		fprintf( fp, "Tree = UPGMA (average).\n" );
	else if( treemethod == 'q' )
		fprintf( fp, "Tree = Minimum linkage.\n" );
	else
		fprintf( fp, "Unknown tree.\n" );

    if( use_fft )
    {
        fprintf( fp, "FFT on\n" );
        if( dorp == 'd' )
            fprintf( fp, "Basis : 4 nucleotides\n" );
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
	else
        fprintf( fp, "FFT off\n" );
	fflush( fp );
}
	 

int main( int argc, char *argv[] )
{
	static int  *nlen;	
	static float *selfscore;
	int nogaplen;
	static char **name, **seq;
	static char **mseq1, **mseq2;
	static char **bseq;
	static float **iscore, **iscore_kozo;
	static double *eff, *eff_kozo, *eff_kozo_mapped = NULL;
	int i, j, ien, ik, jk;
	static int ***topol, ***topol_kozo;
	static int *addmem;
	static Treedep *dep;
	static float **len, **len_kozo;
	FILE *prep;
	FILE *infp;
	FILE *orderfp;
	FILE *hat2p;
	double unweightedspscore;
	int alignmentlength;
	char *mergeoralign;
	int foundthebranch;
	int nsubalignments, maxmem;
	int **subtable;
	int *insubtable;
	int *preservegaps;
	char ***subalnpt;

	char c;
	int alloclen;
	LocalHom **localhomtable = NULL;
	RNApair ***singlerna = NULL;
	float ssi, ssj, bunbo;
	static char *kozoarivec;
	int nkozo;

	arguments( argc, argv );
#ifndef enablemultithread
	nthread = 0;
#endif

	if( fastathreshold < 0.0001 ) constraint = 0; 

	if( inputfile )
	{
		infp = fopen( inputfile, "r" );
		if( !infp ) 
		{
			fprintf( stderr, "Cannot open %s\n", inputfile );
			exit( 1 );
		}
	}
	else    
		infp = stdin;

	getnumlen( infp );
	rewind( infp );


	nkozo = 0;

	if( njob < 2 )
	{
		fprintf( stderr, "At least 2 sequences should be input!\n"
						 "Only %d sequence found.\n", njob ); 
		exit( 1 );
	}

	if( subalignment )
	{
		readsubalignmentstable( njob, NULL, NULL, &nsubalignments, &maxmem );
		fprintf( stderr, "nsubalignments = %d\n", nsubalignments );
		fprintf( stderr, "maxmem = %d\n", maxmem );
		subtable = AllocateIntMtx( nsubalignments, maxmem+1 );
		insubtable = AllocateIntVec( njob );
		for( i=0; i<njob; i++ ) insubtable[i] = 0;
		preservegaps = AllocateIntVec( njob );
		for( i=0; i<njob; i++ ) preservegaps[i] = 0;
		subalnpt = AllocateCharCub( nsubalignments, maxmem, 0 );
		readsubalignmentstable( njob, subtable, preservegaps, NULL, NULL );
	}

	seq = AllocateCharMtx( njob, nlenmax+1 );
	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );

	name = AllocateCharMtx( njob, B+1 );
	nlen = AllocateIntVec( njob );
	selfscore = AllocateFloatVec( njob );

	topol = AllocateIntCub( njob, 2, 0 );
	len = AllocateFloatMtx( njob, 2 );
	iscore = AllocateFloatHalfMtx( njob );
	eff = AllocateDoubleVec( njob );
	kozoarivec = AllocateCharVec( njob );

	mergeoralign = AllocateCharVec( njob );

	dep = (Treedep *)calloc( njob, sizeof( Treedep ) );
	if( nadd ) addmem = AllocateIntVec( nadd+1 );

	if( constraint )
	{
		localhomtable = (LocalHom **)calloc( njob, sizeof( LocalHom *) );
		for( i=0; i<njob; i++ )
		{
			localhomtable[i] = (LocalHom *)calloc( njob, sizeof( LocalHom ) );
			for( j=0; j<njob; j++ )
			{
				localhomtable[i][j].start1 = -1;
				localhomtable[i][j].end1 = -1;
				localhomtable[i][j].start2 = -1;
				localhomtable[i][j].end2 = -1;
				localhomtable[i][j].overlapaa = -1.0;
				localhomtable[i][j].opt = -1.0;
				localhomtable[i][j].importance = -1.0;
				localhomtable[i][j].next = NULL;
				localhomtable[i][j].korh = 'h';
			}
		}

		fprintf( stderr, "Loading 'hat3' ... " );
		prep = fopen( "hat3", "r" );
		if( prep == NULL ) ErrorExit( "Make hat3." );
		readlocalhomtable( prep, njob, localhomtable, kozoarivec );
		fclose( prep );
		fprintf( stderr, "\ndone.\n" );


		nkozo = 0;
		for( i=0; i<njob; i++ ) 
		{
//			fprintf( stderr, "kozoarivec[%d] = %d\n", i, kozoarivec[i] );
			if( kozoarivec[i] ) nkozo++;
		}
		if( nkozo )
		{
			topol_kozo = AllocateIntCub( nkozo, 2, 0 );
			len_kozo = AllocateFloatMtx( nkozo, 2 );
			iscore_kozo = AllocateFloatHalfMtx( nkozo );
			eff_kozo = AllocateDoubleVec( nkozo );
			eff_kozo_mapped = AllocateDoubleVec( njob );
		}


//		outlocalhom( localhomtable, njob );

#if 0
		fprintf( stderr, "Extending localhom ... " );
		extendlocalhom2( njob, localhomtable );
		fprintf( stderr, "done.\n" );
#endif
	}

#if 0
	readData( infp, name, nlen, seq );
#else
	readData_pointer( infp, name, nlen, seq );
	fclose( infp );
#endif

	constants( njob, seq );

#if 0
	fprintf( stderr, "params = %d, %d, %d\n", penalty, penalty_ex, offset );
#endif

	initSignalSM();

	initFiles();

	WriteOptions( trap_g );

	c = seqcheck( seq );
	if( c )
	{
		fprintf( stderr, "Illegal character %c\n", c );
		exit( 1 );
	}

//	writePre( njob, name, nlen, seq, 0 );

	if( treein )
	{
#if 0
		if( nkozo )
		{
			fprintf( stderr, "Both structure and user tree have been given. Not yet supported!\n" );
			exit( 1 );
		}
#endif
		fprintf( stderr, "Loading a tree ... " );
		loadtree( njob, topol, len, name, nlen, dep );
		fprintf( stderr, "\ndone.\n\n" );
	}
	else
	{
		if( tbutree == 0 )
		{
			for( i=1; i<njob; i++ ) 
			{
				if( nlen[i] != nlen[0] ) 
				{
					fprintf( stderr, "Input pre-aligned seqences or make hat2.\n" );
					exit( 1 );
				}
			}
	
			fprintf( stderr, "Making a distance matrix .. \n" );
			fflush( stderr );
			ien = njob-1;
			for( i=0; i<njob; i++ ) 
			{
				selfscore[i] = naivepairscore11( seq[i], seq[i], penalty );
			}
#ifdef enablemultithread
			if( nthread > 0 )
			{
				distancematrixthread_arg_t *targ;
				Jobtable jobpos;
				pthread_t *handle;
				pthread_mutex_t mutex;

				jobpos.i = 0;
				jobpos.j = 0;

				targ = calloc( nthread, sizeof( distancematrixthread_arg_t ) );
				handle = calloc( nthread, sizeof( pthread_t ) );
				pthread_mutex_init( &mutex, NULL );

				for( i=0; i<nthread; i++ )
				{
					targ[i].thread_no = i;
					targ[i].njob = njob;
					targ[i].selfscore = selfscore;
					targ[i].iscore = iscore;
					targ[i].seq = seq;
					targ[i].jobpospt = &jobpos;
					targ[i].mutex = &mutex;

					pthread_create( handle+i, NULL, distancematrixthread, (void *)(targ+i) );
				}

				for( i=0; i<nthread; i++ )
				{
					pthread_join( handle[i], NULL );
				}
				pthread_mutex_destroy( &mutex );
				free( handle );
				free( targ );
			}
			else
#endif
			{
				for( i=0; i<ien; i++ ) 
				{
					if( i % 10 == 0 )
					{
						fprintf( stderr, "\r% 5d / %d", i, ien );
						fflush( stderr );
					}
					ssi = selfscore[i];
					for( j=i+1; j<njob; j++ ) 
					{
						ssj = selfscore[j];
						bunbo = MIN( ssi, ssj );
						if( bunbo == 0.0 )
							iscore[i][j-i] = 2.0; // 2013/Oct/17 2bai
						else
//							iscore[i][j-i] = 1.0 - naivepairscore11( seq[i], seq[j], penalty ) / MIN( selfscore[i], selfscore[j] );
							iscore[i][j-i] = ( 1.0 - naivepairscore11( seq[i], seq[j], penalty ) / bunbo ) * 2.0; // 2013/Oct/17 2bai
		
#if 0
						fprintf( stderr, "### ssj = %f\n", ssj );
						fprintf( stderr, "### selfscore[i] = %f\n", selfscore[i] );
						fprintf( stderr, "### selfscore[j] = %f\n", selfscore[j] );
						fprintf( stderr, "### rawscore = %f\n", naivepairscore11( seq[i], seq[j], penalty ) );
#endif
					}
				}
			}
			fprintf( stderr, "\ndone.\n\n" );
			fflush( stderr );
		}
		else
		{
			if( multidist )
			{
				fprintf( stderr, "Loading 'hat2n' (aligned sequences - new sequences) ... " );
				prep = fopen( "hat2n", "r" );
				if( prep == NULL ) ErrorExit( "Make hat2." );
				readhat2_floathalf_pointer( prep, njob, name, iscore );
				fclose( prep );
				fprintf( stderr, "done.\n" );
			
				fprintf( stderr, "Loading 'hat2i' (aligned sequences) ... " );
				prep = fopen( "hat2i", "r" );
				if( prep == NULL ) ErrorExit( "Make hat2i." );
				readhat2_floathalf_pointer( prep, njob-nadd, name, iscore );
				fclose( prep );
				fprintf( stderr, "done.\n" );
			}
			else
			{
				fprintf( stderr, "Loading 'hat2' ... " );
				prep = fopen( "hat2", "r" );
				if( prep == NULL ) ErrorExit( "Make hat2." );
				readhat2_floathalf_pointer( prep, njob, name, iscore );
				fclose( prep );
				fprintf( stderr, "done.\n" );
			}
		}
#if 1
		if( distout )
		{
			hat2p = fopen( "hat2", "w" );
			WriteFloatHat2_pointer_halfmtx( hat2p, njob, name, iscore );
			fclose( hat2p );
		}
#endif
		if( nkozo )
		{
			ien = njob-1;
			ik = 0;
			for( i=0; i<ien; i++ )
			{
				jk = ik+1;
				for( j=i+1; j<njob; j++ ) 
				{
					if( kozoarivec[i] && kozoarivec[j] )
					{
						iscore_kozo[ik][jk-ik] = iscore[i][j-i];
					}
					if( kozoarivec[j] ) jk++;
				}
				if( kozoarivec[i] ) ik++;
			}
		}

		fprintf( stderr, "Constructing a UPGMA tree ... " );
		fflush( stderr );
		if( topin )
		{
			fprintf( stderr, "--topin has been disabled\n" );
			exit( 1 );
//			fprintf( stderr, "Loading a topology ... " );
//			loadtop( njob, iscore, topol, len );
//			fprintf( stderr, "\ndone.\n\n" );
		}
		else if( subalignment ) // merge error no tame
		{
			fixed_supg_float_realloc_nobk_halfmtx_treeout_constrained( njob, iscore, topol, len, name, nlen, dep, nsubalignments, subtable );
		}
		else if( treeout ) // merge error no tame
		{
			fixed_musclesupg_float_realloc_nobk_halfmtx_treeout( njob, iscore, topol, len, name, nlen, dep );
		}
		else
		{
			fixed_musclesupg_float_realloc_nobk_halfmtx( njob, iscore, topol, len, dep, 1 );
		}
//		else 
//			ErrorExit( "Incorrect tree\n" );

		if( nkozo )
		{
//			for( i=0; i<nkozo-1; i++ )
//				for( j=i+1; j<nkozo; j++ )
//					fprintf( stderr, "iscore_kozo[%d][%d] =~ %f\n", i, j, iscore_kozo[i][j-i] );
			fixed_musclesupg_float_realloc_nobk_halfmtx( nkozo, iscore_kozo, topol_kozo, len_kozo, NULL, 1 );
		}
		fprintf( stderr, "\ndone.\n\n" );
		fflush( stderr );
	}


	orderfp = fopen( "order", "w" );
	if( !orderfp )
	{
		fprintf( stderr, "Cannot open 'order'\n" );
		exit( 1 );
	}
	for( i=0; (j=topol[njob-2][0][i])!=-1; i++ )
	{
		fprintf( orderfp, "%d\n", j );
	}
	for( i=0; (j=topol[njob-2][1][i])!=-1; i++ )
	{
		fprintf( orderfp, "%d\n", j );
	}
	fclose( orderfp );

	if( treeout && noalign ) 
	{
		writeData_pointer( prep_g, njob, name, nlen, seq );
		fprintf( stderr, "\n" ); 
		SHOWVERSION;
		return( 0 );
	}

//	countnode( njob, topol, node0 );
	if( tbrweight )
	{
		weight = 3; 
#if 0
		utree = 0; counteff( njob, topol, len, eff ); utree = 1;
#else
		counteff_simple_float_nostatic( njob, topol, len, eff );
		for( i=njob-nadd; i<njob; i++ ) eff[i] /= (double)100;
#if 0
		fprintf( stderr, "######  WEIGHT = \n" );
		for( i=0; i<njob; i++ )
		{
			fprintf( stderr, "w[%d] = %f\n", i, eff[i] );
		}
		exit( 1 );
#endif
		if( nkozo )
		{
//			counteff_simple_float( nkozo, topol_kozo, len_kozo, eff_kozo ); // single weight nanode iranai
			for( i=0,j=0; i<njob; i++ )
			{
				if( kozoarivec[i] )
				{
//					eff_kozo_mapped[i] = eff_kozo[j]; //
					eff_kozo_mapped[i] = eff[i];      // single weight
					j++;
				}
				else
					eff_kozo_mapped[i] = 0.0;
//				fprintf( stderr, "eff_kozo_mapped[%d] = %f\n", i, eff_kozo_mapped[i] );
//				fprintf( stderr, "            eff[%d] = %f\n", i, eff[i] );
			}
		}


#endif
	}
	else
	{
		for( i=0; i<njob; i++ ) eff[i] = 1.0;
		if( nkozo ) 
		{
			for( i=0; i<njob; i++ ) 
			{
				if( kozoarivec[i] ) 
					eff_kozo_mapped[i] = 1.0;
				else
					eff_kozo_mapped[i] = 0.0;
			}
		}
	}

	FreeFloatHalfMtx( iscore, njob );
	FreeFloatMtx( len );

	alloclen = nlenmax*2+1; //chuui!
	bseq = AllocateCharMtx( njob, alloclen );


	if( nadd )
	{
		alignmentlength = strlen( seq[0] );
		for( i=0; i<njob-nadd; i++ )
		{
			if( alignmentlength != strlen( seq[i] ) )
			{
				fprintf( stderr, "#################################################################################\n" );
				fprintf( stderr, "# ERROR!                                                                        #\n" );
				fprintf( stderr, "# The original%4d sequences must be aligned                                    #\n", njob-nadd );
				fprintf( stderr, "#################################################################################\n" );
				exit( 1 );
			}
		}
		if( addprofile )
		{
			alignmentlength = strlen( seq[njob-nadd] );
			for( i=njob-nadd; i<njob; i++ )
			{
				if( alignmentlength != strlen( seq[i] ) )
				{
					fprintf( stderr, "###############################################################################\n" );
					fprintf( stderr, "# ERROR!                                                                      #\n" );
					fprintf( stderr, "# The%4d additional sequences must be aligned                                #\n", nadd );
					fprintf( stderr, "# Otherwise, try the '--add' option, instead of '--addprofile' option.        #\n" );
					fprintf( stderr, "###############################################################################\n" );
					exit( 1 );
				}
			}
			for( i=0; i<nadd; i++ ) addmem[i] = njob-nadd+i;
			addmem[nadd] = -1;
			foundthebranch = 0;
			for( i=0; i<njob-1; i++ )
			{
				if( samemember( topol[i][0], addmem ) ) // jissainiha nai
				{
					mergeoralign[i] = '1';
					foundthebranch = 1;
				}
				else if( samemember( topol[i][1], addmem ) )
				{
					mergeoralign[i] = '2';
					foundthebranch = 1;
				}
				else
				{
					mergeoralign[i] = 'n';
				}
			}
			if( !foundthebranch )
			{
				fprintf( stderr, "###############################################################################\n" );
				fprintf( stderr, "# ERROR!                                                                      #\n" );
				fprintf( stderr, "# There is no appropriate position to add the%4d sequences in the guide tree.#\n", nadd );
				fprintf( stderr, "# Check whether the%4d sequences form a monophyletic cluster.                #\n", nadd );
				fprintf( stderr, "# If not, try the '--add' option, instead of the '--addprofile' option.       #\n" );
				fprintf( stderr, "############################################################################### \n" );
				exit( 1 );
			}
			commongappick( nadd, seq+njob-nadd );
			for( i=njob-nadd; i<njob; i++ ) strcpy( bseq[i], seq[i] );
		}
		else
		{
			for( i=0; i<njob-1; i++ ) mergeoralign[i] = 'n';
			for( j=njob-nadd; j<njob; j++ )
			{
				addmem[0] = j;
				addmem[1] = -1;
				for( i=0; i<njob-1; i++ )
				{
					if( samemember( topol[i][0], addmem ) ) // arieru
					{
//						fprintf( stderr, "HIT!\n" );
						if( mergeoralign[i] != 'n' ) mergeoralign[i] = 'w';
						else mergeoralign[i] = '1';
					}
					else if( samemember( topol[i][1], addmem ) )
					{
//						fprintf( stderr, "HIT!\n" );
						if( mergeoralign[i] != 'n' ) mergeoralign[i] = 'w';
						else mergeoralign[i] = '2';
					}
				}
			}
	
			for( i=0; i<nadd; i++ ) addmem[i] = njob-nadd+i;
			addmem[nadd] = -1;
			for( i=0; i<njob-1; i++ )
			{
				if( includemember( topol[i][0], addmem ) && includemember( topol[i][1], addmem ) )
				{
					mergeoralign[i] = 'w';
				}
				else if( includemember( topol[i][0], addmem ) )
				{
					mergeoralign[i] = '1';
				}
				else if( includemember( topol[i][1], addmem ) )
				{
					mergeoralign[i] = '2';
				}
			}
#if 0
			for( i=0; i<njob-1; i++ )
			{
				fprintf( stderr, "mem0 = " );
				for( j=0; topol[i][0][j]>-1; j++ )	fprintf( stderr, "%d ", topol[i][0][j] );
				fprintf( stderr, "\n" );
				fprintf( stderr, "mem1 = " );
				for( j=0; topol[i][1][j]>-1; j++ )	fprintf( stderr, "%d ", topol[i][1][j] );
				fprintf( stderr, "\n" );
				fprintf( stderr, "i=%d, mergeoralign[] = %c\n", i, mergeoralign[i] );
			}
#endif
			for( i=njob-nadd; i<njob; i++ ) gappick0( bseq[i], seq[i] );
		}

		commongappick( njob-nadd, seq );
		for( i=0; i<njob-nadd; i++ ) strcpy( bseq[i], seq[i] );
	}
//--------------- kokokara ----
	else if( subalignment )
	{
		for( i=0; i<njob-1; i++ ) mergeoralign[i] = 'a';
		for( i=0; i<nsubalignments; i++ )
		{
			fprintf( stderr, "Checking subalignment %d:\n", i+1 );
			alignmentlength = strlen( seq[subtable[i][0]] );
//			for( j=0; subtable[i][j]!=-1; j++ )
//				fprintf( stderr, " %d. %-30.30s\n", subtable[i][j]+1, name[subtable[i][j]]+1 );
			for( j=0; subtable[i][j]!=-1; j++ )
			{
				if( subtable[i][j] >= njob )
				{
					fprintf( stderr, "No such sequence, %d.\n", subtable[i][j]+1 );
					exit( 1 );
				}
				if( alignmentlength != strlen( seq[subtable[i][j]] ) )
				{
					fprintf( stderr, "\n" );
					fprintf( stderr, "###############################################################################\n" );
					fprintf( stderr, "# ERROR!\n" );
					fprintf( stderr, "# Subalignment %d must be aligned.\n", i+1 );
					fprintf( stderr, "# Please check the alignment lengths of following sequences.\n" );
					fprintf( stderr, "#\n" );
					fprintf( stderr, "# %d. %-10.10s -> %d letters (including gaps)\n", subtable[i][0]+1, name[subtable[i][0]]+1, alignmentlength );
					fprintf( stderr, "# %d. %-10.10s -> %d letters (including gaps)\n", subtable[i][j]+1, name[subtable[i][j]]+1, (int)strlen( seq[subtable[i][j]] ) );
					fprintf( stderr, "#\n" );
					fprintf( stderr, "# See http://mafft.cbrc.jp/alignment/software/merge.html for details.\n" );
					if( subalignmentoffset )
					{
						fprintf( stderr, "#\n" );
						fprintf( stderr, "# You specified seed alignment(s) consisting of %d sequences.\n", subalignmentoffset );
						fprintf( stderr, "# In this case, the rule of numbering is:\n" );
						fprintf( stderr, "#   The aligned seed sequences are numbered as 1 .. %d\n", subalignmentoffset );
						fprintf( stderr, "#   The input sequences to be aligned are numbered as %d .. %d\n", subalignmentoffset+1, subalignmentoffset+njob );
					}
					fprintf( stderr, "###############################################################################\n" );
					fprintf( stderr, "\n" );
					exit( 1 );
				}
				insubtable[subtable[i][j]] = 1;
			}
			for( j=0; j<njob-1; j++ )
			{
				if( includemember( topol[j][0], subtable[i] ) && includemember( topol[j][1], subtable[i] ) )
				{
					mergeoralign[j] = 'n';
				}
			}
			foundthebranch = 0;
			for( j=0; j<njob-1; j++ )
			{
				if( samemember( topol[j][0], subtable[i] ) || samemember( topol[j][1], subtable[i] ) )
				{
					foundthebranch = 1;
					fprintf( stderr, " -> OK\n" );
					break;
				}
			}
			if( !foundthebranch )
			{
				system( "cp infile.tree GuideTree" ); // tekitou
				fprintf( stderr, "\n" );
				fprintf( stderr, "###############################################################################\n" );
				fprintf( stderr, "# ERROR!\n" );
				fprintf( stderr, "# Subalignment %d does not form a monophyletic cluster\n", i+1 );
				fprintf( stderr, "# in the guide tree ('GuideTree' in this directory) internally computed.\n" );
				fprintf( stderr, "# If you really want to use this subalignment, pelase give a tree with --treein \n" );
				fprintf( stderr, "# http://mafft.cbrc.jp/alignment/software/treein.html\n" );
				fprintf( stderr, "# http://mafft.cbrc.jp/alignment/software/merge.html\n" );
				if( subalignmentoffset )
				{
					fprintf( stderr, "#\n" );
					fprintf( stderr, "# You specified seed alignment(s) consisting of %d sequences.\n", subalignmentoffset );
					fprintf( stderr, "# In this case, the rule of numbering is:\n" );
					fprintf( stderr, "#   The aligned seed sequences are numbered as 1 .. %d\n", subalignmentoffset );
					fprintf( stderr, "#   The input sequences to be aligned are numbered as %d .. %d\n", subalignmentoffset+1, subalignmentoffset+njob );
				}
				fprintf( stderr, "############################################################################### \n" );
				fprintf( stderr, "\n" );
				exit( 1 );
			}
//			commongappick( seq[subtable[i]], subalignment[i] ); // irukamo
		}
#if 0
		for( i=0; i<njob-1; i++ )
		{
			fprintf( stderr, "STEP %d\n", i+1 );
			fprintf( stderr, "group1 = " );
			for( j=0; topol[i][0][j] != -1; j++ )
				fprintf( stderr, "%d ", topol[i][0][j]+1 );
			fprintf( stderr, "\n" );
			fprintf( stderr, "group2 = " );
			for( j=0; topol[i][1][j] != -1; j++ )
				fprintf( stderr, "%d ", topol[i][1][j]+1 );
			fprintf( stderr, "\n" );
			fprintf( stderr, "%d -> %c\n\n", i, mergeoralign[i] );
		}
#endif

		for( i=0; i<njob; i++ ) 
		{
			if( insubtable[i] ) strcpy( bseq[i], seq[i] );
			else gappick0( bseq[i], seq[i] );
		}

		for( i=0; i<nsubalignments; i++ ) 
		{
			for( j=0; subtable[i][j]!=-1; j++ ) subalnpt[i][j] = bseq[subtable[i][j]];
			if( !preservegaps[i] ) commongappick( j, subalnpt[i] );
		}

		FreeIntMtx( subtable );
		free( insubtable );
		for( i=0; i<nsubalignments; i++ ) free( subalnpt[i] );
		free( subalnpt );
		free( preservegaps );
	}
//--------------- kokomade ----
	else
	{
		for( i=0; i<njob; i++ ) gappick0( bseq[i], seq[i] );
		for( i=0; i<njob-1; i++ ) mergeoralign[i] = 'a';
	}

	if( rnakozo && rnaprediction == 'm' )
	{
		singlerna = (RNApair ***)calloc( njob, sizeof( RNApair ** ) );
		prep = fopen( "hat4", "r" );
		if( prep == NULL ) ErrorExit( "Make hat4 using mccaskill." );
		fprintf( stderr, "Loading 'hat4' ... " );
		for( i=0; i<njob; i++ )
		{
			nogaplen = strlen( bseq[i] );
			singlerna[i] = (RNApair **)calloc( nogaplen+1, sizeof( RNApair * ) );
			for( j=0; j<nogaplen; j++ )
			{
				singlerna[i][j] = (RNApair *)calloc( 1, sizeof( RNApair ) );
				singlerna[i][j][0].bestpos = -1;
				singlerna[i][j][0].bestscore = -1.0;
			}
			singlerna[i][nogaplen] =  NULL;
//			fprintf( stderr, "### reading bpp %d ...\n", i );
			readmccaskill( prep, singlerna[i], nogaplen );
		}
		fclose( prep );
		fprintf( stderr, "\ndone.\n" );
	}
	else
		singlerna = NULL;


	fprintf( stderr, "Progressive alignment ... \n" );

#ifdef enablemultithread
	if( nthread > 0 && nadd == 0 )
	{
		treebasethread_arg_t *targ;	
		int jobpos;
		pthread_t *handle;
		pthread_mutex_t mutex;
		pthread_cond_t treecond;
		int *fftlog;
		int nrun;
		int nthread_yoyu;

		nthread_yoyu = nthread * 1;
		nrun = 0;
		jobpos = 0;
		targ = calloc( nthread_yoyu, sizeof( treebasethread_arg_t ) );
		fftlog = AllocateIntVec( njob );
		handle = calloc( nthread_yoyu, sizeof( pthread_t ) );
		pthread_mutex_init( &mutex, NULL );
		pthread_cond_init( &treecond, NULL );

		for( i=0; i<njob; i++ ) dep[i].done = 0;
		for( i=0; i<njob; i++ ) fftlog[i] = 1;

		if( constraint )
			calcimportance( njob, eff, bseq, localhomtable );
//			dontcalcimportance( njob, eff, bseq, localhomtable ); // CHUUUUIIII!!!

		for( i=0; i<nthread_yoyu; i++ )
		{
			targ[i].thread_no = i;
			targ[i].nrunpt = &nrun;
			targ[i].njob = njob;
			targ[i].nlen = nlen;
			targ[i].jobpospt = &jobpos;
			targ[i].topol = topol;
			targ[i].dep = dep;
			targ[i].aseq = bseq;
			targ[i].effarr = eff;
			targ[i].alloclenpt = &alloclen;
			targ[i].localhomtable = localhomtable;
			targ[i].singlerna = singlerna;
			targ[i].effarr_kozo = eff_kozo_mapped;
			targ[i].fftlog = fftlog;
			targ[i].mergeoralign = mergeoralign;
			targ[i].mutex = &mutex;
			targ[i].treecond = &treecond;

			pthread_create( handle+i, NULL, treebasethread, (void *)(targ+i) );
		}

		for( i=0; i<nthread_yoyu; i++ )
		{
			pthread_join( handle[i], NULL );
		}
		pthread_mutex_destroy( &mutex );
		pthread_cond_destroy( &treecond );
		free( handle );
		free( targ );
		free( fftlog );
	}
	else
#endif

		treebase( nlen, bseq, nadd, mergeoralign, mseq1, mseq2, topol, dep, eff, &alloclen, localhomtable, singlerna, eff_kozo_mapped );
	fprintf( stderr, "\ndone.\n" );
	if( scoreout )
	{
		unweightedspscore = plainscore( njob, bseq );
		fprintf( stderr, "\nSCORE %s = %.0f, ", "(treebase)", unweightedspscore );
		fprintf( stderr, "SCORE / residue = %f", unweightedspscore / ( njob * strlen( bseq[0] ) ) );
		fprintf( stderr, "\n\n" );
	}

#if 0
	if( constraint )
	{
		LocalHom *tmppt1, *tmppt2;
		for( i=0; i<njob; i++ )
		{
			for( j=0; j<njob; j++ )
			{
				tmppt1 = localhomtable[i]+j;
				while( tmppt2 = tmppt1->next )
				{
					free( (void *)tmppt1 );
					tmppt1 = tmppt2;
				}
				free( (void *)tmppt1 );
			}
			free( (void *)(localhomtable[i]+j) );
		}
		free( (void *)localhomtable );
	}
#endif

	fprintf( trap_g, "done.\n" );
	fclose( trap_g );
	free( mergeoralign );
	if( rnakozo && rnaprediction == 'm' ) 
	{
		if( singlerna ) // nen no tame
		{
			for( i=0; i<njob; i++ ) 
			{
				for( j=0; singlerna[i][j]!=NULL; j++ )
				{
					if( singlerna[i][j] ) free( singlerna[i][j] );
				}
				if( singlerna[i] ) free( singlerna[i] );
			}
			free( singlerna );
			singlerna = NULL;
		}
	}

	writeData_pointer( prep_g, njob, name, nlen, bseq );
#if 0
	writeData( stdout, njob, name, nlen, bseq );
	writePre( njob, name, nlen, bseq, !contin );
	writeData_pointer( prep_g, njob, name, nlen, aseq );
#endif
#if IODEBUG
	fprintf( stderr, "OSHIMAI\n" );
#endif

	if( constraint ) FreeLocalHomTable( localhomtable, njob );

	SHOWVERSION;
	return( 0 );
}
