#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0

#define END_OF_VEC -1

int nadd;
float thresholdtorev;
int dodp;
int addfragment;

typedef struct _thread_arg
{
	int iend; 
	char **seq;
	char *tmpseq;
	int *res;
	int **spointt;
	short *table1;
	int iq;
#ifdef enablemultithread
	int *jshare;
	int thread_no;
	pthread_mutex_t *mutex_counter;
#endif
} thread_arg_t;



void arguments( int argc, char *argv[] )
{
    int c;

	nthread = 1;
	inputfile = NULL;
	nadd = 0;
	dodp = 0;
	alg = 'a';
	alg = 'm';
	dorp = NOTSPECIFIED;
	fmodel = 0;
//	ppenalty = (int)( -2.0 * 1000 - 0.5 );
//	ppenalty_ex = (int)( -0.1 * 1000 - 0.5 );
//	poffset = (int)( 0.1 * 1000 - 0.5 ); 
	ppenalty = NOTSPECIFIED;
	ppenalty_ex = NOTSPECIFIED;
	poffset = NOTSPECIFIED;
	kimuraR = 2;
	pamN = 200;
	thresholdtorev = 0.1;
	addfragment = 0;


    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( (c = *++argv[0]) )
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
				case 'C':
					nthread = myatoi( *++argv );
					fprintf( stderr, "nthread = %d\n", nthread );
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
				case 'j':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = JTT;
					fprintf( stderr, "jtt/kimura %d\n", pamN );
					--argc;
					goto nextoption;
				case 't':
					thresholdtorev = atof( *++argv );
					fprintf( stderr, "thresholdtorev = %f\n", thresholdtorev );
					--argc; 
					goto nextoption;
				case 'd':
					dodp = 1;
					break;
				case 'F':
					addfragment = 1;
					break;
#if 1
				case 'a':
					fmodel = 1;
					break;
#endif
				case 'S':
					alg = 'S';
					break;
				case 'M':
					alg = 'M';
					break;
				case 'm':
					alg = 'm';
					break;
				case 'G':
					alg = 'G';
					break;
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
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
}




static int maxl;
static int tsize;

void seq_grp_nuc( int *grp, char *seq )
{
	int tmp;
	int *grpbk = grp;
	while( *seq )
	{
		tmp = amino_grp[(int)*seq++];
		if( tmp < 4 )
			*grp++ = tmp;
		else
//			fprintf( stderr, "WARNING : Unknown character %c\r", *(seq-1) );
			;
	}
	*grp = END_OF_VEC;
	if( grp - grpbk < 6 )
	{
//		fprintf( stderr, "\n\nWARNING: Too short.\nPlease also consider use mafft-ginsi, mafft-linsi or mafft-ginsi.\n\n\n" );
//		exit( 1 );
		*grpbk = -1;
	}
}

void seq_grp( int *grp, char *seq )
{
	int tmp;
	int *grpbk = grp;
	while( *seq )
	{
		tmp = amino_grp[(int)*seq++];
		if( tmp < 6 )
			*grp++ = tmp;
		else
//			fprintf( stderr, "WARNING : Unknown character %c\r", *(seq-1) );
			;
	}
	*grp = END_OF_VEC;
	if( grp - grpbk < 6 )
	{
//		fprintf( stderr, "\n\nWARNING: Too short.\nPlease also consider use mafft-ginsi, mafft-linsi or mafft-ginsi.\n\n\n" );
//		exit( 1 );
		*grpbk = -1;
	}
}

void makecompositiontable_p( short *table, int *pointt )
{
	int point;

	while( ( point = *pointt++ ) != END_OF_VEC )
		table[point]++;
}


void makepointtable_nuc( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ *  1024;
	point += *n++ *   256;
	point += *n++ *    64;
	point += *n++ *    16;
	point += *n++ *     4;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 1024;
		point *= 4;
		point += *n++;
		*pointt++ = point;
	}
	*pointt = END_OF_VEC;
}

void makepointtable( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ *  7776;
	point += *n++ *  1296;
	point += *n++ *   216;
	point += *n++ *    36;
	point += *n++ *     6;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 7776;
		point *= 6;
		point += *n++;
		*pointt++ = point;
	}
	*pointt = END_OF_VEC;
}

static int commonsextet_p2( short *table, int *pointt )
{
	int value = 0;
	short tmp;
	int point;
	short *memo;
	int *ct;
	int *cp;

	if( *pointt == -1 )
		return( 0 );

	memo = (short *)calloc( tsize, sizeof( short ) );
	if( !memo ) ErrorExit( "Cannot allocate memo\n" );
	ct = (int *)calloc( MIN( maxl, tsize )+1, sizeof( int ) ); // chuui!!
	if( !ct ) ErrorExit( "Cannot allocate memo\n" );

	cp = ct;
	while( ( point = *pointt++ ) != END_OF_VEC )
	{
		tmp = memo[point]++;
		if( tmp < table[point] )
			value++;
		if( tmp == 0 ) *cp++ = point;
	}
	*cp = END_OF_VEC;
	
	cp =  ct;
	while( *cp != END_OF_VEC )
		memo[*cp++] = 0;

	free( memo );
	free( ct );
	return( value );
}

static void	*directionthread( void *arg )
{
	thread_arg_t *targ = (thread_arg_t *)arg;
	int iend = targ->iend;
	char **seq = targ->seq;
	char *tmpseq = targ->tmpseq;
	int *res = targ->res;
	int **spointt = targ->spointt;
	short *table1 = targ->table1;
	int iq = targ->iq;
#ifdef enablemultithread
	int thread_no = targ->thread_no;
	int *jshare = targ->jshare; 
#endif
	int j;
	char **mseq1, **mseq2;


	if( dodp ) // nakuserukamo
	{
		mseq1 = AllocateCharMtx( 1, 0 );
		mseq2 = AllocateCharMtx( 1, 0 );
	}

	j = -1;
	while( 1 )
	{
#ifdef enablemultithread
		if( nthread )
		{
			pthread_mutex_lock( targ->mutex_counter );
			j = *jshare;
			if( j == iend )
			{
				fprintf( stderr, "\r %d / %d (thread %d)   \r", iq, njob, thread_no );
				pthread_mutex_unlock( targ->mutex_counter );
				break;
			}
			++(*jshare);
			pthread_mutex_unlock( targ->mutex_counter );
		}
		else
#endif
		{
			j++;
			if( j == iend ) 
			{
				fprintf( stderr, "\r %d / %d  \r", iq, njob );
				break;
			}
		}


		if( dodp )
		{
//			strcpy( mseq1[0], tmpseq );
//			strcpy( mseq2[0], seq[j] );
			mseq1[0] = tmpseq;
			mseq2[0] = seq[j];
//			res[j] = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, mseq1, mseq2, 0 );
			res[j] = L__align11_noalign( mseq1, mseq2 );
		}
		else
			res[j] = commonsextet_p2( table1, spointt[j] );
	}
	if( dodp ) // nakuserukamo
	{
		free( mseq1 );
		free( mseq2 );
//		G__align11_noalign( NULL, 0, 0, NULL, NULL, 0 );
		L__align11_noalign( NULL, NULL );
	}
//	else
//		if( nthread )  // inthread == 0 no toki free suru to, error. nazeda
//			commonsextet_p( NULL, NULL );
	return( NULL );
}

int main( int argc, char *argv[] )
{
	static int  *nlen;	
	static int  *nogaplen;	
	static char **name, **seq;
	int i, j, istart, iend;
	FILE *infp;
//	FILE *adfp;
	char c;

	int *grpseq;
	char *tmpseq, *revseq;
	int  **pointt, **pointt_rev, **spointt;
	float res_forward, res_reverse, res_max;
	int ires, mres, mres2;
	int *res;
	static short *table1, *table1_rev;
	static char **mseq1f, **mseq1r, **mseq2;

	arguments( argc, argv );
#ifndef enablemultithread
	nthread = 0;
#endif

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

	if( alg == 'a' )
	{
		if( nlenmax < 10000 )
			alg = 'G';
		else
			alg = 'S';
	}

	seq = AllocateCharMtx( njob, nlenmax*1+1 );

#if 0
	Read( name, nlen, seq );
	readData( infp, name, nlen, seq );
#else
    name = AllocateCharMtx( njob, B+1 );
    nlen = AllocateIntVec( njob ); 
    nogaplen = AllocateIntVec( njob ); 
	readData_pointer( infp, name, nlen, seq );
	fclose( infp );

	if( dorp != 'd' )
	{
		fprintf( stderr, "Not necessary!\n" );
		for( i=0; i<njob; i++ ) 
			fprintf( stdout, "_F_%-10.10s\n", name[i]+1 );
		exit( 1 );
	}
#endif

	constants( njob, seq );


#if 0
	fprintf( stderr, "params = %d, %d, %d\n", penalty, penalty_ex, offset );
#endif

	initSignalSM();

	initFiles();

	c = seqcheck( seq );
	if( c )
	{
		fprintf( stderr, "Illegal character %c\n", c );
		exit( 1 );
	}

	fprintf( stderr, "\n" );
	if( alg == 'G' ) // do to the first sequence
	{
		mseq1f = AllocateCharMtx( 1, nlenmax+nlenmax );
		mseq1r = AllocateCharMtx( 1, nlenmax+nlenmax );
		mseq2 = AllocateCharMtx( 1, nlenmax+nlenmax );
	    tmpseq = AllocateCharVec( MAX( nlenmax, B ) +1 );

		gappick0( mseq1f[0], seq[0] );
		sreverse( mseq1r[0], mseq1f[0] );
		strcpy( seq[0], mseq1f[0] );

		if( nadd )
			istart = njob - nadd;
		else
			istart = 1;

		fprintf( stderr, "\n" );

		for( i=0; i<istart; i++ )
		{
			gappick0( tmpseq, seq[i] );
			strcpy( seq[i], tmpseq );
			strcpy( tmpseq, name[i] );
			strcpy( name[i], "_F_" );
			strncpy( name[i]+3, tmpseq+1, 10 );
			name[i][13] = 0;
		}
		for( i=istart; i<njob; i++ ) 
		{
			gappick0( mseq2[0], seq[i] );

//			res_forward = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, mseq1f, mseq2, 0 );
//			res_reverse = G__align11_noalign( n_dis_consweight_multi, penalty, penalty_ex, mseq1r, mseq2, 0 );

			res_forward = L__align11_noalign( mseq1f, mseq2 );
			res_reverse = L__align11_noalign( mseq1r, mseq2 );
#if 0

			strcpy( mseq2[0], seq[i] );
			strcpy( mseq1f[0], seq[0] );
			res_forward = G__align11( n_dis_consweight_multi, mseq1f, mseq2, nlenmax*2, 0, 0 );
			fprintf( stdout, "%s\n", mseq1f[0] );
			fprintf( stdout, "%s\n", mseq2[0] );

			strcpy( mseq2[0], seq[i] );
			sreverse( mseq1r[0], seq[0] );
			res_reverse = G__align11( n_dis_consweight_multi, mseq1r, mseq2, nlenmax*2, 0, 0 );
			fprintf( stdout, "%s\n", mseq1r[0] );
			fprintf( stdout, "%s\n", mseq2[0] );
#endif

//			fprintf( stdout, "\nscore_for(%d,%d) = %f\n", 0, i, res_forward );
//			fprintf( stdout, "score_rev(%d,%d) = %f\n", 0, i, res_reverse );
			res_max = MAX(res_reverse,res_forward);
			if( (res_reverse-res_forward)/res_max > thresholdtorev ) // tekitou
			{
//				fprintf( stderr, "REVERSE!!!\n" );
				sreverse( seq[i], mseq2[0] );

				strcpy( tmpseq, name[i] );
				strcpy( name[i], "_R_" );
				strncpy( name[i]+3, tmpseq+1, 10 );
				name[i][13] = 0;
			}
			else
			{
				strcpy( seq[i], mseq2[0] );

				strcpy( tmpseq, name[i] );
				strcpy( name[i], "_F_" );
				strncpy( name[i]+3, tmpseq+1, 10 );
				name[i][13] = 0;
			}
		}
		FreeCharMtx( mseq1f );
		FreeCharMtx( mseq1r );
		FreeCharMtx( mseq2 );
		free( tmpseq );
	}
	else if( alg == 'm' )
	{
		if( dodp ) // nakuserukamo
		{
			mseq1f = AllocateCharMtx( 1, nlenmax+1);
			mseq1r = AllocateCharMtx( 1, nlenmax+1 );
			mseq2 = AllocateCharMtx( 1, nlenmax+1 );
		}
		else
		{
			spointt = AllocateIntMtx( njob, 0 ); 
			pointt = AllocateIntMtx( njob, nlenmax+1 );
			pointt_rev = AllocateIntMtx( njob, nlenmax+1 );
		}
	    tmpseq = AllocateCharVec( MAX( nlenmax, B ) +1 );
	    revseq = AllocateCharVec( nlenmax+1 );
		grpseq = AllocateIntVec( nlenmax+1 );
		res = AllocateIntVec( njob );
		if( dorp == 'd' ) tsize = (int)pow( 4, 6 );
		else              tsize = (int)pow( 6, 6 ); // iranai

		maxl = 0;
		for( i=0; i<njob; i++ )
		{
			gappick0( tmpseq, seq[i] );
			nogaplen[i] = strlen( tmpseq );
			if( nogaplen[i] > maxl ) maxl = nogaplen[i];
		}

		if( nadd )
			iend = njob - nadd;
		else
			iend = 1;

		for( i=0; i<iend; i++ )
		{
//			fprintf( stdout, "%d, SKIP\n", i );
			gappick0( tmpseq, seq[i] );
			strcpy( seq[i], tmpseq );
//			if( !nadd ) strcpy( seq[i], tmpseq ); // seq ha tsukawanaikara ii.

			if( !dodp )
			{
				seq_grp_nuc( grpseq, tmpseq );
				makepointtable_nuc( pointt[i], grpseq );
				spointt[i] = pointt[i];
			}

			strcpy( tmpseq, name[i] );
			strcpy( name[i], "_F_" );
			strncpy( name[i]+3, tmpseq+1, 10 );
			name[i][13] = 0;
		}

		if( nadd )
			istart = njob - nadd;
		else
			istart = 1;

		fprintf( stderr, "\n" );
		for( i=istart; i<njob; i++ ) 
		{
//			fprintf( stderr, "\r %d / %d ", i, njob );
			gappick0( tmpseq, seq[i] );
			strcpy( seq[i], tmpseq );
			sreverse( revseq, tmpseq );

			if( !dodp )
			{
				table1 = (short *)calloc( tsize, sizeof( short ) );
				if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
				table1_rev = (short *)calloc( tsize, sizeof( short ) );
				if( !table1_rev ) ErrorExit( "Cannot allocate table1_rev\n" );
				seq_grp_nuc( grpseq, tmpseq );
				makepointtable_nuc( pointt[i], grpseq );
				makecompositiontable_p( table1, pointt[i] );
				seq_grp_nuc( grpseq, revseq );
				makepointtable_nuc( pointt_rev[i], grpseq );
				makecompositiontable_p( table1_rev, pointt_rev[i] );
			}

			if( nadd && addfragment )
				iend = njob-nadd;
			else
				iend = i;


#ifdef enablemultithread
			if( nthread )
			{
				pthread_t *handle;
				pthread_mutex_t mutex_counter;
				thread_arg_t *targ;
				int *jsharept;
		
				targ = calloc( nthread, sizeof( thread_arg_t ) );
				handle = calloc( nthread, sizeof( pthread_t ) );
				pthread_mutex_init( &mutex_counter, NULL );
				jsharept = calloc( 1, sizeof(int) );
				*jsharept = 0;
		
				for( j=0; j<nthread; j++ )
				{
					targ[j].iend = iend;
					targ[j].seq = seq;
					targ[j].tmpseq = tmpseq; 
					targ[j].res = res; 
					targ[j].spointt = spointt; 
					targ[j].table1 = table1; 
					targ[j].jshare = jsharept;
					targ[j].iq = i;
					targ[j].mutex_counter = &mutex_counter;
					targ[j].thread_no = j;
					pthread_create( handle+j, NULL, directionthread, (void *)(targ+j) );
				}
				for( j=0; j<nthread; j++ ) pthread_join( handle[j], NULL );
				pthread_mutex_destroy( &mutex_counter );
				free( handle );
				free( targ );
				free( jsharept );
			}
			else
#endif
			{
				thread_arg_t *targ;
				targ = calloc( 1, sizeof( thread_arg_t ) );
				targ[0].iend = iend;
				targ[0].seq = seq;
				targ[0].tmpseq = tmpseq; 
				targ[0].res = res; 
				targ[0].spointt = spointt; 
				targ[0].table1 = table1; 
				targ[0].iq = i; 
				directionthread( targ );
				free( targ );
			}


			mres = mres2 = 0;
			for( j=0; j<iend; j++ )
			{
				ires = res[j];
//				fprintf( stdout, "ires (%d,%d) = %d\n", i, j, ires );
//				fflush( stdout );
				if( ires>mres2 ) 
				{
					if( ires>mres ) 
					{
						mres2 = mres;
						mres = ires;
					}
					else
						mres2 = ires;
				}
			}
			res_forward = (float)( mres + mres2 ) / 2;

#ifdef enablemultithread
			if( nthread )
			{
				pthread_t *handle;
				pthread_mutex_t mutex_counter;
				thread_arg_t *targ;
				int *jsharept;
		
				targ = calloc( nthread, sizeof( thread_arg_t ) );
				handle = calloc( nthread, sizeof( pthread_t ) );
				pthread_mutex_init( &mutex_counter, NULL );
				jsharept = calloc( 1, sizeof(int) );
				*jsharept = 0;
		
				for( j=0; j<nthread; j++ )
				{
					targ[j].iend = iend;
					targ[j].seq = seq;
					targ[j].tmpseq = revseq; 
					targ[j].res = res; 
					targ[j].spointt = spointt; 
					targ[j].table1 = table1_rev; 
					targ[j].jshare = jsharept;
					targ[j].iq = i;
					targ[j].mutex_counter = &mutex_counter;
					targ[j].thread_no = j;
					pthread_create( handle+j, NULL, directionthread, (void *)(targ+j) );
				}
				for( j=0; j<nthread; j++ ) pthread_join( handle[j], NULL );
				pthread_mutex_destroy( &mutex_counter );
				free( handle );
				free( targ );
				free( jsharept );
			}
			else
#endif
			{
				thread_arg_t *targ;
				targ = calloc( 1, sizeof( thread_arg_t ) );
				targ[0].iend = iend;
				targ[0].seq = seq;
				targ[0].tmpseq = revseq; 
				targ[0].res = res; 
				targ[0].spointt = spointt;
				targ[0].table1 = table1_rev; 
				targ[0].iq = i; 
				directionthread( targ );
				free( targ );
			}

			mres = mres2 = 0;
			for( j=0; j<iend; j++ )
			{
				ires = res[j];
				if( ires>mres2 )
				{
					if( ires>mres ) 
					{
						mres2 = mres;
						mres = ires;
					}
					else
						mres2 = ires;
				}
			}
			res_reverse = (float)( mres + mres2 ) / 2;

//			fprintf( stdout, "\n" );
//			fprintf( stdout, "score_for(%d,%d) = %f\n", 0, i, res_forward );
//			fprintf( stdout, "score_rev(%d,%d) = %f\n", 0, i, res_reverse );
//			fflush( stdout );
			res_max = MAX(res_reverse,res_forward);
			if( (res_reverse-res_forward)/res_max > thresholdtorev ) // tekitou

			{
				strcpy( seq[i], revseq );

				strcpy( tmpseq, name[i] );
				strcpy( name[i], "_R_" );
				strncpy( name[i]+3, tmpseq+1, 10 );
				name[i][13] = 0;
				if( !dodp ) spointt[i] = pointt_rev[i];
			}
			else
			{
				strcpy( tmpseq, name[i] );
				strcpy( name[i], "_F_" );
				strncpy( name[i]+3, tmpseq+1, 10 );
				name[i][13] = 0;
				if( !dodp ) spointt[i] = pointt[i];
			}

			if( !dodp )
			{
				free( table1 );
				free( table1_rev );
			}
		}

		free( grpseq );
		free( tmpseq );
		free( revseq );
		free( res );
		if( dodp )
		{
			FreeCharMtx( mseq1f );
			FreeCharMtx( mseq1r );
			FreeCharMtx( mseq2 );
		}
		else
		{
			FreeIntMtx( pointt );
			FreeIntMtx( pointt_rev );
			free( spointt );
		}
	}
	else
	{
		fprintf( stderr, "Unknown alg %c\n", alg );
		exit( 1 );
	}
//	writeData_pointer( stdout, njob, name, nlen, seq );
	for( i=0; i<njob; i++ ) 
	{
//		fprintf( stdout, ">%s\n", name[i] );
//		fprintf( stdout, "%s\n", seq[i] );
		fprintf( stdout, "%s\n", name[i] );
	}

	fprintf( stderr, "\n" );
	SHOWVERSION;
	return( 0 );
}

