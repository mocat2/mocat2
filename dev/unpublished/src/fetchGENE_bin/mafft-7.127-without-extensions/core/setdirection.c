#include "mltaln.h"

#define DEBUG 0

char *directionfile;

void arguments( int argc, char *argv[] )
{
    int c;

	inputfile = NULL;
	directionfile = NULL;

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( (c = *++argv[0]) )
		{
            switch( c )
            {
				case 'd':
					directionfile = *++argv;
					fprintf( stderr, "directionfile = %s\n", directionfile );
					--argc;
					goto nextoption;
				case 'i':
					inputfile = *++argv;
					fprintf( stderr, "inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
                default:
                    fprintf( stderr, "illegal option %c\n", c );
                    argc = 0;
                    break;
            }
		}
		nextoption:
			;
	}
    if( argc != 0 ) 
    {
        fprintf( stderr, "options: Check source file !\n" );
        exit( 1 );
    }
}



int main( int argc, char *argv[] )
{
	FILE *infp;
	FILE *difp;
	int nlenmin;
	char **name;
	char **seq;
	char *tmpseq;
	char line[100];
	int *nlen;
	int i;

	arguments( argc, argv );

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

	if( directionfile )
	{
		difp = fopen( directionfile, "r" );
		if( !difp )
		{
			fprintf( stderr, "Cannot open %s\n", directionfile );
			exit( 1 );
		}
	}
	else
	{
		fprintf( stderr, "Give directionfile!\n" );
	}


	dorp = NOTSPECIFIED;
	getnumlen_casepreserve( infp, &nlenmin );

	fprintf( stderr, "%d x %d - %d %c\n", njob, nlenmax, nlenmin, dorp );

	seq = AllocateCharMtx( njob, nlenmax+1 );
	tmpseq = AllocateCharVec( MAX( B, nlenmax )+1 );
	name = AllocateCharMtx( njob, B+1 );
	nlen = AllocateIntVec( njob );

	readData_pointer_casepreserve( infp, name, nlen, seq );

	for( i=0; i<njob; i++ )
	{
		fgets( line, 99, difp );
		if( line[0] != '_' )
		{
			fprintf( stderr, "Format error!\n" );
			exit( 1 );
		}
		if( line[1] == 'R' )
		{
			sreverse( tmpseq, seq[i] );
			strcpy( seq[i], tmpseq );

			strncpy( tmpseq, name[i]+1, B-3 );
			tmpseq[B-3] = 0;
			strcpy( name[i]+1, "_R_" );
			strcpy( name[i]+4, tmpseq );
		}
		else if( line[1] == 'F' )
		{
			;
		}
		else
		{
			fprintf( stderr, "Format error!\n" );
			exit( 1 );
		}
	}


	for( i=0; i<njob; i++ )
	{
		fprintf( stdout, ">%s\n", name[i]+1 );
		fprintf( stdout, "%s\n", seq[i] );
	}

	free( nlen );
	FreeCharMtx( seq );
	FreeCharMtx( name );
	free( tmpseq );

	return( 0 );
}
