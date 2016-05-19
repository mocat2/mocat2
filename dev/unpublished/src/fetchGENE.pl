#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Cwd;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2013
# This code is released under GNU GPL v3.

# VARIABLES
my ( $noalign, @align, $manual, @input_files, %align );
chomp( my $date = `date +\%Y\%b\%d_\%H\%M\%S` );
my $AQUA = "AQUA.1.3.tcl";
my $cwd  = getcwd;
my $bin  = "$cwd/fetchGENE_bin";

# USAGE
my $usage = "
 fetchGENE extracts selected genes from genomes and metagenomes in an easy and accurate manner.

 Help & Manual
   fetchGENE.pl -h|help

 Usage
   fetchGENE.pl -align [clustalo mafft muscle] -input [-input files-]
  
   Description:
     -align   : specify one or multiple sequence alignment programs to use
                muscle may take very long for alignments with 1000s of sequences
                clustalo is generally producing the best results
                see Muller et al (2010) Bioinformatics for complete details 
     -input   : specify one or multiple input files, these can be either fasta
                files with a set of sequences, or already multiple alignments
                in fasta format
  
   Optional parameters:
     -noalign : specify this instead of align, if sequences are already aligned
     -tmp     : by default, temporary files are stored in the currewnt working directory,
                this setting stores temporary files in the specified location instead
==============================================================================================
";

# OPTIONS
GetOptions(
	'noalign'    => \$noalign,
	'align:s{,}' => \@align,
	'h|help'     => \$manual,
	'bin:s'      => \$bin,
	'input:s{,}' => \@input_files,
	'tmp:s'      => \$cwd,
);

# STARTUP SCREEN
print "
==============================================================================================
                 fetchGENE v2.0 - extract selected gene groups from sequences                 
              Copyright (c) 2014 J. Roat Kultima, S. Sunagawa, D. R. Mende & EMBL             
==============================================================================================
";

# CHECKS
if ($manual) {
	pod2usage( -exitstatus => 0, -verbose => 2 );
	exit 0;
}
die " $usage ERROR & EXIT: Please specify either -noalign or -align [mafft muscle clustalo]" if ( !( $noalign || $align[0] ) || ( $noalign && $align[0] ) );
die " $usage ERROR & EXIT: Please specify input file(s) using -input" if !( $input_files[0] );

open my $JOB, ">", "$cwd/job.$date" or die " $usage ERROR & EXIT: Cannot write to $cwd/job.$date $!";

if ( $align[0] ) {

	# Set aligners
	foreach my $align (@align) {
		( $align eq 'mafft' || $align eq 'muscle' || $align eq 'clustalo' ) or die " $usage ERROR & EXIT: Specify -align [clustalo mafft muscle]";
		$align{$align} = 1;
	}
	$align{clustalo} = 0 unless $align{clustalo};
	$align{mafft} = 0 unless $align{mafft};
	$align{muscle} = 0 unless $align{muscle};

	foreach my $input (@input_files) {
		my ( $filename, $dir ) = fileparse($input);
		mkdir "$cwd/tmp.$date.$filename" or die " $usage ERROR & EXIT: Cannot create temporary directory $cwd/tmp.$date.$filename $!";
		print $JOB "$bin/$AQUA $input $cwd/tmp.$date.$filename $bin $align{clustalo} $align{mafft} $align{muscle}; ";
	}

}
else {
	print $usage;
	exit 0;
}

# EXIT
exit 0;

# MANUAL
__END__


=head1 Name

B<fetchGENE>

=head1 Version

2.0

=head1 Software description

fetchGENE extracts selected genes from genomes and metagenomes in an easy and accurate manner.
