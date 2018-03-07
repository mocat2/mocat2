#!/usr/bin/env perl
#usage: getFastaRecordLengths.pl --infile=[fasta file] --outfile=[location and name of outfile]
use strict;
use warnings;
use Getopt::Long;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

my $infile  = "";
my $outfile = "";
my $zip;

GetOptions(
	"infile=s"  => \$infile,
	"outfile=s" => \$outfile,
	"zip"       => \$zip
);

if ($zip) {
	open( INFILE, "gunzip -dc $infile | " ) or die "Cannot open the data file for data\n";
}
else {
	open( INFILE, "<$infile" ) or die "Cannot open the data file for data\n";
}
open( OUTFILE, ">$outfile" ) or die "Cannot open the data file for data output\n";

my $seqName;
my $seqLen = 0;

while ( my $currentFastaLine = <INFILE> ) {
	chomp($currentFastaLine);

	if ( $currentFastaLine =~ /^>/ ) {
		if ( $seqLen > 0 ) {
			print OUTFILE "$seqName\t$seqLen\n";
			$seqLen = 0;
		}
		$seqName = substr( $currentFastaLine, 1, length($currentFastaLine) - 1 );
		my @seqNameArray = split( /\s+/, $seqName );
		$seqName = $seqNameArray[0];
	}
	else {
		$seqLen += length($currentFastaLine);
	}
}
print OUTFILE "$seqName\t$seqLen\n";

