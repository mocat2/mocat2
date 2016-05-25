#!/usr/bin/env perl

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

# Takes an abundance tabke as input and produces fractions file as output

# Usage: MOCATFraction.pl -in INFILE -out OUTFILE [-colnames 1 -rownames 1 -sep '\t']
# NOTE: colnames and rownames are number of header rows and number of rownames rows (rownames is always 1 I suppose), and this number
# excludes the number of comment lines starting with #. This means, even though MOCAT tables has 4 comment lines, colnames is still 1.

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my ( $in, $out, @sum );
my $rownames = 1;
my $colnames = 1;
my $sep      = "\t";

GetOptions(
	'in:s'       => \$in,
	'out:s'      => \$out,
	'colnames:s' => \$colnames,
	'rownames:s' => \$rownames,
	'sep:s'      => \$sep
);

open IN, "<$in" or die "ERROR & EXIT: Cannot read $in";
my $counter = 0;
while (<IN>) {
	chomp;
	unless (m/^#/) {
		$counter++;
		if ( $counter > $colnames ) {
			my @line = split $sep, $_;
			for my $i ( $rownames .. scalar @line - 1 ) {
				unless ( $line[$i] eq 'NA' ) {
					if ( looks_like_number( $line[$i] ) ) {
						$sum[$i] += $line[$i];
					}
				}
			}
		}
	}
}
close IN;

open IN,  "<$in"  or die "ERROR & EXIT: Cannot read $in";
open OUT, ">$out" or die "ERROR & EXIT: Cannot write $out";
$counter = 0;
while (<IN>) {
	chomp;
	my $final_line = "";
	if (m/^#/) {
		$final_line = $_;
		if ( $final_line =~ m/^# MOCAT/ ) {
			$final_line = $final_line . " (fractions)";
		}
	}
	elsif ( $counter < $colnames ) {
		$counter++;
		$final_line = $_;
	}
	else {
		my @line = split $sep, $_;
		my @tokens = ();
		for my $i ( 0 .. $rownames - 1 ) {
			push @tokens, $line[$i];
		}
		for my $i ( $rownames .. scalar @line - 1 ) {
			if ( ( $line[$i] eq 'NA' ) or ( $line[$i] eq 'not_calculated' ) ) {
				push @tokens, 'NA';
			}
			else {
				if ($sum[$i] > 0) {
					push @tokens, ( $line[$i] / $sum[$i] );
				} else {
					push @tokens, ( 0 );
				}
			}
		}
		$final_line = join $sep, @tokens;
	}
	print OUT "$final_line\n";
}
close IN;
close OUT;

exit 0;
