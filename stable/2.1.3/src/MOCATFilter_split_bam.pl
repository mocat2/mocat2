#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

# this script takes an input pipe and splits it into multiple gzipped sam files

my $INSERTS = 300000000; # among largest TARA samples has 215042748 inserts, we just don't want to split into more than specified files
my $ZIP = "pigz -p 8 -4";
my $OUT = "OUT";
my %IDs;
my $files = 4;

GetOptions(
	'zip=s' => \$ZIP,
	'inserts=s' => \$INSERTS, 
	'out=s' => \$OUT,
	'files=s' => \$files
);

print STDERR "MOCATFilter_split_bam.pl (".getLoggingTime().") : starting [zip='$ZIP' inserts=$INSERTS out=$OUT files=$files]\n";


my $line;
my $prev_insert = '';
my $fileNumber = 1;
my $split_at = $INSERTS/$files * $fileNumber;
my $insertCounter = 0;
open my $OUTFILE, " | $ZIP > $OUT.part.$fileNumber.sam.gz" or die "ERROR & EXIT: Cannot write to: | $ZIP > $OUT.part.$fileNumber.sam.gz";
while (<STDIN>) {
	$line = $_;
	my @line = split "\t", $line;
	$line[0] =~ s/\/[0-9]$//;
	if ($insertCounter <= $split_at || ($insertCounter > $split_at && $line[0] eq $prev_insert) ) {
		$IDs{$line[2]} = undef;
		print $OUTFILE $line;
	} else {
		close $OUTFILE;
		open my $OUTFILEr, " | $ZIP > $OUT.part.$fileNumber.sam.gz.references.gz" or die "ERROR & EXIT: Cannot write to: | $ZIP > $OUT.part.$fileNumber.sam.gz.references.gz";
		while ((my $key, my $value) = each %IDs) {
			print $OUTFILEr "$key\n";
		}
		close $OUTFILEr;
		%IDs = ();
		$fileNumber++;
		open $OUTFILE, " | $ZIP > $OUT.part.$fileNumber.sam.gz" or die "ERROR & EXIT: Cannot write to: | $ZIP > $OUT.part.$fileNumber.sam.gz";
		print $OUTFILE $line;
		$split_at = $INSERTS/$files * $fileNumber;
	}
	if ($prev_insert ne $line[0]) {
		$insertCounter++;
	}
	$prev_insert = $line[0];
}
open my $OUTFILEr, " | $ZIP > $OUT.part.$fileNumber.sam.gz.references.gz" or die "ERROR & EXIT: Cannot write to: | $ZIP > $OUT.part.$fileNumber.sam.gz.references.gz";
while ((my $key, my $value) = each %IDs) {
	print $OUTFILEr "$key\n";
}
close $OUTFILEr;
close $OUTFILE;

print STDERR "MOCATFilter_split_bam.pl (".getLoggingTime().") : done \n";

sub getLoggingTime {

	my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) = localtime(time);
	my $nice_timestamp = sprintf( "%04d%02d%02d %02d:%02d:%02d", $year + 1900, $mon + 1, $mday, $hour, $min, $sec );
	return $nice_timestamp;
}

exit 0;
