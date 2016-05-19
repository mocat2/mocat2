#!/usr/bin/env perl

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

# this script is intended to take reads originating from converting colorspace (e.g. SOLiD) reads to base space
# and filter them on account of the fact that the converted sequence is invalid at any position following an 
# ambiguous assignment; this may shorten the read, and if shortened too much, remove it entirely

use Getopt::Long;

# DIRS

my $indir; # input dir
my $originaldir; # store originals dir
my $outdir; # output dir - might be same as input dir
my $scoreThreshold; # what score is considered unassigned?
my $lengthThreshold; # how short can the trimmed read be for us to still keep it?

# "perl filterConvertedColorspaceReads.pl -indir=<some dir> -outdir=<some dir> -originaldir=<some dir> -scoreThreshold=<some value> -lengthThreshold=<some value>"

GetOptions (
                "indir=s" => \$indir, 
                "originaldir=s" => \$originaldir,
                "outdir=s" => \$outdir,
                "scoreThreshold=s" => \$scoreThreshold,
                "lengthThreshold=s" => \$lengthThreshold
);

# code for checking whether a character has quality enough for inclusions

my $alphabet = "!\"#\$%&'()*+,-./0123456789:;<=>?\@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

sub checkAllowed {

	my $char = $_[0];

	if (index ($alphabet, $char) > index ($alphabet, $scoreThreshold)) {

		return 1;

	}

	else {

		return 0;

	}

};

# make sure directories exist

if (! (-e $originaldir && -d $originaldir)) {

	system ("mkdir $originaldir");

}

if (! (-e $outdir && -d $outdir)) {

	system ("mkdir $outdir");

}

# backup files

system ("mv $indir/*.fq $originaldir/ 2>/dev/null");
system ("mv $indir/*.fq.gz $originaldir/ 2>/dev/null");

# process all zipped files

foreach $file (<$originaldir/*.fq.gz>) {

	my $name = "";
	my $newname = "";

	if ($file =~ /([^\/]+\.fq)\.gz$/) {
		
		$name = $1;
		$newname = $name;

		if ($name =~ /^(\S+)\.([^\.]+)\.fq$/) {

			$newname = "$1.SOLiD_trimmed.$2.fq";
	
		}

		elsif ($name =~ /^(\S+)\.fq$/) {

			$newname = "$1.SOLiD_trimmed.fq";
	
		}

	}

	else {

		die ("Issue parsing file name from $file\n");

	}

	open (FH, "gzip -dc $file |");
	open (FH2, " | gzip -2 -c > $outdir/$newname.gz");

	my $line = 1;

	my $id, $seq, $tag, $qual;

	while (<FH>) {

		$aLine = $_;
		chomp ($aLine);

		if ($aLine =~ /^\#/) {

			next;

		}

		if ($line % 4 == 1) {

			# line 1 of post

			$id = $aLine;

		}

		if ($line % 4 == 2) {

			# line 2 of post

			$seq = $aLine;

		}

		if ($line % 4 == 3) {

			# line 3 of post

			$tag = $aLine;

		}

		if ($line % 4 == 0) {

			# line 4 of post

			$qual = $aLine;

			# then, process!

			my $i;

			for ($i = 0; $i < length ($qual); $i++) {

				if (! checkAllowed (substr ($qual, $i, 1))) {

					last;

				}

			}

			# i is now position of first non-OK value

			$qual = substr ($qual, 0, $i);
			$seq = substr ($seq, 0, $i);

			if (length ($qual) >= $lengthThreshold) {

				print FH2 "$id\n$seq\n$tag\n$qual\n";

			}

		}

		$line++;

	}

	close (FH2);
	close (FH);

}

# done

# process all non-zipped files

foreach $file (<$originaldir/*.fq>) {

	my $name = "";
	my $newname = "";

	if ($file =~ /([^\/]+\.fq)$/) {
		
		$name = $1;
		$newname = $name;

		if ($name =~ /^(\S+)\.([^\.]+)\.fq$/) {

			$newname = "$1.SOLiD_trimmed.$2.fq";
	
		}

		elsif ($name =~ /^(\S+)\.fq$/) {

			$newname = "$1.SOLiD_trimmed.fq";
	
		}

	}

	else {

		die ("Issue parsing file name from $file\n");

	}

	open (FH, $file);
	open (FH2, " | gzip -2 -c > $outdir/$newname.gz");

	my $line = 1;

	my $id, $seq, $tag, $qual;

	while (<FH>) {

		$aLine = $_;
		chomp ($aLine);

		if ($aLine =~ /^\#/) {

			next;

		}

		if ($line % 4 == 1) {

			# line 1 of post

			$id = $aLine;

		}

		if ($line % 4 == 2) {

			# line 2 of post

			$seq = $aLine;

		}

		if ($line % 4 == 3) {

			# line 3 of post

			$tag = $aLine;

		}

		if ($line % 4 == 0) {

			# line 4 of post

			$qual = $aLine;

			# then, process!

			my $i;

			for ($i = 0; $i < length ($qual); $i++) {

				if (! checkAllowed (substr ($qual, $i, 1))) {

					last;

				}

			}

			# i is now position of first non-OK value

			$qual = substr ($qual, 0, $i);
			$seq = substr ($seq, 0, $i);

			if (length ($qual) >= $lengthThreshold) {

				print FH2 "$id\n$seq\n$tag\n$qual\n";

			}

		}

		$line++;

	}

	close (FH2);
	close (FH);

}
