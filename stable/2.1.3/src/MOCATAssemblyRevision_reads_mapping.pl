#!/usr/bin/env perl
use warnings;

#
# Run reads mapping automatically
#

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

use strict;
use Getopt::Long;
use File::Basename;

#sub psystem { print ">> @_\n"; system("@_"); }
sub psystem {
	my $cmd = shift;
	print STDERR "SYSTEM CALL (reads_mapping): $cmd";
	system($cmd) == 0 or die("\nERROR & EXIT: system($cmd) failed: $!\n");
}

our ( $ref, $qs, $q1, $q2, $prog, $ins, $mis, $indel, $hits, $outdir, $Bin, $processors ) = ( "", "", "", "", "bwa", 500, 3, 2, 3, "bwa.aln", "", "8" );
my $cmdline = `date && pwd` . "\$ perl $0 \\\n @ARGV\n";

GetOptions(
	"ref:s"        => \$ref,
	"q1:s"         => \$q1,
	"q2:s"         => \$q2,
	"qs:s"         => \$qs,
	"prog:s"       => \$prog,
	"ins:i"        => \$ins,
	"mis:i"        => \$mis,
	"indel:i"      => \$indel,
	"hits:i"       => \$hits,
	"outdir:s"     => \$outdir,
	"bindir:s"     => \$Bin,
	"processors:s" => \$processors

);

our $bwa = "$Bin/bwa";
chomp( our $systemType = `uname -s` );
if ( $systemType =~ m/Darwin/ ) {
	$bwa = "$Bin/bwa_OSX";
}
our $bwt         = "$Bin/2bwt-builder";
our $soap        = "/$Bin/soap2.21";
our $bwt_builder = "$Bin/2bwt-builder";
our $samtools    = "$Bin/samtools";

die "
reads_mapping.pl - run reads mapping automatically, soap or bwa.

Usage: perl $0 \\
        --ref <ref.fsa> --q1 </path/to/q1> --q2 </path/to/q2> --qs </path/to/qs> [OPTIONS]
    
Resource Requirement: 2 cpu, 5G memory.
Options:
    --prog   STR - set program to run, soap or bwa [$prog]
    --ins    INT - set insert size, needed when running soap [$ins]
    --mis    INT - set allow mis-matches [$mis]
    --indel  INT - set allow indels [$indel]
    --hits   INT - set maximum hits to output [$hits]
    --outdir STR - set output directory [$outdir]
    --bindir
    \n" unless ( -f $ref and -f $q1 and -f $q2 );

psystem("mkdir $outdir") unless ( -d $outdir );
system("echo '$cmdline' >>$outdir/command_line.save");

&bwt  if ( $prog eq "bwa" );
&soap if ( $prog eq "soap" );

##
# run bwt
sub bwt {
	my ( $rq1, $rq2, $rqs );
	$rq1 = basename $q1;
	$rq2 = basename $q2;
	$rqs = basename $qs;
	$rq1 =~ s/\.gz$//;
	$rq2 =~ s/\.gz$//;
	$rqs =~ s/\.gz$//;

	psystem("$bwa index $ref") unless ( -e "$ref.bwt" );
	my $aln = "$bwa aln -e 63 -L -i 15 -o 1 -I";

	my $sampe = "$outdir/" . substr( $rq1, 0, length($rq1) - 5 ) . ".sampe";
	my $samse = "$outdir/" . substr( $rqs, 0, length($rqs) - 3 ) . ".samse";

	psystem("$aln -t $processors $ref $q1 > $outdir/$rq1.sai");
	psystem("$aln -t $processors $ref $q2 > $outdir/$rq2.sai");
	psystem("$bwa sampe -n $hits $ref $outdir/$rq1.sai $outdir/$rq2.sai $q1 $q2 | $samtools view -b -S -T $ref - | $samtools sort -m 3500000000 - $sampe");

	#psystem("$samtools view -b -S -T $ref $sampe | $samtools sort -m 3500000000 - $sampe");

	psystem("$aln -t $processors $ref $qs | $bwa samse -n $hits $ref - $qs | $samtools view -b -S -T $ref - | $samtools sort -m 3500000000 - $samse");

	#    psystem("$bwa samse -n $hits $ref $outdir/$rqs.sai $qs > $samse");
	#   psystem("$samtools view -b -S -T $ref $samse | $samtools sort -m 3500000000 - $samse");

	#   psystem("$aln -t $processors $ref $q1 | $bwa samse -n $hits $ref - $q1 | $samtools view -b -S -T $ref -o - - | $samtools sort -m 3500000000 - $result1");
	#   psystem("$aln -t $processors $ref $q2 | $bwa samse -n $hits $ref - $q2 | $samtools view -b -S -T $ref -o - - | $samtools sort -m 3500000000 - $result2");
	#   psystem("$aln -t $processors $ref $qs | $bwa samse -n $hits $ref - $qs | $samtools view -b -S -T $ref -o - - | $samtools sort -m 3500000000 - $resultS");
}

# run soap
sub soap {

	#die "Not implemented yet!\n";
	my ( $rq1, $rq2 );
	$rq1 = basename $q1;
	$rq2 = basename $q2;
	$rq1 =~ s/\.gz$//;
	$rq2 =~ s/\.gz$//;

	psystem("$bwt_builder $ref") unless ( -e "$ref.index.bwt" );

	my $min_ins = $ins - 0.2 * $ins;
	my $max_ins = $ins + 0.2 * $ins;

	my $varg;
	if ( $hits == 0 ) { $varg = 0; }
	if ( $hits == 1 ) { $varg = 1; }
	if ( $hits >= 2 ) { $varg = 2; }

	my $result = "$outdir/" . substr( $rq1, 0, length($rq1) - 5 );
	psystem( $soap, "-a $q1 -b $q2 -D $ref.index -o $result.pe -2 $result.se -m $min_ins -x $max_ins -v $mis -g $indel -r $hits" );
}

