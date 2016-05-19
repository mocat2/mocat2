#!/usr/bin/perl
use strict;
use warnings;

# Converts old i8 assembly format into new i8 assembly format.
# Run from the current folder containing sample folders.
# Usage: X.pl <SAMPLE FILE> <DATA TYPE: solexaqa, fastx>

unless (exists $ARGV[0] && exists $ARGV[1]) {
  die "Missing sample file and data type.\n";
}

print "Loading samples...\n";
my @samples;
my $mode = $ARGV[1];
open IN, '<', $ARGV[0];
while (<IN>) {
  chomp;
  push @samples, $_;
}
close IN;

print "Processing samples";
foreach my $sample (@samples) {
  unless ($sample eq "") {
    print ".";

    # Get Kmer
    my $statsFile;
    $statsFile = "$sample/$sample.stats_read_trim_filter_$mode";
    if (-s "$statsFile" ) {
      open IL, '<', "$statsFile";
    } else {
      die "\nERROR & EXIT: Missing insert file $statsFile\n";	
    }      
    my $temp;
    while (<IL>) {
      chomp; $temp = $_;
    }
    close IL;
    my @temp = split(/\t/, $temp);
    my $kmer = $temp[4];

    chomp(my $file = `ls -1 $sample/assembly.$mode.K$kmer/$sample.scaftig`);
    $file =~ s/.*$sample\/assembly.$mode/\.\.\/assembly.$mode/;
    system "mkdir -p $sample/assembly.FROMOLDi8.$mode.K$kmer/";
    system "mkdir -p $sample/stats";
    system "mkdir -p $sample/reads.screened.FROMOLDi8.$mode";
    system "ln -s $file $sample/assembly.FROMOLDi8.$mode.K$kmer/$sample.assembly.FROMOLDi8.$mode.K$kmer.scaftig.500";

    chomp($file = `ls -1 $sample/assembly.$mode.K$kmer/$sample.scafSeq`);
    system "perl /g/bork4/kultima/GIT/MOCAT/public/src/MOCATAssembly_getScaf.pl 60 $file > $sample/assembly.FROMOLDi8.$mode.K$kmer/$sample.assembly.FROMOLDi8.$mode.K$kmer.scaftig";

    open IN, '<', "$sample/$sample.stats_assembly_insert_sizes";
    open OUT, '>', "$sample/stats/$sample.assembly.FROMOLDi8.$mode.K$kmer.inserts.stats";
    while (<IN>) {
      chomp $_;
      my @line = split ("\t", $_);
      print OUT "$line[0].screened.FROMOLDi8\t$line[1]\n";
    }
    close IN;
    close OUT;

    system "ln -s ../$sample.stats_read_trim_filter_$mode $sample/stats/$sample.screen.FROMOLDi8.$mode.stats";
    system " for f in $sample/reads.processed.$mode/*gz; do h=`echo \$f | sed 's/.*reads.processed.$mode\\///'`; g=`echo \$f | sed 's/.*reads.processed.$mode\\///' | sed 's/all.fq.gz/screened.FROMOLDi8.all.fq.gz/' | sed 's/.single/.screened.FROMOLDi8.single/' | sed 's/.pair/.screened.FROMOLDi8.pair/'`; ln -s ../reads.processed.$mode/\$h $sample/reads.screened.FROMOLDi8.$mode/\$g; done";
  }
}
print " Done.\n";

exit 0;
