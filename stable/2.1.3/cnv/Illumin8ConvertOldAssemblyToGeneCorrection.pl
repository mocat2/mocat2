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
  print ".";
  chomp(my $file = `ls -1 $sample/assembly.$mode.K*/$sample.scaftig`);
  $file =~ s/.*$sample\/assembly.$mode/\.\.\/assembly.$mode/;
  system "mkdir -p $sample/assembly.FROMOLDi8.$mode.K0/";
  system "mkdir -p $sample/stats";
  system "ln -s $file $sample/assembly.FROMOLDi8.$mode.K0/$sample.assembly.FROMOLDi8.$mode.K0.scaftig";
  system "echo -e \"Reads\\tBases\\tMax\\tAvg\\tKmer\\n0\\t0\\t0\\t0\\t0\" > $sample/stats/$sample.screen.FROMOLDi8.$mode.stats";
}
print " Done.\n";

exit 0;
