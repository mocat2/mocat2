#!/usr/bin/env perl

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

use strict;
use warnings;

print STDERR "GP:Prodigal:aux: starting script\n";

my $nt_tmp = $ARGV[0];
my $aa_tmp = $ARGV[1];
my $sum    = $ARGV[2];

my $aa = $aa_tmp;
$aa =~ s/.tmp$//;
my $nt = $nt_tmp;
$nt =~ s/.tmp$//;
open my $AAt, "<$aa_tmp" or die "ERROR & EXIT: Missing $aa_tmp\nDid the Gene prediction step complete?";
open my $NTt, "<$nt_tmp" or die "ERROR & EXIT: Missing $nt_tmp\nDid the Gene prediction step complete?";
open my $SUM, '>', $sum or die "ERROR & EXIT: Missing $sum";
open my $AA, '>', $aa or die "ERROR & EXIT: Missing $aa\nDid the Gene prediction step complete?";
open my $NT, '>', $nt or die "ERROR & EXIT: Missing $nt\nDid the Gene prediction step complete?";

print STDERR "GP:Prodigal:aux: opened files, processing...\n";

for my $var ($AAt, $NTt) {
  my $varout;
  if ($var == $AAt) {
    $varout = $AA;
  } else {
    $varout = $NT;
  }
  while (<$var>) {
    if (m/^>/) {
      chomp;
      my @line = split(/ /, $_);
      $line[0] =~ s/\_(\d+)\_(\d+)$/_gene$1_$2/;
      my $strand;
      if ($line[6] eq '1') {
	$strand = '+';
      } elsif ($line[6] eq '-1') {
	$strand = '-';
      } else {
	die "STRAND ERROR...\n";
      }
      my $start = $line[2];
      my $stop = $line[4];
      my $length = $stop - $start + 1;
      if ($length <= 0) {
	die "LENGTH PARSING ERROR...\n";
      }
      my $start_codon;
      my $stop_codon;
      my $complete;
      if ($line[8] =~ m/partial=00/) {
	$start_codon = 'yes';
	$stop_codon = 'yes';
	$complete = "complete";
      } elsif ($line[8] =~ m/partial=10/) {
	$start_codon = 'no';
	$stop_codon = 'yes';
	$complete = "incomplete";
      } elsif ($line[8] =~ m/partial=01/) {
	$start_codon = 'yes';
	$stop_codon = 'no';
	$complete = "incomplete";
      } elsif ($line[8] =~ m/partial=11/) {
	$start_codon = 'no';
	$stop_codon = 'no';
	$complete = "incomplete";
      } else {
	die "GENE COMPLETENESS PARSING ERROR...\n";
      }
      print $varout "$line[0] strand:$strand start:$start stop:$stop length:$length start_codon:$start_codon stop_codon:$stop_codon gene_type:$complete\n";
      if ($varout == $NT) {
	$line[0] =~ s/^>//;
	print $SUM "$line[0]\t$strand\t$start\t$stop\t$length\t$start_codon\t$stop_codon\t$complete\n";
      }
    } else {
      print $varout $_;
    }    
  }
}
    
close $AA;
close $NT;
close $AAt;
close $NTt;
close $SUM;

print STDERR "GP:Prodigal:aux: processed and completed files. deleting $aa.tmp and $nt.tmp\n";
system "rm -f $nt.tmp $aa.tmp";

exit 0;
