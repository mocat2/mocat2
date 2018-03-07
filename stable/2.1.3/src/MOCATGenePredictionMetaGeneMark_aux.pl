#!/usr/bin/env perl

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

use strict;
use warnings;
use MOCATCore;

my %conf;
my $sample                      = $ARGV[0];
$conf{MOCAT_data_type}          = $ARGV[1];
my $assembly    = $ARGV[2];
my $reads       = $ARGV[3];
my $cwd                         = $ARGV[4];
$conf{gene_prediction_input}    = $ARGV[5];
my $kmer                        = $ARGV[6];
my $GP_min_len                  = $ARGV[7];

# Create job
my $lst = "$cwd/$sample/gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$GP_min_len/$sample.gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.MetaGeneMark.lst";
my $nt = "$cwd/$sample/gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$GP_min_len/$sample.gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.MetaGeneMark.fna";
my $aa = "$cwd/$sample/gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$GP_min_len/$sample.gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.MetaGeneMark.faa";
my $scaftig = "$cwd/$sample/$assembly.$reads.$conf{MOCAT_data_type}.K$kmer/$sample.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$conf{gene_prediction_input}";
my $sum = "$cwd/$sample/gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.$GP_min_len/$sample.gene.prediction.$assembly.$reads.$conf{MOCAT_data_type}.K$kmer.MetaGeneMark.tab";

open IN, '<', $lst or die "ERROR & EXIT: Missing $lst\nDid the Gene prediction step complete?";
open AA, '>', $aa or die "ERROR & EXIT: Cannot open output file $aa for gene prediction";
open NT, '>', $nt or die "ERROR & EXIT: Cannot open output file $nt for gene prediction";
open SUM, '>', $sum or die "ERROR & EXIT: Cannot open output file $sum for gene prediction";
my %genes;
my %genesSUM;
    
my $section = -1;
while (<IN>) {
  chomp;
  if ($_ =~ /^FASTA definition line/) {
    $section = -1;
    next;
  }
  if ($_ =~ /^Predicted genes/) {
    $section = 0;
    %genes = ();
    next;
  }
  if ($_ =~ /^Predicted proteins:/) {
    $section = 1;
    next;
  }
  if ($_ =~ /^Nucleotide sequence of predicted genes:/) {
    $section = 2;
    next;
  }
  if ($section == 0) {
    if ($_ =~ /^\s*(\d+)\s+(\+|-)\s+(<*)(\d+)\s+(>*)(\d+)\s+(\d+)\s+/) {
      my ($gene_id, $strand, $langle, $start, $rangle, $end, $length) = ($1, $2, $3, $4, $5, $6, $7);
      my $incomplete5  = 0;
      my $incomplete3  = 0;
      my $incomplete = 0;
      my ($start_codon, $stop_codon);
      if ($langle eq "<") {
	$incomplete = 1;
	$incomplete5 = 1;
      }
      if ($rangle eq ">") {
	$incomplete = 1;
	$incomplete3 = 1;
      }
      if ($strand eq "+") {
	$start_codon = ($incomplete5 == 1)?"no":"yes";
	$stop_codon  = ($incomplete3 == 1)?"no":"yes";
      } else {
	$stop_codon  = ($incomplete5 == 1)?"no":"yes";
	$start_codon = ($incomplete3 == 1)?"no":"yes";
      }
      my $gene;
      my $complete = ($incomplete == 1)?"incomplete":"complete";
      $gene = "PoPIiI_gene$1 strand:$2 start:$4 stop:$6 length:$7 start_codon:$start_codon stop_codon:$stop_codon gene_type:$complete";
      my $geneSUM = "PoPIiI_gene$1\t$2\t$4\t$6\t$7\t$start_codon\t$stop_codon\t$complete";
      $genes{$1} = $gene;
      $genesSUM{$1} = $geneSUM;
    }
  }
  if ($section == 1) {
    if ($_ =~ /^\>gene\_(\d+)\|.*>(\S+)\s+length.*/) {
      my $n2 = $2;
      my $n1 = $1;  
      $genes{$n1} =~ s/PoPIiI/$n2/;
      $genesSUM{$n1} =~ s/PoPIiI/$n2/;
      print SUM "$genesSUM{$n1}\n";
      print AA ">$genes{$n1}\n";
    } else {
      if ($_ =~ /^\w+$/) {
	print AA "$_\n";
      }	  
    }
  }
  if ($section == 2) {
    if ($_ =~ /^\>gene\_(\d+)\|.*>(\S+)\s+length.*/) {
      $genes{$1} =~ s/PoPIiI/$2/;
      print NT ">$genes{$1}\n";
    } else {
      if ($_ =~ /^\w+$/) {
	print NT "$_\n";
      }	  
    }
  }
}
  
close AA;
close NT;
close SUM;
close IN;

exit 0;
