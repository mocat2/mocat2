#!/usr/bin/env perl

package Smash::Utils::Codon;
use strict;
use warnings;
use base "Exporter";

our @EXPORT_OK = qw($CodonTable $SynSites $NonSynSites $StartCodonTable are_synonymous);
our %EXPORT_TAGS = (all => [@EXPORT_OK]);

my $StartCodonTable  = {};
my $CodonTable  = {};
my $SynCodons   = {};
our $SynSites    = {};
our $NonSynSites = {};

sub init {
	my $translation_table = shift;

	$CodonTable = 
		{ 
		TTT=>"F", TTC=>"F", TTA=>"L", TTG=>"L", TTN=>"X",
		TCT=>"S", TCC=>"S", TCA=>"S", TCG=>"S", TCN=>"S",
		TAT=>"Y", TAC=>"Y", TAA=>"*", TAG=>"*", TAN=>"X",
		TGT=>"C", TGC=>"C", TGA=>"*", TGG=>"W", TGN=>"X",
		CTT=>"L", CTC=>"L", CTA=>"L", CTG=>"L", CTN=>"L",
		CCT=>"P", CCC=>"P", CCA=>"P", CCG=>"P", CCN=>"P",
		CAT=>"H", CAC=>"H", CAA=>"Q", CAG=>"Q", CAN=>"X",
		CGT=>"R", CGC=>"R", CGA=>"R", CGG=>"R", CGN=>"R",
		ATT=>"I", ATC=>"I", ATA=>"I", ATG=>"M", ATN=>"X",
		ACT=>"T", ACC=>"T", ACA=>"T", ACG=>"T", ACN=>"T",
		AAT=>"N", AAC=>"N", AAA=>"K", AAG=>"K", AAN=>"X",
		AGT=>"S", AGC=>"S", AGA=>"R", AGG=>"R", AGN=>"X",
		GTT=>"V", GTC=>"V", GTA=>"V", GTG=>"V", GTN=>"V",
		GCT=>"A", GCC=>"A", GCA=>"A", GCG=>"A", GCN=>"A",
		GAT=>"D", GAC=>"D", GAA=>"E", GAG=>"E", GAN=>"X",
		GGT=>"G", GGC=>"G", GGA=>"G", GGG=>"G", GGN=>"G"
		};

	map {$StartCodonTable->{$_} = $CodonTable->{$_}} keys %$CodonTable;

	if ($translation_table == 1) {
		map {$StartCodonTable->{$_} = "M"} qw(ATG);
	} elsif ($translation_table == 4) {
		map {$StartCodonTable->{$_} = "M"} qw(TTA TTG CTG ATT ATC ATA ATG GTG);
		$CodonTable->{TGA} = "W";
	} elsif ($translation_table == 11) {
		map {$StartCodonTable->{$_} = "M"} qw(TTG CTG ATT ATC ATA ATG GTG);
	}

	my @codons = keys %$CodonTable;
	foreach my $c1 (@codons) {
		foreach my $c2 (@codons) {
			if ($CodonTable->{$c1} eq $CodonTable->{$c2}) {
				$SynCodons->{$c1}->{$c2} = 1;
			} else {
				$SynCodons->{$c1}->{$c2} = 0;
			}
		}
	}

	# Count synonymous and nonsynonymous sites

	my @bases = qw(A C G T);
	foreach my $b1 (0..3) {
		foreach my $b2 (0..3) {
			foreach my $b3 (0..3) {
				my $codon = join("", @bases[$b1,$b2,$b3]);
				next if $CodonTable->{$codon} eq "*";
				my $syn_sites = 0;
				my $nonsyn_sites = 0;
				my $sub;
				foreach my $s1 (0..3) {
					next if $s1 == $b1;
					$sub = join("", @bases[$s1,$b2,$b3]);
					next if $CodonTable->{$sub} eq "*";
					if ($SynCodons->{$codon}->{$sub} == 1) {
						$syn_sites++;
					} else {
						$nonsyn_sites++;
					}
				}
				foreach my $s2 (0..3) {
					next if $s2 == $b2;
					$sub = join("", @bases[$b1,$s2,$b3]);
					next if $CodonTable->{$sub} eq "*";
					if ($SynCodons->{$codon}->{$sub} == 1) {
						$syn_sites++;
					} else {
						$nonsyn_sites++;
					}
				}
				foreach my $s3 (0..3) {
					next if $s3 == $b3;
					$sub = join("", @bases[$b1,$b2,$s3]);
					next if $CodonTable->{$sub} eq "*";
					if ($SynCodons->{$codon}->{$sub} == 1) {
						$syn_sites++;
					} else {
						$nonsyn_sites++;
					}
				}
				my $total = $syn_sites + $nonsyn_sites;
				$SynSites->{$codon}    = 3*$syn_sites/$total;
				$NonSynSites->{$codon} = 3*$nonsyn_sites/$total;
				#printf "%s: SYN:%2.2f,NON:%2.2f\n", $codon, $SynSites->{$codon}, $NonSynSites->{$codon};
			}
		}
	}
}

sub are_synonymous {
	my ($c1, $c2, $start_codon) = @_;
	($c1, $c2) = map {uc($_)} ($c1, $c2);
	if ($start_codon) {
		return ($StartCodonTable->{$c1} eq $StartCodonTable->{$c2});
	} else {
		return $SynCodons->{$c1}->{$c2};
	}
}

sub count_sites {
	my $gene = shift;

	die "Gene should be in codons for counting sites" if length($gene)%3;

	my ($S, $N) = (0, 0);
	$gene = uc($gene);
	while ($gene ne "") {
		my $codon = substr($gene, 0, 3);
		$S += ($SynSites->{$codon} || 0);
		$N += ($NonSynSites->{$codon} || 0);
		substr($gene, 0, 3) = "";
	}

	return ($S, $N);
}

1;
