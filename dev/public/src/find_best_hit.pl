#! /usr/bin/env perl

use strict;
use warnings;
use Sys::Hostname;
use POSIX;

#print STDERR "## $0 @ARGV ##\n";

print "##\n";
print "## $0 @ARGV\n"; 
print "##\n";
print "## "; print scalar localtime; print "\n";
print "## "; print hostname; print " - "; print `pwd`; 
print "##\n";

######################################################################################################################
## Original author: Christian von Mering
## Modifications  : Jeroen Raes, Mani Arumugam
##
## This script will reduce the BLAST summary to one best hit per 'region'. 
##
## parameters for this script:
##
## -i input_file       a file from which to parse the information [mandatory].
## -excl_sh            setting this flag will not allow self-hits to be reported as best hit [optional].
## -excl_100h          setting this flag will not allow 100% identical hits to be reported as best hit [optional]. 
## -maxic fraction     'maximum identity cutoff' specify here exactly how low-identity a hit has to be before it is
##                     reported. 
## -match word         if provided, select the subset of sequences from the input set which maps the word [optional].
## -f blast-flavor     flavor of BLASTP used to generate the alignment; must be WU or NCBI [mandatory].
######################################################################################################################

my $param_input_file = undef;
my $param_exclude_self_hit = 0;
my $param_exclude_100_hit = 0;
my $param_match_word = undef;
my $param_maxic_fraction = 10000;
my $param_blast_flavor = undef;

my $param = shift @ARGV;
while ($param) {
    unless ($param =~ /\A-\w/) { die "illegal parameter '$param'\n"; }
    if ($param eq "-i") {
	$param_input_file = shift @ARGV or die "  error:  -i parameter needs value!\n";
	$param = shift @ARGV;
	next;
    }
    if ($param eq "-excl_sh") {
	$param_exclude_self_hit = 1;
	$param = shift @ARGV;
	next;
    }
    if ($param eq "-excl_100h") {
	$param_exclude_100_hit = 1;
	$param = shift @ARGV;
	next;
    }
    if ($param eq "-match") {
	$param_match_word = shift @ARGV or die "  error; -match parameter needs value!\n";
	$param = shift @ARGV;
	next;
    }
    if ($param eq "-maxic") {
	$param_maxic_fraction = shift @ARGV or die "  error; -maxic parameter needs value!\n";
	$param = shift @ARGV;
	next;
    }
    if ($param eq "-f") {
	$param_blast_flavor = shift @ARGV or die "  error; -f parameter needs value!\n";
	$param = shift @ARGV;
	next;
    }

    die "illegal parameter '$param'\n";
}

unless ($param_input_file) {
    print STDERR "error: need '-i' parameter. Read the Code.!\n";
    exit;
}

unless ($param_blast_flavor) {
    print STDERR "error: need '-f' parameter. Read the Code.!\n";
    exit;
}

my %best_hit_per_region_partner = ();
my %best_hit_per_region_score = ();
my %best_hit_per_region_rest_info = ();

open (FH, "$param_input_file") or die "cannot read file '$param_input_file'!\n";

while (<FH>) {

    chomp; next if /\A\#/;

    my ($region, $partner, $e_value, $bitscore, $identity);
    my @rest;
    if (uc($param_blast_flavor) eq "WU") {
        ($region, $partner, $e_value, undef, $bitscore, @rest) = split;
        $identity = $rest[5];
    } elsif (uc($param_blast_flavor) eq "NCBI") {
        ($region, $partner, $identity, @rest) = split;
        $e_value = $rest[7];
        $bitscore = $rest[8];
    } else {
        die "Unknown BLAST flavor: $param_blast_flavor\n";
    }

    next unless $identity <= $param_maxic_fraction;

    if ($param_exclude_self_hit) {
	next if $region eq $partner;
    }

    if ($param_exclude_100_hit) {
	next if $identity == 1.0;
    }

    if ($param_match_word) {
	next unless $region =~ /$param_match_word/;
    }

    my $rest_info = join " ", @rest;

    my $best_hit_partner = "unassigned";
    my $best_hit_rest_info = "unassigned";
    my $best_hit_score = 0;

    $best_hit_partner = $best_hit_per_region_partner{$region} if exists $best_hit_per_region_partner{$region};
    $best_hit_score = $best_hit_per_region_score{$region} if exists $best_hit_per_region_score{$region};
    $best_hit_rest_info = $best_hit_per_region_rest_info{$region} if exists $best_hit_per_region_rest_info{$region};

    if ($bitscore > $best_hit_score) {
	
	$best_hit_partner = $partner;
	$best_hit_score = $bitscore;
	$best_hit_rest_info = $rest_info;
    }

    $best_hit_per_region_partner{$region} = $best_hit_partner;
    $best_hit_per_region_score{$region} = $best_hit_score;
    $best_hit_per_region_rest_info{$region} = $best_hit_rest_info;
}

foreach my $region (sort {$best_hit_per_region_score{$a} <=> $best_hit_per_region_score{$b}} keys %best_hit_per_region_partner) {
    
    my $best_hit = $best_hit_per_region_partner{$region};
    my $best_score = $best_hit_per_region_score{$region};
    my $best_rest_info = $best_hit_per_region_rest_info{$region};

    print "$region $best_hit $best_score $best_rest_info\n";
}

