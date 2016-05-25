#!/usr/bin/env perl 
# 
# Reformat a fasta file 
#

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

use strict;
use Getopt::Long;

our $nucmer = "/share/raid1/genome/bin/nucmer";
our $mummerplot = "/share/raid1/genome/bin/mummerplot";

sub psystem { print ">> @_\n"; system "@_"; }

my ($locus, $length, $ref, $same);
GetOptions(
    "locus:s"  => \$locus,  # Output locus, like locus_1
    "length:s" => \$length,  # Output read length
    "ref:s"    => \$ref,
    "same"     => \$same,
);
$same = 1 unless $locus;
$length = 70 unless $length;

die "
  Usage: perl $0 OPTIONS <input.fa> <output.fa>
  OPTIONS:
    --locus  <locus>        - new locus, default is not set
    --length <rds_length>   - new read length, default 70
    --ref    <ref.fsa>      - set reference sequence, 
                              used to sort scaffolds, default is not set
    --same                  - reserve original title, default is on
\n" unless ($#ARGV == 1 and -e $ARGV[0]);

my ($infile, $outfile) = @ARGV;

open IN, "$infile" or die "open error! $infile $!\n";
open OUT, ">$outfile" or die "open error! $outfile $!\n";

print "## Reformating scaffold file $infile\n";
local $/ = '>';
my $count = 1; my %scafs;
while (<IN>) 
{
    chomp; s/\r//g;
    next unless $_;
    my ($t, $seq) = split(/\n/, $_, 2);
    my ($title, $t) = split(/\s+/, $t, 2);
    $seq =~ s/\n//g;
    
    if ( $same ){
        print OUT ">$title\n";
        $scafs{$title} = $seq;
    } else {
        print OUT ">${locus}${count}\n";
        $scafs{"${locus}${count}"} = $seq;
        $count ++;
    }

    while ($seq) {
        my $sub_seq = substr $seq, 0, $length;
        print OUT "$sub_seq\n";
        $seq = substr $seq, $length;
    }
}
close IN; close OUT;

# die "C22399 --> $scafs{C22399}\n";

local $/ = "\n";
if ($ref and -e $ref ) {
    my $ref_short = "ref.tmp.$$.fsa";
    my $prefix_out = "nucmer.tmp.$$";
    my @sort;
    psystem("echo '>Refseqs' >$ref_short");
    psystem("grep -v '>' $ref >>$ref_short");
    psystem("$nucmer -p $prefix_out $ref_short $outfile");
    psystem("$mummerplot -t png -l -p $prefix_out $prefix_out.delta >/dev/null 2>/dev/null");
    open IN, "$prefix_out.gp" or die "open $prefix_out.gp error, $!\n";
    open OUT, ">$outfile.sort";
    my ($t, $s);
    while ( <IN> ) {
        # push(@sort, $1) if m/^\s"(\S+)"\s\d+/;
        if ( m/^\s"(\S+)"\s\d+/ ) {
            $t = $1; 
            if ( $t =~ m/^\*/ ) {
                $t = substr($t, 1);
                $s = reverse $scafs{$t};
                $s =~ tr/atcgATCG/TAGCTAGC/;
                $t .= "_rc";
            } else {
                $s = $scafs{$t};
            }
            print OUT ">$t\n";
            while ($s) {
                my $ss = substr $s, 0, $length;
                print OUT "$ss\n";
                $s = substr $s, $length;
            }
        } 
    }
    # print "@sort \n";
    psystem("rm $prefix_out.{c,f,r}* $ref_short");
    close IN; close OUT;
} else {
    print "## Reference not set, or $ref not exist, job finished.\n\n";
}

