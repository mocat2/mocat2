#!/usr/bin/env perl
use warnings;
use strict;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

unless(@ARGV == 2){
    &usage();
    exit;
}

my($contig,$lcutoff_list) = @ARGV;
my @lcutoff = split /:/,$lcutoff_list;

my(@seqlen,@sum,@num,@max,@min,@n50,@n90);
$#seqlen = $#sum = $#num = $#max = $#min = $#n50 = $#n90 = $#lcutoff;
for(my $i=0;$i<@lcutoff;$i++){
    $sum[$i] = 0;
    $num[$i] = 0;
    $max[$i] = 0;
    $min[$i] = 100000000;
    $n50[$i] = 0;
    $n90[$i] = 0;
}

my $seq = '';
open CFA, $contig or die "$!\n";
while (<CFA>) {
    chomp;
    if (/^>/) {
	&foo if ($seq ne '');
	$seq = '';
    } else {
	$seq .= $_;
    }
}
&foo if ($seq ne '');

for(my $i=0;$i<@seqlen;$i++){
    my @temp = reverse sort {$a <=> $b} @{$seqlen[$i]};
    my $tsum=0;
    foreach my $l(@temp){
        $tsum += $l;
        if($tsum >= $sum[$i] / 2){
            $n50[$i] = $l;
            last;
        }
    }
    $tsum = 0;
    foreach my $l(@temp){
        $tsum += $l;
        if($tsum >= $sum[$i]*0.9){
            $n90[$i] = $l;
            last;
        }
    }
}

print STDOUT "ctg_num\tctg_length_all\tctg_n50\tctg_n90\tctg_max\tctg_min\n";
for(my $i=0;$i<@num;$i++){
    print STDOUT "$num[$i]\t$sum[$i]\t$n50[$i]\t$n90[$i]\t$max[$i]\t$min[$i]\n";
}

sub foo {
    my $len = length $seq;
    for(my $i=0;$i<@lcutoff;$i++){
        next if($len < $lcutoff[$i]);
        push @{$seqlen[$i]},$len;
        $sum[$i] += $len;
        $num[$i]++;
        $max[$i] = $len if($max[$i] < $len);
        $min[$i] = $len if($min[$i] > $len);
    }
}


sub usage{
    print <<EOD
    usage: perl $0 contig.fa lcutoff_list
    eg: perl $0 K23_L90.scafSeq 0:200:500:1000:2000:4000
EOD
}
