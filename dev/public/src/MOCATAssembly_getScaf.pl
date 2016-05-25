#!/usr/bin/env perl
use warnings;
use strict;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

unless (scalar @ARGV == 3) {
  die "SCAF ERROR: Incorrect arguments. Check you code!\n";
}

my $min_length = $ARGV[0];
my $file = $ARGV[1];
my $sample = $ARGV[2];
my $name = '';
my $seq = '';

open IN, '<', $file;

while(<IN>){
    if(/^>(\S+)/){
	&print_scafftig($name, $seq) if($seq);
	$name = $1;
	$seq  = '';
    } else {
	chomp;
	$seq .= $_;
    }
}
&print_scafftig($name, $seq) if($seq);

1;

sub print_scafftig {
    my $name = shift;
    my $seq  = shift;
    my $id = 1;
    while($seq=~/([ATGCatgc]+)/g){
my $s = $1;
next if(length($s) < $min_length);
print ">$sample\_$name\_$id length=".length($s)."\n";
while($s=~/(.{1,60})/g){
print "$1\n";
}
$id ++;
    }
}
