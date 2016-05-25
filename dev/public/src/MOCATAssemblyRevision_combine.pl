#!/usr/bin/env perl
use warnings;
use strict;

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

die "Usage $0 [in1] [in2] [out]" unless @ARGV==3;
open IN,$ARGV[0] or die "$!";
open IN2,$ARGV[1] or die "$!";
open OT,">$ARGV[2]" or die "$!";
print OT "Sample\t#Single base error\t#Small insertion\t#Total length of insertion\t#Small deletion\t#Total length of deletion\t#Ridiculous region\t#Total length of ridiculous region\n";
<IN>;
<IN2>;
while(<IN>){
    chomp;
    chomp (my $info=<IN2>);
    my @info=(split /\s+/,$_)[0,1,2,3,4,5];
    my $info2=join "\t",@info;
    my @info1=(split /\s+/,$info)[1,2];
    my $info3=join "\t",@info1;
    print OT $info2,"\t",$info3,"\n";
}
close IN;
close IN2;
close OT;

