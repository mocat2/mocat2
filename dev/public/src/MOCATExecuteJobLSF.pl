#!/usr/bin/env perl

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

use strict;
use warnings;
my $i = 1;

chomp (my $cwd=`pwd`);

while (<>){
    my $command = $_;
    chomp $command;
    open (OUT,">RJSS");
    print OUT     "#\!/bin/sh \n";    
    print OUT     "#BSUB -o $cwd/$i.log \n";
    print OUT     "#BSUB -e $cwd/$i.error \n";
    print OUT     "#BSUB -J MOCATJob_$i \n";
    print OUT     "$command \n";
    close OUT;
    $i++;
}
