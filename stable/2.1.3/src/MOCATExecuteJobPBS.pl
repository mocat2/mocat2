#!/usr/bin/env perl

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL, 2012-2016
# This code is released under GNU GPL v3.

use strict;
use warnings;
my $i = 1;

while (<>){
    my $command = $_;
    chomp $command;
    open (OUT,">RJSS");
    print OUT     "#\!/bin/sh \n";    
    print OUT     "#PBS -cwd \n";
    print OUT     "#\$ -o $i.log \n";
    print OUT     "#\$ -e $i.error \n";
    print OUT     "#\$ -N MOCATJob_$i \n";
    print OUT     "$command \n";
    close OUT;
    $i++;
}
