#!/usr/bin/env perl

# This code is part of the MOCAT analysis pipeline
# This code is released under GNU GPL v3.

use strict;
use warnings;
my $i = 1;

chomp (my $cwd=`pwd`);


while (<>){
    my $command = $_;
    chomp $command;
    open (OUT,">RJSS");
    print OUT "#\!/bin/bash\n";
    print OUT "$command\n";
    close OUT;
    $i++;
}

